"""Strict own64/400 B1 pole allocation for separately authorized DZ or TZ data."""

from __future__ import annotations

import argparse
import json
import multiprocessing
import os
import signal
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from qscat.units import hartree_to_ev

from projects.ukrmol_co import compressed_pole_reference as core
from projects.ukrmol_co.compressed_a2_contrast_run import save_json
from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.outgoing_map import outgoing_map, physical_s
from projects.ukrmol_co.paired_pole_character import real_checks
from projects.ukrmol_co.pole_character import canonical_labels, screen_witnesses
from projects.ukrmol_co.pole_reference import (
    bordered_search,
    contour_witness,
    derivative_residue,
    matching,
    pole_difference,
)
from projects.ukrmol_co.saved_matrices import read_interface, scattering_matrix
from projects.ukrmol_co.stretched_real_guard import own_model

PROPOSAL = "projects/ukrmol_co/stretched-dz-pole-proposal-20261102.json"


def numerical_proposal(source, execution):
    path = source / execution.get("proposal_path", PROPOSAL)
    if (
        not path.resolve().is_relative_to(source.resolve())
        or core.digest(path) != execution["proposal_sha256"]
    ):
        raise ValueError("Changed frozen pole proposal")
    proposal = json.loads(path.read_text())
    if "base_dz_proposal" in proposal:
        base = source / proposal["base_dz_proposal"]
        if (
            not base.resolve().is_relative_to(source.resolve())
            or core.digest(base) != proposal["base_dz_proposal_sha256"]
        ):
            raise ValueError("Changed shared numerical allocation")
        proposal = json.loads(base.read_text()) | proposal
    return proposal


def resolve(source, execution):
    if "base_execution_contract" in execution:
        path = source / execution["base_execution_contract"]
        if (
            not path.resolve().is_relative_to(source.resolve())
            or core.digest(path) != execution["base_execution_sha256"]
        ):
            raise ValueError("Changed original pole execution contract")
        base = json.loads(path.read_text())
        inputs = base["input_sha256"] | execution["input_sha256"]
        execution = base | execution
        execution["input_sha256"] = inputs
    return numerical_proposal(source, execution) | execution


def authenticate(root, source, c):
    for base, key in [(root, "input_sha256"), (source, "source_sha256")]:
        for name, want in c[key].items():
            path = base / name
            if not path.resolve().is_relative_to(base.resolve()) or core.digest(path) != want:
                raise ValueError("Changed/out-of-root pole prerequisite: " + name)
    proposal = numerical_proposal(source, c)
    if not c["execution_authorized"] or any(
        c[k] != v
        for k, v in proposal.items()
        if k not in ["status", "execution_authorized", "pole_search_released"]
    ):
        raise ValueError("Outside the explicitly authorized own-model proposal")
    case = {"stretched/cas12/dz/B1": "dz", "stretched/cas12/tz/B1": "tz"}.get(c["case"])
    if case is None:
        raise ValueError("Only the two declared finite own64/400 B1 models are supported")
    receipt = json.loads((root / c["model"]["receipt"]).read_text())
    if not (
        receipt["case"] == case
        and receipt["inner_export_verdict"] is True
        and receipt["public_reader_exit"] == 0
        and receipt["original_exits"] == [0] * 6
        and receipt["contract_sha256"] == c["model"]["inner_contract_sha256"]
        and (receipt["roots"], receipt["target_count"], receipt["boundary_channels"])
        == (2048, 64, 400)
        and receipt["spectral_truncation_verdict"] is None
        and not receipt["full_sweep_released"]
    ):
        raise ValueError("Matching own-model public native receipt required")
    reading = json.loads((root / c["completion_reading"]).read_text())
    r = reading["result"]
    if not (
        reading["reader_exit"] == 0
        and r["combined_20degree_real_prerequisite_verdict"] is True
        and r["original_exits"] == [0] * 6
        and r["zero_map_failure_original_exits"] == [2] * 6
        and r["original40degree_scientific_exit"] == 2
        and not r["complex_use_released"]
    ):
        raise ValueError("Independent combined real prerequisite required")
    seed = json.loads((root / c["seed_reading"]).read_text())["result"]
    if not (
        seed["real_prerequisite_verdict"]
        and not seed["pole_measurement"]
        and seed["cases"][case]["index"] == c["near_feature_real_index"]
    ):
        raise ValueError("Own native seed provenance differs")
    if "failed_stop_reading" in c:
        stop = json.loads((root / c["failed_stop_reading"]).read_text())
        if not (
            stop["reader_exit"] == 0
            and stop["result"]["original_exits"] == [2] * 6
            and stop["result"]["real_maps"] == stop["result"]["complex_maps"] == 0
            and stop["result"]["error"] == "KeyError('interface')"
        ):
            raise ValueError("Original zero-map adapter stop receipt required")
    if case == "tz":
        dz = json.loads((root / c["dz_pole_reading"]).read_text())
        if not (
            dz["reader_exit"] == 0
            and dz["result"]["original_exits"] == [0] * 6
            and dz["result"]["source_commit"] == "fbd16ee059f0c5fa682444a0c3af2dba4fe090fc"
            and dz["result"]["finite_reference_pole_verdict"] is True
            and dz["result"]["physical_pi_support_verdict"] is True
        ):
            raise ValueError("Independent qualified own-DZ pole required before own-TZ attempt")
    return receipt


def load_own(root, c, receipt):
    model = own_model(root / c["model"]["native_directory"], receipt)
    if (len(model["poles"]), model["amplitudes"].shape, len(model["thresholds"])) != (
        2048,
        (2048, 400),
        400,
    ):
        raise ValueError("Strict own64/400/2048 model required; no slicing or fallback")
    directory = root / c["model"]["native_directory"]
    labels = read_interface(directory / "fort.10", directory / "fort.21")["channel_labels"][:6]
    canonical_labels(labels)
    if (
        labels != c["open_labels"]
        or np.count_nonzero(model["thresholds"] < c["reference_energy_rydberg"]) != 6
    ):
        raise ValueError("Own open-channel identities differ")
    return model


def initialize(root, output, c, receipt):
    # Only populate the state consumed by the existing calculate function.
    # Never call its compressed loader, initializer, or experiment owner.
    core.MODEL = load_own(Path(root), c, receipt)
    core.OUTPUT = Path(output)
    core.CONTRACT = c
    with np.load(core.OUTPUT / "probes.npz", allow_pickle=False) as z:
        core.PROBES = z["x"]


def require_box(z, c):
    left, right, bottom, top = c["search_box_incident_hartree"]
    if not (left <= z.real <= right and bottom <= z.imag <= top):
        raise ValueError("Evaluation leaves frozen search box")


def run(root, source, output, c):
    started = time.monotonic()
    report = dict(
        case=c["case"],
        scientific_exit=2,
        disposition="inconclusive",
        finite_reference_pole_verdict=None,
        physical_pi_support_verdict=None,
        real_rows=[],
        roots=[],
        evaluations=[],
        native_calls=0,
        route_holdout_s_matrices_analyzed=0,
        tz_started=c["case"] == "stretched/cas12/tz/B1",
        automatic_successor=False,
        original_failures_requalified=False,
        full_sweep_released=False,
    )

    def save():
        report["wall_seconds"] = time.monotonic() - started
        for name in ["memory.peak", "cpu.stat"]:
            path = Path("/sys/fs/cgroup") / name
            if path.is_file():
                report["cgroup_" + name] = path.read_text()
        save_json(output / "result.json", report)

    def reserve(z, setting, amplitudes, kind):
        if time.monotonic() - started > c["resources"]["science_cap_seconds"]:
            raise TimeoutError("Combined science cap")
        if len(report["real_rows"]) + len(report["evaluations"]) >= c["resources"]["total_map_cap"]:
            raise ValueError("Frozen total map cap exhausted")
        require_box(z, c)
        row = dict(
            index=len(report["evaluations"]),
            z_hartree=[float(z.real), float(z.imag)],
            setting=setting,
            amplitudes=amplitudes,
            kind=kind,
            file=f"evaluation-{len(report['evaluations']):03d}.npz",
            completed=False,
        )
        report["evaluations"].append(row)
        save()
        return row

    cache = {}

    def evaluate(z, setting, amplitudes=False, kind="root"):
        key = (complex(z), tuple(setting), amplitudes)
        if key not in cache:
            row = reserve(z, setting, amplitudes, kind)
            row.update(core.calculate(row))
            save()
            with np.load(output / row["file"], allow_pickle=False) as data:
                cache[key] = {n: data[n] for n in ["q", "y", "a"]}
        return cache[key]

    try:
        receipt = authenticate(root, source, c)
        model = load_own(root, c, receipt)
        index = c["near_feature_real_index"]
        energy = model["energies"][index]
        for setting in c["near_feature_real_settings"]:
            row = dict(index=index, setting=setting, completed=False)
            report["real_rows"].append(row)
            save()
            endpoint, angle, order, rtol = setting
            q = surface_matrix(model["poles"], model["amplitudes"], model["ground"] + energy / 2)
            y, a, d = outgoing_map(
                model["thresholds"],
                model["angular"],
                model["c2"],
                model["c3"],
                energy,
                18,
                endpoint,
                angle,
                order,
                rtol=rtol,
                amplitudes=True,
                timeout=c["resources"]["per_map_seconds"],
                max_nfev=c["resources"]["rhs_cap"],
            )
            actual, invariants = physical_s(q, y, a)
            reference = scattering_matrix(model["saved"][index])
            np.savez(
                output / f"real-{angle}.npz",
                q=q,
                y=y,
                a=a,
                s=actual,
                reference=reference,
                energy_rydberg=energy,
            )
            row.update(
                completed=True,
                diagnostics=d,
                invariants=invariants,
                **real_checks(actual, reference, invariants, c),
            )
            save()
            if not all(row["checks"].values()):
                raise ValueError("Near-feature real gate failed before complex use")
        x, _ = np.linalg.qr(
            np.random.default_rng(c["projection_seed"]).normal(size=(2048, c["projection_columns"]))
        )
        np.savez(output / "probes.npz", x=x)
        initialize(root, output, c, receipt)
        candidate = None
        for stage, (setting, limit) in enumerate(
            zip(c["settings"], c["root_evaluation_caps"], strict=True)
        ):
            seeds = (
                [complex(*z) for z in c["seeds_incident_hartree"]]
                if candidate is None
                else [candidate, candidate + c["root_seed_refinement_delta_hartree"]]
            )
            row = dict(stage=stage, setting=setting, completed=False)
            report["roots"].append(row)
            save()

            def observe(trace, row=row):
                row["search"] = trace
                save()

            candidate, trace = bordered_search(
                lambda z, setting=setting: evaluate(z, setting),
                seeds,
                c["search_box_incident_hartree"],
                limit,
                observer=observe,
            )
            row.update(
                completed=True,
                z_hartree=[float(candidate.real), float(candidate.imag)],
                position_ev=float(hartree_to_ev(candidate.real)),
                full_width_ev=float(hartree_to_ev(-2 * candidate.imag)),
                search=trace,
            )
            save()
        errors = [pole_difference(complex(*row["z_hartree"]), candidate) for row in report["roots"]]
        report["root_setting_differences"] = errors
        if max(max(e.values()) for e in errors) > c["same_root_position_full_width_budget_ev"]:
            raise ValueError("Same-root numerical gate failed")
        radius = c["contour_radius_hartree"]
        for delta in [radius, -radius, 1j * radius, -1j * radius]:
            require_box(candidate + delta, c)
        center = evaluate(candidate, c["settings"][-1], True, "residue-center")
        derivatives = []
        for ordinal, h in enumerate(c["derivative_steps_hartree"]):
            plus = evaluate(candidate + h, c["settings"][-1], kind="derivative")
            minus = evaluate(candidate - h, c["settings"][-1], kind="derivative")
            derivative = (matching(plus["q"], plus["y"]) - matching(minus["q"], minus["y"])) / (
                2 * h
            )
            residue, d = derivative_residue(center["q"], center["y"], center["a"], derivative)
            np.savez(output / f"derivative-{ordinal}.npz", derivative=derivative, residue=residue)
            derivatives.append((residue, d))
        points = candidate + radius * np.exp(2j * np.pi * np.arange(32) / 32)
        rows = [reserve(z, c["settings"][-1], True, "contour") for z in points]
        with ProcessPoolExecutor(
            max_workers=8,
            mp_context=multiprocessing.get_context("spawn"),
            initializer=initialize,
            initargs=(root, output, c, receipt),
        ) as pool:
            futures = {pool.submit(core.calculate, row): row for row in rows}
            first_error = None
            for future in as_completed(futures):
                row = futures[future]
                if future.cancelled():
                    row["cancelled_before_execution"] = True
                    save()
                    continue
                try:
                    row.update(future.result())
                except Exception as error:
                    row["error"] = repr(error)
                    if first_error is None:
                        first_error = error
                        for other in futures:
                            other.cancel()
                save()
            if first_error is not None:
                raise first_error
        interior, outgoing = [], []
        for row in rows:
            with np.load(output / row["file"], allow_pickle=False) as data:
                interior.append(data["interior"])
                outgoing.append(data["outgoing"])
        witness, moments = contour_witness(
            candidate,
            points,
            interior,
            outgoing,
            [row["det_phase"] for row in rows],
            c["same_root_position_full_width_budget_ev"],
        )
        report["contour_witness"] = witness
        for ordinal, moment in enumerate(moments):
            np.savez(
                output / f"moments-{ordinal}.npz",
                c0=moment["c0"],
                c1=moment["c1"],
                residue=moment["residue"],
            )
        reference = moments[-1]["residue"]
        norm = np.linalg.norm(reference)
        comparison = []
        for residue, d in derivatives:
            error = float(np.linalg.norm(residue - reference) / norm) if norm else 1e300
            comparison.append(
                dict(
                    **d,
                    relative_contour_residue_error=error,
                    checks=dict(
                        center_null=d["sigma_min"] <= 1e-8,
                        residue=error <= c["derivative_residue_relative_budget"],
                        transpose=d["boundary_transpose_relative_error"]
                        <= c["derivative_residue_relative_budget"],
                    ),
                )
            )
        report["derivative_witnesses"] = comparison
        if not all(witness["checks"].values()) or not all(
            all(r["checks"].values()) for r in comparison
        ):
            raise ValueError("Contour/residue/isolation gate failed")
        report["finite_reference_pole_verdict"] = True
        report["pole_incident_hartree"] = [float(candidate.real), float(candidate.imag)]
        screen = screen_witnesses(
            [m["residue"] for m in moments] + [r for r, _ in derivatives],
            c["open_labels"],
            -2 * candidate.imag,
            c["character"],
        )
        disposition = screen["disposition"]
        report.update(
            screen=screen,
            disposition=disposition,
            physical_pi_support_verdict=True
            if disposition == "pass"
            else False
            if disposition == "reject"
            else None,
            scientific_exit=0 if disposition == "pass" else 1 if disposition == "reject" else 2,
        )
    except Exception as error:
        report.update(error=repr(error), scientific_exit=2)
    finally:
        try:
            authenticate(root, source, c)
            report["inputs_unchanged"] = True
        except Exception as error:
            report.update(
                integrity_error=repr(error), scientific_exit=2, disposition="inconclusive"
            )
        report["new_maps_reserved"] = len(report["real_rows"]) + len(report["evaluations"])
        save()
    return report["scientific_exit"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["root", "source", "output", "contract"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--profiled", action="store_true")
    args = parser.parse_args()
    c = resolve(args.source, json.loads(args.contract.read_text()))
    if args.profiled:

        def timeout(*_):
            raise TimeoutError("Authorized pole science wall cap")

        signal.signal(signal.SIGALRM, timeout)
        signal.alarm(c["resources"]["science_cap_seconds"])
        return run(args.root, args.source, args.output, c)
    from projects.ukrmol_co.run import run_profiled

    args.output.mkdir(exist_ok=False)
    (args.output / "scratch").mkdir()
    return run_profiled(
        args.output,
        os.environ | dict(TMPDIR=str(args.output / "scratch")),
        [sys.executable, "-m", "projects.ukrmol_co.stretched_dz_pole", *sys.argv[1:], "--profiled"],
    )


if __name__ == "__main__":
    raise SystemExit(main())
