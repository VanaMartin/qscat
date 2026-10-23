"""One frozen compressed-anchor boundary-root / resolvent-contour experiment."""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
import resource
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from qscat.units import hartree_to_ev

from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.native_full_reconstruction import load_model
from projects.ukrmol_co.outgoing_map import outgoing_map
from projects.ukrmol_co.pole_reference import (
    bordered_search,
    contour_witness,
    derivative_residue,
    matching,
    pole_difference,
    projected_resolvent,
)


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def initialize(source, output, contract):
    global MODEL, OUTPUT, CONTRACT, PROBES
    MODEL = load_model(Path(source))
    OUTPUT = Path(output)
    CONTRACT = contract
    with np.load(OUTPUT / "probes.npz") as data:
        PROBES = data["x"]


def calculate(row):
    z = complex(*row["z_hartree"])
    endpoint, angle, order, rtol = row["setting"]
    model = MODEL
    y, a, diagnostics = outgoing_map(
        model["thresholds"],
        model["angular"],
        model["c2"],
        model["c3"],
        2 * z,
        18,
        endpoint,
        angle,
        order,
        reference=CONTRACT["reference_energy_rydberg"],
        rtol=rtol,
        amplitudes=row["amplitudes"],
        continued_amplitudes=row["amplitudes"],
        timeout=CONTRACT["resources"]["per_map_seconds"],
        max_nfev=CONTRACT["resources"]["rhs_cap"],
    )
    total = model["ground"] + z
    q = surface_matrix(model["poles"], model["amplitudes"], total)
    sign, logabs = np.linalg.slogdet(matching(q, y))
    capture = dict(q=q, y=y, a=a, z_hartree=np.array(z), total_hartree=np.array(total))
    if row["kind"] == "contour":
        r, t = projected_resolvent(model["poles"], model["amplitudes"], total, y, a, PROBES, q)
        capture.update(interior=r, outgoing=t)
    np.savez(OUTPUT / row["file"], **capture)
    return dict(
        completed=True,
        diagnostics=diagnostics,
        det_phase=float(np.angle(sign)),
        log_abs_det_m=float(logabs) if np.isfinite(logabs) else None,
    )


def run(root, output, contract):
    output.mkdir(exist_ok=False)
    started = time.monotonic()
    raw = root / contract["raw_directory"]
    report = dict(
        case=contract["case"],
        evaluations=[],
        roots=[],
        finite_reference_pole_verdict=None,
        independent_exterior_verdict=None,
        physical_co_accuracy_verdict=None,
        width_curve_verdict=None,
        scientific_exit=2,
        native_calls=0,
        route_holdout_s_matrices_analyzed=0,
        full_sweep_released=False,
        scope=contract["claim_boundary"],
    )

    def save():
        usage = resource.getrusage(resource.RUSAGE_SELF)
        report.update(
            wall_seconds=time.monotonic() - started,
            process_cpu_seconds=usage.ru_utime + usage.ru_stime,
            peak_parent_rss_platform_units=usage.ru_maxrss,
        )
        for name in ["memory.peak", "cpu.stat"]:
            path = Path("/sys/fs/cgroup") / name
            if path.is_file():
                report["cgroup_" + name] = path.read_text()
        temporary = output / "result.partial"
        temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        temporary.replace(output / "result.json")

    def verify():
        for name, want in contract["input_sha256"].items():
            if digest(root / name) != want:
                raise ValueError(f"Changed prerequisite: {name}")
        manifest = json.loads((raw / "input-manifest.json").read_text())
        if len(manifest) != 123:
            raise ValueError("Changed original input inventory")
        for name, want in manifest.items():
            if digest(raw / name) != want:
                raise ValueError(f"Changed original input: {name}")
        return len(manifest)

    def reserve(z, setting, amplitudes, kind):
        if len(report["evaluations"]) >= contract["resources"]["total_map_cap"]:
            raise ValueError("Total frozen map cap exhausted")
        box = contract["search_box_incident_hartree"]
        if not (box[0] <= z.real <= box[1] and box[2] <= z.imag <= box[3]):
            raise ValueError("Evaluation leaves frozen search box")
        index = len(report["evaluations"])
        row = dict(
            index=index,
            z_hartree=[float(z.real), float(z.imag)],
            setting=setting,
            amplitudes=amplitudes,
            kind=kind,
            file=f"evaluation-{index:03d}.npz",
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
            row.update(calculate(row))
            save()
            with np.load(output / row["file"]) as data:
                cache[key] = {name: data[name] for name in ["q", "y", "a"]}
        return cache[key]

    def reject(reason):
        report.update(finite_reference_pole_verdict=False, scientific_exit=1, stopped_at=reason)
        return 1

    try:
        report["inputs_verified_before"] = verify()
        canary = json.loads((root / contract["canary_result"]).read_text())
        if not (
            canary["canary_verdict"] is True
            and len(canary["rows"]) == 11
            and all(all(row["checks"].values()) for row in canary["rows"])
        ):
            raise ValueError("Finite outgoing-map prerequisite is not qualified")
        for module in ["test_outgoing_map", "test_pole_reference"]:
            with (output / (module + ".log")).open("x") as log:
                result = subprocess.run(
                    [sys.executable, "-m", "projects.ukrmol_co." + module],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    timeout=120,
                )
            report[module + "_exit"] = result.returncode
            save()
            if result.returncode:
                raise ValueError("Controls failed before CO pole use")
        if contract["case"] not in ["compressed/nine/B1", "compressed/nine/B2"]:
            raise ValueError("Only frozen development anchor allowed")
        model = load_model(raw / contract["case"])
        rng = np.random.default_rng(contract["projection_seed"])
        x, _ = np.linalg.qr(rng.normal(size=(len(model["poles"]), contract["projection_columns"])))
        np.savez(output / "probes.npz", x=x)
        initialize(raw / contract["case"], output, contract)
        candidate = None
        for stage, (setting, limit) in enumerate(
            zip(contract["settings"], contract["root_evaluation_caps"], strict=True)
        ):
            if candidate is None:
                seeds = [complex(*z) for z in contract["seeds_incident_hartree"]]
            else:
                seeds = [candidate, candidate + contract["root_seed_refinement_delta_hartree"]]
            record = dict(stage=stage, setting=setting, completed=False)
            report["roots"].append(record)
            save()

            def observe(value, record=record):
                record["search"] = value
                save()

            candidate, trace = bordered_search(
                lambda z, setting=setting: evaluate(z, setting),
                seeds,
                contract["search_box_incident_hartree"],
                limit,
                observer=observe,
            )
            record.update(
                completed=True,
                z_hartree=[float(candidate.real), float(candidate.imag)],
                position_ev=float(hartree_to_ev(candidate.real)),
                full_width_ev=float(hartree_to_ev(-2 * candidate.imag)),
                search=trace,
            )
            save()
        budget = contract["same_root_position_full_width_budget_ev"]
        errors = [pole_difference(complex(*row["z_hartree"]), candidate) for row in report["roots"]]
        report["root_setting_differences"] = errors
        if not all(max(error.values()) <= budget for error in errors):
            return reject("root-setting error gate")
        box = contract["search_box_incident_hartree"]
        radius = contract["contour_radius_hartree"]
        if (
            min(
                candidate.real - box[0],
                box[1] - candidate.real,
                candidate.imag - box[2],
                box[3] - candidate.imag,
            )
            <= radius
        ):
            raise ValueError("Candidate lacks frozen contour clearance")
        setting = contract["settings"][-1]
        center = evaluate(candidate, setting, True, "residue-center")
        derivatives = []
        for h in contract["derivative_steps_hartree"]:
            plus = evaluate(candidate + h, setting, kind="derivative")
            minus = evaluate(candidate - h, setting, kind="derivative")
            derivative = (matching(plus["q"], plus["y"]) - matching(minus["q"], minus["y"])) / (
                2 * h
            )
            residue, diagnostics = derivative_residue(
                center["q"], center["y"], center["a"], derivative
            )
            np.savez(
                output / f"derivative-{len(derivatives)}.npz",
                derivative=derivative,
                residue=residue,
            )
            derivatives.append((residue, diagnostics))
            save()
        points = candidate + radius * np.exp(
            2j * np.pi * np.arange(contract["contour_nodes"]) / contract["contour_nodes"]
        )
        rows = [reserve(z, setting, True, "contour") for z in points]
        with ProcessPoolExecutor(
            max_workers=contract["resources"]["contour_workers"],
            mp_context=multiprocessing.get_context("spawn"),
            initializer=initialize,
            initargs=(raw / contract["case"], output, contract),
        ) as pool:
            futures = {pool.submit(calculate, row): row for row in rows}
            first_error = None
            for future in as_completed(futures):
                row = futures[future]
                if future.cancelled():
                    row["cancelled_before_execution"] = True
                    save()
                    continue
                try:
                    row.update(future.result())
                except BaseException as error:
                    row["error"] = repr(error)
                    if first_error is None:
                        first_error = error
                        for other in futures:
                            other.cancel()
                save()
            if first_error is not None:
                raise first_error
        interior = []
        outgoing = []
        for row in rows:
            with np.load(output / row["file"]) as data:
                interior.append(data["interior"])
                outgoing.append(data["outgoing"])
        witness, moments = contour_witness(
            candidate, points, interior, outgoing, [r["det_phase"] for r in rows], budget
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
        norm = float(np.linalg.norm(reference))
        comparison = []
        for residue, diagnostics in derivatives:
            error = float(np.linalg.norm(residue - reference) / norm) if norm else 1e300
            comparison.append(
                dict(
                    **diagnostics,
                    relative_contour_residue_error=error,
                    checks=dict(
                        center_null=diagnostics["sigma_min"] <= 1e-8,
                        residue=error <= contract["derivative_residue_relative_budget"],
                        transpose=diagnostics["boundary_transpose_relative_error"]
                        <= contract["derivative_residue_relative_budget"],
                    ),
                )
            )
        report["derivative_witnesses"] = comparison
        passed = all(witness["checks"].values()) and all(
            all(v["checks"].values()) for v in comparison
        )
        if not passed:
            return reject("contour/residue/isolation gate")
        report.update(
            finite_reference_pole_verdict=True,
            scientific_exit=0,
            pole_incident_hartree=[float(candidate.real), float(candidate.imag)],
        )
        return 0
    except BaseException as error:
        report.update(error=repr(error), scientific_exit=2)
        return 2
    finally:
        try:
            report["inputs_verified_after"] = verify()
        except BaseException as error:
            report.update(integrity_error=repr(error), scientific_exit=2)
        save()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["root", "output", "contract"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    run(args.root, args.output, json.loads(args.contract.read_text()))
    raise SystemExit(json.loads((args.output / "result.json").read_text())["scientific_exit"])
