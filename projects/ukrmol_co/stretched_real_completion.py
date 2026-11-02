"""Authorized missing20degree real checks only; reuses rather than repeats DZ80."""

from __future__ import annotations

import argparse
import json
import os
import signal
import sys
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.compressed_a2_contrast_run import save_json
from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.outgoing_map import outgoing_map, physical_s
from projects.ukrmol_co.paired_pole_character import real_checks
from projects.ukrmol_co.saved_matrices import scattering_matrix
from projects.ukrmol_co.stretched_basis_properties import digest
from projects.ukrmol_co.stretched_real_guard import authenticate as guard_authenticate
from projects.ukrmol_co.stretched_real_guard import own_model

MAPS = [("dz", 60, 120), ("tz", 60, 120), ("tz", 80, 60), ("tz", 80, 90), ("tz", 80, 120)]


def resolve(source, execution):
    base = source / execution["base_guard_contract"]
    if (
        not base.resolve().is_relative_to(source.resolve())
        or digest(base) != execution["base_guard_sha256"]
    ):
        raise ValueError("Changed original method/input contract")
    original = json.loads(base.read_text())
    inherited = [
        "models",
        "input_sha256",
        "source_sha256",
        "real_s_absolute_budget",
        "real_invariant_absolute_budget",
        "real_phase_sum_budget_rad",
        "real_flux_absolute_budget",
    ]
    result = {key: original[key] for key in inherited} | execution
    for key in ["input_sha256", "source_sha256"]:
        result[key] = original[key] | execution[key]
    return result


def authenticate(root, source, contract):
    if (
        not contract["execution_authorized"]
        or contract["maps"] != [list(row) for row in MAPS]
        or (
            contract["angle_degrees"],
            contract["boundary_bohr"],
            contract["order"],
            contract["rtol"],
        )
        != (20, 18, 10, 1e-12)
    ):
        raise ValueError("Outside the explicitly authorized five-map completion")
    receipts = guard_authenticate(root, source, contract)
    previous = json.loads((root / contract["previous_public_reading"]).read_text())
    prior = previous["result"]
    if not (
        previous["reader_exit"] == 0
        and prior["contour_study_verdict"] is True
        and prior["original_exits"] == [0] * 6
        and prior["original_guard_scientific_exit"] == 2
        and prior["original_guard_requalified"] is False
        and not prior["tz_started"]
    ):
        raise ValueError("Independent public actual-DZ80 success required")
    original = json.loads((root / contract["original_guard_result"]).read_text())
    if original["scientific_exit"] != 2 or original["real_guard_verdict"] is not None:
        raise ValueError("Original40degree stop must remain unchanged")
    return receipts


def compute(root, source, output, c):
    started = time.monotonic()
    result = dict(
        scientific_exit=2,
        real_completion_verdict=None,
        real_rows=[],
        real_maps_reserved=0,
        previous_dz80_repeated=False,
        native_calls=0,
        new_complex_maps=0,
        original_guard_requalified=False,
        route_holdout_s_matrices_analyzed=0,
        automatic_successor=False,
        full_sweep_released=False,
    )

    def save():
        result["science_wall_seconds"] = time.monotonic() - started
        for name in ["memory.peak", "cpu.stat"]:
            path = Path("/sys/fs/cgroup") / name
            if path.exists():
                result["cgroup_" + name] = path.read_text()
        save_json(output / "result.json", result)

    try:
        receipts = authenticate(root, source, c)
        models = {}
        captures = []
        for case, index, endpoint in MAPS:
            if time.monotonic() - started > c["science_timeout_seconds"]:
                raise TimeoutError("Completion science cap")
            if case not in models:
                models[case] = own_model(
                    root / c["models"][case]["native_directory"], receipts[case]
                )
            model = models[case]
            energy = model["energies"][index]
            q = surface_matrix(model["poles"], model["amplitudes"], model["ground"] + energy / 2)
            reference = scattering_matrix(model["saved"][index])
            row = dict(case=case, index=index, endpoint=endpoint, completed=False)
            result["real_rows"].append(row)
            result["real_maps_reserved"] += 1
            save()
            y, a, d = outgoing_map(
                model["thresholds"],
                model["angular"],
                model["c2"],
                model["c3"],
                energy,
                c["boundary_bohr"],
                endpoint,
                c["angle_degrees"],
                c["order"],
                rtol=c["rtol"],
                amplitudes=True,
                timeout=c["per_map_seconds"],
                max_nfev=c["rhs_cap"],
            )
            actual, invariants = physical_s(q, y, a)
            np.savez(
                output / f"{case}-{index}-{endpoint}.npz",
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
            if case == "tz" and index == 80:
                captures.append((y, a))
            save()
            if (index == 60 or endpoint == 120) and not all(row["checks"].values()):
                raise ValueError("Standalone or final TZ80 real physical gate failed")
        stability = dict(
            y_absolute=float(abs(captures[-1][0] - captures[-2][0]).max()),
            a_absolute=float(abs(captures[-1][1] - captures[-2][1]).max()),
        )
        result["tz80_last_two_stability"] = stability
        if max(stability.values()) > c["tz80_last_two_y_and_a_stability_budget"]:
            raise ValueError("TZ80 endpoint stability failed")
        authenticate(root, source, c)
        result.update(scientific_exit=0, real_completion_verdict=True, inputs_unchanged=True)
    except Exception as error:
        result["error"] = repr(error)
    finally:
        save()
    return result["scientific_exit"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["root", "source", "output", "contract"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--profiled", action="store_true")
    args = parser.parse_args()
    c = resolve(args.source, json.loads(args.contract.read_text()))
    if args.profiled:

        def timeout(*_):
            raise TimeoutError("Authorized completion science wall cap")

        signal.signal(signal.SIGALRM, timeout)
        signal.alarm(c["science_timeout_seconds"])
        return compute(args.root, args.source, args.output, c)
    from projects.ukrmol_co.run import run_profiled

    args.output.mkdir(exist_ok=False)
    (args.output / "scratch").mkdir()
    return run_profiled(
        args.output,
        os.environ | dict(TMPDIR=str(args.output / "scratch")),
        [
            sys.executable,
            "-m",
            "projects.ukrmol_co.stretched_real_completion",
            *sys.argv[1:],
            "--profiled",
        ],
    )


if __name__ == "__main__":
    raise SystemExit(main())
