"""Finite energy-corner and finest-step endpoint qualification of saved CO models."""

from __future__ import annotations

import argparse
import json
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from itertools import pairwise
from pathlib import Path

import numpy as np
from scipy.linalg import LinAlgWarning

from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.native_exterior_attribution import algebra, differences
from projects.ukrmol_co.native_exterior_study import cayley
from projects.ukrmol_co.native_full_reconstruction import digest, gated_pool, load_model, save
from projects.ukrmol_co.numerov_increment import increment_numerov_k
from projects.ukrmol_co.saved_matrices import read_k


def configurations(steps):
    """Three spacing rows, two extra endpoints and one independent tail-order row."""
    return [(h, 400.0, 8) for h in steps] + [
        (steps[-1], 200.0, 8),
        (steps[-1], 800.0, 8),
        (steps[-1], 400.0, 6),
    ]


def evaluate(rows, reference, steps, budget):
    expected = set(configurations(steps))
    matrices = {
        tuple(row[k] for k in ["step_bohr", "endpoint_bohr", "tail_order"]): np.asarray(
            row["k"], float
        )
        for row in rows
    }
    if len(rows) != 6 or len(matrices) != 6 or set(matrices) != expected:
        raise ValueError("Missing, duplicate or undeclared configuration")
    if any(k.shape != reference.shape or not np.isfinite(k).all() for k in matrices.values()):
        raise ValueError("Nonfinite or incompatible K matrix")
    fine = steps[-1]
    step_changes = [
        differences(matrices[a, 400, 8], matrices[b, 400, 8]) for a, b in pairwise(steps)
    ]
    endpoint_changes = [
        differences(matrices[fine, a, 8], matrices[fine, b, 8]) for a, b in [(200, 400), (400, 800)]
    ]
    tail_change = differences(matrices[fine, 400, 6], matrices[fine, 400, 8])
    native = {str(int(r)): differences(matrices[fine, r, 8], reference) for r in [200, 400, 800]}
    checks = {
        family + "_" + key: value <= budget
        for family, metrics in [
            ("step", step_changes[-1]),
            ("endpoint", endpoint_changes[-1]),
            ("order", tail_change),
            *[("native_" + radius, metrics) for radius, metrics in native.items()],
        ]
        for key, value in metrics.items()
    }
    checks["spatial_trend"] = step_changes[-1]["s"] <= max(step_changes[0]["s"] / 8, 1e-8)
    for key in ["unitarity", "reciprocity"]:
        checks[key] = all(algebra(k)[key] <= budget for k in matrices.values())
    return dict(
        verdict=all(checks.values()),
        checks=checks,
        metrics=dict(
            step_changes=step_changes,
            endpoint_changes=endpoint_changes,
            order_change=tail_change,
            finest_endpoint_native=native,
        ),
        all_configuration_native=[
            dict(configuration=list(key), **differences(k, reference))
            for key, k in matrices.items()
        ],
    )


def identity(task):
    args, item, contract = task
    case, index, steps = item["case"], item["energy_index"], item["steps_bohr"]
    destination = args.output / case / f"{index:03d}"
    destination.mkdir(parents=True)
    result = dict(case=case, energy_index=index, verdict=None, rows=[], checks={})
    save(destination / "report.json", result)
    begin = time.monotonic()
    try:
        model = load_model(args.raw / case)
        energy = float(model["energies"][index])
        if energy.hex() != item["energy_rydberg_hex"]:
            raise ValueError("Saved energy bits differ from contract")
        folder = args.prior if item["reference"] == "64" else args.saved
        energies, native, header = read_k(
            folder / "result/native" / case / item["reference"] / model["unit"]
        )
        if any(
            header[key] != model["header"][key]
            for key in [
                "target_components",
                "boundary_channels",
                "symmetry_charge",
                "geometry_mass",
                "open_channels",
                "stored_triangle_only",
            ]
        ):
            raise ValueError("Native reference header differs")
        if item["reference"] == "64":
            if native.shape != (1, 6, 6) or energies[0] != energy:
                raise ValueError("Single-energy native64 identity differs")
            reference = native[0]
        else:
            if not np.array_equal(energies, model["energies"]):
                raise ValueError("Full-grid native32 energy bits differ")
            e16, n16, h16 = read_k(args.saved / "result/native" / case / "16" / model["unit"])
            if not np.array_equal(e16, energies) or h16 != header:
                raise ValueError("Native16 control identity differs")
            result["native_reference_change"] = differences(n16[index], native[index])
            result["reference_checks"] = {
                k: v <= contract["absolute_budget"]
                for k, v in result["native_reference_change"].items()
            }
            if not all(result["reference_checks"].values()):
                result["verdict"] = False
                return result
            reference = native[index]
        result["energy_rydberg"] = energy
        q = surface_matrix(model["poles"], model["amplitudes"], model["ground"] + energy / 2)
        reused = []
        if item["reference"] == "64":
            old = json.loads((args.increment / "result/report.json").read_text())
            reused = [
                row
                for row in old["rows"]
                if row["case"] == case and row["endpoint_bohr"] == 400 and row["tail_order"] == 8
            ]
            if len(reused) != 3 or any(row["energy_rydberg"] != energy for row in reused):
                raise ValueError("Reused step rows differ")
        with (destination / "rows.jsonl").open("x") as stream:
            for step, endpoint, order in configurations(steps):
                previous = next(
                    (
                        row
                        for row in reused
                        if row["step_bohr"] == step
                        and row["endpoint_bohr"] == endpoint
                        and row["tail_order"] == order
                    ),
                    None,
                )
                if previous is not None:
                    row = dict(previous, reused=True)
                else:
                    start = time.monotonic()
                    with warnings.catch_warnings():
                        warnings.simplefilter("error", LinAlgWarning)
                        k, diagnostics = increment_numerov_k(
                            q,
                            model["thresholds"],
                            model["angular"],
                            model["c2"],
                            model["c3"],
                            energy,
                            18,
                            endpoint,
                            step,
                            order,
                            1e-10,
                        )
                    s = cayley(k)
                    row = dict(
                        step_bohr=step,
                        endpoint_bohr=endpoint,
                        tail_order=order,
                        energy_rydberg=energy,
                        k=k.tolist(),
                        s_real=s.real.tolist(),
                        s_imag=s.imag.tolist(),
                        invariants=algebra(k),
                        diagnostics=diagnostics,
                        wall_seconds=time.monotonic() - start,
                        reused=False,
                    )
                result["rows"].append(row)
                stream.write(json.dumps(row, allow_nan=False) + "\n")
                stream.flush()
                save(destination / "report.json", result)
        result.update(evaluate(result["rows"], reference, steps, contract["absolute_budget"]))
    except (
        OSError,
        ValueError,
        RuntimeError,
        AssertionError,
        np.linalg.LinAlgError,
        LinAlgWarning,
    ) as error:
        result.update(verdict=None, failure=f"{type(error).__name__}: {error}")
    finally:
        result["wall_seconds"] = time.monotonic() - begin
        save(destination / "report.json", result)
    return result


def verify(args, contract):
    for folder, key in [
        (args.raw, "input_sha256"),
        (args.saved, "saved_sha256"),
        (args.prior, "prior_sha256"),
        (args.increment, "increment_sha256"),
    ]:
        for name, want in contract[key].items():
            if digest(folder / name) != want:
                raise ValueError(f"Pinned input differs: {folder / name}")
    for name, want in contract["prerequisite_sha256"].items():
        path = Path(__file__).with_name(name)
        if digest(path) != want:
            raise ValueError(f"Prerequisite differs: {name}")
    prerequisite = json.loads(
        Path(__file__).with_name("verification-numerov-endpoint-analytic-20261019.json").read_text()
    )
    if prerequisite["analytic_endpoint_verdict"] is not True:
        raise ValueError("Analytic endpoint prerequisite did not pass")
    public = json.loads(
        Path(__file__).with_name("verification-numerov-increment-public-20261019.json").read_text()
    )
    if not (
        public["experiment_verdict"] is True
        and public["public_numerical_reexecution_passed"] is True
        and public["all_final_payloads_fresh_downloaded"] is True
    ):
        raise ValueError("Previous increment public qualification is incomplete")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ["raw", "saved", "prior", "increment", "output", "contract"]:
        parser.add_argument("--" + key, type=Path, required=True)
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text())
    args.output.mkdir()
    report = dict(
        contract=contract,
        endpoint=[],
        corners=[],
        stage_verdicts={},
        experiment_verdict=None,
        original_two_point_verdict=True,
        original_attribution_verdict=False,
        original_attribution_exits=[1, 1],
        original_reconstruction_verdict=False,
        original_reconstruction_exits=[1, 1],
        native_full_k_reconstruction_verdict=False,
        co_complex_pole_verdict=None,
        native_invocations=0,
        new_inner_solves=0,
        full_sweep_released=False,
    )
    save(args.output / "report.json", report)
    begin = time.monotonic()
    try:
        verify(args, contract)

        def receive(stage, row):
            report[stage].append(row)
            save(args.output / "report.json", report)

        with ProcessPoolExecutor(max_workers=8) as pool:
            for stage, workers in [("endpoint", 2), ("corners", 8)]:
                tasks = [
                    (args, item, contract)
                    for item in contract["identities"]
                    if item["stage"] == stage
                ]
                verdict = gated_pool(
                    pool, identity, tasks, workers, lambda row, stage=stage: receive(stage, row)
                )
                report["stage_verdicts"][stage] = verdict
                if verdict is not True:
                    report["experiment_verdict"] = verdict
                    break
            else:
                expected = {(x["case"], x["energy_index"]) for x in contract["identities"]}
                found = [
                    (x["case"], x["energy_index"]) for x in report["endpoint"] + report["corners"]
                ]
                report["experiment_verdict"] = (
                    len(found) == len(expected) and set(found) == expected
                )
    except (ValueError, OSError, RuntimeError, AssertionError) as error:
        report.update(experiment_verdict=None, failure=f"{type(error).__name__}: {error}")
    finally:
        try:
            verify(args, contract)
            report["all_input_hashes_verified_after"] = True
        except (ValueError, OSError) as error:
            report.update(experiment_verdict=None, integrity_failure=str(error))
        report["wall_seconds"] = time.monotonic() - begin
        save(args.output / "report.json", report)
    if report["experiment_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
