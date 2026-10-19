"""Prospective fixed-grid cancellation-reduced control of two saved CO identities."""

from __future__ import annotations

import argparse
import json
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.linalg import LinAlgWarning

from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.native_exterior_attribution import algebra, differences
from projects.ukrmol_co.native_exterior_study import cayley
from projects.ukrmol_co.native_full_reconstruction import digest, gated_pool, load_model, save
from projects.ukrmol_co.numerov_increment import increment_numerov_k
from projects.ukrmol_co.saved_matrices import read_k, unit


def control(task):
    raw, output, case, configuration, contract = task
    step, endpoint, order = configuration
    result = dict(case=case, step_bohr=step, endpoint_bohr=endpoint, tail_order=order, verdict=None)
    destination = output / "independent" / case / f"h{step}-r{endpoint}-p{order}.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    save(destination, result)
    try:
        model = load_model(raw / case)
        energy = model["energies"][contract["energy_index"]]
        if energy != contract["energy_rydberg"]:
            raise ValueError("Requested energy identity differs")
        q = surface_matrix(model["poles"], model["amplitudes"], model["ground"] + energy / 2)
        begin = time.monotonic()
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
        invariants = algebra(k)
        result.update(
            k=k.tolist(),
            s_real=s.real.tolist(),
            s_imag=s.imag.tolist(),
            diagnostics=diagnostics,
            invariants=invariants,
            energy_rydberg=float(energy),
            wall_seconds=time.monotonic() - begin,
            verdict=all(x <= contract["absolute_budget"] for x in invariants.values()),
        )
    except (ValueError, OSError, RuntimeError, np.linalg.LinAlgError, LinAlgWarning) as error:
        result["failure"] = f"{type(error).__name__}: {error}"
    save(destination, result)
    return result


def evaluate(rows, reference, budget):
    """Fixed families, with independent K/S/phase vetoes for each comparison."""
    expected = {(h, 400.0, 8) for h in [0.025, 0.0125, 0.00625]}
    expected |= {(0.0125, r, 8) for r in [200.0, 800.0]}
    matrices = {
        (r["step_bohr"], r["endpoint_bohr"], r["tail_order"]): np.asarray(r["k"]) for r in rows
    }
    if len(rows) != 5 or set(matrices) != expected:
        raise ValueError("Missing, duplicate or undeclared candidate control")
    step_changes = [
        differences(matrices[a, 400.0, 8], matrices[b, 400.0, 8])
        for a, b in [(0.025, 0.0125), (0.0125, 0.00625)]
    ]
    endpoint_changes = [
        differences(matrices[0.0125, a, 8], matrices[0.0125, b, 8])
        for a, b in [(200.0, 400.0), (400.0, 800.0)]
    ]
    native = differences(matrices[0.00625, 400.0, 8], reference)
    checks = {
        family + "_" + k: v <= budget
        for family, metrics in [
            ("step", step_changes[-1]),
            ("endpoint", endpoint_changes[-1]),
            ("native", native),
        ]
        for k, v in metrics.items()
    }
    checks["spatial_trend"] = step_changes[-1]["s"] <= max(step_changes[0]["s"] / 8, 1e-8)
    for key in ["unitarity", "reciprocity"]:
        checks[key] = all(algebra(k)[key] <= budget for k in matrices.values())
    return dict(
        verdict=all(checks.values()),
        checks=checks,
        metrics=dict(
            step_changes=step_changes, endpoint_changes=endpoint_changes, native_comparison=native
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ["raw", "prior", "saved", "contract", "output"]:
        parser.add_argument("--" + key, type=Path, required=True)
    args = parser.parse_args()
    c = json.loads(args.contract.read_text())
    for directory, pins in [
        (args.raw, c["input_sha256"]),
        (args.prior, c["prior_sha256"]),
        (args.saved, c["saved_sha256"]),
    ]:
        for name, want in pins.items():
            if digest(directory / name) != want:
                raise ValueError(f"Input digest differs: {name}")
    qualification = Path(__file__).with_name(
        "verification-numerov-increment-analytic-20261019.json"
    )
    if digest(qualification) != c["analytic_qualification_sha256"]:
        raise ValueError("Analytic prerequisite differs")
    if json.loads(qualification.read_text())["analytic_qualification_verdict"] is not True:
        raise ValueError("Analytic qualification did not pass")
    args.output.mkdir()
    report = dict(
        contract=c,
        rows=[],
        cases=[],
        experiment_verdict=None,
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

    def receive(row):
        report["rows"].append(row)
        save(args.output / "report.json", report)

    with ProcessPoolExecutor(max_workers=4) as pool:
        tasks = [
            (args.raw, args.output, case, config, c)
            for config in c["configurations"]
            for case in c["cases"]
        ]
        verdict = gated_pool(pool, control, tasks, 4, receive)
    report["invariant_stage_verdict"] = verdict
    if verdict is True:
        prior = json.loads((args.prior / "result/report.json").read_text())
        for case in c["cases"]:
            rows = [row for row in report["rows"] if row["case"] == case]
            _, native, _ = read_k(
                args.prior
                / "result/native"
                / case
                / "64"
                / unit(args.raw / case / "rsolve.inp", "lukmt")
            )
            result = evaluate(rows, native[0], c["absolute_budget"])
            originals = json.loads(
                (args.saved / "result/independent" / case / "018.json").read_text()
            )["rows"]
            originals += [x for x in prior["independent"] if x["case"] == case]
            result["same_grid_original_differences"] = []
            for row in rows:
                old = next(
                    x
                    for x in originals
                    if all(x[k] == row[k] for k in ["step_bohr", "endpoint_bohr", "tail_order"])
                )
                result["same_grid_original_differences"].append(
                    dict(
                        configuration=[
                            row[k] for k in ["step_bohr", "endpoint_bohr", "tail_order"]
                        ],
                        **differences(np.asarray(row["k"]), np.asarray(old["k"])),
                    )
                )
            report["cases"].append(dict(case=case, **result))
        report["experiment_verdict"] = all(x["verdict"] for x in report["cases"])
    else:
        report["experiment_verdict"] = verdict
    for directory, pins in [
        (args.raw, c["input_sha256"]),
        (args.prior, c["prior_sha256"]),
        (args.saved, c["saved_sha256"]),
    ]:
        for name, want in pins.items():
            if digest(directory / name) != want:
                raise ValueError(f"Post-run digest differs: {name}")
    report["all_input_hashes_verified_after"] = True
    save(args.output / "report.json", report)
    if report["experiment_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
