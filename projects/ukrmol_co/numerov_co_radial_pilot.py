"""Fixed long-range saved-CO qualification around the unchanged coupled solver."""

from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.linalg import LinAlgWarning

from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.native_exterior_attribution import algebra, differences
from projects.ukrmol_co.native_full_reconstruction import digest, load_model, save, verify_inputs
from projects.ukrmol_co.numerov_co_canary import condition_estimate, oracle_checks
from projects.ukrmol_co.numerov_co_canary_condition_replay import condition_drift
from projects.ukrmol_co.numerov_native_coupled import (
    REFERENCE_SOLVE,
    ResidualAudit,
    installed_solver,
    native_refined_solve,
    reference,
)
from projects.ukrmol_co.numerov_native_residual_control import canonical_hash
from projects.ukrmol_co.saved_matrices import read_k


def sample_indices(endpoint, step):
    """Both recurring systems at fixed radii, plus initialization/end/matching."""
    n = round((endpoint - 18) / step)
    indices = {0, 2 * n - 2, 2 * n - 1}
    for radius in [25, 50, 100, 150, 200, 300]:
        i = round((radius - 18) / step)
        if 1 <= i < n:
            indices.update([2 * i - 1, 2 * i])
    return indices


def solve_role(index, n):
    if index == 0:
        return "initialization"
    if index == 2 * n - 1:
        return "matching"
    return "potential" if index % 2 else "ratio"


def propagate(args, energy, config, directory, *, dtype=np.longdouble):
    """Instrument original full propagation; do not replace recurrence/tail algebra."""
    directory.mkdir(parents=True)
    endpoint, step, backend = (config[k] for k in ["endpoint_bohr", "step_bohr", "backend"])
    n = round((endpoint - 18) / step)
    sampled = sample_indices(endpoint, step)
    captured = {0, 2 * n - 2, 2 * n - 1}
    row = dict(config, samples=[], captures=[], solve_calls=0, experiment_verdict=None)
    audit = ResidualAudit()
    role_times = {key: 0.0 for key in ["initialization", "potential", "ratio", "matching"]}
    bins = {}
    started = time.perf_counter()

    def measured(a, b, stats):
        index = row["solve_calls"]
        if index >= 2 * n:
            raise ValueError("Undeclared full-propagation solve")
        role = solve_role(index, n)
        radius = 18 + ((index + 1) // 2) * step if role == "potential" else 18 + index // 2 * step
        if role == "matching":
            radius = endpoint
        before = time.perf_counter()
        x = (
            native_refined_solve(a, b, stats, observer=audit if index in sampled else None)
            if backend == "native"
            else REFERENCE_SOLVE(a, b, stats)
        )
        elapsed = time.perf_counter() - before
        role_times[role] += elapsed
        # Audited solve time is excluded from recurring throughput observations.
        if role in ["potential", "ratio"] and index not in sampled:
            label = str(int(radius // 25) * 25)
            bucket = bins.setdefault(label, dict(solves=0, seconds=0.0))
            bucket["solves"] += 1
            bucket["seconds"] += elapsed
        if index in sampled:
            rcond = condition_estimate(a)
            if rcond <= 0:
                raise ValueError("Nonpositive sampled reciprocal condition estimate")
            row["samples"].append(
                dict(
                    index=index,
                    role=role,
                    radius_bohr=radius,
                    reciprocal_condition_1=rcond,
                    a_sha256=canonical_hash(a),
                    b_sha256=canonical_hash(b),
                    x_sha256=canonical_hash(x),
                )
            )
        if index in captured and backend == "native":
            path = directory / f"solve-{index:06d}.npz"
            np.savez(path, a=a, b=b, x=x)
            row["captures"].append(dict(file=path.name, sha256=digest(path), index=index))
        row["solve_calls"] += 1
        if row["solve_calls"] % 1024 == 0 or row["solve_calls"] == 2 * n:
            save(
                directory / "progress.json",
                dict(
                    config=config,
                    completed_solves=row["solve_calls"],
                    expected_solves=2 * n,
                    last_radius_bohr=radius,
                    elapsed_seconds=time.perf_counter() - started,
                    diagnostics=stats,
                    samples=row["samples"],
                ),
            )
        return x

    with installed_solver(measured):
        k, stats = reference.typed_coupled_k(
            *args, energy, 18, endpoint, step, 8, 1e-10, dtype=dtype
        )
    elapsed = time.perf_counter() - started
    row.update(
        k=np.asarray(k, float).tolist(),
        extended_k_sha256=canonical_hash(k),
        diagnostics=stats,
        invariants=algebra(np.asarray(k, float)),
        residual_audit=audit.summary() if backend == "native" else None,
        role_solve_seconds=role_times,
        unaudited_recurring_bins=bins,
        kernel_wall_seconds=elapsed,
    )
    return row


def row_checks(row, budget):
    n = round((row["endpoint_bohr"] - 18) / row["step_bohr"])
    checks = dict(
        solve_count=row["solve_calls"] == row["diagnostics"]["large_system_solves"] == 2 * n,
        radial_intervals=row["diagnostics"]["radial_intervals"] == n,
        real_and_complex_precision=row["diagnostics"]["significand_bits"]
        == row["diagnostics"]["complex_significand_bits"]
        == 64,
        backward_error=row["diagnostics"]["maximum_solve_backward_error"] <= 32 * 2.0**-63,
        corrections=row["diagnostics"]["maximum_refinement_corrections"] <= 3,
        sampled_indices={x["index"] for x in row["samples"]}
        == sample_indices(row["endpoint_bohr"], row["step_bohr"]),
        condition_observations=all(x["reciprocal_condition_1"] > 0 for x in row["samples"]),
        unitarity=row["invariants"]["unitarity"] <= budget,
        reciprocity=row["invariants"]["reciprocity"] <= budget,
    )
    if row["backend"] == "native":
        a = row["residual_audit"]
        checks.update(
            residuals_exact=a["evaluations"] == a["exact_residuals"],
            gate_decisions=a["gate_disagreements"] == 0,
            sampled_solves_audited=a["accepted_solves"] == len(row["samples"]),
            independent_backward_error=a["maximum_accepted_numpy_error"] <= 32 * 2.0**-63,
            precision=a["precision_preserved"],
            capture_count=len(row["captures"]) == 3,
        )
    return {key: bool(value) for key, value in checks.items()}


def comparisons(rows, native, old, budget):
    """Separate finite-step, backend, saved-increment and native-reference gates."""
    by_id = {r["id"]: r for r in rows}
    metrics, checks = {}, {}

    def compare(label, a, b):
        metrics[label] = differences(np.asarray(a), np.asarray(b))
        checks.update({f"{label}_{k}": v <= budget for k, v in metrics[label].items()})

    for identifier in ["native-200-fine", "reference-200-fine", "native-400-fine"]:
        if identifier in by_id:
            row = by_id[identifier]
            compare(identifier + "_native64", row["k"], native)
            previous = next(
                r
                for r in old["rows"]
                if (r["step_bohr"], r["endpoint_bohr"], r["tail_order"])
                == (0.00625, row["endpoint_bohr"], 8)
            )
            compare(identifier + "_increment", row["k"], previous["k"])
    if all(key in by_id for key in ["native-200-coarse", "native-200-medium", "native-200-fine"]):
        coarse, medium, fine = [
            by_id["native-200-" + key]["k"] for key in ["coarse", "medium", "fine"]
        ]
        metrics["coarse_step"] = differences(np.asarray(coarse), np.asarray(medium))
        compare("fine_step", medium, fine)
        checks["s_step_trend"] = metrics["fine_step"]["s"] <= max(
            metrics["coarse_step"]["s"] / 8, 1e-8
        )
    if all(key in by_id for key in ["native-200-fine", "reference-200-fine"]):
        left, right = [by_id[key] for key in ["native-200-fine", "reference-200-fine"]]
        compare("backend", left["k"], right["k"])
        checks["extended_k_exact"] = left["extended_k_sha256"] == right["extended_k_sha256"]
        checks["diagnostics_exact"] = left["diagnostics"] == right["diagnostics"]
        strip = lambda r: [  # noqa: E731
            {k: v for k, v in sample.items() if k != "reciprocal_condition_1"}
            for sample in r["samples"]
        ]
        checks["sampled_operands_exact"] = strip(left) == strip(right)
        metrics["auxiliary_rcond_relative_drift"] = condition_drift(left, right)
        checks["auxiliary_rcond_rounding"] = (
            metrics["auxiliary_rcond_relative_drift"] <= 32 * 2**-52
        )
    if all(key in by_id for key in ["native-200-fine", "native-400-fine"]):
        compare("endpoint", by_id["native-200-fine"]["k"], by_id["native-400-fine"]["k"])
    return {key: bool(v) for key, v in checks.items()}, metrics


def worker(task):
    raw, output, contract, config = task
    directory = output / config["id"]
    try:
        model = load_model(raw / contract["case"])
        energy = model["energies"][contract["energy_index"]]
        if float(energy).hex() != contract["energy_rydberg_hex"]:
            raise ValueError("Archived energy identity differs")
        q = surface_matrix(model["poles"], model["amplitudes"], model["ground"] + energy / 2)
        args = [q] + [model[key] for key in ["thresholds", "angular", "c2", "c3"]]
        row = propagate(args, energy, config, directory)
        row.update(energy_rydberg_hex=float(energy).hex(), case=contract["case"])
        row["checks"] = row_checks(row, contract["absolute_budget"])
        row["oracles"] = oracle_checks(directory) if config["backend"] == "native" else []
        if config["backend"] == "native":
            row["checks"]["oracle_count"] = len(row["oracles"]) == 3
        for i, oracle in enumerate(row["oracles"]):
            row["checks"].update({f"oracle_{i}_{k}": bool(v) for k, v in oracle["checks"].items()})
        row["experiment_verdict"] = all(row["checks"].values())
    except (
        OSError,
        ValueError,
        RuntimeError,
        AssertionError,
        LinAlgWarning,
        np.linalg.LinAlgError,
    ) as e:
        row = dict(config, experiment_verdict=None, failure=f"{type(e).__name__}: {e}")
    directory.mkdir(parents=True, exist_ok=True)
    save(directory / "report.json", row)
    return row


def verify(raw, prior, contract):
    raw_contract = Path(__file__).with_name("numerov-corner-contract.json")
    if digest(raw_contract) != contract["raw_contract_sha256"]:
        raise ValueError("Raw contract changed")
    verify_inputs(raw, json.loads(raw_contract.read_text()))
    for name, want in contract["prior_sha256"].items():
        if digest(prior / name) != want:
            raise ValueError(f"Prior evidence changed: {name}")
    for name, want in contract["prerequisite_sha256"].items():
        if digest(Path(__file__).with_name(name)) != want:
            raise ValueError(f"Qualification prerequisite changed: {name}")


def study(raw, prior, output, contract):
    output.mkdir()
    report = dict(
        contract=contract,
        rows=[],
        checks={},
        metrics={},
        experiment_verdict=None,
        original_interior_analytic_verdict=False,
        original_interior_scientific_exits=[1, 1],
        historical_unique_energy_coverage=132,
        new_qualified_co_identities=0,
        native_full_k_reconstruction_verdict=False,
        co_successor_submitted=False,
        native_ukrmol_invocations=0,
        new_inner_solves=0,
    )
    save(output / "report.json", report)
    try:
        import qscat_extended_kernels as native

        if (
            tuple(native.capability()) != (16, 64, 16384, True)
            or np.finfo(np.longdouble).nmant != 63
        ):
            raise ValueError("Unqualified scalar capability")
        verify(raw, prior, contract)
        old = json.loads((prior / "corner.json").read_text())
        if (
            old["verdict"] is not True
            or old["case"] != contract["case"]
            or old["energy_index"] != contract["energy_index"]
            or float(old["energy_rydberg"]).hex() != contract["energy_rydberg_hex"]
        ):
            raise ValueError("Prior corner qualification differs")
        energies, native_k, header = read_k(prior / "native64.bin")
        model = load_model(raw / contract["case"])
        if (
            native_k.shape != (1, 6, 6)
            or float(energies[0]).hex() != contract["energy_rydberg_hex"]
            or any(
                header[key] != value
                for key, value in model["header"].items()
                if key != "actual_energy_records"
            )
        ):
            raise ValueError("Native reference identity/header differs")
        with ProcessPoolExecutor(max_workers=contract["resources"]["workers"]) as pool:
            for stage in ["initial", "endpoint"]:
                configs = [r for r in contract["configurations"] if r["stage"] == stage]
                futures = [pool.submit(worker, (raw, output, contract, r)) for r in configs]
                for future in as_completed(futures):
                    row = future.result()
                    report["rows"].append(row)
                    save(output / "report.json", report)
                if any(r["experiment_verdict"] is not True for r in report["rows"]):
                    report["experiment_verdict"] = (
                        False
                        if any(r["experiment_verdict"] is False for r in report["rows"])
                        else None
                    )
                    return report
                checks, metrics = comparisons(
                    report["rows"], native_k[0], old, contract["absolute_budget"]
                )
                report.update(checks=checks, metrics=metrics)
                save(output / "report.json", report)
                if not all(checks.values()):
                    report["experiment_verdict"] = False
                    return report
        report["checks"].update(
            configuration_coverage={r["id"] for r in report["rows"]}
            == {r["id"] for r in contract["configurations"]},
            solve_count=sum(r["solve_calls"] for r in report["rows"]) == contract["mixed_solves"],
        )
        report["experiment_verdict"] = all(report["checks"].values())
    except (OSError, ValueError, RuntimeError, AssertionError, np.linalg.LinAlgError) as e:
        report["failure"] = f"{type(e).__name__}: {e}"
    finally:
        try:
            verify(raw, prior, contract)
            report["input_hashes_verified_after"] = True
        except (OSError, ValueError) as e:
            report.update(experiment_verdict=None, integrity_failure=str(e))
        save(output / "report.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["raw", "prior", "output", "contract"]:
        parser.add_argument("--" + name, type=Path, required=True)
    a = parser.parse_args()
    contract = json.loads(a.contract.read_text())
    if study(a.raw, a.prior, a.output, contract)["experiment_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
