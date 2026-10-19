"""Finite native-residual integration controls using retained synthetic references."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from functools import partial
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import LinAlgWarning

from projects.ukrmol_co.coupled_tail import exact_tail, tail_coefficient
from projects.ukrmol_co.native_exterior_attribution import algebra, differences
from projects.ukrmol_co.native_exterior_study import cayley
from projects.ukrmol_co.numerov_arithmetic_control import decimal_matrix
from projects.ukrmol_co.numerov_coupled_precision import mp_error
from projects.ukrmol_co.numerov_dense_scaling import solve_roles, timed_propagation
from projects.ukrmol_co.numerov_native_coupled import (
    REFERENCE_SOLVE,
    ResidualAudit,
    installed_solver,
    native_refined_solve,
)
from projects.ukrmol_co.numerov_typed_coupled import typed_coupled_k


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def embed(q, thresholds, angular, c2, c3, energy):
    """Append uncoupled closed spectators to force the four-channel mixed path."""
    n = len(thresholds)
    if not 1 <= n <= 3:
        raise ValueError("Only the declared small reference embeddings are supported")
    matrices = []
    for source, diagonal in [(q, 0.3), (c2, 0.0), (c3, 0.0)]:
        out = np.eye(4) * diagonal
        out[:n, :n] = source
        matrices.append(out)
    ts = np.r_[thresholds, [energy + 2 + i for i in range(4 - n)]]
    return matrices[0], ts, np.r_[angular, np.zeros(4 - n, int)], *matrices[1:]


def propagation(args, coordinate, backend, audit_enabled):
    audit = ResidualAudit() if audit_enabled else None
    solver = (
        partial(native_refined_solve, observer=audit) if backend == "native" else REFERENCE_SOLVE
    )
    begin = time.perf_counter()
    with installed_solver(solver):
        k, diagnostics = typed_coupled_k(*args, *coordinate, dtype=np.longdouble)
    return dict(
        backend=backend,
        k_decimal=decimal_matrix(k),
        k=np.asarray(k, float).tolist(),
        diagnostics=diagnostics,
        invariants=algebra(np.asarray(k, float)),
        wall_seconds=time.perf_counter() - begin,
        residual_audit=audit.summary() if audit else None,
    )


def pair_checks(row, baseline, expected_solves):
    diagnostics = row["diagnostics"]
    checks = dict(
        mixed_solve_count=diagnostics["large_system_solves"] == expected_solves,
        backward_error=diagnostics["maximum_solve_backward_error"]
        <= 32 * np.finfo(np.longdouble).eps,
        corrections=diagnostics["maximum_refinement_corrections"] <= 3,
        precision=diagnostics["significand_bits"] == diagnostics["complex_significand_bits"] == 64,
    )
    if baseline is not None:
        checks["paired_k_exact"] = row["k_decimal"] == baseline["k_decimal"]
        checks["paired_diagnostics_exact"] = diagnostics == baseline["diagnostics"]
        checks["paired_invariants_exact"] = row["invariants"] == baseline["invariants"]
    if row["residual_audit"] is not None:
        audit = row["residual_audit"]
        checks.update(
            residuals_exact=audit["evaluations"] == audit["exact_residuals"],
            residual_decisions=audit["gate_disagreements"] == 0,
            accepted_numpy_error=audit["maximum_accepted_numpy_error"]
            <= 32 * np.finfo(np.longdouble).eps,
            audited_solves=audit["accepted_solves"] == expected_solves,
            audit_precision=audit["precision_preserved"],
        )
    return checks


def timing_summary(rows, minimum):
    summaries = []
    for model in ["free", "varying"]:
        means = {}
        for backend in ["reference", "native"]:
            for n in [32, 128]:
                selected = [
                    r
                    for r in rows
                    if (r["model"], r["backend"], r["intervals"]) == (model, backend, n)
                ]
                if len(selected) != 2:
                    raise ValueError("Incomplete fixed timing family")
                means[backend, n] = float(np.mean([r["kernel_wall_seconds"] for r in selected]))
        marginal = {b: (means[b, 128] - means[b, 32]) / 96 for b in ["reference", "native"]}
        speedup = marginal["reference"] / marginal["native"] if marginal["native"] > 0 else None
        summaries.append(
            dict(
                model=model,
                mean_kernel_seconds={
                    b: [means[b, n] for n in [32, 128]] for b in ["reference", "native"]
                },
                marginal_seconds_per_interval=marginal,
                marginal_speedup=speedup,
                n128_total_speedup=means["reference", 128] / means["native", 128],
                hypothesis_native_marginal_speedup_at_least_minimum=speedup is not None
                and speedup >= minimum,
            )
        )
    return summaries


def study(c, prior_directory, output):
    report = dict(
        contract=c,
        experiment_verdict=None,
        checks={},
        rows=[],
        stages=[],
        timing_analysis=[],
        original_interior_analytic_verdict=False,
        original_interior_scientific_exits=[1, 1],
        co_successor_submitted=False,
        new_co_independent_solves=0,
        new_inner_solves=0,
        native_full_k_reconstruction_verdict=False,
        co_complex_pole_verdict=None,
        unique_energy_coverage=132,
        uncomputed_identities=462,
    )
    begin = time.monotonic()

    def save():
        temporary = output.with_suffix(".partial")
        temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        temporary.replace(output)

    def append(row, checks):
        row["checks"] = {k: bool(v) for k, v in checks.items()}
        prefix = str(len(report["rows"]))
        report["checks"].update({prefix + "_" + k: v for k, v in row["checks"].items()})
        report["rows"].append(row)
        save()

    def paired(stage, label, args, coordinate, old, reference=None):
        baseline = None
        n = round((coordinate[2] - coordinate[1]) / coordinate[3])
        roles = solve_roles("free" if not np.any(args[3]) and not np.any(args[4]) else "varying", n)
        for backend in ["reference", "native"]:
            row = propagation(args, coordinate, backend, backend == "native")
            row.update(stage=stage, model=label, intervals=n, coordinate=list(coordinate))
            checks = pair_checks(row, baseline, len(roles))
            with mp.workdps(80):
                measured = mp.matrix(row["k_decimal"])
                row["archived_k_difference"] = mp_error(measured, mp.matrix(old["k_decimal"]))
                if "discrete_k_decimal" in old:
                    row["discrete_arithmetic_error"] = mp_error(
                        measured, mp.matrix(old["discrete_k_decimal"])
                    )
                    checks["discrete_arithmetic"] = (
                        row["discrete_arithmetic_error"] <= c["budgets"][stage + "_discrete_K"]
                    )
            checks["archived_k"] = (
                row["archived_k_difference"] <= c["budgets"][stage + "_archived_K"]
            )
            if reference is not None:
                row["continuous_errors"] = differences(np.asarray(row["k"]), reference)
                if (
                    stage == "stress"
                    or coordinate[3] == coupled["contract"]["short"]["steps_bohr"][-1]
                ):
                    checks["continuous"] = all(v <= 1e-6 for v in row["continuous_errors"].values())
                    expected_invariants = algebra(reference)
                    checks["reference_invariants"] = all(
                        abs(row["invariants"][k] - v) <= (1e-9 if stage == "short" else 1e-6)
                        for k, v in expected_invariants.items()
                    )
            if stage == "physical":
                with mp.workdps(80):
                    expected = np.array(
                        exact_tail(1.0, coupled["contract"]["physical"]["dipole"]).tolist(), complex
                    )
                row["exact_S_error"] = float(np.max(abs(cayley(np.asarray(row["k"])) - expected)))
                if coordinate[3] == coupled["contract"]["physical"]["steps_bohr"][-1]:
                    checks["physical_S"] = row["exact_S_error"] <= 1e-6
                    checks["physical_invariants"] = all(
                        v <= 1e-6 for v in row["invariants"].values()
                    )
            append(row, checks)
            baseline = row
        return all(r["checks"][k] for r in report["rows"][-2:] for k in r["checks"])

    save()
    try:
        import qscat_extended_kernels as native

        if (
            tuple(native.capability()) != (16, 64, 16384, True)
            or np.finfo(np.longdouble).nmant + 1 != 64
        ):
            raise ValueError("Unsupported native/NumPy extended scalar")
        priors = {}
        for name, want in c["prior_report_sha256"].items():
            path = prior_directory / (name + ".json")
            if digest(path) != want:
                raise ValueError("Prior report hash differs: " + name)
            priors[name] = json.loads(path.read_text())
            if priors[name]["experiment_verdict"] is not True:
                raise ValueError("Prior report did not pass: " + name)
        coupled, dense, arithmetic = [priors[k] for k in ["coupled", "dense", "arithmetic"]]
        for stage in ["short", "physical", "dense_accuracy", "stress", "timing"]:
            first = len(report["rows"])
            passed = True
            if stage == "short":
                cfg = coupled["contract"]["short"]
                for model in cfg["models"]:
                    args = embed(
                        *(model[k] for k in ["q", "thresholds", "angular", "c2", "c3"]),
                        cfg["energy_rydberg"],
                    )
                    ref = np.asarray(
                        next(r for r in coupled["references"] if r["model"] == model["name"])["k"]
                    )
                    for h in cfg["steps_bohr"]:
                        old = next(
                            r
                            for r in coupled["rows"]
                            if r.get("stage") == "short"
                            and r["backend"] == "longdouble"
                            and r["model"] == model["name"]
                            and r["step_bohr"] == h
                        )
                        passed &= paired(
                            stage,
                            model["name"],
                            args,
                            [
                                cfg["energy_rydberg"],
                                cfg["radius_bohr"],
                                cfg["endpoint_bohr"],
                                h,
                                cfg["tail_order"],
                            ],
                            old,
                            ref,
                        )
            elif stage == "physical":
                cfg = coupled["contract"]["physical"]
                zero = np.zeros((2, 2))
                args = embed(
                    zero,
                    np.zeros(2),
                    cfg["angular"],
                    tail_coefficient(cfg["dipole"]),
                    zero,
                    cfg["energy_rydberg"],
                )
                for h in cfg["steps_bohr"]:
                    old = next(
                        r
                        for r in coupled["rows"]
                        if r.get("stage") == "physical"
                        and r["backend"] == "longdouble"
                        and r["step_bohr"] == h
                    )
                    passed &= paired(
                        stage,
                        "inverse_square",
                        args,
                        [
                            cfg["energy_rydberg"],
                            cfg["radius_bohr"],
                            cfg["endpoint_bohr"],
                            h,
                            cfg["tail_order"],
                        ],
                        old,
                    )
            elif stage == "stress":
                cfg = coupled["contract"]["stress"]
                energy = float.fromhex(cfg["energy_rydberg_hex"])
                for endpoint in cfg["endpoints_bohr"]:
                    source = next(
                        r
                        for r in arithmetic["oracle_controls"]
                        if r["coordinate"]["coupled_closed"]
                        and r["coordinate"]["step_bohr"] == cfg["step_bohr"]
                        and r["coordinate"]["endpoint_bohr"] == endpoint
                    )
                    zero = np.zeros_like(source["q"])
                    args = embed(
                        source["q"],
                        source["thresholds"],
                        np.zeros(len(zero), int),
                        zero,
                        zero,
                        energy,
                    )
                    old = next(
                        r
                        for r in coupled["rows"]
                        if r.get("stage") == "stress"
                        and r["backend"] == "longdouble"
                        and r["endpoint_bohr"] == endpoint
                    )
                    ref = np.array(old["continuous_k_decimal"], float)
                    passed &= paired(
                        stage,
                        str(endpoint),
                        args,
                        [energy, cfg["radius_bohr"], endpoint, cfg["step_bohr"], cfg["tail_order"]],
                        old,
                        ref,
                    )
            else:
                cfg = dense["contract"]
                for repetition in range(1 if stage == "dense_accuracy" else 2):
                    models = ["free", "varying"] if repetition == 0 else ["varying", "free"]
                    for model in models:
                        for n in (
                            [8]
                            if stage == "dense_accuracy"
                            else ([32, 128] if repetition == 0 else [128, 32])
                        ):
                            baseline = None
                            backends = (
                                ["reference", "native"]
                                if repetition == 0
                                else ["native", "reference"]
                            )
                            pair = []
                            for backend in backends:
                                audit = (
                                    ResidualAudit()
                                    if stage == "dense_accuracy" and backend == "native"
                                    else None
                                )
                                solver = (
                                    partial(native_refined_solve, observer=audit)
                                    if backend == "native"
                                    else REFERENCE_SOLVE
                                )
                                with installed_solver(solver):
                                    row = timed_propagation(cfg, model, n, "longdouble")
                                row.update(
                                    stage=stage,
                                    backend=backend,
                                    repetition=repetition,
                                    residual_audit=audit.summary() if audit else None,
                                )
                                old = next(
                                    r
                                    for r in dense["rows"]
                                    if r["model"] == model
                                    and r["intervals"] == n
                                    and r["backend"] == "longdouble"
                                    and r["repetition"] == 0
                                )
                                checks = pair_checks(row, baseline, len(solve_roles(model, n)))
                                checks["archived_k_exact"] = row["k_decimal"] == old["k_decimal"]
                                checks["archived_diagnostics_exact"] = (
                                    row["diagnostics"] == old["diagnostics"]
                                )
                                checks["role_sequence"] = [
                                    r["role"]
                                    for r in row["calls"]
                                    if r["function"] == "refined_solve"
                                ] == solve_roles(model, n)
                                checks["timer_partition"] = (
                                    row["timed_function_seconds"]
                                    <= row["kernel_wall_seconds"] + 1e-6
                                )
                                if stage == "dense_accuracy":
                                    oracle = next(
                                        r
                                        for r in dense["oracles"]
                                        if r["model"] == model and r["intervals"] == n
                                    )
                                    with mp.workdps(80):
                                        row["discrete_arithmetic_error"] = mp_error(
                                            mp.matrix(row["k_decimal"]),
                                            mp.matrix(oracle["k_decimal"]),
                                        )
                                    checks["discrete_arithmetic"] = (
                                        row["discrete_arithmetic_error"] <= 1e-12
                                    )
                                append(row, checks)
                                baseline = row
                                pair.append(row)
                            passed &= all(v for r in pair for v in r["checks"].values())
            if stage in ["short", "physical"]:
                for model in sorted({r["model"] for r in report["rows"][first:]}):
                    for backend in ["reference", "native"]:
                        family = [
                            r
                            for r in report["rows"][first:]
                            if r["model"] == model and r["backend"] == backend
                        ]
                        errors = [
                            r["continuous_errors"]["k"] if stage == "short" else r["exact_S_error"]
                            for r in family
                        ]
                        trend = errors[-1] <= errors[-2] / 8 or max(errors[-2:]) <= 1e-10
                        report["checks"][f"{stage}_{model}_{backend}_trend"] = bool(trend)
                        passed &= trend
            report["stages"].append(
                dict(
                    stage=stage,
                    first_row=first,
                    rows=len(report["rows"]) - first,
                    verdict=bool(passed),
                )
            )
            save()
            if not passed:
                report["experiment_verdict"] = False
                return report
        report["timing_analysis"] = timing_summary(
            [r for r in report["rows"] if r["stage"] == "timing"],
            c["timing_speedup_hypothesis_minimum"],
        )
        report["checks"]["declared_row_count"] = len(report["rows"]) == c["new_propagations"]
        report["checks"]["declared_interval_count"] = (
            sum(r["intervals"] for r in report["rows"]) == c["new_radial_intervals"]
        )
        report["checks"]["declared_large_solve_count"] = (
            sum(r["diagnostics"]["large_system_solves"] for r in report["rows"])
            == c["expected_large_system_solves"]
        )
        report["checks"]["declared_audited_solve_count"] = (
            sum(
                r["residual_audit"]["accepted_solves"]
                for r in report["rows"]
                if r["residual_audit"] is not None
            )
            == c["accuracy_native_audited_solves"]
        )
        report["experiment_verdict"] = all(report["checks"].values())
    except (
        OSError,
        ValueError,
        RuntimeError,
        AssertionError,
        np.linalg.LinAlgError,
        ZeroDivisionError,
        ImportError,
        LinAlgWarning,
    ) as error:
        report.update(experiment_verdict=None, failure=f"{type(error).__name__}: {error}")
    finally:
        report["wall_seconds"] = time.monotonic() - begin
        save()
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ["contract", "priors", "output"]:
        parser.add_argument("--" + key, type=Path, required=True)
    args = parser.parse_args()
    report = study(json.loads(args.contract.read_text()), args.priors, args.output)
    if report["experiment_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
