"""Bounded synthetic dense timing attribution using unchanged numerical modules."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import LinAlgWarning

import projects.ukrmol_co.numerov_typed_coupled as kernel
from projects.ukrmol_co.native_exterior_attribution import algebra
from projects.ukrmol_co.numerov_arithmetic_control import decimal_matrix, precision_capability
from projects.ukrmol_co.numerov_coupled_oracle import direct_discrete_k
from projects.ukrmol_co.numerov_coupled_precision import decimal_mp, mp_error
from projects.ukrmol_co.numerov_discrete_free import discrete_free_k


def inputs(c, model):
    n = c["channels"]
    q = np.eye(n) * c["diagonal_q_bohr"]
    active = c["active_indices"]
    q[np.ix_(active, active)] = c["active_q"]
    for i, value in enumerate(c["other_open_q"], 2):
        q[i, i] = value
    thresholds = np.array(
        [0.0] * c["open_channels"] + [c["closed_threshold_rydberg"]] * (n - c["open_channels"])
    )
    c2, c3 = np.zeros((n, n)), np.zeros((n, n))
    if model == "varying":
        c2[np.ix_(active, active)] = c["varying_active_c2"]
        c3[np.ix_(active, active)] = c["varying_active_c3"]
    return q, thresholds, np.zeros(n, int), c2, c3


def oracle(c, model, intervals, digits):
    endpoint = c["radius_bohr"] + intervals * c["step_bohr"]
    with mp.workdps(digits):
        expected = mp.zeros(c["open_channels"])
        ts = [0.0, 0.0, c["closed_threshold_rydberg"]]
        if model == "free":
            block = discrete_free_k(
                c["active_q"],
                ts,
                c["energy_rydberg"],
                c["radius_bohr"],
                endpoint,
                c["step_bohr"],
                digits,
            )
        else:
            block = direct_discrete_k(
                c["active_q"],
                ts,
                [0, 0, 0],
                c["varying_active_c2"],
                c["varying_active_c3"],
                c["energy_rydberg"],
                c["radius_bohr"],
                endpoint,
                c["step_bohr"],
                c["tail_orders"][model],
                digits,
            )
        for i in range(2):
            for j in range(2):
                expected[i, j] = block[i, j]
        for i, q in enumerate(c["other_open_q"], 2):
            expected[i, i] = discrete_free_k(
                [[q]],
                [0.0],
                c["energy_rydberg"],
                c["radius_bohr"],
                endpoint,
                c["step_bohr"],
                digits,
            )[0, 0]
        return +expected


def solve_roles(model, intervals):
    if model == "free":
        return (
            ["initialization", "cached_potential"]
            + ["recurring_defect"] * (intervals - 1)
            + ["matching"]
        )
    return (
        ["initialization"]
        + ["recurring_potential", "recurring_defect"] * (intervals - 1)
        + ["matching"]
    )


def timed_propagation(c, model, intervals, backend):
    dtype = np.float64 if backend == "binary64" else np.longdouble
    original = {
        name: getattr(kernel, name) for name in ["refined_solve", "taylor_change", "typed_tail"]
    }
    roles = solve_roles(model, intervals)
    calls = []
    solve_index = 0

    def wrapped(name):
        def measured(*args, **kwargs):
            nonlocal solve_index
            role = name
            if name == "refined_solve":
                if solve_index >= len(roles):
                    raise ValueError("Numerical kernel made an undeclared solve")
                role = roles[solve_index]
                solve_index += 1
            begin = time.perf_counter()
            result = original[name](*args, **kwargs)
            elapsed = time.perf_counter() - begin
            array = result[0][0] if name == "typed_tail" else result
            calls.append(
                dict(
                    function=name,
                    role=role,
                    wall_seconds=elapsed,
                    result_dtype=np.asarray(array).dtype.name,
                )
            )
            return result

        return measured

    for name in original:
        setattr(kernel, name, wrapped(name))
    args = inputs(c, model)
    begin = time.perf_counter()
    try:
        k, diagnostics = kernel.typed_coupled_k(
            *args,
            c["energy_rydberg"],
            c["radius_bohr"],
            c["radius_bohr"] + intervals * c["step_bohr"],
            c["step_bohr"],
            c["tail_orders"][model],
            dtype=dtype,
        )
    finally:
        elapsed = time.perf_counter() - begin
        for name, value in original.items():
            setattr(kernel, name, value)
    totals = {}
    for call in calls:
        totals[call["role"]] = totals.get(call["role"], 0.0) + call["wall_seconds"]
    recurring = sum(v for key, v in totals.items() if key.startswith("recurring_"))
    return dict(
        model=model,
        intervals=intervals,
        backend=backend,
        k_decimal=decimal_matrix(k),
        k=np.asarray(k, float).tolist(),
        invariants=algebra(np.asarray(k, float)),
        diagnostics=diagnostics,
        kernel_wall_seconds=elapsed,
        timed_function_seconds=sum(totals.values()),
        timing_by_role=totals,
        recurring_solve_seconds=recurring,
        setup_solve_seconds=sum(
            v
            for key, v in totals.items()
            if not key.startswith("recurring_") and key not in ["taylor_change", "typed_tail"]
        ),
        untimed_remainder_seconds=elapsed - sum(totals.values()),
        calls=calls,
    )


def timing_analysis(rows, c):
    fits = []
    for model in c["models"]:
        for backend in c["backends"]:
            means, spreads = [], []
            counts = c["interval_counts"]
            for n in counts:
                pair = [
                    r["kernel_wall_seconds"]
                    for r in rows
                    if (r["model"], r["backend"], r["intervals"]) == (model, backend, n)
                ]
                if len(pair) != c["repetitions"]:
                    raise ValueError("Timing fit lacks a declared repetition")
                means.append(float(np.mean(pair)))
                spreads.append(max(pair) - min(pair))
            x, y = np.asarray(counts, float), np.asarray(means)
            beta = float(np.sum((x - x.mean()) * (y - y.mean())) / np.sum((x - x.mean()) ** 2))
            alpha = float(y.mean() - beta * x.mean())
            residual = float(np.max(abs(alpha + beta * x - y) / y))
            marginal = (means[-1] - means[-2]) / (counts[-1] - counts[-2])
            long_rows = [
                r
                for r in rows
                if (r["model"], r["backend"], r["intervals"]) == (model, backend, counts[-1])
            ]
            recurring = float(
                np.mean([r["recurring_solve_seconds"] / (counts[-1] - 1) for r in long_rows])
            )
            fits.append(
                dict(
                    model=model,
                    backend=backend,
                    mean_kernel_wall_seconds=means,
                    repetition_spread_seconds=spreads,
                    affine_intercept_seconds=alpha,
                    affine_slope_seconds_per_interval=beta,
                    maximum_relative_affine_residual=residual,
                    marginal_32_to_128_seconds_per_interval=marginal,
                    direct_recurring_solve_seconds_per_interval=recurring,
                    synthetic_60_finest400_affine_days=60
                    * (alpha + beta * (400 - 18) / 0.0015625)
                    / 86400,
                    synthetic_60_finest400_marginal_days=60
                    * marginal
                    * (400 - 18)
                    / 0.0015625
                    / 86400,
                    synthetic_60_finest400_direct_recurring_days=60
                    * recurring
                    * (400 - 18)
                    / 0.0015625
                    / 86400,
                )
            )
    extended = {r["model"]: r for r in fits if r["backend"] == "longdouble"}
    ratio = (
        extended["varying"]["marginal_32_to_128_seconds_per_interval"]
        / extended["free"]["marginal_32_to_128_seconds_per_interval"]
    )
    free = [
        r for r in rows if (r["model"], r["backend"], r["intervals"]) == ("free", "longdouble", 128)
    ]
    recurring = float(np.mean([r["recurring_solve_seconds"] for r in free]))
    fixed = float(
        np.mean(
            [r["timing_by_role"]["taylor_change"] + r["timing_by_role"]["typed_tail"] for r in free]
        )
    )
    hypotheses = dict(
        all_affine_fits_positive_and_within_20percent=all(
            r["affine_slope_seconds_per_interval"] > 0
            and r["maximum_relative_affine_residual"] <= 0.2
            for r in fits
        ),
        free_extended_N128_recurring_solves_exceed_Taylor_plus_tail=recurring > fixed,
        varying_to_free_extended_marginal_rate_between_1_and_3=1 <= ratio <= 3,
    )
    return dict(
        fits=fits,
        hypotheses=hypotheses,
        varying_to_free_extended_marginal_ratio=ratio,
        free_extended_N128_recurring_solve_seconds=recurring,
        free_extended_N128_Taylor_plus_tail_seconds=fixed,
        scope=(
            "Synthetic298 matrices and exercised conditioning, one worker. "
            "Setup/interval separation does not establish native CO throughput, "
            "conditioning or an absolute-error bound."
        ),
    )


def study(prior_path, c, output):
    report = dict(
        contract=c,
        experiment_verdict=None,
        checks={},
        oracles=[],
        rows=[],
        timing_analysis={},
        original_interior_analytic_verdict=False,
        original_interior_scientific_exits=[1, 1],
        co_successor_submitted=False,
        new_co_independent_solves=0,
        new_unique_energy_identities=0,
        candidate_qualified_identities=14,
        unique_energy_coverage=132,
        uncomputed_identities=462,
        native_invocations=0,
        new_inner_solves=0,
        native_full_k_reconstruction_verdict=False,
        co_complex_pole_verdict=None,
        full_sweep_released=False,
    )
    begin = time.monotonic()

    def save():
        (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    save()
    try:
        report["precision_bits"] = precision_capability()
        if hashlib.sha256(prior_path.read_bytes()).hexdigest() != c["prior_report_sha256"]:
            raise ValueError("Prior coupled report hash differs")
        prior = json.loads(prior_path.read_text())
        if prior["experiment_verdict"] is not True:
            raise ValueError("Coupled prerequisite did not pass")
        refs = {}
        for model in c["models"]:
            for n in c["interval_counts"]:
                low, high = [oracle(c, model, n, d) for d in c["oracle_digits"]]
                drift = mp_error(low, high)
                report["checks"][f"oracle_{model}_{n}"] = drift < c["oracle_precision_budget"]
                report["oracles"].append(
                    dict(
                        model=model, intervals=n, k_decimal=decimal_mp(high), precision_drift=drift
                    )
                )
                refs[model, n] = high
                save()
                if not report["checks"][f"oracle_{model}_{n}"]:
                    report["experiment_verdict"] = False
                    return report
        for repetition in range(c["repetitions"]):
            models = c["models"] if repetition == 0 else c["models"][::-1]
            counts = c["interval_counts"] if repetition == 0 else c["interval_counts"][::-1]
            backends = c["backends"] if repetition == 0 else c["backends"][::-1]
            for model in models:
                for n in counts:
                    coordinate_checks = []
                    for backend in backends:
                        result = timed_propagation(c, model, n, backend)
                        result["repetition"] = repetition
                        with mp.workdps(80):
                            error = mp_error(mp.matrix(result["k_decimal"]), refs[model, n])
                        result["discrete_arithmetic_error"] = error
                        prefix = f"{repetition}_{model}_{n}_{backend}"
                        roles = solve_roles(model, n)
                        diagnostics = result["diagnostics"]
                        checks = dict(
                            arithmetic=error <= c["arithmetic_K_budgets"][backend],
                            invariants=all(
                                x <= c["invariant_budget"] for x in result["invariants"].values()
                            ),
                            solve_count=diagnostics["large_system_solves"] == len(roles)
                            and len(
                                [x for x in result["calls"] if x["function"] == "refined_solve"]
                            )
                            == len(roles),
                            fixed_call_counts=sum(
                                x["function"] == "taylor_change" for x in result["calls"]
                            )
                            == 3
                            and sum(x["function"] == "typed_tail" for x in result["calls"]) == 1,
                            timer_partition=result["timed_function_seconds"]
                            <= result["kernel_wall_seconds"]
                            + c["disjoint_timer_sum_slack_seconds"],
                            backward_error=diagnostics["maximum_solve_backward_error"]
                            <= 32 * 2.0 ** (1 - report["precision_bits"][backend]),
                        )
                        if model == "free" and n == 8:
                            old = next(
                                x
                                for x in prior["rows"]
                                if x["stage"] == "dense" and x["backend"] == backend
                            )
                            checks["old_baseline_exact"] = result["k_decimal"] == old["k_decimal"]
                        if repetition:
                            old = next(
                                x
                                for x in report["rows"]
                                if (x["model"], x["intervals"], x["backend"], x["repetition"])
                                == (model, n, backend, 0)
                            )
                            checks["repetition_exact"] = result["k_decimal"] == old["k_decimal"]
                        report["checks"].update(
                            {prefix + "_" + key: bool(value) for key, value in checks.items()}
                        )
                        coordinate_checks.extend(checks.values())
                        report["rows"].append(result)
                        save()
                    if not all(coordinate_checks):
                        report["experiment_verdict"] = False
                        return report
        report["timing_analysis"] = timing_analysis(report["rows"], c)
        report["experiment_verdict"] = all(report["checks"].values())
    except (
        OSError,
        ValueError,
        RuntimeError,
        AssertionError,
        np.linalg.LinAlgError,
        ZeroDivisionError,
        LinAlgWarning,
    ) as error:
        report.update(experiment_verdict=None, failure=f"{type(error).__name__}: {error}")
    finally:
        report["wall_seconds"] = time.monotonic() - begin
        save()
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ["prior", "contract", "output"]:
        parser.add_argument("--" + key, type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir()
    report = study(args.prior, json.loads(args.contract.read_text()), args.output)
    if report["experiment_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
