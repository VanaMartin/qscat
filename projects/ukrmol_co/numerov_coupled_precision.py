"""Finite coupled-precision qualification and synthetic298-channel cost control."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import LinAlgWarning

from projects.ukrmol_co.coupled_tail import exact_tail, tail_coefficient, tail_waves
from projects.ukrmol_co.native_exterior_attribution import algebra, differences
from projects.ukrmol_co.native_exterior_study import cayley
from projects.ukrmol_co.numerov_arithmetic_control import (
    decimal_matrix,
    precision_capability,
)
from projects.ukrmol_co.numerov_coupled_oracle import continuous_k, direct_discrete_k, mp_tail
from projects.ukrmol_co.numerov_discrete_free import discrete_free_k
from projects.ukrmol_co.numerov_increment_study import exact_free_k
from projects.ukrmol_co.numerov_typed_coupled import refined_solve, typed_coupled_k, typed_tail

BACKENDS = [("binary64", np.float64), ("longdouble", np.longdouble)]


def mp_error(a, b):
    with mp.workdps(80):
        return float(max(abs(x - y) for x, y in zip(a, b, strict=True)))


def decimal_mp(matrix):
    with mp.workdps(80):
        return [[mp.nstr(x, 80) for x in row] for row in matrix.tolist()]


def row(q, thresholds, angular, c2, c3, energy, radius, endpoint, step, order, backend, dtype):
    begin = time.monotonic()
    k, diagnostics = typed_coupled_k(
        q, thresholds, angular, c2, c3, energy, radius, endpoint, step, order, dtype=dtype
    )
    decimals = decimal_matrix(k)
    return dict(
        backend=backend,
        k_decimal=decimals,
        k=np.asarray(k, float).tolist(),
        diagnostics=diagnostics,
        wall_seconds=time.monotonic() - begin,
        invariants=algebra(np.asarray(k, float)),
    )


def preflight(checks, report):
    report["precision_bits"] = precision_capability()
    for backend, dtype in BACKENDS:
        for n in [12, 298]:
            v = np.arange(1, n + 1, dtype=dtype) / 512
            a = np.eye(n, dtype=dtype) + np.outer(v, v) / 512
            rhs = np.column_stack([v, np.ones(n, dtype=dtype)])
            stats = dict(
                maximum_solve_backward_error=0.0,
                maximum_refinement_corrections=0,
                large_system_solves=0,
            )
            answer = refined_solve(a, rhs, stats)
            with mp.workdps(80):
                vm = mp.matrix([mp.mpf(int(i)) / 512 for i in range(1, n + 1)])
                bm = mp.matrix([[vm[i], 1] for i in range(n)])
                expected = bm - vm * (vm.T * bm) / (512 + (vm.T * vm)[0])
                error = mp_error(mp.matrix(decimal_matrix(answer)), expected)
            checks[f"algebra_{backend}_{n}"] = error <= (2e-15 if backend == "binary64" else 2e-17)
            checks[f"algebra_{backend}_{n}_residual"] = (
                stats["maximum_solve_backward_error"] <= 32 * np.finfo(dtype).eps
            )
            report["preflight"].append(
                dict(
                    kind="rank_one",
                    backend=backend,
                    channels=n,
                    forward_error=error,
                    diagnostics=stats,
                )
            )
        ts = np.array([0.0, 0.3, 1.8], dtype=dtype)
        ang = np.array([0, 1, 2])
        c2 = np.array([[2.0, -0.4, 0.1], [-0.4, 6.0, 0.2], [0.1, 0.2, 1.0]], dtype=dtype)
        c3 = np.array([[0.1, 0.3, -0.08], [0.3, -0.2, 0.12], [-0.08, 0.12, 0.05]], dtype=dtype)
        values, derivatives, diagnostics = typed_tail(
            ts, ang, dtype(1), c2, c3, 8, dtype(20), dtype(1e-10)
        )
        with mp.workdps(80):
            # The physical input is the original binary64 coefficient value.
            expected, prime = mp_tail(
                [mp.mpf(float(t)) for t in ts],
                list(ang),
                mp.mpf(1),
                mp.matrix(np.asarray(c2, float).tolist()),
                mp.matrix(np.asarray(c3, float).tolist()),
                8,
                mp.mpf(20),
            )
            errors = [
                mp_error(mp.matrix(decimal_matrix(a)), b)
                for a, b in zip((*values, *derivatives), (*expected, *prime), strict=True)
            ]
        checks[f"tail_{backend}_mp"] = max(errors) <= (2e-13 if backend == "binary64" else 2e-15)
        checks[f"tail_{backend}_bits"] = (
            diagnostics["complex_significand_bits"] == report["precision_bits"][backend]
        )
        dipole = np.asarray(tail_coefficient(1.2), dtype=dtype)
        val, der, _ = typed_tail(
            np.zeros(2, dtype=dtype),
            np.array([1, 2]),
            dtype(1),
            dipole,
            np.zeros((2, 2), dtype=dtype),
            12,
            dtype(20),
            dtype(1e-10),
        )
        (out, prime), _ = tail_waves(1.0, 20.0, 1.2)
        error = max(
            np.max(abs(np.asarray(val[1] + 1j * val[0], complex) - out)),
            np.max(abs(np.asarray(der[1] + 1j * der[0], complex) - prime)),
        )
        checks[f"tail_{backend}_hankel"] = error <= 1e-11
        report["preflight"].append(
            dict(
                kind="tail",
                backend=backend,
                mp_errors=errors,
                hankel_error=float(error),
                diagnostics=diagnostics,
            )
        )


def short_stage(c, checks, report, save):
    for model in c["models"]:
        q, ts, c2, c3 = (np.asarray(model[k]) for k in ["q", "thresholds", "c2", "c3"])
        ang = model["angular"]
        energy, radius, endpoint, order = (
            c[k] for k in ["energy_rydberg", "radius_bohr", "endpoint_bohr", "tail_order"]
        )
        refs = []
        for rtol, atol in c["continuous_tolerances"]:
            reference, _stats = continuous_k(
                q, ts, ang, c2, c3, energy, radius, endpoint, order, rtol, atol
            )
            refs.append(reference)
        stable = float(np.max(abs(refs[0] - refs[1])))
        key = model["name"]
        checks[f"{key}_continuous_reference"] = stable <= c["continuous_reference_stability_budget"]
        report["references"].append(dict(model=key, k=refs[1].tolist(), stability=stable))
        if not checks[f"{key}_continuous_reference"]:
            return
        for step in c["steps_bohr"]:
            low, high = [
                direct_discrete_k(q, ts, ang, c2, c3, energy, radius, endpoint, step, order, digits)
                for digits in c["oracle_digits"]
            ]
            precision = mp_error(low, high)
            checks[f"{key}_{step}_oracle_precision"] = precision < c["oracle_precision_budget"]
            if not checks[f"{key}_{step}_oracle_precision"]:
                return
            for backend, dtype in BACKENDS:
                result = row(
                    q, ts, ang, c2, c3, energy, radius, endpoint, step, order, backend, dtype
                )
                with mp.workdps(80):
                    error = mp_error(mp.matrix(result["k_decimal"]), high)
                metrics = differences(np.asarray(result["k"]), refs[1])
                prefix = f"{key}_{step}_{backend}"
                checks[prefix + "_arithmetic"] = error <= c["discrete_arithmetic_budgets"][backend]
                if step == c["steps_bohr"][-1]:
                    checks[prefix + "_continuous"] = all(
                        x <= c["finest_continuous_K_S_phase_budget"] for x in metrics.values()
                    )
                    invariant_reference = algebra(refs[1])
                    checks[prefix + "_invariant_reference"] = all(
                        abs(result["invariants"][k] - invariant_reference[k])
                        <= c["finest_invariant_reference_difference_budget"]
                        for k in invariant_reference
                    )
                result.update(
                    stage="short",
                    model=key,
                    step_bohr=step,
                    discrete_arithmetic_error=error,
                    discrete_k_decimal=decimal_mp(high),
                    discrete_oracle_precision=precision,
                    continuous_errors=metrics,
                )
                report["rows"].append(result)
                save()
        for backend, _ in BACKENDS:
            family = [
                r for r in report["rows"] if r.get("model") == key and r["backend"] == backend
            ]
            first, last = [r["continuous_errors"]["k"] for r in family[-2:]]
            checks[f"{key}_{backend}_trend"] = (
                last <= first / c["raw_K_error_reduction_minimum"]
                or max(first, last) <= c["both_errors_trend_floor"]
            )


def physical_stage(c, checks, report, save):
    c2 = tail_coefficient(c["dipole"])
    zero = np.zeros((2, 2))
    with mp.workdps(80):
        expected = np.array(exact_tail(1.0, c["dipole"]).tolist(), complex)
    for backend, dtype in BACKENDS:
        family = []
        for step in c["steps_bohr"]:
            result = row(
                zero,
                np.zeros(2),
                c["angular"],
                c2,
                zero,
                c["energy_rydberg"],
                c["radius_bohr"],
                c["endpoint_bohr"],
                step,
                c["tail_order"],
                backend,
                dtype,
            )
            error = float(np.max(abs(cayley(np.asarray(result["k"])) - expected)))
            if step == c["steps_bohr"][-1]:
                checks[backend + "_physical_S"] = error <= c["finest_S_and_invariant_budget"]
                checks[backend + "_physical_invariants"] = all(
                    x <= c["finest_S_and_invariant_budget"] for x in result["invariants"].values()
                )
            family.append(error)
            result.update(stage="physical", step_bohr=step, exact_S_error=error)
            report["rows"].append(result)
            save()
        checks[backend + "_physical_trend"] = (
            family[-1] <= family[-2] / c["raw_S_error_reduction_minimum"]
            or max(family[-2:]) <= c["both_errors_trend_floor"]
        )


def stress_stage(c, prior, checks, report, save):
    energy = float.fromhex(c["energy_rydberg_hex"])
    ratios = []
    for endpoint in c["endpoints_bohr"]:
        source = next(
            x
            for x in prior["oracle_controls"]
            if x["coordinate"]["coupled_closed"]
            and (x["coordinate"]["step_bohr"], x["coordinate"]["endpoint_bohr"])
            == (c["step_bohr"], endpoint)
        )
        q, ts = np.asarray(source["q"]), np.asarray(source["thresholds"])
        zero = np.zeros_like(q)
        discrete = discrete_free_k(q, ts, energy, c["radius_bohr"], endpoint, c["step_bohr"], 100)
        continuous = exact_free_k(q, ts, energy, c["radius_bohr"], 100)
        old = next(
            x
            for x in prior["rows"]
            if x["backend"] == "longdouble"
            and x["coordinate"]["coupled_closed"]
            and x["coordinate"]["endpoint_bohr"] == endpoint
            and x["coordinate"]["step_bohr"] == c["step_bohr"]
        )
        errors = {}
        for backend, dtype in BACKENDS:
            result = row(
                q,
                ts,
                [0] * len(q),
                zero,
                zero,
                energy,
                c["radius_bohr"],
                endpoint,
                c["step_bohr"],
                c["tail_order"],
                backend,
                dtype,
            )
            with mp.workdps(80):
                measured = mp.matrix(result["k_decimal"])
                error = mp_error(measured, discrete)
                continuous_error = mp_error(measured, continuous)
                old_error = mp_error(measured, mp.matrix(old["k_decimal"]))
            errors[backend] = error
            metrics = differences(np.asarray(result["k"]), np.array(continuous.tolist(), float))
            checks[f"{endpoint}_{backend}_invariants"] = all(
                x <= c["continuous_K_S_phase_and_invariant_budget"]
                for x in result["invariants"].values()
            )
            if backend == "longdouble":
                checks[f"{endpoint}_extended_arithmetic"] = error <= c["extended_roundoff_K_budget"]
                checks[f"{endpoint}_extended_continuous"] = continuous_error <= c[
                    "continuous_K_S_phase_and_invariant_budget"
                ] and all(
                    x <= c["continuous_K_S_phase_and_invariant_budget"] for x in metrics.values()
                )
                checks[f"{endpoint}_previous_extended"] = (
                    old_error <= c["previous_extended_K_difference_budget"]
                )
            result.update(
                stage="stress",
                endpoint_bohr=endpoint,
                discrete_arithmetic_error=error,
                continuous_error=continuous_error,
                previous_extended_difference=old_error,
                continuous_errors=metrics,
                discrete_k_decimal=decimal_mp(discrete),
                continuous_k_decimal=decimal_mp(continuous),
            )
            report["rows"].append(result)
            save()
        ratios.append(errors["binary64"] / errors["longdouble"])
    report["hypotheses"] = dict(
        all_three_precision_gain_at_least_100fold=all(
            x >= c["hypothesis_all_three_precision_gain_minimum"] for x in ratios
        ),
        roundoff_reduction_ratios=ratios,
    )


def dense_stage(c, checks, report, save):
    n = c["channels"]
    q = np.eye(n) * c["diagonal_q_bohr"]
    active = c["active_indices"]
    q[np.ix_(active, active)] = c["active_q"]
    for i, value in enumerate(c["other_open_q"], 2):
        q[i, i] = value
    ts = np.array(
        [0.0] * c["open_channels"] + [c["closed_threshold_rydberg"]] * (n - c["open_channels"])
    )
    zero = np.zeros((n, n))
    energy, radius, endpoint, step = (
        c[k] for k in ["energy_rydberg", "radius_bohr", "endpoint_bohr", "step_bohr"]
    )
    with mp.workdps(100):
        expected = mp.zeros(c["open_channels"])
        block = discrete_free_k(
            np.asarray(c["active_q"]), ts[active], energy, radius, endpoint, step, 100
        )
        for i in range(2):
            for j in range(2):
                expected[i, j] = block[i, j]
        for i, value in enumerate(c["other_open_q"], 2):
            expected[i, i] = discrete_free_k([[value]], [0.0], energy, radius, endpoint, step, 100)[
                0, 0
            ]
    costs = {}
    for backend, dtype in BACKENDS:
        result = row(
            q,
            ts,
            np.zeros(n, int),
            zero,
            zero,
            energy,
            radius,
            endpoint,
            step,
            c["tail_order"],
            backend,
            dtype,
        )
        with mp.workdps(80):
            error = mp_error(mp.matrix(result["k_decimal"]), expected)
        checks[backend + "_dense_arithmetic"] = error <= c["discrete_K_budgets"][backend]
        checks[backend + "_dense_invariants"] = all(
            x <= c["invariant_budget"] for x in result["invariants"].values()
        )
        checks[backend + "_dense_backward_residual"] = (
            result["diagnostics"]["maximum_solve_backward_error"] <= 32 * np.finfo(dtype).eps
        )
        result.update(
            stage="dense", discrete_arithmetic_error=error, discrete_k_decimal=decimal_mp(expected)
        )
        report["rows"].append(result)
        costs[backend] = result["wall_seconds"]
        save()
    report["synthetic_cost"] = dict(
        extended_to_binary64_wall_ratio=costs["longdouble"] / costs["binary64"],
        extended_wall_seconds=costs["longdouble"],
        binary64_wall_seconds=costs["binary64"],
        intervals=8,
        extended_seconds_per_interval=costs["longdouble"] / 8,
        rough_60_solve_400bohr_finest_seconds=60 * (400 - 18) / 0.0015625 * costs["longdouble"] / 8,
        scope=(
            "Synthetic298 free matrices, eight intervals, single worker. "
            "Includes startup/tail/matching; not a measured CO cost or successor authorization."
        ),
    )


def study(prior_path, contract, output):
    report = dict(
        contract=contract,
        experiment_verdict=None,
        checks={},
        stages=[],
        preflight=[],
        references=[],
        rows=[],
        hypotheses={},
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
        if hashlib.sha256(prior_path.read_bytes()).hexdigest() != contract["prior_report_sha256"]:
            raise ValueError("Archived arithmetic report hash differs")
        prior = json.loads(prior_path.read_text())
        if prior["experiment_verdict"] is not True:
            raise ValueError("Arithmetic prerequisite did not pass")
        for stage, run in [
            ("preflight", lambda checks: preflight(checks, report)),
            ("short", lambda checks: short_stage(contract["short"], checks, report, save)),
            ("physical", lambda checks: physical_stage(contract["physical"], checks, report, save)),
            (
                "stress",
                lambda checks: stress_stage(contract["stress"], prior, checks, report, save),
            ),
            ("dense", lambda checks: dense_stage(contract["dense"], checks, report, save)),
        ]:
            checks = {}
            run(checks)
            checks = {key: bool(value) for key, value in checks.items()}
            report["checks"].update({stage + "_" + k: v for k, v in checks.items()})
            verdict = all(checks.values())
            report["stages"].append(
                dict(
                    stage=stage,
                    verdict=verdict,
                    checks_passed=sum(checks.values()),
                    checks_total=len(checks),
                )
            )
            save()
            if not verdict:
                report["experiment_verdict"] = False
                return report
        report["experiment_verdict"] = True
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
