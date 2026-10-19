"""Bounded saved-CO recurrence prefixes through the frozen coupled propagator."""

from __future__ import annotations

import argparse
import json
import time
from functools import partial
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import LinAlgWarning, get_lapack_funcs, lu_factor

import projects.ukrmol_co.numerov_typed_coupled as kernel
from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.native_full_reconstruction import digest, load_model, save, verify_inputs
from projects.ukrmol_co.numerov_dense_scaling import solve_roles
from projects.ukrmol_co.numerov_native_coupled import (
    REFERENCE_SOLVE,
    ResidualAudit,
    installed_solver,
    native_refined_solve,
)
from projects.ukrmol_co.numerov_native_residual_control import canonical_hash, exact, reference


class _PrefixComplete(Exception):
    """Deliberate stop after the declared recurrence, before endpoint matching."""


def condition_estimate(a):
    """Binary64 LU one-norm reciprocal-condition estimate, not an error bound."""
    a64 = np.asarray(a, float)
    factors, _ = lu_factor(a64)
    rcond, info = get_lapack_funcs("gecon", (a64,))(factors, np.linalg.norm(a64, 1))
    if info or not np.isfinite(rcond) or rcond < 0:
        raise ValueError("Invalid LAPACK condition estimate")
    return float(rcond)


def prefix(args, energy, n, h, backend, *, audited=False, captures=None, dtype=np.longdouble):
    """Execute original initialization and N-1 recurrences, restoring all hooks."""
    roles = solve_roles("varying", n)[:-1]
    audit = ResidualAudit() if audited and backend == "native" else None
    solver = (
        partial(native_refined_solve, observer=audit) if backend == "native" else REFERENCE_SOLVE
    )
    calls, samples = [], []
    stats, last = None, None
    original_taylor = kernel.taylor_change
    taylor_calls = 0

    def measured(a, b, diagnostics):
        nonlocal stats, last
        index = len(calls)
        if index >= len(roles):
            raise ValueError("Undeclared prefix solve")
        started = time.perf_counter()
        last = solver(a, b, diagnostics)
        calls.append(dict(role=roles[index], wall_seconds=time.perf_counter() - started))
        stats = diagnostics
        if audited:
            samples.append(
                dict(
                    index=index,
                    role=roles[index],
                    reciprocal_condition_1=condition_estimate(a),
                    a_sha256=canonical_hash(a),
                    b_sha256=canonical_hash(b),
                    x_sha256=canonical_hash(last),
                )
            )
            if captures is not None and index in [0, len(roles) - 1]:
                np.savez(captures / f"solve-{index:03d}.npz", a=a, b=b, x=last)
        return last

    def bounded_taylor(*values):
        nonlocal taylor_calls
        taylor_calls += 1
        if taylor_calls == 2:
            if len(calls) != len(roles):
                raise ValueError("Prefix ended before the declared recurrence")
            raise _PrefixComplete
        return original_taylor(*values)

    started = time.perf_counter()
    try:
        with installed_solver(measured):
            kernel.taylor_change = bounded_taylor
            try:
                kernel.typed_coupled_k(*args, energy, 18.0, 18.0 + n * h, h, 8, dtype=dtype)
            except _PrefixComplete:
                pass
            else:
                raise ValueError("Prefix did not stop before matching")
    finally:
        kernel.taylor_change = original_taylor
    elapsed = time.perf_counter() - started
    return dict(
        intervals=n,
        backend=backend,
        audited=audited,
        final_delta_sha256=canonical_hash(last),
        final_delta_maximum=float(np.max(abs(last))),
        diagnostics=stats,
        kernel_wall_seconds=elapsed,
        calls=calls,
        samples=samples,
        residual_audit=audit.summary() if audit else None,
        recurring_seconds=sum(x["wall_seconds"] for x in calls[1:]),
        matching_evaluations=0,
        last_recurrence_radius_bohr=18 + (n - 1) * h,
    )


def oracle_checks(directory):
    """Independent high-precision residual entries of actual captured accepted X."""
    rows = []
    for path in sorted(directory.glob("solve-*.npz")):
        with np.load(path) as data:
            a, b, x = [data[k] for k in ["a", "b", "x"]]
        positions = [(0, 0), (len(a) // 2, len(b[0]) // 2), (len(a) - 1, len(b[0]) - 1)]
        low, high = [reference(a, b, x, positions, digits) for digits in [80, 110]]
        residual = b - a @ x
        with mp.workdps(110):
            norm = lambda value: max(mp.fsum(abs(exact(v)) for v in row) for row in value)  # noqa: E731
            denominator = norm(a) * norm(x) + norm(b)
            precision = (
                max(abs(left - right) for left, right in zip(low, high, strict=True)) / denominator
            )
            error = (
                max(
                    abs(exact(residual[i, j]) - want)
                    for (i, j), want in zip(positions, high, strict=True)
                )
                / denominator
            )
        rows.append(
            dict(
                file=path.name,
                sha256=digest(path),
                positions=positions,
                references={
                    str(digits): [mp.nstr(x, digits) for x in values]
                    for digits, values in [(80, low), (110, high)]
                },
                reference_scaled_stability=float(precision),
                residual_scaled_error=float(error),
                checks=dict(
                    precision=precision <= mp.mpf("1e-60"), residual=error <= 32 * mp.mpf(2) ** -63
                ),
            )
        )
    return rows


def decisions(row, baseline=None):
    checks = dict(
        solve_count=row["diagnostics"]["large_system_solves"] == 2 * row["intervals"] - 1,
        backward_error=row["diagnostics"]["maximum_solve_backward_error"] <= 32 * 2.0**-63,
        corrections=row["diagnostics"]["maximum_refinement_corrections"] <= 3,
        role_sequence=[x["role"] for x in row["calls"]]
        == solve_roles("varying", row["intervals"])[:-1],
        stopped_before_matching=row["matching_evaluations"] == 0,
    )
    if baseline is not None:
        checks.update(
            final_delta_exact=row["final_delta_sha256"] == baseline["final_delta_sha256"],
            diagnostics_exact=row["diagnostics"] == baseline["diagnostics"],
            sampled_operands_exact=row["samples"] == baseline["samples"],
        )
    if row["audited"]:
        checks["condition_estimates"] = len(row["samples"]) == 2 * row["intervals"] - 1 and all(
            x["reciprocal_condition_1"] > 0 for x in row["samples"]
        )
    if row["residual_audit"] is not None:
        a = row["residual_audit"]
        checks.update(
            residuals_exact=a["evaluations"] == a["exact_residuals"],
            gate_decisions=a["gate_disagreements"] == 0,
            accepted_solves=a["accepted_solves"] == 2 * row["intervals"] - 1,
            independent_backward_error=a["maximum_accepted_numpy_error"] <= 32 * 2.0**-63,
            precision=a["precision_preserved"],
        )
    return {key: bool(value) for key, value in checks.items()}


def timing_summary(rows):
    result = []
    for case in dict.fromkeys(x["case"] for x in rows):
        means = {
            (b, n): np.mean(
                [
                    x["kernel_wall_seconds"]
                    for x in rows
                    if (x["case"], x["backend"], x["intervals"]) == (case, b, n)
                ]
            )
            for b in ["reference", "native"]
            for n in [32, 128]
        }
        rates = {b: float((means[b, 128] - means[b, 32]) / 96) for b in ["reference", "native"]}
        result.append(
            dict(
                case=case,
                marginal_seconds_per_interval=rates,
                marginal_speedup=rates["reference"] / rates["native"]
                if rates["native"] > 0
                else None,
                full_prefix_n128_speedup=float(means["reference", 128] / means["native", 128]),
                native_recurring_seconds_per_interval=float(
                    np.mean(
                        [
                            x["recurring_seconds"] / 127
                            for x in rows
                            if (x["case"], x["backend"], x["intervals"]) == (case, "native", 128)
                        ]
                    )
                ),
            )
        )
    return result


def study(c, raw, output):
    output.mkdir()
    raw_contract_path = Path(__file__).with_name("numerov-corner-contract.json")
    if digest(raw_contract_path) != c["raw_contract_sha256"]:
        raise ValueError("Raw input contract changed")
    raw_contract = json.loads(raw_contract_path.read_text())
    report = dict(
        contract=c,
        experiment_verdict=None,
        rows=[],
        checks={},
        oracles=[],
        timing_analysis=[],
        original_interior_analytic_verdict=False,
        original_interior_scientific_exits=[1, 1],
        historical_unique_energy_coverage=132,
        uncomputed_identities=462,
        new_qualified_co_identities=0,
        co_successor_submitted=False,
        native_full_k_reconstruction_verdict=False,
        native_ukrmol_invocations=0,
        new_inner_solves=0,
    )

    def checkpoint():
        save(output / "report.json", report)

    checkpoint()
    try:
        import qscat_extended_kernels as native

        if (
            tuple(native.capability()) != (16, 64, 16384, True)
            or np.finfo(np.longdouble).nmant + 1 != 64
        ):
            raise ValueError("Unqualified scalar capability")
        verify_inputs(raw, raw_contract)
        for stage in ["accuracy", "timing"]:
            for repetition in range(1 if stage == "accuracy" else 2):
                items = (
                    c["identities"]
                    if stage == "accuracy"
                    else [x for x in c["identities"] if x["timing"]]
                )
                for item in items if repetition == 0 else list(reversed(items)):
                    model = load_model(raw / item["case"])
                    energy = model["energies"][item["energy_index"]]
                    if float(energy).hex() != item["energy_rydberg_hex"]:
                        raise ValueError("Archived energy identity differs")
                    q = surface_matrix(
                        model["poles"], model["amplitudes"], model["ground"] + energy / 2
                    )
                    args = [q] + [model[k] for k in ["thresholds", "angular", "c2", "c3"]]
                    for n in (
                        [8]
                        if stage == "accuracy"
                        else ([32, 128] if repetition == 0 else [128, 32])
                    ):
                        baseline = None
                        for backend in (
                            ["reference", "native"] if repetition == 0 else ["native", "reference"]
                        ):
                            captures = output / "operands" / item["case"]
                            if stage == "accuracy" and backend == "native":
                                captures.mkdir(parents=True)
                            row = prefix(
                                args,
                                energy,
                                n,
                                c["step_bohr"],
                                backend,
                                audited=stage == "accuracy",
                                captures=captures
                                if stage == "accuracy" and backend == "native"
                                else None,
                            )
                            row.update(
                                case=item["case"],
                                energy_index=item["energy_index"],
                                energy_rydberg_hex=float(energy).hex(),
                                stage=stage,
                                repetition=repetition,
                            )
                            row["checks"] = decisions(row, baseline)
                            report["checks"].update(
                                {f"{len(report['rows'])}_{k}": v for k, v in row["checks"].items()}
                            )
                            report["rows"].append(row)
                            checkpoint()
                            baseline = row
                        if not all(v for r in report["rows"][-2:] for v in r["checks"].values()):
                            report["experiment_verdict"] = False
                            return report
            if stage == "accuracy":
                for item in c["identities"]:
                    for oracle in oracle_checks(output / "operands" / item["case"]):
                        oracle["case"] = item["case"]
                        report["checks"].update(
                            {
                                f"oracle_{len(report['oracles'])}_{k}": bool(v)
                                for k, v in oracle["checks"].items()
                            }
                        )
                        report["oracles"].append(oracle)
                if not all(report["checks"].values()):
                    report["experiment_verdict"] = False
                    return report
        report["timing_analysis"] = timing_summary(
            [x for x in report["rows"] if x["stage"] == "timing"]
        )
        report["checks"].update(
            row_count=len(report["rows"]) == c["prefix_executions"],
            solve_count=sum(x["diagnostics"]["large_system_solves"] for x in report["rows"])
            == c["mixed_solves"],
            oracle_count=len(report["oracles"]) == 12,
        )
        report["experiment_verdict"] = all(report["checks"].values())
    except (
        OSError,
        ValueError,
        RuntimeError,
        AssertionError,
        ImportError,
        LinAlgWarning,
        np.linalg.LinAlgError,
    ) as error:
        report["failure"] = f"{type(error).__name__}: {error}"
    finally:
        verify_inputs(raw, raw_contract)
        checkpoint()
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["contract", "raw", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    a = parser.parse_args()
    if study(json.loads(a.contract.read_text()), a.raw, a.output)["experiment_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
