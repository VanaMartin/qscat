"""Independent saved numerical decisions, call partitions and timing-fit audit."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from projects.ukrmol_co.native_reconstruction_replay import scattering
from projects.ukrmol_co.numerov_coupled_oracle import direct_discrete_k
from projects.ukrmol_co.numerov_discrete_free import discrete_free_k


def reference(c, model, intervals, digits):
    e, r, h = c["energy_rydberg"], c["radius_bohr"], c["step_bohr"]
    endpoint = r + h * intervals
    with mp.workdps(digits):
        thresholds = [0.0, 0.0, c["closed_threshold_rydberg"]]
        if model == "free":
            active = discrete_free_k(c["active_q"], thresholds, e, r, endpoint, h, digits)
        else:
            active = direct_discrete_k(
                c["active_q"],
                thresholds,
                [0, 0, 0],
                c["varying_active_c2"],
                c["varying_active_c3"],
                e,
                r,
                endpoint,
                h,
                c["tail_orders"][model],
                digits,
            )
        full = mp.zeros(6)
        for i in range(2):
            for j in range(2):
                full[i, j] = active[i, j]
        for i, q in enumerate(c["other_open_q"], 2):
            full[i, i] = discrete_free_k([[q]], [0.0], e, r, endpoint, h, digits)[0, 0]
        return +full


def difference(a, b):
    with mp.workdps(100):
        return float(max(abs(x - y) for x, y in zip(a, b, strict=True)))


def fits_and_hypotheses(rows, c):
    fits = []
    for model in c["models"]:
        for backend in c["backends"]:
            counts = np.asarray(c["interval_counts"], float)
            means = np.array(
                [
                    np.mean(
                        [
                            r["kernel_wall_seconds"]
                            for r in rows
                            if (r["model"], r["backend"], r["intervals"])
                            == (model, backend, int(n))
                        ]
                    )
                    for n in counts
                ]
            )
            intercept, slope = np.linalg.lstsq(
                np.column_stack([np.ones(3), counts]), means, rcond=None
            )[0]
            residual = float(np.max(abs(intercept + slope * counts - means) / means))
            marginal = float((means[2] - means[1]) / (counts[2] - counts[1]))
            fits.append(
                dict(
                    model=model,
                    backend=backend,
                    intercept=float(intercept),
                    slope=float(slope),
                    relative_residual=residual,
                    marginal=marginal,
                )
            )
    extended = {r["model"]: r for r in fits if r["backend"] == "longdouble"}
    ratio = extended["varying"]["marginal"] / extended["free"]["marginal"]
    free = [
        r for r in rows if (r["model"], r["backend"], r["intervals"]) == ("free", "longdouble", 128)
    ]
    recurring = np.mean(
        [
            sum(
                call["wall_seconds"] for call in r["calls"] if call["role"].startswith("recurring_")
            )
            for r in free
        ]
    )
    fixed = np.mean(
        [
            sum(
                call["wall_seconds"]
                for call in r["calls"]
                if call["function"] in ["taylor_change", "typed_tail"]
            )
            for r in free
        ]
    )
    hypotheses = dict(
        all_affine_fits_positive_and_within_20percent=all(
            f["slope"] > 0 and f["relative_residual"] <= 0.2 for f in fits
        ),
        free_extended_N128_recurring_solves_exceed_Taylor_plus_tail=bool(recurring > fixed),
        varying_to_free_extended_marginal_rate_between_1_and_3=1 <= ratio <= 3,
    )
    return fits, hypotheses


def verify(root):
    r = json.loads((root / "collected/result/report.json").read_text())
    c = r["contract"]
    e = json.loads((root / "collected/execution.json").read_text())
    if (
        hashlib.sha256((root / "prior-report.json").read_bytes()).hexdigest()
        != c["prior_report_sha256"]
    ):
        raise ValueError("Prior coupled report hash differs")
    for n, w in json.loads((root / "collected/collection-manifest.json").read_text()).items():
        if hashlib.sha256((root / "collected" / n).read_bytes()).hexdigest() != w:
            raise ValueError("Collected payload hash differs")
    prior = json.loads((root / "prior-report.json").read_text())
    checks = {}
    refs = {}
    for model in c["models"]:
        for n in c["interval_counts"]:
            low, high = [reference(c, model, n, d) for d in c["oracle_digits"]]
            checks[f"oracle_{model}_{n}"] = difference(low, high) < c["oracle_precision_budget"]
            saved = next(x for x in r["oracles"] if (x["model"], x["intervals"]) == (model, n))
            with mp.workdps(100):
                saved_k = mp.matrix(saved["k_decimal"])
            if difference(saved_k, high) > 1e-60:
                raise ValueError("Saved discrete oracle differs")
            refs[model, n] = high
    seen = {}
    roles_audited = 0
    for row in r["rows"]:
        rep, model, n, backend = (row[k] for k in ["repetition", "model", "intervals", "backend"])
        identity = (rep, model, n, backend)
        if (
            identity in seen
            or rep not in [0, 1]
            or model not in c["models"]
            or n not in c["interval_counts"]
            or backend not in c["backends"]
        ):
            raise ValueError("Duplicate or undeclared benchmark row")
        seen[identity] = row
        with mp.workdps(100):
            measured = mp.matrix(row["k_decimal"])
        array = np.array(row["k_decimal"], float)
        if not np.array_equal(array, np.asarray(row["k"])):
            raise ValueError("Decimal/binary64 presentation differs")
        error = difference(measured, refs[model, n])
        if abs(error - row["discrete_arithmetic_error"]) > 1e-20:
            raise ValueError("Discrete arithmetic error differs")
        s = scattering(array)
        invariant = max(
            float(np.max(abs(s.conj().T @ s - np.eye(len(s))))), float(np.max(abs(s - s.T)))
        )
        roles = (
            ["initialization", "cached_potential"] + ["recurring_defect"] * (n - 1) + ["matching"]
            if model == "free"
            else ["initialization"]
            + ["recurring_potential", "recurring_defect"] * (n - 1)
            + ["matching"]
        )
        calls = row["calls"]
        actual_roles = [x["role"] for x in calls if x["function"] == "refined_solve"]
        if actual_roles != roles:
            raise ValueError("Frozen solve-role order differs")
        roles_audited += len(roles)
        if any(x["wall_seconds"] < 0 or not np.isfinite(x["wall_seconds"]) for x in calls):
            raise ValueError("Invalid elapsed time")
        totals = {}
        for call in calls:
            totals[call["role"]] = totals.get(call["role"], 0.0) + call["wall_seconds"]
        if (
            totals != row["timing_by_role"]
            or abs(sum(totals.values()) - row["timed_function_seconds"]) > 1e-12
        ):
            raise ValueError("Disjoint function-time partition differs")
        dtype = row["diagnostics"]["dtype"]
        if any(x["result_dtype"] != dtype for x in calls):
            raise ValueError("A wrapped function lost dtype")
        prefix = f"{rep}_{model}_{n}_{backend}"
        local = dict(
            arithmetic=error <= c["arithmetic_K_budgets"][backend],
            invariants=invariant <= c["invariant_budget"],
            solve_count=row["diagnostics"]["large_system_solves"] == len(roles)
            and len(actual_roles) == len(roles),
            fixed_call_counts=sum(x["function"] == "taylor_change" for x in calls) == 3
            and sum(x["function"] == "typed_tail" for x in calls) == 1,
            timer_partition=row["timed_function_seconds"]
            <= row["kernel_wall_seconds"] + c["disjoint_timer_sum_slack_seconds"],
            backward_error=row["diagnostics"]["maximum_solve_backward_error"]
            <= 32 * 2.0 ** (1 - r["precision_bits"][backend]),
        )
        if model == "free" and n == 8:
            old = next(
                x for x in prior["rows"] if x["stage"] == "dense" and x["backend"] == backend
            )
            local["old_baseline_exact"] = row["k_decimal"] == old["k_decimal"]
        if rep:
            local["repetition_exact"] = row["k_decimal"] == seen[0, model, n, backend]["k_decimal"]
        checks.update({prefix + "_" + key: bool(value) for key, value in local.items()})
    if checks != r["checks"]:
        raise ValueError("Independently reconstructed benchmark decisions differ")
    expected = {
        (rep, model, n, backend)
        for rep in range(2)
        for model in c["models"]
        for n in c["interval_counts"]
        for backend in c["backends"]
    }
    verdict = all(checks.values()) and set(seen) == expected
    if (
        verdict is not r["experiment_verdict"]
        or e["original_command_exit"] != (0 if verdict else 1)
        or e["original_container_exit"] != e["original_command_exit"]
    ):
        raise ValueError("Aggregate benchmark verdict or original exits differ")
    fits, hypotheses = fits_and_hypotheses(r["rows"], c)
    if hypotheses != r["timing_analysis"]["hypotheses"]:
        raise ValueError("Independent timing hypotheses differ")
    for fit in fits:
        stored = next(
            x
            for x in r["timing_analysis"]["fits"]
            if (x["model"], x["backend"]) == (fit["model"], fit["backend"])
        )
        for key, saved_key in [
            ("intercept", "affine_intercept_seconds"),
            ("slope", "affine_slope_seconds_per_interval"),
            ("relative_residual", "maximum_relative_affine_residual"),
            ("marginal", "marginal_32_to_128_seconds_per_interval"),
        ]:
            if abs(fit[key] - stored[saved_key]) > 1e-10:
                raise ValueError("Independent affine fit differs")
    return dict(
        saved_decision_replay_passed=True,
        experiment_verdict=verdict,
        checks_passed=sum(checks.values()),
        checks_total=len(checks),
        propagation_rows_verified=len(r["rows"]),
        solve_call_roles_audited=roles_audited,
        timer_partitions_audited=len(r["rows"]),
        output_dtype_audits_passed=True,
        affine_fits_verified=4,
        timing_hypotheses=hypotheses,
        original_exits=[e["original_command_exit"], e["original_container_exit"]],
        shape_provenance=(
            "Original per-call logs store role, elapsed time and output dtype, not array shapes; "
            "shapes follow fixed numerical source and are not independently measured here."
        ),
        original_interior_analytic_verdict=False,
        original_interior_scientific_exits=[1, 1],
        co_successor_submitted=False,
        new_co_independent_solves=0,
        new_unique_energy_identities=0,
        candidate_qualified_identities=14,
        unique_energy_coverage=132,
        uncomputed_identities=462,
        native_full_k_reconstruction_verdict=False,
        co_complex_pole_verdict=None,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(verify(args.root), indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
