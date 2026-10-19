"""Independent saved-decision and source-line timing audit for mixed solves."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def verify(root):
    report = json.loads((root / "collected/result/report.json").read_text())
    contract = report["contract"]
    prior_path = root / "prior-report.json"
    if hashlib.sha256(prior_path.read_bytes()).hexdigest() != contract["prior_report_sha256"]:
        raise ValueError("Prior dense report differs")
    prior = json.loads(prior_path.read_text())
    execution = json.loads((root / "collected/execution.json").read_text())
    for name, want in json.loads((root / "collected/collection-manifest.json").read_text()).items():
        if hashlib.sha256((root / "collected" / name).read_bytes()).hexdigest() != want:
            raise ValueError("Collected payload differs")
    c = prior["contract"]
    checks = {}
    summaries = []
    identities = []
    capture_count = 0
    solve_count = 0
    statements = {
        "residual = b - a @ x": "residual_evaluation",
        "factors = lu_factor(np.asarray(a, float), check_finite=True)": "factorization_and_cast",
        "x = np.asarray(lu_solve(factors, np.asarray(b, float)), dtype=dtype)": (
            "initial_solve_and_cast"
        ),
        "x += np.asarray(lu_solve(factors, np.asarray(residual, float)), dtype=dtype)": (
            "correction_and_cast"
        ),
        "denominator = norm(a) * norm(x) + norm(b)": "error_norms",
        "error = norm(residual) / denominator if denominator else dtype(0)": "error_norms",
    }
    for row in report["rows"]:
        model, backend, traced = (row[k] for k in ["model", "backend", "traced"])
        identities.append((model, backend, traced))
        old = next(
            x
            for x in prior["rows"]
            if (x["model"], x["backend"], x["intervals"], x["repetition"])
            == (model, backend, contract["intervals"], 0)
        )
        n = contract["intervals"]
        roles = (
            ["initialization", "cached_potential"] + ["recurring_defect"] * (n - 1) + ["matching"]
            if model == "free"
            else ["initialization"]
            + ["recurring_potential", "recurring_defect"] * (n - 1)
            + ["matching"]
        )
        shape, dtype_ok, timer_ok = True, True, True
        residual_counts = []
        recurring = 0.0
        categories = {}
        per_role = {}
        for call in row["calls"]:
            solve_count += 1
            expected_b = [
                c["channels"],
                c["open_channels"] if call["role"] == "matching" else c["channels"],
            ]
            shape = (
                shape
                and call["a_shape"] == [c["channels"]] * 2
                and call["b_shape"] == call["x_shape"] == expected_b
            )
            expected_dtype = old["diagnostics"]["dtype"]
            dtype_ok = (
                dtype_ok and call["dtype"] == call["b_dtype"] == call["x_dtype"] == expected_dtype
            )
            rebuilt = {}
            residual_count = 0
            for line in call["lines"]:
                category = statements.get(line["source"], "other")
                if category != line["category"] or line["wall_seconds"] < 0 or line["calls"] < 1:
                    raise ValueError("Invalid source-line category/timing")
                rebuilt[category] = rebuilt.get(category, 0.0) + line["wall_seconds"]
                if category == "residual_evaluation":
                    residual_count += line["calls"]
            if rebuilt != call["category_seconds"]:
                raise ValueError("Disjoint source-line partition differs")
            residual_counts.append(residual_count)
            timer_ok = (
                timer_ok
                and 0
                <= sum(rebuilt.values())
                <= call["wall_seconds"] + contract["timer_partition_slack_seconds"]
            )
            if call["role"].startswith("recurring_"):
                recurring += call["wall_seconds"]
                value = per_role.setdefault(
                    call["role"], dict(wall_seconds=0.0, category_seconds={})
                )
                value["wall_seconds"] += call["wall_seconds"]
                for key, seconds in rebuilt.items():
                    categories[key] = categories.get(key, 0.0) + seconds
                    value["category_seconds"][key] = (
                        value["category_seconds"].get(key, 0.0) + seconds
                    )
        prefix = f"{traced}_{model}_{backend}_"
        local = dict(
            archived_K_exact=row["k_decimal"] == old["k_decimal"] and row["k"] == old["k"],
            archived_diagnostics_exact=row["diagnostics"] == old["diagnostics"],
            roles=[x["role"] for x in row["calls"]] == roles,
            shapes=shape,
            dtypes=dtype_ok,
            timer_partitions=timer_ok,
            backward_error=row["diagnostics"]["maximum_solve_backward_error"]
            <= 32 * 2.0 ** (1 - report["precision_bits"][backend]),
            bounded_corrections=row["diagnostics"]["maximum_refinement_corrections"]
            <= contract["maximum_refinement_corrections"],
            captures=len(row["captures"]) == ((1 if model == "free" else 2) if traced else 0),
        )
        if traced:
            local["residual_counts"] = all(1 <= x <= 4 for x in residual_counts)
        checks.update({prefix + key: bool(value) for key, value in local.items()})
        for capture in row["captures"]:
            if (
                hashlib.sha256(
                    (root / "collected/result" / capture["file"]).read_bytes()
                ).hexdigest()
                != capture["sha256"]
            ):
                raise ValueError("Captured operand hash differs")
            if (
                capture["a_shape"] != capture["b_shape"]
                or capture["b_shape"] != capture["x_shape"]
                or capture["interval"] != contract["capture_interval"]
            ):
                raise ValueError("Captured operand metadata differs")
            capture_count += 1
        summaries.append(
            dict(
                model=model,
                backend=backend,
                traced=traced,
                recurring_wall_seconds=recurring,
                category_seconds=categories,
                roles=per_role,
            )
        )
    expected = [(m, b, False) for m in contract["models"] for b in contract["backends"]]
    expected += [(m, b, True) for m, b, _ in expected[::-1]]
    if identities != expected or checks != report["checks"]:
        raise ValueError("Frozen row order or qualification decisions differ")
    hypotheses = {}
    h = contract["hypotheses"]
    for model in contract["models"]:
        pair = {
            x["traced"]: x for x in summaries if (x["model"], x["backend"]) == (model, "longdouble")
        }
        measured = pair[True]
        values = measured["category_seconds"]
        fraction = values["residual_evaluation"] / measured["recurring_wall_seconds"]
        lu = sum(
            values.get(k, 0.0)
            for k in ["factorization_and_cast", "initial_solve_and_cast", "correction_and_cast"]
        )
        ratio = measured["recurring_wall_seconds"] / pair[False]["recurring_wall_seconds"]
        measured.update(
            residual_fraction=fraction,
            residual_to_lu_ratio=values["residual_evaluation"] / lu,
            traced_to_untraced_recurring_ratio=ratio,
        )
        hypotheses[model] = dict(
            residual_fraction_at_least_80percent=fraction
            >= h["each_extended_model_residual_fraction_minimum"],
            residual_exceeds_lu_at_least_10fold=values["residual_evaluation"] / lu
            >= h["each_extended_model_residual_to_factor_initial_correction_ratio_minimum"],
            trace_overhead_within_declared_range=h[
                "each_extended_model_traced_to_untraced_recurring_ratio_range"
            ][0]
            <= ratio
            <= h["each_extended_model_traced_to_untraced_recurring_ratio_range"][1],
        )
    if dict(rows=summaries, hypotheses=hypotheses) != report["summary"]:
        raise ValueError("Independent source-line timing summary differs")
    verdict = all(checks.values())
    if (
        verdict is not report["experiment_verdict"]
        or [execution["original_command_exit"], execution["original_container_exit"]]
        != [0 if verdict else 1] * 2
    ):
        raise ValueError("Original verdict/exits differ")
    return dict(
        saved_decision_replay_passed=True,
        experiment_verdict=verdict,
        checks_passed=sum(checks.values()),
        checks_total=len(checks),
        propagation_rows_verified=len(report["rows"]),
        solve_calls_audited=solve_count,
        capture_hashes_verified=capture_count,
        capture_numerical_values_reexecuted=False,
        capture_scope=(
            "Saved local verifier audits hashes and metadata. Wider-host numerical replay "
            "compares loaded operand values and accepted solve backward errors."
        ),
        timing_hypotheses=hypotheses,
        source_statement_partition_verified=True,
        original_exits=[0, 0],
        original_interior_analytic_verdict=False,
        original_interior_scientific_exits=[1, 1],
        co_successor_submitted=False,
        new_co_independent_solves=0,
        new_unique_energy_identities=0,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(verify(args.root), indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
