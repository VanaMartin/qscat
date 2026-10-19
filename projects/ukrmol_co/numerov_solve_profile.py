"""Source-line timing of unchanged mixed solves and bounded operand capture."""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import LinAlgWarning

import projects.ukrmol_co.numerov_typed_coupled as kernel
from projects.ukrmol_co.numerov_arithmetic_control import decimal_matrix, precision_capability
from projects.ukrmol_co.numerov_dense_scaling import inputs, solve_roles


def line_categories(function):
    lines, start = inspect.getsourcelines(function)
    expressions = {
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
    source = {start + i: line.strip() for i, line in enumerate(lines)}
    if not set(expressions).issubset(source.values()):
        raise ValueError("Mixed-solve source statements differ from the timing contract")
    return source, {n: expressions.get(line, "other") for n, line in source.items()}


def trace_solve(function, a, b, stats):
    """Account parent-frame line intervals; child time is included once."""
    if sys.gettrace() is not None:
        raise ValueError("Existing Python tracer prevents isolated line attribution")
    source, categories = line_categories(function)
    code = function.__code__
    previous, started = None, None
    seconds, counts = {}, {}

    def trace(frame, event, _argument):
        nonlocal previous, started
        if frame.f_code is not code:
            return None
        now = time.perf_counter()
        if event in ("line", "return"):
            if previous is not None:
                seconds[previous] = seconds.get(previous, 0.0) + now - started
            if event == "line":
                previous = frame.f_lineno
                counts[previous] = counts.get(previous, 0) + 1
                started = time.perf_counter()
            else:
                previous = None
        return trace

    begin = time.perf_counter()
    try:
        sys.settrace(trace)
        solution = function(a, b, stats)
    finally:
        sys.settrace(None)
    elapsed = time.perf_counter() - begin
    lines = [
        dict(
            line=n,
            source=source[n],
            category=categories[n],
            calls=counts[n],
            wall_seconds=seconds.get(n, 0.0),
        )
        for n in sorted(counts)
    ]
    totals = {}
    for line in lines:
        totals[line["category"]] = totals.get(line["category"], 0.0) + line["wall_seconds"]
    return solution, dict(wall_seconds=elapsed, lines=lines, category_seconds=totals)


def profiled_row(c, model, backend, traced, output, capture_interval):
    function = kernel.refined_solve
    roles = solve_roles(model, c["profile_intervals"])
    calls, captures, occurrences = [], [], {}
    dtype = np.float64 if backend == "binary64" else np.longdouble

    def measured(a, b, stats):
        index = len(calls)
        if index >= len(roles):
            raise ValueError("Undeclared mixed solve")
        role = roles[index]
        occurrences[role] = occurrences.get(role, 0) + 1
        begin = time.perf_counter()
        if traced:
            solution, timing = trace_solve(function, a, b, stats)
        else:
            solution = function(a, b, stats)
            timing = dict(wall_seconds=time.perf_counter() - begin, lines=[], category_seconds={})
        calls.append(
            timing
            | dict(
                role=role,
                a_shape=list(a.shape),
                b_shape=list(b.shape),
                x_shape=list(solution.shape),
                dtype=a.dtype.name,
                b_dtype=b.dtype.name,
                x_dtype=solution.dtype.name,
            )
        )
        if traced and role.startswith("recurring_") and occurrences[role] == capture_interval:
            captures.append((role, a.copy(), b.copy(), solution.copy()))
        return solution

    kernel.refined_solve = measured
    n = c["profile_intervals"]
    begin = time.perf_counter()
    try:
        k, diagnostics = kernel.typed_coupled_k(
            *inputs(c, model),
            c["energy_rydberg"],
            c["radius_bohr"],
            c["radius_bohr"] + n * c["step_bohr"],
            c["step_bohr"],
            c["tail_orders"][model],
            dtype=dtype,
        )
    finally:
        kernel.refined_solve = function
    wall = time.perf_counter() - begin
    saved = []
    for role, a, b, x in captures:
        name = f"operands-{model}-{backend}-{role}.npz"
        path = output / name
        np.savez(path, a=a, b=b, x=x)
        saved.append(
            dict(
                file=name,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                role=role,
                interval=capture_interval,
                dtype=a.dtype.name,
                significand_bits=np.finfo(a.dtype).nmant + 1,
                a_shape=list(a.shape),
                b_shape=list(b.shape),
                x_shape=list(x.shape),
            )
        )
    return dict(
        model=model,
        backend=backend,
        traced=traced,
        kernel_wall_seconds=wall,
        k_decimal=decimal_matrix(k),
        k=np.asarray(k, float).tolist(),
        diagnostics=diagnostics,
        calls=calls,
        captures=saved,
    )


def summarize(rows, contract):
    summaries = []
    hypotheses = {}
    for row in rows:
        recurring = [x for x in row["calls"] if x["role"].startswith("recurring_")]
        totals = {}
        roles = {}
        for call in recurring:
            per_role = roles.setdefault(call["role"], dict(wall_seconds=0.0, category_seconds={}))
            per_role["wall_seconds"] += call["wall_seconds"]
            for key, seconds in call["category_seconds"].items():
                totals[key] = totals.get(key, 0.0) + seconds
                per_role["category_seconds"][key] = (
                    per_role["category_seconds"].get(key, 0.0) + seconds
                )
        total = sum(x["wall_seconds"] for x in recurring)
        summaries.append(
            dict(
                model=row["model"],
                backend=row["backend"],
                traced=row["traced"],
                recurring_wall_seconds=total,
                category_seconds=totals,
                roles=roles,
            )
        )
    for model in contract["models"]:
        pair = {
            x["traced"]: x for x in summaries if (x["model"], x["backend"]) == (model, "longdouble")
        }
        timed = pair[True]
        values = timed["category_seconds"]
        residual = values["residual_evaluation"]
        lu = sum(
            values.get(key, 0.0)
            for key in ["factorization_and_cast", "initial_solve_and_cast", "correction_and_cast"]
        )
        fraction = residual / timed["recurring_wall_seconds"]
        ratio = timed["recurring_wall_seconds"] / pair[False]["recurring_wall_seconds"]
        timed.update(
            residual_fraction=fraction,
            residual_to_lu_ratio=residual / lu,
            traced_to_untraced_recurring_ratio=ratio,
        )
        h = contract["hypotheses"]
        hypotheses[model] = dict(
            residual_fraction_at_least_80percent=fraction
            >= h["each_extended_model_residual_fraction_minimum"],
            residual_exceeds_lu_at_least_10fold=residual / lu
            >= h["each_extended_model_residual_to_factor_initial_correction_ratio_minimum"],
            trace_overhead_within_declared_range=h[
                "each_extended_model_traced_to_untraced_recurring_ratio_range"
            ][0]
            <= ratio
            <= h["each_extended_model_traced_to_untraced_recurring_ratio_range"][1],
        )
    return dict(rows=summaries, hypotheses=hypotheses)


def study(prior_path, contract, output):
    report = dict(contract=contract, experiment_verdict=None, checks={}, rows=[], summary={})
    begin = time.monotonic()

    def save():
        (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    try:
        if (
            hashlib.sha256(Path(kernel.__file__).read_bytes()).hexdigest()
            != contract["numerical_module_sha256"]
        ):
            raise ValueError("Pinned coupled numerical module differs")
        if hashlib.sha256(prior_path.read_bytes()).hexdigest() != contract["prior_report_sha256"]:
            raise ValueError("Pinned dense scaling report differs")
        prior = json.loads(prior_path.read_text())
        if prior["experiment_verdict"] is not True:
            raise ValueError("Dense scaling prerequisite did not pass")
        report["precision_bits"] = precision_capability()
        c = prior["contract"] | dict(profile_intervals=contract["intervals"])
        coordinates = [(m, b) for m in contract["models"] for b in contract["backends"]]
        for traced in [False, True]:
            for model, backend in coordinates if not traced else coordinates[::-1]:
                row = profiled_row(c, model, backend, traced, output, contract["capture_interval"])
                old = next(
                    x
                    for x in prior["rows"]
                    if (x["model"], x["backend"], x["intervals"], x["repetition"])
                    == (model, backend, contract["intervals"], 0)
                )
                calls = row["calls"]
                dtype = "float64" if backend == "binary64" else np.dtype(np.longdouble).name
                checks = dict(
                    archived_K_exact=row["k_decimal"] == old["k_decimal"] and row["k"] == old["k"],
                    archived_diagnostics_exact=row["diagnostics"] == old["diagnostics"],
                    roles=[x["role"] for x in calls] == solve_roles(model, contract["intervals"]),
                    shapes=all(
                        x["a_shape"] == [c["channels"]] * 2
                        and x["b_shape"]
                        == x["x_shape"]
                        == [
                            c["channels"],
                            c["open_channels"] if x["role"] == "matching" else c["channels"],
                        ]
                        for x in calls
                    ),
                    dtypes=all(x["dtype"] == x["b_dtype"] == x["x_dtype"] == dtype for x in calls),
                    timer_partitions=all(
                        0
                        <= sum(x["category_seconds"].values())
                        <= x["wall_seconds"] + contract["timer_partition_slack_seconds"]
                        for x in calls
                    ),
                    backward_error=row["diagnostics"]["maximum_solve_backward_error"]
                    <= 32 * 2.0 ** (1 - report["precision_bits"][backend]),
                    bounded_corrections=row["diagnostics"]["maximum_refinement_corrections"]
                    <= contract["maximum_refinement_corrections"],
                    captures=len(row["captures"])
                    == ((1 if model == "free" else 2) if traced else 0),
                )
                if traced:
                    checks["residual_counts"] = all(
                        sum(
                            x["calls"]
                            for x in call["lines"]
                            if x["category"] == "residual_evaluation"
                        )
                        in range(1, 5)
                        for call in calls
                    )
                report["checks"].update(
                    {
                        f"{traced}_{model}_{backend}_{key}": bool(value)
                        for key, value in checks.items()
                    }
                )
                report["rows"].append(row)
                save()
                if not all(checks.values()):
                    report["experiment_verdict"] = False
                    return report
        report["summary"] = summarize(report["rows"], contract)
        report["experiment_verdict"] = all(report["checks"].values())
    except (
        OSError,
        ValueError,
        RuntimeError,
        AssertionError,
        np.linalg.LinAlgError,
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
