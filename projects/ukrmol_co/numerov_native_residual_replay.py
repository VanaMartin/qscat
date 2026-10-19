"""Independent saved-metric decision audit and exact primitive numerical replay."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


def saved_decisions(report):
    c = report["contract"]
    if (
        report["native_capability"] != [16, 64, 16384, True]
        or report["numpy_significand_bits"] != 64
    ):
        raise ValueError("Precision capability differs")
    checks = {}
    references = 0
    for row in report["rows"]:
        native = "native" in row["output_numeric_hashes"]
        bits = 64 if native else 53
        eps = 2.0 ** (1 - bits)
        if any(not math.isfinite(v) for v in row["oracle_scaled_errors"].values()):
            raise ValueError("Nonfinite saved oracle error")
        local = dict(
            inputs_unchanged=row["checks"]["inputs_unchanged"],
            shape_dtype_owned=row["checks"]["shape_dtype_owned"],
            independent_oracle_pair=row["oracle_pair_scaled_error"]
            <= c["oracle_pair_scaled_error_maximum"],
            oracle_accuracy=all(
                v <= c["scaled_error_eps_factor"] * eps
                for v in row["oracle_scaled_errors"].values()
            ),
        )
        if native:
            local["full_matrix_differential"] = (
                row["full_matrix_differential_scaled_error"] <= c["scaled_error_eps_factor"] * eps
            )
        if row["name"].startswith("operands-"):
            local["accepted_X_backward_gate"] = all(
                v <= 32 * eps for v in row["backward_errors"].values()
            )
            expected = [
                [i, j] for i in c["capture_reference_rows"] for j in c["capture_reference_columns"]
            ]
            if row["reference_positions"] != expected or row["shapes"] != [[298, 298]] * 3:
                raise ValueError("Captured reference selection or shape differs")
        if row["name"] == "low_bit_signal":
            # Scalar output identity is checked by the complete numerical replay.
            local["signal_not_downcast"] = row["checks"]["signal_not_downcast"]
        if local != row["checks"]:
            raise ValueError("Saved row decisions differ from recorded metrics")
        references += len(row["reference_positions"])
        checks.update({row["name"] + "_" + key: value for key, value in local.items()})
    if checks != report["checks"] or all(checks.values()) is not report["experiment_verdict"]:
        raise ValueError("Overall numerical decisions differ")
    timing = []
    for row in report["timings"]:
        actual_order = [(item["repetition"], item["backend"]) for item in row["rows"]]
        expected_order = [
            (i, backend) for i, order in enumerate(c["timing"]["order"]) for backend in order
        ]
        if actual_order != expected_order:
            raise ValueError("Timing order differs")
        means = {
            backend: sum(item["wall_seconds"] for item in row["rows"] if item["backend"] == backend)
            / 2
            for backend in ("numpy", "native")
        }
        speedup = means["numpy"] / means["native"]
        hypothesis = speedup >= c["timing"]["minimum_speedup_each_extended_capture"]
        numerical = next(item for item in report["rows"] if item["name"] == row["name"])
        if any(
            item["output_numeric_hash"] != numerical["output_numeric_hashes"][item["backend"]]
            for item in row["rows"]
        ):
            raise ValueError("Timing result differs from qualified output")
        if (
            means != row["mean_seconds"]
            or speedup != row["speedup"]
            or hypothesis is not row["performance_hypothesis_passed"]
        ):
            raise ValueError("Saved timing calculations differ")
        timing.append(hypothesis)
    if (
        all(timing) is not report["all_timing_hypotheses_passed"]
        or len(timing) != c["native_capture_cases"]
    ):
        raise ValueError("Timing hypotheses differ")
    if len(report["rows"]) != c["captured_triples"] + len(c["small_cases"]):
        raise ValueError("Finite case scope differs")
    return dict(
        saved_decisions_verified=True,
        experiment_verdict=report["experiment_verdict"],
        checks_passed=sum(checks.values()),
        checks_total=len(checks),
        rows_verified=len(report["rows"]),
        oracle_pairs_verified=len(report["rows"]),
        reference_entries_verified=references,
        all_timing_hypotheses_passed=all(timing),
        saved_metadata_observations_retained=[
            "inputs_unchanged",
            "shape_dtype_owned",
            "signal_not_downcast",
        ],
        numeric_observations_independently_reexecuted=False,
    )


def complete_replay(inputs, old, output):
    from projects.ukrmol_co.numerov_native_residual_control import study

    study(inputs, old["contract"], output / "report.json")
    fresh = json.loads((output / "report.json").read_text())
    exact_rows = fresh["rows"] == old["rows"]
    checks = fresh["checks"] == old["checks"]
    passed = (
        exact_rows and checks and fresh["experiment_verdict"] is old["experiment_verdict"] is True
    )
    result = dict(
        public_numerical_reexecution_passed=passed,
        all_rows_exact=exact_rows,
        all_qualification_decisions_exact=checks,
        extension_sha256=fresh["extension_sha256"],
        original_extension_sha256=old["extension_sha256"],
        timings=fresh["timings"],
        all_timing_hypotheses_passed=fresh["all_timing_hypotheses_passed"],
        timing_identity_is_not_a_numerical_gate=True,
        saved_audit=saved_decisions(fresh),
        new_propagations=0,
        new_co_independent_solves=0,
        new_unique_energy_identities=0,
    )
    (output / "replay.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--inputs", type=Path)
    args = parser.parse_args()
    old = json.loads(args.report.read_text())
    if args.inputs:
        args.output.mkdir()
        result = complete_replay(args.inputs, old, args.output)
        if not result["public_numerical_reexecution_passed"]:
            raise SystemExit(1)
    else:
        result = saved_decisions(old)
        result["report_sha256"] = hashlib.sha256(args.report.read_bytes()).hexdigest()
        args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
