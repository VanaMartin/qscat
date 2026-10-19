"""Reconstruct integrated native-residual decisions and optionally rerun the scope."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from projects.ukrmol_co.coupled_tail import exact_tail
from projects.ukrmol_co.native_exterior_attribution import algebra, differences
from projects.ukrmol_co.native_exterior_study import cayley
from projects.ukrmol_co.numerov_dense_scaling import solve_roles
from projects.ukrmol_co.numerov_native_integration import digest, study, timing_summary


def maximum_decimal_difference(a, b):
    with mp.workdps(80):
        aa, bb = mp.matrix(a), mp.matrix(b)
        if (aa.rows, aa.cols) != (bb.rows, bb.cols):
            raise ValueError("Reference shape differs")
        return float(max(abs(x - y) for x, y in zip(aa, bb, strict=True)))


def expected_coordinates(coupled, dense):
    cfg = coupled["contract"]
    coordinates = []
    for model in cfg["short"]["models"]:
        for h in cfg["short"]["steps_bohr"]:
            coordinates += [("short", model["name"], h, b, 0) for b in ["reference", "native"]]
    for h in cfg["physical"]["steps_bohr"]:
        coordinates += [("physical", "inverse_square", h, b, 0) for b in ["reference", "native"]]
    for model in ["free", "varying"]:
        coordinates += [("dense_accuracy", model, 8, b, 0) for b in ["reference", "native"]]
    for endpoint in cfg["stress"]["endpoints_bohr"]:
        coordinates += [("stress", str(endpoint), endpoint, b, 0) for b in ["reference", "native"]]
    for repetition in range(2):
        for model in ["free", "varying"] if repetition == 0 else ["varying", "free"]:
            for n in [32, 128] if repetition == 0 else [128, 32]:
                for b in ["reference", "native"] if repetition == 0 else ["native", "reference"]:
                    coordinates.append(("timing", model, n, b, repetition))
    return coordinates


def audit(report, prior_directory):
    c = report["contract"]
    priors = {}
    for name, want in c["prior_report_sha256"].items():
        path = prior_directory / (name + ".json")
        if digest(path) != want:
            raise ValueError("Prior checksum differs: " + name)
        priors[name] = json.loads(path.read_text())
        if priors[name]["experiment_verdict"] is not True:
            raise ValueError("Rejected prerequisite")
    coupled, dense = priors["coupled"], priors["dense"]
    coordinates = expected_coordinates(coupled, dense)
    reconstructed = {}
    previous = None
    for index, row in enumerate(report["rows"]):
        stage, model, parameter, backend, _repetition = coordinates[index]
        if stage in ["short", "physical"]:
            actual_parameter = row["coordinate"][3]
        elif stage == "stress":
            actual_parameter = row["coordinate"][2]
        else:
            actual_parameter = row["intervals"]
        if (
            row["stage"],
            row["model"],
            actual_parameter,
            row["backend"],
            row.get("repetition", 0),
        ) != coordinates[index]:
            raise ValueError("Undeclared row ordering or coordinate")
        k = np.array(row["k_decimal"], float)
        if (
            not np.isfinite(k).all()
            or not np.array_equal(k, row["k"])
            or algebra(k) != row["invariants"]
        ):
            raise ValueError("Saved K or invariant arithmetic differs")
        if stage in ["short", "physical", "stress"]:
            family = coupled["contract"][stage]
            h = parameter if stage != "stress" else family["step_bohr"]
            endpoint = parameter if stage == "stress" else family["endpoint_bohr"]
            energy = (
                float.fromhex(family["energy_rydberg_hex"])
                if stage == "stress"
                else family["energy_rydberg"]
            )
            coordinate = [energy, family["radius_bohr"], endpoint, h, family["tail_order"]]
            if coordinate != row["coordinate"]:
                raise ValueError("Changed physical coordinate")
            n = round((endpoint - family["radius_bohr"]) / h)
            if row["intervals"] != n:
                raise ValueError("Interval count differs")
            roles = solve_roles("free" if stage == "stress" else "varying", n)
        else:
            roles = solve_roles(model, row["intervals"])
        d = row["diagnostics"]
        checks = dict(
            mixed_solve_count=d["large_system_solves"] == len(roles),
            backward_error=d["maximum_solve_backward_error"] <= 32 * 2.0**-63,
            corrections=d["maximum_refinement_corrections"] <= 3,
            precision=d["significand_bits"] == d["complex_significand_bits"] == 64,
        )
        if index % 2:
            checks.update(
                paired_k_exact=row["k_decimal"] == previous["k_decimal"],
                paired_diagnostics_exact=d == previous["diagnostics"],
                paired_invariants_exact=row["invariants"] == previous["invariants"],
            )
        observed = row["residual_audit"]
        audited = backend == "native" and stage != "timing"
        if audited != (observed is not None):
            raise ValueError("Residual audit/timing separation differs")
        if audited:
            checks.update(
                residuals_exact=observed["evaluations"] == observed["exact_residuals"],
                residual_decisions=observed["gate_disagreements"] == 0,
                accepted_numpy_error=observed["maximum_accepted_numpy_error"] <= 32 * 2.0**-63,
                audited_solves=observed["accepted_solves"] == len(roles),
                audit_precision=observed["precision_preserved"],
            )
            if (
                observed["evaluations"] < observed["accepted_solves"]
                or observed["maximum_iteration"] > 3
            ):
                raise ValueError("Impossible residual observation counts")
            if (
                observed["evaluations"] == observed["exact_residuals"]
                and observed["maximum_scaled_difference"] != 0
            ):
                raise ValueError("Exact residual count contradicts nonzero differential")
            channels = 298 if stage == "dense_accuracy" else 4
            if observed["last"]["a_shape"] != [channels, channels] or observed["last"][
                "b_shape"
            ] != [channels, len(k)]:
                raise ValueError("Last audited matching layout differs")
        if stage in ["short", "physical", "stress"]:
            candidates = [
                r
                for r in coupled["rows"]
                if r.get("stage") == stage and r["backend"] == "longdouble"
            ]
            if stage == "short":
                old = next(
                    r for r in candidates if r["model"] == model and r["step_bohr"] == parameter
                )
            elif stage == "physical":
                old = next(r for r in candidates if r["step_bohr"] == parameter)
            else:
                old = next(r for r in candidates if r["endpoint_bohr"] == parameter)
            drift = maximum_decimal_difference(row["k_decimal"], old["k_decimal"])
            if drift != row["archived_k_difference"]:
                raise ValueError("Archived difference arithmetic differs")
            checks["archived_k"] = drift <= c["budgets"][stage + "_archived_K"]
            if stage != "physical":
                error = maximum_decimal_difference(row["k_decimal"], old["discrete_k_decimal"])
                if error != row["discrete_arithmetic_error"]:
                    raise ValueError("Discrete reference error differs")
                checks["discrete_arithmetic"] = error <= c["budgets"][stage + "_discrete_K"]
                ref = (
                    np.array(old["continuous_k_decimal"], float)
                    if stage == "stress"
                    else np.asarray(
                        next(r for r in coupled["references"] if r["model"] == model)["k"]
                    )
                )
                metrics = differences(k, ref)
                if metrics != row["continuous_errors"]:
                    raise ValueError("Continuous error arithmetic differs")
                if stage == "stress" or parameter == family["steps_bohr"][-1]:
                    checks["continuous"] = all(v <= 1e-6 for v in metrics.values())
                    checks["reference_invariants"] = all(
                        abs(row["invariants"][key] - value) <= (1e-9 if stage == "short" else 1e-6)
                        for key, value in algebra(ref).items()
                    )
            else:
                with mp.workdps(80):
                    expected = np.array(exact_tail(1.0, family["dipole"]).tolist(), complex)
                error = float(np.max(abs(cayley(k) - expected)))
                if error != row["exact_S_error"]:
                    raise ValueError("Physical S error arithmetic differs")
                if parameter == family["steps_bohr"][-1]:
                    checks["physical_S"] = error <= 1e-6
                    checks["physical_invariants"] = all(
                        v <= 1e-6 for v in row["invariants"].values()
                    )
        else:
            old = next(
                r
                for r in dense["rows"]
                if (r["model"], r["intervals"], r["backend"], r["repetition"])
                == (model, parameter, "longdouble", 0)
            )
            checks["archived_k_exact"] = row["k_decimal"] == old["k_decimal"]
            checks["archived_diagnostics_exact"] = d == old["diagnostics"]
            calls = row["calls"]
            checks["role_sequence"] = [
                r["role"] for r in calls if r["function"] == "refined_solve"
            ] == roles
            total = sum(r["wall_seconds"] for r in calls)
            if abs(total - row["timed_function_seconds"]) > 1e-10:
                raise ValueError("Timer arithmetic differs")
            checks["timer_partition"] = (
                row["timed_function_seconds"] <= row["kernel_wall_seconds"] + 1e-6
            )
            if stage == "dense_accuracy":
                oracle = next(
                    r
                    for r in dense["oracles"]
                    if (r["model"], r["intervals"]) == (model, parameter)
                )
                error = maximum_decimal_difference(row["k_decimal"], oracle["k_decimal"])
                if error != row["discrete_arithmetic_error"]:
                    raise ValueError("Dense oracle error differs")
                checks["discrete_arithmetic"] = error <= 1e-12
        if checks != row["checks"]:
            raise ValueError("Row decisions differ at " + str(index))
        reconstructed.update({str(index) + "_" + key: bool(value) for key, value in checks.items()})
        previous = row
    for stage in ["short", "physical"]:
        if not any(s["stage"] == stage for s in report["stages"]):
            continue
        for model in sorted({r["model"] for r in report["rows"] if r["stage"] == stage}):
            for backend in ["reference", "native"]:
                family_rows = [
                    r
                    for r in report["rows"]
                    if (r["stage"], r["model"], r["backend"]) == (stage, model, backend)
                ]
                errors = [
                    r["continuous_errors"]["k"] if stage == "short" else r["exact_S_error"]
                    for r in family_rows
                ]
                reconstructed[f"{stage}_{model}_{backend}_trend"] = (
                    errors[-1] <= errors[-2] / 8 or max(errors[-2:]) <= 1e-10
                )
    complete = len(report["rows"]) == len(coordinates)
    if "declared_row_count" in report["checks"]:
        if not complete:
            raise ValueError("Final counts precede complete coordinate coverage")
        reconstructed.update(
            declared_row_count=len(report["rows"]) == c["new_propagations"],
            declared_interval_count=sum(r["intervals"] for r in report["rows"])
            == c["new_radial_intervals"],
            declared_large_solve_count=sum(
                r["diagnostics"]["large_system_solves"] for r in report["rows"]
            )
            == c["expected_large_system_solves"],
            declared_audited_solve_count=sum(
                r["residual_audit"]["accepted_solves"]
                for r in report["rows"]
                if r["residual_audit"] is not None
            )
            == c["accuracy_native_audited_solves"],
        )
        if (
            timing_summary(
                [r for r in report["rows"] if r["stage"] == "timing"],
                c["timing_speedup_hypothesis_minimum"],
            )
            != report["timing_analysis"]
        ):
            raise ValueError("Timing summary differs")
    if reconstructed != report["checks"]:
        raise ValueError("Full decision map differs")
    first = 0
    for stage in report["stages"]:
        rows = [r for r in report["rows"] if r["stage"] == stage["stage"]]
        gate_values = [value for r in rows for value in r["checks"].values()]
        gate_values += [
            value for key, value in reconstructed.items() if key.startswith(stage["stage"] + "_")
        ]
        if (
            stage["first_row"] != first
            or stage["rows"] != len(rows)
            or stage["verdict"] != all(gate_values)
        ):
            raise ValueError("Stage summary differs")
        first += len(rows)
    if report["experiment_verdict"] is True and not (
        complete and "declared_row_count" in reconstructed and all(reconstructed.values())
    ):
        raise ValueError("Scientific verdict differs")
    return dict(
        saved_decisions_reconstructed=True,
        rows=len(report["rows"]),
        checks=len(reconstructed),
        checks_passed=sum(reconstructed.values()),
        experiment_verdict=report["experiment_verdict"],
        timing_hypotheses_are_observations=True,
        residual_audit_scope=(
            "Saved aggregate observations; numerical reexecution remeasures each residual"
        ),
    )


def numerical_records(report):
    """Timing partitions are new observations; numeric rows/audits are invariant."""
    rows = copy.deepcopy(report["rows"])
    for row in rows:
        for key in [
            "wall_seconds",
            "kernel_wall_seconds",
            "timed_function_seconds",
            "timing_by_role",
            "recurring_solve_seconds",
            "setup_solve_seconds",
            "untimed_remainder_seconds",
        ]:
            row.pop(key, None)
        for call in row.get("calls", []):
            call.pop("wall_seconds")
        row["checks"].pop("timer_partition", None)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ["report", "priors", "output"]:
        parser.add_argument("--" + key, type=Path, required=True)
    parser.add_argument("--reexecute", action="store_true")
    args = parser.parse_args()
    report = json.loads(args.report.read_text())
    summary = audit(report, args.priors)
    if args.reexecute:
        args.output.mkdir()
        fresh = study(report["contract"], args.priors, args.output / "report.json")
        summary["fresh_saved_audit"] = audit(fresh, args.priors)
        summary["numeric_rows_exact"] = numerical_records(fresh) == numerical_records(report)
        summary["fresh_verdict_matches"] = (
            fresh["experiment_verdict"] is report["experiment_verdict"]
        )
        summary["fresh_scientific_exit"] = 0 if fresh["experiment_verdict"] is True else 1
        (args.output / "verification.json").write_text(json.dumps(summary, indent=2) + "\n")
        if not summary["numeric_rows_exact"] or not summary["fresh_verdict_matches"]:
            raise SystemExit(1)
    else:
        args.output.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
