"""Independent saved-CO prefix decision audit and fixed six-prefix replay."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np

from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.native_full_reconstruction import digest, load_model, verify_inputs
from projects.ukrmol_co.numerov_co_canary import oracle_checks, prefix, timing_summary
from projects.ukrmol_co.numerov_native_residual import native_residual
from projects.ukrmol_co.numerov_native_residual_control import canonical_hash


def coordinates(c):
    result = []
    for item in c["identities"]:
        for backend in ["reference", "native"]:
            result.append(("accuracy", item["case"], 8, backend, 0))
    items = [x for x in c["identities"] if x["timing"]]
    for rep in range(2):
        for item in items if rep == 0 else list(reversed(items)):
            for n in [32, 128] if rep == 0 else [128, 32]:
                for backend in ["reference", "native"] if rep == 0 else ["native", "reference"]:
                    result.append(("timing", item["case"], n, backend, rep))
    return result


def numerical_row(row):
    result = copy.deepcopy(row)
    for name in ["kernel_wall_seconds", "recurring_seconds", "checks"]:
        result.pop(name, None)
    for call in result["calls"]:
        call.pop("wall_seconds")
    return result


def audit(report, operands):
    c = report["contract"]
    expected = coordinates(c)
    checks = {}
    previous = None
    for i, row in enumerate(report["rows"]):
        if (
            row["stage"],
            row["case"],
            row["intervals"],
            row["backend"],
            row["repetition"],
        ) != expected[i]:
            raise ValueError("Unexpected row coordinate/order")
        item = next(x for x in c["identities"] if x["case"] == row["case"])
        if (item["energy_index"], item["energy_rydberg_hex"]) != (
            row["energy_index"],
            row["energy_rydberg_hex"],
        ):
            raise ValueError("Energy identity differs")
        n = row["intervals"]
        roles = ["initialization"] + ["recurring_potential", "recurring_defect"] * (n - 1)
        stats = row["diagnostics"]
        decisions = dict(
            solve_count=stats["large_system_solves"] == 2 * n - 1,
            backward_error=stats["maximum_solve_backward_error"] <= 32 * 2.0**-63,
            corrections=stats["maximum_refinement_corrections"] <= 3,
            role_sequence=[x["role"] for x in row["calls"]] == roles,
            stopped_before_matching=row["matching_evaluations"] == 0,
        )
        if row["last_recurrence_radius_bohr"] != 18 + (n - 1) * c["step_bohr"]:
            raise ValueError("Recurrence radius differs")
        if i % 2:
            decisions.update(
                final_delta_exact=row["final_delta_sha256"] == previous["final_delta_sha256"],
                diagnostics_exact=stats == previous["diagnostics"],
                sampled_operands_exact=row["samples"] == previous["samples"],
            )
        if row["audited"] != (row["stage"] == "accuracy"):
            raise ValueError("Accuracy/timing instrumentation differs")
        if row["audited"]:
            decisions["condition_estimates"] = len(row["samples"]) == 2 * n - 1 and all(
                x["reciprocal_condition_1"] > 0 for x in row["samples"]
            )
        else:
            if row["samples"] or row["residual_audit"] is not None:
                raise ValueError("Timing has accuracy instrumentation")
            if sum(x["wall_seconds"] for x in row["calls"]) > row["kernel_wall_seconds"] + 1e-6:
                raise ValueError("Disjoint solve times exceed prefix time")
        observed = row["residual_audit"]
        if (observed is not None) != (row["stage"] == "accuracy" and row["backend"] == "native"):
            raise ValueError("Unexpected residual observer")
        if observed is not None:
            decisions.update(
                residuals_exact=observed["evaluations"] == observed["exact_residuals"],
                gate_decisions=observed["gate_disagreements"] == 0,
                accepted_solves=observed["accepted_solves"] == 2 * n - 1,
                independent_backward_error=observed["maximum_accepted_numpy_error"]
                <= 32 * 2.0**-63,
                precision=observed["precision_preserved"],
            )
        if decisions != row["checks"]:
            raise ValueError("Saved row decisions differ: " + str(i))
        checks.update({f"{i}_{name}": bool(value) for name, value in decisions.items()})
        previous = row
    fresh_oracles = []
    for item in c["identities"]:
        directory = operands / item["case"]
        for oracle in oracle_checks(directory):
            oracle["positions"] = [list(position) for position in oracle["positions"]]
            oracle["case"] = item["case"]
            with np.load(directory / oracle["file"]) as data:
                a, b, x = [data[k] for k in ["a", "b", "x"]]
                if any(v.shape != (298, 298) for v in [a, b, x]):
                    raise ValueError("Captured operand shape differs")
                if not np.array_equal(native_residual(a, b, x), b - a @ x):
                    raise ValueError("Captured full native residual differs")
                index = int(oracle["file"][6:9])
                row = next(
                    r
                    for r in report["rows"]
                    if r["case"] == item["case"]
                    and r["stage"] == "accuracy"
                    and r["backend"] == "native"
                )
                sample = row["samples"][index]
                for key, array in zip(["a", "b", "x"], [a, b, x], strict=True):
                    if sample[key + "_sha256"] != canonical_hash(array):
                        raise ValueError("Captured operand does not match actual solve")
            checks.update(
                {f"oracle_{len(fresh_oracles)}_{k}": bool(v) for k, v in oracle["checks"].items()}
            )
            fresh_oracles.append(oracle)
    if fresh_oracles != report["oracles"]:
        raise ValueError("Independent captured references differ")
    if len(report["rows"]) == len(expected):
        checks.update(
            row_count=len(report["rows"]) == c["prefix_executions"],
            solve_count=sum(r["diagnostics"]["large_system_solves"] for r in report["rows"])
            == c["mixed_solves"],
            oracle_count=len(report["oracles"]) == 12,
        )
        if (
            timing_summary([r for r in report["rows"] if r["stage"] == "timing"])
            != report["timing_analysis"]
        ):
            raise ValueError("Timing arithmetic differs")
    if checks != report["checks"]:
        raise ValueError("Complete decision map differs")
    if report["experiment_verdict"] is not True or not all(checks.values()):
        raise ValueError("Canary qualification did not pass")
    return dict(
        experiment_verdict=True,
        saved_decisions_reconstructed=True,
        rows=len(report["rows"]),
        checks_passed=sum(checks.values()),
        checks_total=len(checks),
        full_captured_operand_triples=12,
        paired_oracle_entries=36,
        all_full_captured_residuals_exact=True,
        timing_identity_is_not_a_gate=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["report", "operands", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--raw", type=Path)
    a = parser.parse_args()
    report = json.loads(a.report.read_text())
    summary = audit(report, a.operands)
    if a.raw is not None:
        raw_contract_path = Path(__file__).with_name("numerov-corner-contract.json")
        if digest(raw_contract_path) != report["contract"]["raw_contract_sha256"]:
            raise ValueError("Raw input contract changed")
        raw_contract = json.loads(raw_contract_path.read_text())
        verify_inputs(a.raw, raw_contract)
        rows = []
        for item in report["contract"]["identities"]:
            model = load_model(a.raw / item["case"])
            energy = model["energies"][item["energy_index"]]
            if float(energy).hex() != item["energy_rydberg_hex"]:
                raise ValueError("Replay energy identity differs")
            q = surface_matrix(model["poles"], model["amplitudes"], model["ground"] + energy / 2)
            args = [q] + [model[key] for key in ["thresholds", "angular", "c2", "c3"]]
            fresh = prefix(args, energy, 8, report["contract"]["step_bohr"], "native", audited=True)
            fresh.update(
                case=item["case"],
                energy_index=item["energy_index"],
                energy_rydberg_hex=float(energy).hex(),
                stage="accuracy",
                repetition=0,
            )
            old = next(
                r
                for r in report["rows"]
                if r["case"] == item["case"]
                and r["stage"] == "accuracy"
                and r["backend"] == "native"
            )
            if numerical_row(fresh) != numerical_row(old):
                raise ValueError("Fresh audited CO prefix differs")
            rows.append(fresh)
        verify_inputs(a.raw, raw_contract)
        summary.update(
            fresh_native_prefixes=6,
            fresh_native_solves=90,
            all_fresh_numerical_records_exact=True,
            rows_reexecuted=rows,
        )
    a.output.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
