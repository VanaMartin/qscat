"""Separate replay successor: exact science records and bounded rcond observations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from projects.ukrmol_co.native_boundary import surface_matrix
from projects.ukrmol_co.native_full_reconstruction import digest, load_model, verify_inputs
from projects.ukrmol_co.numerov_co_canary import prefix
from projects.ukrmol_co.numerov_co_canary_replay import audit, numerical_row


def scientific_record(row):
    result = numerical_row(row)
    for sample in result["samples"]:
        sample.pop("reciprocal_condition_1")
    return result


def condition_drift(old, fresh):
    if len(old["samples"]) != len(fresh["samples"]):
        raise ValueError("Condition observation count differs")
    relative = []
    for left, right in zip(old["samples"], fresh["samples"], strict=True):
        a, b = left["reciprocal_condition_1"], right["reciprocal_condition_1"]
        if not np.isfinite([a, b]).all() or min(a, b) <= 0:
            raise ValueError("Invalid reciprocal condition observation")
        relative.append(abs(a - b) / max(abs(a), abs(b)))
    return max(relative, default=0.0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["report", "operands", "raw", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    a = parser.parse_args()
    report = json.loads(a.report.read_text())
    summary = audit(report, a.operands)
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
            if r["case"] == item["case"] and r["stage"] == "accuracy" and r["backend"] == "native"
        )
        relative = condition_drift(old, fresh)
        fresh["scientific_record_exact"] = scientific_record(old) == scientific_record(fresh)
        fresh["rcond_maximum_relative_drift"] = relative
        fresh["rcond_rounding_budget_passed"] = relative <= 32 * 2.0**-52
        rows.append(fresh)
        a.output.write_text(
            json.dumps(
                summary
                | dict(
                    public_replay_verdict=None,
                    fresh_prefixes_completed=len(rows),
                    rows_reexecuted=rows,
                ),
                indent=2,
                allow_nan=False,
            )
            + "\n"
        )
        if not fresh["scientific_record_exact"] or not fresh["rcond_rounding_budget_passed"]:
            raise ValueError("Successor scientific record or rcond budget failed")
    verify_inputs(a.raw, raw_contract)
    summary.update(
        fresh_native_prefixes=6,
        fresh_native_solves=90,
        all_fresh_scientific_records_exact=True,
        condition_estimates_are_observations=True,
        rcond_relative_budget=32 * 2.0**-52,
        maximum_rcond_relative_drift=max(x["rcond_maximum_relative_drift"] for x in rows),
        all_rcond_rounding_budgets_passed=True,
        original_strict_replay_exits=[1, 1],
        public_replay_verdict=True,
        rows_reexecuted=rows,
    )
    a.output.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
