"""Known-pole calibration for the prospective CO phase-extraction rule."""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.phase_extraction import estimate_phase


def analytic_phase(energy: np.ndarray, case: dict) -> tuple[np.ndarray, float]:
    """Generate phases from unitary complex S factors, independently of the fitter."""
    center, width = case["position_hartree"], case["full_width_hartree"]
    x = (energy - center) / width
    background = sum(c * x**i for i, c in enumerate(case["background"]))
    background += case.get("exponential_background_amplitude", 0) * np.exp(x / 3)
    perturbation = case.get("phase_perturbation_rad", 0) * np.sin(17 * x + 0.3)
    s = np.exp(2j * (background + perturbation))
    if not case.get("omit_resonance", False):
        s *= (center - energy + 0.5j * width) / (center - energy - 0.5j * width)
    if "neighbor_position_hartree" in case:
        neighbor, gamma = case["neighbor_position_hartree"], case["neighbor_full_width_hartree"]
        s *= (neighbor - energy + 0.5j * gamma) / (neighbor - energy - 0.5j * gamma)
    return np.angle(s) / 2, float(max(abs(abs(s) - 1)))


def calibrate(contract: dict) -> dict:
    """Retain every frozen case; stop native qualification on any calibration failure."""
    start = time.monotonic()
    settings = contract["estimator"] | {"hartree_to_ev": contract["hartree_to_ev"]}
    calibration, rows = contract["calibration"], []
    for case in calibration["cases"]:
        for points in calibration["grid_points"]:
            energy = np.linspace(*case["range_hartree"], points)
            phase, unitarity = analytic_phase(energy, case)
            extraction = estimate_phase(energy, phase, case["windows_hartree"], settings)
            verdict = extraction["estimator_verdict"]
            record = {
                "case": case["name"],
                "grid_points": points,
                "energy_hartree": energy.tolist(),
                "phase_rad": phase.tolist(),
                "maximum_unitarity_error": unitarity,
                "extraction": extraction,
                "known_pole_truth_passed": None,
                "calibration_requirement_passed": False,
            }
            truth_passed = None
            if not case.get("omit_resonance") and "neighbor_position_hartree" not in case:
                conversion = contract["hartree_to_ev"]
                if "position_range_ev" in extraction:
                    center_error = max(
                        abs(v - case["position_hartree"] * conversion)
                        for v in extraction["position_range_ev"]
                    )
                    width_error = max(
                        abs(v - case["full_width_hartree"] * conversion)
                        for v in extraction["full_width_range_ev"]
                    )
                    position_budget = calibration["truth_position_budget_ev"]
                    width_budget = max(
                        calibration["truth_full_width_floor_ev"],
                        calibration["truth_full_width_relative_budget"]
                        * case["full_width_hartree"]
                        * conversion,
                    )
                    if case["required"] == "accept_with_0.1mev_truth":
                        position_budget = width_budget = calibration["resolved_truth_budget_ev"]
                    truth_passed = bool(
                        center_error <= position_budget and width_error <= width_budget
                    )
                    record.update(
                        maximum_position_truth_error_ev=center_error,
                        maximum_full_width_truth_error_ev=width_error,
                    )
            record["known_pole_truth_passed"] = truth_passed
            requirement = case["required"]
            if requirement == "accept_with_0.1mev_truth":
                passed = verdict is True and truth_passed is True
            elif requirement == "truth_if_accepted":
                passed = verdict is not True or truth_passed is True
            elif requirement == "not_accept":
                passed = verdict is not True
            elif requirement == "stress_diagnostic":
                passed = True
            else:
                raise ValueError(f"Unknown calibration requirement: {requirement}")
            record["calibration_requirement_passed"] = bool(passed and unitarity <= 3e-15)
            rows.append(record)
    return {
        "verified_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "contract": contract,
        "source_sha256": {
            str(p.relative_to(Path(__file__).resolve().parents[2])): hashlib.sha256(
                p.read_bytes()
            ).hexdigest()
            for p in (
                Path(__file__),
                Path(__file__).with_name("phase_extraction.py"),
                Path(__file__).with_name("phase_fit.py"),
            )
        },
        "cases": rows,
        "case_count": len(rows),
        "calibration_verdict": all(row["calibration_requirement_passed"] for row in rows),
        "native_application_released": all(row["calibration_requirement_passed"] for row in rows),
        "accepted_case_count": sum(row["extraction"]["estimator_verdict"] is True for row in rows),
        "reader_wall_seconds": time.monotonic() - start,
        "native_solves": 0,
        "full_sweep_released": False,
        "scope": "Known-pole finite-matrix calibration; acceptance is provisional phase "
        "extraction, not native S-matrix pole certification. All original native extraction "
        "rejections remain preserved.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = calibrate(json.loads(args.contract.read_text()))
    with args.output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if result["calibration_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
