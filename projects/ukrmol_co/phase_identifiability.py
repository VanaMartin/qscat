"""Known wrong poles test the limits of local unitary phase representations."""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import time
from pathlib import Path

import numpy as np


def oracle_s(energy: np.ndarray, case: dict) -> np.ndarray:
    """Generate the independent nonpolynomial-background S factor on real energies."""
    center, width = case["position_hartree"], case["full_width_hartree"]
    x = (np.asarray(energy) - center) / width
    background = sum(c * x**i for i, c in enumerate(case["background"]))
    background += case["exponential_background_amplitude"] * np.exp(x / 3)
    return np.exp(2j * background) * (
        (center - energy + 0.5j * width) / (center - energy - 0.5j * width)
    )


def candidate_s(energy: np.ndarray, fit: dict) -> np.ndarray:
    """Evaluate an entire nonzero background times an exactly known pole factor."""
    x = (np.asarray(energy) - fit["midpoint_hartree"]) / fit["scale_hartree"]
    background = np.polynomial.chebyshev.chebval(x, fit["background_chebyshev_coefficients_rad"])
    displacement = fit["position_hartree"] - energy
    half = fit["full_width_hartree"] / 2
    return np.exp(2j * background) * (displacement + 1j * half) / (displacement - 1j * half)


def fit_background(
    energy: np.ndarray, values: np.ndarray, window: list, pole: list, degree: int
) -> dict:
    """Linear projection at a fixed pole, with no access to audit phases or grids."""
    midpoint, scale = (window[0] + window[1]) / 2, (window[1] - window[0]) / 2
    matrix = np.polynomial.chebyshev.chebvander((energy - midpoint) / scale, degree)
    # Unwrap only development phases; no full-grid branch information is supplied.
    phase = np.unwrap(np.angle(values)) / 2
    resonance = np.arctan2(pole[1] / 2, pole[0] - energy)
    coefficients, _, rank, singular = np.linalg.lstsq(matrix, phase - resonance, rcond=None)
    return {
        "position_hartree": float(pole[0]),
        "full_width_hartree": float(pole[1]),
        "midpoint_hartree": midpoint,
        "scale_hartree": scale,
        "background_degree": degree,
        "background_chebyshev_coefficients_rad": coefficients.tolist(),
        "background_matrix_rank": int(rank),
        "background_matrix_singular_values": singular.tolist(),
        "background_matrix_condition_number": float(singular[0] / singular[-1]),
        "full_rank_and_finite": bool(rank == degree + 1 and np.all(np.isfinite(coefficients))),
    }


def _score(energy, case, fit):
    truth, candidate = oracle_s(energy, case), candidate_s(energy, fit)
    errors = abs(np.angle(candidate * truth.conj())) / 2
    worst = int(np.argmax(errors))
    return {
        "points": len(energy),
        "maximum_phase_error_rad_mod_pi": float(errors[worst]),
        "maximum_error_energy_hartree": float(energy[worst]),
        "maximum_unitarity_error": float(
            max(np.max(abs(abs(truth) - 1)), np.max(abs(abs(candidate) - 1)))
        ),
    }


def audit(contract: dict) -> dict:
    """Retain the complete fixed matrix and gate mandatory true-pole controls."""
    start, rows, grids = time.monotonic(), [], []
    for case in contract["cases"]:
        for points in contract["grid_points"]:
            grid = np.linspace(*case["range_hartree"], points)
            reserve = np.arange(points) % contract["audit_modulus"] == contract["audit_remainder"]
            for window in case["windows_hartree"]:
                mask = (grid >= window[0] - 1e-14) & (grid <= window[1] + 1e-14)
                develop, held = grid[mask & ~reserve], grid[mask & reserve]
                low, high = window
                exclusion = contract["development_exclusion_fraction"] * (high - low)
                develop = develop[
                    np.all(abs(develop[:, None] - np.array(window)) > exclusion, axis=1)
                ]
                dense = low + (high - low) * (
                    (np.arange(contract["dense_midpoint_grid_points"]) + 0.5)
                    / contract["dense_midpoint_grid_points"]
                )
                dense = dense[np.all(abs(dense[:, None] - develop) > exclusion, axis=1)]
                index = len(grids)
                grids.append(
                    {
                        "case": case["name"],
                        "original_grid_points": points,
                        "window_hartree": window,
                        "development_energy_hartree": develop.tolist(),
                        "reserved_energy_hartree": held.tolist(),
                        "dense_audit_energy_hartree": dense.tolist(),
                        "endpoint_energy_hartree": window,
                    }
                )
                for degree in contract["background_degrees"]:
                    for factor in contract["exact_pole_width_factors"]:
                        width = factor * case["full_width_hartree"]
                        fit = fit_background(
                            develop,
                            oracle_s(develop, case),
                            window,
                            [case["position_hartree"], width],
                            degree,
                        )
                        scores = {
                            name: _score(energy, case, fit)
                            for name, energy in (
                                ("development", develop),
                                ("reserved", held),
                                ("dense", dense),
                                ("endpoints", np.array(window)),
                            )
                        }
                        audits_pass = (
                            all(
                                v["maximum_phase_error_rad_mod_pi"]
                                <= contract["phase_fidelity_rad_mod_pi"]
                                and v["maximum_unitarity_error"] <= contract["unitarity_budget"]
                                for name, v in scores.items()
                                if name != "development"
                            )
                            if fit["full_rank_and_finite"]
                            else None
                        )
                        beyond = abs(width - case["full_width_hartree"]) > max(
                            contract["full_width_floor_ev"] / contract["hartree_to_ev"],
                            contract["full_width_relative_budget"] * case["full_width_hartree"],
                        )
                        rows.append(
                            {
                                "case": case["name"],
                                "original_grid_points": points,
                                "window_hartree": window,
                                "grid_record": index,
                                "exact_pole_width_factor": factor,
                                "fit": fit,
                                "scores": scores,
                                "all_audit_sets_passed": audits_pass,
                                "known_width_error_exceeds_scientific_budget": bool(beyond),
                                "width_counterexample": bool(audits_pass is True and beyond),
                                "mandatory_true_pole_control": factor == 1 and degree == 11,
                            }
                        )
    required = [row for row in rows if row["mandatory_true_pole_control"]]
    passed = (
        len(rows) == contract["expected_cells"]
        and len(required) == 45
        and all(row["all_audit_sets_passed"] is True for row in required)
    )
    return {
        "verified_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "contract": contract,
        "source_sha256": {
            str(Path(__file__).relative_to(Path(__file__).resolve().parents[2])): hashlib.sha256(
                Path(__file__).read_bytes()
            ).hexdigest()
        },
        "grids": grids,
        "cells": rows,
        "cell_count": len(rows),
        "mandatory_control_count": len(required),
        "mandatory_controls_passed": sum(row["all_audit_sets_passed"] is True for row in required),
        "representation_qualification_verdict": bool(passed),
        "width_counterexample_count": sum(row["width_counterexample"] for row in rows),
        "reader_wall_seconds": time.monotonic() - start,
        "native_application_released": False,
        "native_solves": 0,
        "full_sweep_released": False,
        "scope": "Finite unused-grid local analytic/unitary counterexamples. Exponential "
        "polynomial backgrounds are not qualified causal finite-range CO models. No "
        "native fit, order selection, interval supremum bound or uniqueness certification.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(json.loads(args.contract.read_text()))
    with args.output.open("x") as stream:
        stream.write(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if report["representation_qualification_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
