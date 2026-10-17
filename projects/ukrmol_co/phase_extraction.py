"""Prospective phase-model selection with untouched audits and width counterexamples."""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from projects.ukrmol_co.phase_fit import fit_phase, predict_phase


def _phase_error(prediction, observed):
    return float(np.max(abs((prediction - observed + np.pi / 2) % np.pi - np.pi / 2)))


def _valid(fit):
    return fit.get("optimizer_success", False) and not fit.get("at_parameter_bound", True)


def _starts(energy, phase, order, settings):
    fits = []
    for center_fraction, width_fraction in settings["grid_only_starts"]:
        try:
            span = np.ptp(energy)
            fit = fit_phase(
                energy, phase, order, energy[0] + center_fraction * span, width_fraction * span
            )
        except (ValueError, np.linalg.LinAlgError) as error:
            fit = {"optimizer_success": False, "at_parameter_bound": None, "error": str(error)}
        fits.append(fit | {"grid_only_seed": [center_fraction, width_fraction]})
    return fits


def width_alternatives(energy: np.ndarray, phase: np.ndarray, fit: dict, settings: dict) -> list:
    """Construct beyond-budget widths indistinguishable at the declared phase tolerance.

    Finding one is a counterexample; not finding one is not a uniqueness proof.
    All inputs are development data, never the untouched audit phases.
    """
    midpoint, scale = fit["background_midpoint_hartree"], fit["background_scale_hartree"]
    x = (energy - midpoint) / scale
    background = np.polynomial.polynomial.polyvander(x, fit["nback"] - 1)
    unwrapped = np.unwrap(2 * phase) / 2
    width = fit["full_width_hartree"]
    displacement = settings["alternative_budget_factor"] * max(
        settings["full_width_floor_ev"] / settings["hartree_to_ev"],
        settings["full_width_relative_budget"] * width,
    )
    alternatives = []
    for candidate_width in (width - displacement, width + displacement):
        if candidate_width <= 0:
            continue

        def components(parameter, fixed_width=candidate_width):
            center = midpoint + parameter[0] * scale
            resonant = np.arctan2(fixed_width / 2, center - energy)
            coefficients = np.linalg.lstsq(background, unwrapped - resonant, rcond=None)[0]
            return resonant + background @ coefficients, coefficients

        trials = []
        for center_fraction, _ in settings["grid_only_starts"]:
            center = energy[0] + center_fraction * np.ptp(energy)
            solution = least_squares(
                lambda p: components(p)[0] - unwrapped,
                [(center - midpoint) / scale],
                bounds=([-1], [1]),
                ftol=1e-13,
                xtol=1e-13,
                gtol=1e-13,
                jac="3-point",
                max_nfev=2000,
            )
            prediction, coefficients = components(solution.x)
            error = _phase_error(prediction, phase)
            trials.append(
                {
                    "position_hartree": float(midpoint + solution.x[0] * scale),
                    "full_width_hartree": float(candidate_width),
                    "background_coefficients_radian": coefficients.tolist(),
                    "background_midpoint_hartree": midpoint,
                    "background_scale_hartree": scale,
                    "nback": fit["nback"],
                    "optimizer_success": bool(solution.success),
                    "at_parameter_bound": bool(np.any(solution.active_mask)),
                    "maximum_phase_error_rad_mod_pi": error,
                    "development_prediction_rad": prediction.tolist(),
                    "compatible": bool(
                        solution.success
                        and not np.any(solution.active_mask)
                        and error <= settings["phase_fidelity_rad_mod_pi"]
                    ),
                }
            )
        alternatives.extend(trials)
    return alternatives


def estimate_phase(energy: np.ndarray, phase: np.ndarray, windows: list, settings: dict) -> dict:
    """Select backgrounds prospectively; return a provisional verdict and every fit.

    Hartree energies/full widths and radian phases are used throughout fitting.
    Every window must qualify. No post-audit reselection or native pole claim.
    """
    energy, phase = np.asarray(energy, dtype=float), np.asarray(phase, dtype=float)
    if (
        energy.ndim != 1
        or energy.shape != phase.shape
        or len(energy) < 12
        or not np.all(np.isfinite(energy))
        or not np.all(np.isfinite(phase))
        or not np.all(np.diff(energy) > 0)
    ):
        raise ValueError("Finite phases on a complete increasing energy grid are required")
    indices = np.arange(len(energy))
    audit = indices % settings["audit_modulus"] == settings["audit_remainder"]
    rows, choices = [], []
    for low, high in windows:
        in_window = (energy >= low - 1e-14) & (energy <= high + 1e-14)
        develop = np.flatnonzero(in_window & ~audit)
        held = np.flatnonzero(in_window & audit)
        admissible = []
        for order in settings["background_orders"]:
            row = {
                "window_hartree": [low, high],
                "nback": order,
                "development_indices": develop.tolist(),
                "audit_indices": held.tolist(),
                "starts": [],
                "folds": [],
                "valid": False,
                "predictively_admissible": False,
            }
            rows.append(row)
            if len(develop) < 12 or not len(held):
                row["error"] = "Insufficient development or untouched audit points"
                continue
            row["starts"] = _starts(energy[develop], phase[develop], order, settings)
            for fold in range(settings["inner_folds"]):
                test = np.arange(len(develop)) % settings["inner_folds"] == fold
                train, validation = develop[~test], develop[test]
                fits = _starts(energy[train], phase[train], order, settings)
                valid = all(_valid(f) for f in fits)
                best = min(fits, key=lambda f: f.get("residual_sum_squares_rad2", np.inf))
                row["folds"].append(
                    {
                        "fold": fold,
                        "training_indices": train.tolist(),
                        "validation_indices": validation.tolist(),
                        "starts": fits,
                        "valid": valid,
                        "maximum_phase_error_rad_mod_pi": (
                            _phase_error(predict_phase(energy[validation], best), phase[validation])
                            if valid
                            else None
                        ),
                    }
                )
            row["valid"] = all(_valid(f) for f in row["starts"]) and all(
                f["valid"] for f in row["folds"]
            )
            if not row["valid"]:
                continue
            row["inner_validation_maximum_phase_error_rad_mod_pi"] = max(
                f["maximum_phase_error_rad_mod_pi"] for f in row["folds"]
            )
            row["predictively_admissible"] = bool(
                row["inner_validation_maximum_phase_error_rad_mod_pi"]
                <= settings["phase_fidelity_rad_mod_pi"]
            )
            if row["predictively_admissible"]:
                admissible.append(order)
        choices.append(
            {"window_hartree": [low, high], "selected_order": min(admissible, default=None)}
        )
    result = {
        "rows": rows,
        "frozen_choices_before_audit": choices,
        "estimator_verdict": None,
        "scope": "Provisional phase representation; no native pole or global uniqueness "
        "certification",
    }
    selected_rows = [row for row in rows if row["predictively_admissible"]]
    if not selected_rows:
        result["estimator_verdict"] = False if all(r["valid"] for r in rows) else None
        return result
    envelope, audits, ambiguous = [], [], False
    for row in selected_rows:
        develop, held = row["development_indices"], row["audit_indices"]
        best = min(row["starts"], key=lambda f: f["residual_sum_squares_rad2"])
        envelope.extend(row["starts"])
        envelope.extend(fit for fold in row["folds"] for fit in fold["starts"])
        row["audit_phase_errors_rad_mod_pi"] = [
            _phase_error(predict_phase(energy[held], fit), phase[held]) for fit in row["starts"]
        ]
        audits.extend(row["audit_phase_errors_rad_mod_pi"])
        row["width_alternatives"] = width_alternatives(
            energy[develop], phase[develop], best, settings
        )
        ambiguous |= any(a["compatible"] for a in row["width_alternatives"])
    conversion = settings["hartree_to_ev"]
    centers = np.asarray([f["position_hartree"] for f in envelope]) * conversion
    widths = np.asarray([f["full_width_hartree"] for f in envelope]) * conversion
    width_budget = max(
        settings["full_width_floor_ev"], settings["full_width_relative_budget"] * np.median(widths)
    )
    result.update(
        position_range_ev=[float(min(centers)), float(max(centers))],
        full_width_range_ev=[float(min(widths)), float(max(widths))],
        position_span_ev=float(np.ptp(centers)),
        full_width_span_ev=float(np.ptp(widths)),
        full_width_budget_ev=float(width_budget),
        maximum_audit_phase_error_rad_mod_pi=max(audits),
        beyond_budget_width_counterexample_found=bool(ambiguous),
        reference_fit=min(selected_rows[0]["starts"], key=lambda f: f["residual_sum_squares_rad2"]),
    )
    if all(row["valid"] for row in rows):
        result["estimator_verdict"] = bool(
            all(choice["selected_order"] is not None for choice in choices)
            and np.ptp(centers) <= settings["position_budget_ev"]
            and np.ptp(widths) <= width_budget
            and max(audits) <= settings["phase_fidelity_rad_mod_pi"]
            and not ambiguous
        )
    return result
