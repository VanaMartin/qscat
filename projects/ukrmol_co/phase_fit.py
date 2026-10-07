"""Independent one-candidate eigenphase fits with a polynomial background.

Energies and full widths are Hartree; phases are radians. Fit residuals and
Jacobian conditioning diagnose the representation, not physical uncertainty.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares


def fit_phase(
    energy: np.ndarray,
    phase: np.ndarray,
    nback: int,
    initial_position: float,
    initial_width: float,
) -> dict:
    """Fit a positive-width phase step, eliminating the linear background exactly."""
    energy = np.asarray(energy, dtype=float)
    phase = np.asarray(phase, dtype=float)
    if (
        energy.ndim != 1
        or energy.shape != phase.shape
        or nback not in range(1, 5)
        or energy.size <= nback + 2
        or not np.all(np.isfinite(energy))
        or not np.all(np.isfinite(phase))
        or not np.all(np.diff(energy) > 0)
        or not np.isfinite(initial_position)
        or not energy[0] < initial_position < energy[-1]
        or not np.isfinite(initial_width)
        or initial_width <= 0
    ):
        raise ValueError(
            "A finite increasing grid, interior center and positive width are required"
        )
    phase = np.unwrap(2 * phase) / 2
    midpoint = float((energy[-1] + energy[0]) / 2)
    scale = float((energy[-1] - energy[0]) / 2)
    x = (energy - midpoint) / scale
    background = np.polynomial.polynomial.polyvander(x, nback - 1)

    def components(parameters):
        center = midpoint + parameters[0] * scale
        width = np.exp(parameters[1]) * scale
        resonant = np.arctan2(width / 2, center - energy)
        coefficients = np.linalg.lstsq(background, phase - resonant, rcond=None)[0]
        return resonant + background @ coefficients, coefficients

    def residual(parameters):
        return components(parameters)[0] - phase

    width_ratio = initial_width / scale
    if not 1e-6 < width_ratio < 20:
        raise ValueError("Initial width is outside the diagnostic width bounds")
    solution = least_squares(
        residual,
        [(initial_position - midpoint) / scale, np.log(width_ratio)],
        bounds=([-1, np.log(1e-6)], [1, np.log(20)]),
        ftol=1e-13,
        xtol=1e-13,
        gtol=1e-13,
        max_nfev=2000,
        jac="3-point",
    )
    prediction, coefficients = components(solution.x)
    center = float(midpoint + solution.x[0] * scale)
    width = float(np.exp(solution.x[1]) * scale)
    half_width = width / 2
    denominator = (center - energy) ** 2 + half_width**2
    jacobian = np.column_stack(
        (
            -scale * half_width / denominator,
            (center - energy) * half_width / denominator,
            background,
        )
    )
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    difference = prediction - phase
    return {
        "position_hartree": center,
        "full_width_hartree": width,
        "background_coefficients_radian": coefficients.tolist(),
        "background_midpoint_hartree": midpoint,
        "background_scale_hartree": scale,
        "nback": nback,
        "points": int(energy.size),
        "residual_rms_rad": float(np.sqrt(np.mean(difference**2))),
        "residual_max_abs_rad": float(np.max(np.abs(difference))),
        "residual_sum_squares_rad2": float(difference @ difference),
        "jacobian_singular_values": singular_values.tolist(),
        "jacobian_condition_number": float(singular_values[0] / singular_values[-1]),
        "optimizer_success": bool(solution.success),
        "optimizer_message": solution.message,
        "optimizer_optimality": float(solution.optimality),
        "optimizer_evaluations": int(solution.nfev),
        "at_parameter_bound": bool(np.any(solution.active_mask)),
        "scope": "One-candidate phase representation; no pole or uncertainty certification",
    }


def predict_phase(energy: np.ndarray, fit: dict) -> np.ndarray:
    """Evaluate a recorded phase model on the requested Hartree energies."""
    x = (np.asarray(energy) - fit["background_midpoint_hartree"]) / fit["background_scale_hartree"]
    return np.arctan2(fit["full_width_hartree"] / 2, fit["position_hartree"] - energy) + (
        np.polynomial.polynomial.polyval(x, fit["background_coefficients_radian"])
    )
