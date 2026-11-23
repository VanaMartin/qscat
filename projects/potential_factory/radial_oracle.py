"""Independent scalar radial oracle for synthetic p-wave Gaussian wells."""

from __future__ import annotations

import cmath

import numpy as np
from qscat.exceptions import ConvergenceError
from scipy.integrate import solve_ivp

__all__ = ["gaussian_p_wave_mismatch", "p_wave_boundary", "regular_p_wave"]


def regular_p_wave(
    energy: complex, lam: float, alpha: float, *, radius: float, rtol: float
) -> tuple[complex, complex]:
    """Integrate regular `u=r²f`; return `(f,f')` at a real radius in bohr."""
    origin = 1e-4
    if (
        not np.isfinite([energy, lam, alpha, radius, rtol]).all()
        or lam < 0
        or alpha <= 0
        or radius <= origin
        or rtol <= 0
    ):
        raise ValueError(
            "Finite energy, nonnegative depth and positive range/radius/tolerance required"
        )
    c2 = -(energy + lam) / 5
    c4 = ((energy + lam) ** 2 / 5 + lam * alpha) / 14
    initial = np.array(
        [1 + c2 * origin**2 + c4 * origin**4, 2 * c2 * origin + 4 * c4 * origin**3], complex
    )

    def rhs(r, values):
        f, df = values
        return np.array([df, -4 * df / r - 2 * (energy + lam * np.exp(-alpha * r * r)) * f])

    solution = solve_ivp(
        rhs, (origin, radius), initial, method="DOP853", rtol=rtol, atol=0.01 * rtol
    )
    if not solution.success or not np.isfinite(solution.y[:, -1]).all():
        raise ConvergenceError(f"Regular p-wave integration failed: {solution.message}")
    return complex(solution.y[0, -1]), complex(solution.y[1, -1])


def _free_wave(k: complex, radius: complex, sign: int) -> tuple[complex, complex]:
    phase = cmath.exp(sign * 1j * k * radius)
    factor = 1 + sign * 1j / (k * radius)
    return phase * factor, phase * (sign * 1j * k * factor - sign * 1j / (k * radius**2))


def p_wave_boundary(
    energy: complex, radius: complex, *, endpoint: complex | None = None
) -> tuple[complex, complex]:
    """Free outgoing value/derivative, or a solution vanishing at a finite endpoint."""
    if not np.isfinite([energy, radius]).all() or energy == 0 or radius == 0:
        raise ValueError("Finite nonzero energy and radius required")
    k = cmath.sqrt(2 * energy)
    outgoing, doutgoing = _free_wave(k, radius, 1)
    if endpoint is None:
        return outgoing, doutgoing
    if not np.isfinite(endpoint) or endpoint == 0:
        raise ValueError("Finite nonzero endpoint required")
    incoming, dincoming = _free_wave(k, radius, -1)
    end_out, _ = _free_wave(k, endpoint, 1)
    end_in, _ = _free_wave(k, endpoint, -1)
    return (outgoing * end_in - incoming * end_out, doutgoing * end_in - dincoming * end_out)


def gaussian_p_wave_mismatch(
    energy: complex,
    lam: float,
    alpha: float,
    *,
    radius: float,
    rtol: float,
    endpoint: complex | None = None,
) -> complex:
    """Dimensionless Wronskian residual; neglect the Gaussian beyond `radius`."""
    f, df = regular_p_wave(energy, lam, alpha, radius=radius, rtol=rtol)
    b, db = p_wave_boundary(energy, radius, endpoint=endpoint)
    k = cmath.sqrt(2 * energy)
    du_over_k = (df + 2 * f / radius) / k
    db_over_k = db / k
    scale = np.linalg.norm([f, du_over_k]) * np.linalg.norm([b, db_over_k])
    if not np.isfinite(scale) or scale == 0:
        raise ConvergenceError("Nonfinite or zero radial Wronskian scale")
    return complex((f * db_over_k - du_over_k * b) / scale)
