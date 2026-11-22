"""Six positive local Gaussian shape parameters; no fit, asymptote or state solve."""

from __future__ import annotations

import numpy as np
from numpy.polynomial.chebyshev import chebval


def local_gaussian_parameters(coefficients, radius_bohr) -> tuple[np.ndarray, np.ndarray]:
    """Return lambda(Ha), alpha(a0^-2) on finite real1.9–2.5bohr nuclear R."""
    c, radius = np.asarray(coefficients), np.asarray(radius_bohr)
    if (
        c.shape != (6,)
        or not np.isrealobj(c)
        or not np.isfinite(c).all()
        or not np.isrealobj(radius)
        or radius.dtype.kind not in "fiu"
        or radius.size > 1000
        or not np.isfinite(radius).all()
        or np.any((radius < 1.9) | (radius > 2.5))
    ):
        raise ValueError("Six finite real coefficients and local real nuclear R required")
    x = (radius - 2.2) / 0.3
    logs = np.array([chebval(x, c[:3]), chebval(x, c[3:])])
    if np.max(abs(logs)) > 50:
        raise ValueError("Unsupported exponential parameter scale")
    values = np.exp(logs)
    return values[0], values[1]


def local_gaussian_interaction(
    coefficients, electronic_radius_bohr, nuclear_radius_bohr
) -> np.ndarray:
    """Only -lambda exp(-alpha r²); add V0 and centrifugal terms separately once."""
    lam, alpha = local_gaussian_parameters(coefficients, nuclear_radius_bohr)
    r = np.asarray(electronic_radius_bohr)
    if r.dtype.kind not in "fciu" or r.size > 100_000 or not np.isfinite(r).all():
        raise ValueError("Bounded finite real/complex electronic coordinates required")
    if np.prod(np.broadcast_shapes(r.shape, np.shape(lam))) > 100_000:
        raise ValueError("Broadcast electronic/nuclear point budget exceeded")
    with np.errstate(over="raise", invalid="raise"):
        try:
            value = -lam * np.exp(-alpha * r**2)
        except (FloatingPointError, ValueError) as error:
            raise ValueError("Unsupported electronic contour or broadcast shape") from error
    if not np.isfinite(value).all():
        raise ValueError("Nonfinite local Gaussian interaction")
    return value
