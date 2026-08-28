"""The angular factor of a spherical harmonic, normalized for the one
integral this project performs.

Y_{lm}(theta, phi) = Theta_{lm}(cos theta) * exp(i m phi) / sqrt(2 pi).

The two-centre potential is independent of phi, so every angular matrix
element collapses to a single integral over x = cos(theta) carrying two
Theta factors and no 2 pi. Theta is normalized so that the integral of
Theta_{lm} Theta_{l'm} over x in [-1, 1] is the Kronecker delta -- which
makes an ISOTROPIC potential give a diagonal, unscaled channel matrix, the
identity the whole embedding gate rests on.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from scipy.special import gammaln, lpmv

__all__ = ["theta_lm"]


def theta_lm(l: int, m: int, x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """`Theta_{lm}(x)`, orthonormal in `x = cos(theta)` on `[-1, 1]`.

    `x` is real by construction: exterior complex scaling rotates the RADIAL
    coordinate, never the angular one.
    """
    if l < abs(m):
        raise ValueError(f"theta_lm requires l >= |m|, got l={l}, m={m}")
    norm = np.sqrt((2 * l + 1) / 2 * np.exp(gammaln(l - m + 1) - gammaln(l + m + 1)))
    return np.asarray(norm * lpmv(m, l, np.asarray(x, dtype=np.float64)), dtype=np.float64)
