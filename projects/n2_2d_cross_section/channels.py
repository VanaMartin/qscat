"""Asymptotic channel functions for the exact 2-D e-N2 scattering problem.

The entrance/exit channel of a VE transition is a free electron of momentum
`k` in partial wave `l`, times a neutral vibrational state. The electronic
factor is the ENERGY-NORMALIZED regular free radial solution

    F_{E,l}(r) = sqrt(2/(pi k)) (k r) j_l(k r) = sqrt(2 k / pi) r j_l(k r)

at electron mass 1 -- eMoScat's `sphBesselJEn` (`source/bessel.cpp:50`),
equivalently `sF_en` (`source/coulomb.cpp:75`) with charge 0 since N2 is
neutral. Energy normalization (`<F_E|F_E'> = delta(E-E')`) is what makes the
`sigma = 4 pi^3 |T|^2 / k^2` prefactor correct; getting the constant wrong
rescales every cross section.

`scipy.special.spherical_jn` takes REAL arguments only. That is not a
limitation here: `F` is only ever needed on the UNSCALED region, because a
channel projection that extends onto the exterior-complex-scaled tail is
meaningless and must be masked to zero there anyway (eMoScat zeroes it
explicitly -- `time_independent_model.cpp:149-151`). So `F` is evaluated on
real points and the tail is zeroed; no complex Bessel function is needed.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from scipy.special import spherical_jn

__all__ = ["riccati_bessel_en"]


def riccati_bessel_en(r: npt.NDArray[np.float64], k: float, l: int) -> npt.NDArray[np.float64]:
    """`F_{E,l}(r) = sqrt(2k/pi) r j_l(k r)`, energy-normalized at mass 1.

    `r` must be REAL (see module docstring); `k = sqrt(2E) > 0`.
    """
    if k <= 0.0:
        raise ValueError(f"k must be positive, got {k}")
    rr = np.asarray(r, dtype=np.float64)
    out: npt.NDArray[np.float64] = np.sqrt(2.0 * k / np.pi) * rr * spherical_jn(l, k * rr)
    return out
