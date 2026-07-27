"""Energy-normalized free radial functions for electron-scattering channels.

The entrance/exit channel of a scattering transition is a free electron of
momentum `k` in partial wave `l`. The two energy-normalized free radial
solutions of the free radial Schroedinger equation (electron mass 1) are the
REGULAR Riccati-Bessel function

    F_{E,l}(r) = sqrt(2/(pi k)) (k r) j_l(k r) = sqrt(2 k / pi) r j_l(k r)

(`riccati_bessel_en`, eMoScat's `sphBesselJEn` -- `source/bessel.cpp:50`,
equivalently `sF_en` -- `source/coulomb.cpp:75` with charge 0) and its
OUTGOING Riccati-Hankel sibling

    F^{(1)}_{E,l}(r) = sqrt(2 k / pi) r h_l^{(1)}(k r),  h_l^{(1)} = j_l + i y_l

(`riccati_hankel_en`, eMoScat's `sphHankel1En` -- `TestFunction2d.cpp:207`).
Energy normalization (`<F_E|F_E'> = delta(E-E')`) is what makes the
`sigma = 4 pi^3 |T|^2 / k^2` prefactor correct; getting the constant wrong
rescales every cross section.

`scipy.special.spherical_jn`/`spherical_yn` take REAL arguments only. That is
not a limitation here: these functions are only ever needed on the UNSCALED
region, because a channel projection that extends onto the exterior-complex-
scaled tail is meaningless and must be masked to zero there anyway (eMoScat
zeroes it explicitly -- `time_independent_model.cpp:149-151`). So both
functions are evaluated on real points and the tail is zeroed; no complex
Bessel function is needed.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from scipy.special import spherical_jn, spherical_yn

__all__ = ["riccati_bessel_en", "riccati_hankel_en", "riccati_bessel_en_mass"]


def riccati_bessel_en(r: npt.NDArray[np.float64], k: float, l: int) -> npt.NDArray[np.float64]:
    """`F_{E,l}(r) = sqrt(2k/pi) r j_l(k r)`, energy-normalized at mass 1.

    `r` must be REAL (see module docstring); `k = sqrt(2E) > 0`.
    """
    if k <= 0.0:
        raise ValueError(f"k must be positive, got {k}")
    rr = np.asarray(r, dtype=np.float64)
    out: npt.NDArray[np.float64] = np.sqrt(2.0 * k / np.pi) * rr * spherical_jn(l, k * rr)
    return out


def riccati_hankel_en(r: npt.NDArray[np.float64], k: float, l: int) -> npt.NDArray[np.complex128]:
    """`F^{(1)}_{E,l}(r) = sqrt(2k/pi) r h_l^{(1)}(k r)`, energy-normalized, mass 1.

    `h_l^{(1)} = j_l + i y_l` is the OUTGOING spherical Hankel function, so
    `Re(F^{(1)}) == riccati_bessel_en` and `Im(F^{(1)}) = sqrt(2k/pi) r y_l(k r)`.
    `r` must be REAL (see module docstring); `k = sqrt(2E) > 0`.
    """
    if k <= 0.0:
        raise ValueError(f"k must be positive, got {k}")
    rr = np.asarray(r, dtype=np.float64)
    h1_l = spherical_jn(l, k * rr) + 1j * spherical_yn(l, k * rr)
    out: npt.NDArray[np.complex128] = np.sqrt(2.0 * k / np.pi) * rr.astype(np.complex128) * h1_l
    return out


def riccati_bessel_en_mass(
    r: npt.NDArray[np.float64], k: float, l: int, mu: float
) -> npt.NDArray[np.float64]:
    """`F_{E,l}(r) = sqrt(2 mu k / pi) r j_l(k r)`, energy-normalized at mass `mu`.

    The mass-`mu` generalization of `riccati_bessel_en` (which is the `mu=1`
    case): the energy-normalized (`<F_E|F_E'> = delta(E-E')`) regular radial
    solution for a particle of reduced mass `mu` and momentum `k = sqrt(2 mu E)`.
    Used for the OUTGOING NUCLEAR dissociation wave in the DA/DR exit channel
    (eMoScat `bessel::s_jEn(R, K, mu, l)`). `r` must be REAL (channel
    projections are masked to the unscaled region); `k>0`, `mu>0`.
    """
    if k <= 0.0:
        raise ValueError(f"k must be positive, got {k}")
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu}")
    rr = np.asarray(r, dtype=np.float64)
    out: npt.NDArray[np.float64] = (
        np.sqrt(2.0 * mu * k / np.pi) * rr * spherical_jn(l, k * rr)
    )
    return out
