"""Energy-normalized free radial functions for electron-scattering channels.

The entrance/exit channel of a scattering transition is a free particle of
momentum `k` and reduced mass `mu` in partial wave `l`. The two
energy-normalized free radial solutions of the free radial Schroedinger
equation are the REGULAR Riccati-Bessel function

    F_{E,l}(r) = sqrt(2 mu/(pi k)) (k r) j_l(k r) = sqrt(2 mu k / pi) r j_l(k r)

(`riccati_bessel_en`, eMoScat's `sphBesselJEn` -- `source/bessel.cpp:50`,
equivalently `sF_en` -- `source/coulomb.cpp:75` with charge 0), and its
OUTGOING Riccati-Hankel sibling

    F^{(1)}_{E,l}(r) = sqrt(2 mu k / pi) r h_l^{(1)}(k r),  h_l^{(1)} = j_l + i y_l

(`riccati_hankel_en`, eMoScat's `sphHankel1En` -- `TestFunction2d.cpp:207`).
`mu` enters ONLY the normalization prefactor in both formulas above, never
the momentum argument `k r`; `mu=1.0` (the default on both functions) is the
electron case, bit-for-bit the historical mass-1 functions, while `mu != 1`
is the nuclear case, used for the OUTGOING NUCLEAR dissociation wave in the
DA/DR exit channel.
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

__all__ = [
    "riccati_bessel_en",
    "riccati_bessel_en_mass",
    "riccati_hankel_en",
    "riccati_hankel_en_mass",
]


def riccati_bessel_en(
    r: npt.NDArray[np.float64], k: float, l: int, mu: float = 1.0
) -> npt.NDArray[np.float64]:
    """`F_{E,l}(r) = sqrt(2 mu k / pi) r j_l(k r)`, energy-normalized at mass `mu`.

    The energy-normalized (`<F_E|F_E'> = delta(E-E')`) REGULAR radial
    solution for a particle of reduced mass `mu` and momentum
    `k = sqrt(2 mu E)`. `mu=1.0` (the default) is the electron case --
    bit-for-bit the historical mass-1 function, since `2.0*1.0 == 2.0`
    exactly and `mu` enters ONLY the normalization prefactor, never the
    momentum argument `k r`. `mu != 1` is the nuclear case, used for the
    OUTGOING NUCLEAR dissociation wave in the DA/DR exit channel (eMoScat
    `bessel::s_jEn(R, K, mu, l)`). `r` must be REAL (see module docstring);
    `k > 0`, `mu > 0`.
    """
    if k <= 0.0:
        raise ValueError(f"k must be positive, got {k}")
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu}")
    rr = np.asarray(r, dtype=np.float64)
    out: npt.NDArray[np.float64] = np.sqrt(2.0 * mu * k / np.pi) * rr * spherical_jn(l, k * rr)
    return out


def riccati_hankel_en(
    r: npt.NDArray[np.float64], k: float, l: int, mu: float = 1.0
) -> npt.NDArray[np.complex128]:
    """`F^{(1)}_{E,l}(r) = sqrt(2 mu k / pi) r h_l^{(1)}(k r)`, energy-normalized
    at mass `mu`, `h_l^{(1)} = j_l + i y_l` the OUTGOING spherical Hankel
    function.

    The outgoing sibling of `riccati_bessel_en`, with the same mass
    convention: `mu` enters only the `sqrt(mu)` normalization prefactor
    (`2.0*1.0 == 2.0` exactly, so `mu=1.0` is bit-for-bit the historical
    mass-1 function), never the momentum argument `k r`. So
    `Re(F^{(1)}) == riccati_bessel_en` and `Im(F^{(1)}) = sqrt(2 mu k/pi)
    r y_l(k r)` at the same `(r, k, l, mu)`. `mu != 1` drives the OUTGOING
    NUCLEAR dissociation wave in the DA/DR flux (Wronskian) extractor
    (eMoScat `bessel::sphHankel1En(R, K, mu, l)` -- see
    `qscat.core.td_extractors.Flux`). `r` must be REAL (module docstring);
    `k > 0`, `mu > 0`.
    """
    if k <= 0.0:
        raise ValueError(f"k must be positive, got {k}")
    if mu <= 0.0:
        raise ValueError(f"mu must be positive, got {mu}")
    rr = np.asarray(r, dtype=np.float64)
    h1_l = spherical_jn(l, k * rr) + 1j * spherical_yn(l, k * rr)
    out: npt.NDArray[np.complex128] = (
        np.sqrt(2.0 * mu * k / np.pi) * rr.astype(np.complex128) * h1_l
    )
    return out


def riccati_bessel_en_mass(
    r: npt.NDArray[np.float64], k: float, l: int, mu: float
) -> npt.NDArray[np.float64]:
    """Deprecated alias for `riccati_bessel_en(r, k, l, mu)`: the mass
    generalization lives on the base name now. Kept so existing imports
    keep working; new code should call `riccati_bessel_en` directly."""
    return riccati_bessel_en(r, k, l, mu)


def riccati_hankel_en_mass(
    r: npt.NDArray[np.float64], k: float, l: int, mu: float
) -> npt.NDArray[np.complex128]:
    """Deprecated alias for `riccati_hankel_en(r, k, l, mu)`: the mass
    generalization lives on the base name now. Kept so existing imports
    keep working; new code should call `riccati_hankel_en` directly."""
    return riccati_hankel_en(r, k, l, mu)
