"""Energy-normalized Coulomb functions (the charge-z generalization of the
Riccati-Bessel/Hankel radial functions in `radial.py`).

For a scattering particle of mass `m` in a Coulomb field of charge `z` (the
Sommerfeld parameter `eta = m z / k`, momentum `k`), the energy-normalized
regular / irregular / outgoing radial solutions are

    F_en(x) = sqrt(2 m/(pi k)) F_l(eta, k x)
    G_en(x) = sqrt(2 m/(pi k)) G_l(eta, k x)
    H1_en(x) = sqrt(2 m/(pi k)) (G_l + i F_l)(eta, k x)   [outgoing, H+ = G + iF]

with F_l/G_l the standard regular/irregular Coulomb functions (mpmath, which
accepts COMPLEX arguments -- needed for ECS-rotated x). At z=0 (eta=0),
F_l(0, rho) = rho j_l(rho), so `coulomb_f_en(., ., 0, m, l)` reduces to
`riccati_bessel_en` at mass m -- the differential-test hook. eMoScat's
`sH1` wrapper had a copy-paste bug (returned F, not G+iF); we define H+ = G + iF
correctly. (eMoScat coulomb.cpp / coulcc.f; the DR incident wave uses F_en.)
"""

from __future__ import annotations

import mpmath  # type: ignore[import-untyped]
import numpy as np
import numpy.typing as npt

__all__ = ["coulomb_f_en", "coulomb_g_en", "coulomb_h1_en"]


def _fg(
    x: npt.NDArray[np.complex128], k: float, z: float, m: float, l: int
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
    if k <= 0.0:
        raise ValueError(f"k must be positive, got {k}")
    if m <= 0.0:
        raise ValueError(f"m must be positive, got {m}")
    eta = m * z / k
    xs = np.asarray(x, dtype=np.complex128).ravel()
    pref = np.sqrt(2.0 * m / (np.pi * k))
    f = np.empty(xs.size, dtype=np.complex128)
    g = np.empty(xs.size, dtype=np.complex128)
    for i, xv in enumerate(xs):
        rho = mpmath.mpc(k * xv)
        f[i] = complex(mpmath.coulombf(l, eta, rho))
        g[i] = complex(mpmath.coulombg(l, eta, rho))
    return (pref * f).reshape(np.shape(x)), (pref * g).reshape(np.shape(x))


def coulomb_f_en(
    x: npt.ArrayLike, k: float, z: float, m: float, l: int
) -> npt.NDArray[np.complex128]:
    """`sqrt(2m/pi k) F_l(m z/k, k x)`; reduces to `riccati_bessel_en(x,k,l)` at z=0, m=1."""
    f, _ = _fg(np.asarray(x, dtype=np.complex128), k, z, m, l)
    return np.asarray(f, dtype=np.complex128)


def coulomb_g_en(
    x: npt.ArrayLike, k: float, z: float, m: float, l: int
) -> npt.NDArray[np.complex128]:
    """Energy-normalized irregular Coulomb function `sqrt(2m/pi k) G_l`."""
    _, g = _fg(np.asarray(x, dtype=np.complex128), k, z, m, l)
    return np.asarray(g, dtype=np.complex128)


def coulomb_h1_en(
    x: npt.ArrayLike, k: float, z: float, m: float, l: int
) -> npt.NDArray[np.complex128]:
    """Energy-normalized OUTGOING Coulomb function `sqrt(2m/pi k)(G_l + i F_l)`."""
    f, g = _fg(np.asarray(x, dtype=np.complex128), k, z, m, l)
    return np.asarray(g + 1j * f, dtype=np.complex128)
