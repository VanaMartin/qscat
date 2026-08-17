"""Real-energy fixed-nuclei electronic scattering, the NRM's shared primitive.

Both the physical discrete state (`discrete_state.PhysicalDiscreteState`, the
scattering function at `Re E_res(R)`) and the discrete-continuum coupling
(`coupling.v_dk_plus`, Eq. 18/21) need the scattering solution of a
1-D electronic Hamiltonian at a REAL energy. This module is that solve, in the
driven (scattered-wave) form `qscat.core.driven` uses in 2-D:

    (E I - h) phi_sc = (h - H_free) J_k,   phi+ = J_k + phi_sc

with `J_k` the energy-normalized Riccati-Bessel function masked to the
unscaled region. `h - H_free` is short-ranged for every `h` this package
passes in, so the right-hand side has finite support in the real region.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid, kinetic
from qscat.special import riccati_bessel_en

__all__ = ["free_hamiltonian", "incident_coefficients", "scattering_state"]


def free_hamiltonian(grid: FemDvrEcsGrid, ell: int) -> npt.NDArray[np.complex128]:
    """`H_free = T_r + ell(ell+1)/(2 r^2)` on the electronic grid (mass 1).

    The Hamiltonian whose regular energy-normalized solution is
    `qscat.special.riccati_bessel_en`; the reference against which
    `scattering_state`'s driven source term is formed.
    """
    cent = ell * (ell + 1) / (2.0 * grid.points**2)
    out: npt.NDArray[np.complex128] = kinetic(grid, 1.0) + np.diag(cent)
    return out


def incident_coefficients(grid: FemDvrEcsGrid, k: float, ell: int) -> npt.NDArray[np.complex128]:
    """DVR coefficients of `J_{E,ell}(r)`, masked to the unscaled region.

    `riccati_bessel_en` is a FUNCTION; the DVR coefficient is
    `c_j = J(r_j) sqrt(w_j)`. The ECS tail is zeroed -- the incident wave is
    not defined on the complex contour, and every consumer pairs it with a
    short-ranged operator that vanishes there.
    """
    if k <= 0.0:
        raise ValueError(f"k must be positive, got {k}")
    real = grid.points.imag == 0.0
    out = np.zeros(grid.n, dtype=np.complex128)
    vals = riccati_bessel_en(grid.points[real].real, k, ell)
    out[real] = vals * np.sqrt(grid.weights[real])
    return out


def scattering_state(
    h: npt.NDArray[np.complex128],
    grid: FemDvrEcsGrid,
    energy: float,
    ell: int,
) -> npt.NDArray[np.complex128]:
    """`phi+ = J_k + (E I - h)^-1 (h - H_free) J_k`, DVR coefficients.

    Parameters
    ----------
    h : ndarray
        The electronic Hamiltonian matrix to scatter off -- the full `H_el`
        (Eq. 17) for the physical discrete state, or `P H_el P` (Eq. 18) for
        the background continuum. Must differ from `free_hamiltonian` only by
        a short-ranged operator.
    grid : FemDvrEcsGrid
        The electronic radial grid.
    energy : float
        The REAL electron energy `E = k^2/2` (hartree); must be positive.
    ell : int
        Partial wave.

    Returns
    -------
    ndarray
        The scattering solution as DVR coefficients on the full grid.
    """
    if energy <= 0.0:
        raise ValueError(f"energy must be positive, got {energy}")
    k = float(np.sqrt(2.0 * energy))
    inc = incident_coefficients(grid, k, ell)
    h_free = free_hamiltonian(grid, ell)
    rhs = (h - h_free) @ inc
    a = energy * np.eye(grid.n, dtype=np.complex128) - h
    phi_sc = np.linalg.solve(a, rhs)
    out: npt.NDArray[np.complex128] = inc + phi_sc
    return out
