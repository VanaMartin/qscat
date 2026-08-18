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

__all__ = [
    "free_hamiltonian",
    "incident_coefficients",
    "scattering_state",
    "scattering_state_minus",
]


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


def scattering_state_minus(
    h: npt.NDArray[np.complex128],
    grid: FemDvrEcsGrid,
    energy: float,
    ell: int,
) -> npt.NDArray[np.complex128]:
    """`phi^-`: the INCOMING-boundary partner of `scattering_state`.

    PRA 77 Eq. (34) gives `phi_k^- = (phi_k^+)^*` for a real discrete state in
    the radial case (Eq. 36 is the three-dimensional replacement and does not
    apply here). Conjugation is the physics here -- it is what reverses the
    boundary condition -- so this is the one place in `qscat.core.nrm` where
    `np.conjugate` is correct rather than a c-product violation.

    Under exterior complex scaling conjugation also conjugates the CONTOUR, so
    the identity is used only on the real region and the ECS tail is zeroed.
    That is sufficient for its consumer: Eq. (37)'s integrand carries
    `V_int`, a Gaussian in `r`, which has no support on the tail.
    `test_minus_state_is_purely_incoming_by_hankel_decomposition` gates the
    identity rather than assuming it, by checking that -- beyond the
    potential's support -- the scattered part carries no outgoing component.

    Parameters
    ----------
    h : ndarray
        The electronic Hamiltonian matrix to scatter off -- same argument
        `scattering_state` takes.
    grid : FemDvrEcsGrid
        The electronic radial grid.
    energy : float
        The REAL electron energy `E = k^2/2` (hartree); must be positive.
    ell : int
        Partial wave.

    Returns
    -------
    ndarray
        The incoming-boundary scattering solution as DVR coefficients on the
        full grid, zero on the ECS tail.
    """
    plus = scattering_state(h, grid, energy, ell)
    out = np.conjugate(plus)
    out[grid.points.imag != 0.0] = 0.0
    return np.asarray(out, dtype=np.complex128)
