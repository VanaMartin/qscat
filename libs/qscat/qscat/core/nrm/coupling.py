"""The discrete-continuum coupling `V_dk+(R)` at a real electron energy.

Eq. (21): `V_dk+(R) = <phi_d| H_el |phi_k+>`, with `phi_k+` the P-space
continuum of Eq. (18) at the REAL energy `k^2/2`, obtained under exterior
complex scaling (p. 012710-6).

This is not an optional refinement of the discretized couplings `V_dn` -- it is
the right-hand side of the nuclear equation (Eq. 52). The paper is explicit
that these elements "are at specific real electron energies and hence must be
evaluated directly using Eq. (21)".

Eq. (68), `Gamma(E,R) = 2 pi |V_dk+(R)|^2`, connects this to the width
`qscat.core.lcp` computes from the ECS resonance pole; the two are independent
routes to the same physical quantity and `test_nrm_coupling.py` gates on their
agreement.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from qscat.dvr import FemDvrEcsGrid

from .discrete_state import DiscreteState, electronic_hamiltonian
from .scattering import scattering_state

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["gamma_from_coupling", "v_dk_plus"]


def v_dk_plus(
    grid: FemDvrEcsGrid,
    model: ResonanceModel,
    phi_d: DiscreteState,
    R_values: npt.NDArray[np.float64],
    energy: float,
) -> npt.NDArray[np.complex128]:
    """`V_dk+(R)` (Eq. 21) at real electron `energy`, one entry per `R`.

    Parameters
    ----------
    grid : FemDvrEcsGrid
        The electronic radial grid (exterior-complex-scaled).
    model : ResonanceModel
        Supplies `surface`, `v0` and `ell`.
    phi_d : DiscreteState
        The discrete state; consumed only through `phi_d.phi_d(R)`.
    R_values : ndarray
        Nuclear coordinates (real) at which to evaluate the coupling.
    energy : float
        The real electron energy `E = k^2/2` (hartree), positive.

    Returns
    -------
    ndarray
        Complex `V_dk+(R)`, shape `(R_values.size,)`.
    """
    if energy <= 0.0:
        raise ValueError(f"energy must be positive, got {energy}")
    R = np.asarray(R_values, dtype=np.float64)
    out = np.empty(R.size, dtype=np.complex128)
    ident = np.eye(grid.n, dtype=np.complex128)
    for j in range(R.size):
        d = phi_d.phi_d(float(R[j]))
        h_el = electronic_hamiltonian(grid, model, float(R[j]))
        p = ident - np.outer(d, d)  # Eq. (57)-(58), bilinear -- no conjugation
        php = p @ h_el @ p
        phi_k = scattering_state(php, grid, energy, model.ell)
        # Eq. (21) under the c-product: coefficients pair without conjugation.
        out[j] = d @ (h_el @ phi_k)
    return out


def gamma_from_coupling(
    v_dk: npt.NDArray[np.complex128],
) -> npt.NDArray[np.float64]:
    """`Gamma(E,R) = 2 pi |V_dk+(R)|^2` -- Eq. (68), p. 012710-8."""
    out: npt.NDArray[np.float64] = 2.0 * np.pi * np.abs(v_dk) ** 2
    return out
