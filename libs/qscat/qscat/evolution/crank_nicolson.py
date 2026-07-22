"""Crank-Nicolson time propagator for the time-dependent Schrodinger equation.

Advances psi_{n+1} from psi_n via the Cayley form

    (I + i H dt/2) psi_{n+1} = (I - i H dt/2) psi_n

which is exact to O(dt^3) per step, unconditionally stable, unitary when H
is Hermitian, and norm-decaying when H has a negative-imaginary-part
(absorbing/optical-potential) component. H may be a general complex,
non-Hermitian matrix. The left-hand-side matrix A = I + i H dt/2 is
LU-factored once per call to `make_cn_stepper`, and each returned
`stepper(psi)` reuses that factorization via `scipy.linalg.lu_solve`.

Promoted from `projects/n2_td_cross_section/propagator.py` (Task 1 of the
N2 time-dependent cross-section sub-project) as the general primitive: it
has no dependency on FEM-DVR-ECS or any N2-specific structure, so it lives
here as a reusable `qscat.evolution` propagator. See
`docs/physics/n2-td-cross-section.md` for the physics application (the
N2 resonance's non-Hermitian H_res = T_nuc + diag(V_d - i*Gamma/2)).
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np
import numpy.typing as npt
from scipy.linalg import lu_factor, lu_solve

__all__ = ["make_cn_stepper"]


def make_cn_stepper(
    H: npt.NDArray[np.complexfloating[Any, Any]], dt: float
) -> Callable[[npt.NDArray[np.complexfloating[Any, Any]]], npt.NDArray[np.complex128]]:
    """Build a Crank-Nicolson stepper for Hamiltonian ``H`` and time step ``dt``.

    Parameters
    ----------
    H : ndarray, shape (n, n)
        Hamiltonian matrix (complex, possibly non-Hermitian).
    dt : float
        Time step.

    Returns
    -------
    stepper : Callable[[ndarray], ndarray]
        Advances a state vector ``psi`` by one time step ``dt``.
    """
    n = H.shape[0]
    ident = np.eye(n, dtype=complex)
    A = ident + 0.5j * dt * H  # (I + i H dt/2)
    B = ident - 0.5j * dt * H  # (I - i H dt/2)
    lu = lu_factor(A)

    def stepper(
        psi: npt.NDArray[np.complexfloating[Any, Any]],
    ) -> npt.NDArray[np.complex128]:
        result: npt.NDArray[np.complex128] = lu_solve(lu, B @ psi)
        return result

    return stepper
