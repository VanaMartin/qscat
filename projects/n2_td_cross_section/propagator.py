"""Crank-Nicolson time propagator for the time-dependent Schrodinger equation.

Advances psi_{n+1} from psi_n via the Cayley form

    (I + i H dt/2) psi_{n+1} = (I - i H dt/2) psi_n

which is exact to O(dt^3) per step, unconditionally stable, unitary when H
is Hermitian, and norm-decaying when H has a negative-imaginary-part
(absorbing/optical-potential) component. H may be a general complex,
non-Hermitian matrix. The left-hand-side matrix A = I + i H dt/2 is
LU-factored once per call to `make_cn_stepper`, and each returned
`stepper(psi)` reuses that factorization via `scipy.linalg.lu_solve`.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy.linalg import lu_factor, lu_solve


def make_cn_stepper(H: np.ndarray, dt: float) -> Callable[[np.ndarray], np.ndarray]:
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

    def stepper(psi: np.ndarray) -> np.ndarray:
        return lu_solve(lu, B @ psi)

    return stepper
