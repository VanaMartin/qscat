"""Crank-Nicolson time propagator for the time-dependent Schrodinger equation.

Advances psi_{n+1} from psi_n via the Cayley form

    (I + i H dt/2) psi_{n+1} = (I - i H dt/2) psi_n

which is exact to O(dt^3) per step, unconditionally stable, unitary when H
is Hermitian, and norm-decaying when H has a negative-imaginary-part
(absorbing/optical-potential) component. H may be a general complex,
non-Hermitian matrix. The left-hand-side matrix A = I + i H dt/2 is
LU-factored once per call to `make_cn_stepper`, and each returned
`stepper(psi)` reuses that factorization via `scipy.linalg.lu_solve`.

Promoted from `projects/n2_td_cross_section/propagator.py` as the general
primitive: it has no dependency on FEM-DVR-ECS or any N2-specific structure,
so it lives here as a reusable `qscat.evolution` propagator. See
`docs/physics/n2-td-cross-section.md` for the physics application (the
N2 resonance's non-Hermitian H_res = T_nuc + diag(V_d - i*Gamma/2)).
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from scipy.linalg import lu_factor, lu_solve

from qscat.linalg import SparseLU

__all__ = ["make_cn_stepper", "make_sparse_cn_stepper"]


def make_cn_stepper(
    H: npt.NDArray[np.complexfloating[Any, Any]], dt: float
) -> Callable[[npt.NDArray[np.complexfloating[Any, Any]]], npt.NDArray[np.complex128]]:
    """Build a Crank-Nicolson stepper for Hamiltonian ``H`` and time step ``dt``.

    Parameters
    ----------
    H : ndarray
        ``(n, n)``. Hamiltonian matrix (complex, possibly non-Hermitian).
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


def make_sparse_cn_stepper(
    H: sp.spmatrix, dt: float
) -> Callable[[npt.NDArray[np.complexfloating[Any, Any]]], npt.NDArray[np.complex128]]:
    """Sparse Crank-Nicolson stepper -- the sparse sibling of `make_cn_stepper`.

    Same Cayley form `(I + i H dt/2) psi_{n+1} = (I - i H dt/2) psi_n`, but
    `A = I + i H dt/2` is factored once with `qscat.linalg.SparseLU` and each
    step is a single sparse back-substitution. For the ~1e4-1e5-dimension
    sparse Hamiltonians this targets, dense factorization is infeasible; the
    dense `make_cn_stepper` is retained as this function's differential oracle.

    `H` must be square and sparse. Complex symmetric (ECS) `H` is fine -- no
    Hermiticity is assumed.
    """
    n = H.shape[0]
    # Explicit conversion to a concrete complex128 csc_matrix, mirroring
    # SparseLU's own internal conversion: scipy-stubs' generic `spmatrix`
    # mixin lacks the arithmetic overloads needed for `ident + ... * H`
    # below, so `H` (whatever concrete sparse subtype/dtype it arrives as)
    # is normalized to a type mypy can resolve `+`/`-` against.
    H_csc: sp.csc_matrix[np.complex128] = sp.csc_matrix(H, dtype=np.complex128)
    ident: sp.csc_matrix[np.complex128] = sp.identity(n, format="csc", dtype=np.complex128)
    A = (ident + 0.5j * dt * H_csc).tocsc()
    B = (ident - 0.5j * dt * H_csc).tocsr()
    lu = SparseLU(A)

    def stepper(
        psi: npt.NDArray[np.complexfloating[Any, Any]],
    ) -> npt.NDArray[np.complex128]:
        result: npt.NDArray[np.complex128] = lu.solve(B @ psi)
        return result

    return stepper
