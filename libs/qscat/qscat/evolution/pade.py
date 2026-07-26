"""Order-N diagonal-Pade time propagator for `d/dt psi = -i H psi`.

The Crank-Nicolson stepper (`crank_nicolson.py`) is the order-1 special case of
the diagonal [N,N] Pade approximant of `exp(-i H dt)`. This module generalizes
it to arbitrary order: writing the approximant in product/partial-fraction form
over the Pade denominator roots `{r_i}`,

    exp(-i H dt)  ~  prod_{i=1}^{order} (I - i H dt / r_i) (I + i H dt / r_i)^{-1}

each factor is one sparse mat-vec (numerator) and one back-substitution
(denominator, LU-factored once). The error is `O(dt^(2*order+1))` per step, so
order 3 (`O(dt^7)`) is dramatically more accurate than Crank-Nicolson's
`O(dt^3)` at the same `dt` -- the accuracy needed to make the time-dependent
cross section converge to the time-independent oracle (see
`docs/physics/n2-2d-td-cross-section.md`).

The roots are the roots of the diagonal-Pade denominator of `exp(z)`; the table
below matches eMoScat's `Pade_Roots` (`FemDvrEcs/FemDvrFunctions.cpp`), and
`order=1` (`r=2`) reproduces `make_sparse_cn_stepper` bit-for-bit. Sign
convention: for `exp(-i H dt)` we evaluate the `exp(z)` approximant at
`z = -i H dt`, giving the numerator `(I - i H dt / r_i)` / denominator
`(I + i H dt / r_i)` above; `r=2` -> `(I - iHdt/2)(I + iHdt/2)^{-1}`, the Cayley
form. The `{r_i}` are complex (one real root plus conjugate pairs), so each
factor is complex-symmetric when `H` is (the ECS case) -- `SparseLU` detects
and exploits that per factor.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.linalg import SparseLU

__all__ = ["make_pade_stepper", "pade_roots"]


# Roots `r_i` of the diagonal [N,N] Pade denominator, in the sign convention
# where `prod_i (1 - z/r_i)/(1 + z/r_i) = exp(-z)` -- i.e. `sum_i 1/r_i = 1/2`
# for every order, so the stepper below yields `exp(-i H dt)` uniformly. order 1
# (r=2) is Crank-Nicolson. These are eMoScat's `Pade_Roots` magnitudes with the
# real part taken POSITIVE (eMoScat stores orders >=2 negated, which would flip
# the propagation sign -- verified: sum 1/r_i = 1/2 and the product matches
# `exp(-z)` to O(z^(2N+1)), errors 1e-3/1e-6/1e-9/1e-13 for N=1..4).
_PADE_ROOTS: dict[int, list[complex]] = {
    1: [2.0 + 0.0j],
    2: [3.0 - np.sqrt(3.0) * 1j, 3.0 + np.sqrt(3.0) * 1j],
    3: [
        4.6443707092521712 + 0.0j,
        3.6778146453739144 - 3.5087619195674433j,
        3.6778146453739144 + 3.5087619195674433j,
    ],
    4: [
        5.7924212056407443 - 1.7344682578690075j,
        5.7924212056407443 + 1.7344682578690075j,
        4.2075787943592557 - 5.3148360837135054j,
        4.2075787943592557 + 5.3148360837135054j,
    ],
}


def pade_roots(order: int) -> npt.NDArray[np.complex128]:
    """The `order` diagonal-Pade denominator roots of `exp(z)` (see module docstring).

    `order=1` is Crank-Nicolson (single real root 2). Raises `ValueError` for
    an order not in the table (1-4).
    """
    if order not in _PADE_ROOTS:
        raise ValueError(f"pade order must be one of {sorted(_PADE_ROOTS)}, got {order}")
    return np.array(_PADE_ROOTS[order], dtype=np.complex128)


def make_pade_stepper(
    H: sp.spmatrix, dt: float, order: int = 1
) -> Callable[[npt.NDArray[np.complexfloating[Any, Any]]], npt.NDArray[np.complex128]]:
    """Build an order-`order` diagonal-Pade stepper for `d/dt psi = -i H psi`.

    `H` must be square and sparse (complex; ECS complex-symmetric is fine, no
    Hermiticity assumed). The `order` denominators `(I + i H dt / r_i)` are each
    factored once with `qscat.linalg.SparseLU`; every returned `stepper(psi)`
    applies the product `prod_i (I + iHdt/r_i)^{-1} (I - iHdt/r_i)` in one pass.
    `order=1` is exactly `make_sparse_cn_stepper`. The factors commute (all are
    rational functions of `H`), so their application order is immaterial.
    """
    n = H.shape[0]
    # Concrete complex128 csc, mirroring make_sparse_cn_stepper's conversion so
    # the `ident +/- ... * H` arithmetic resolves (scipy-stubs' generic
    # spmatrix mixin lacks the operator overloads).
    h_csc: sp.csc_matrix[np.complex128] = sp.csc_matrix(H, dtype=np.complex128)
    ident: sp.csc_matrix[np.complex128] = sp.identity(n, format="csc", dtype=np.complex128)

    factors: list[tuple[sp.csr_matrix[np.complex128], SparseLU]] = []
    for root in pade_roots(order):
        coeff = 1j * dt / root
        numerator = (ident - coeff * h_csc).tocsr()  # (I - i H dt / r)
        denominator = (ident + coeff * h_csc).tocsc()  # (I + i H dt / r)
        factors.append((numerator, SparseLU(denominator)))

    def stepper(
        psi: npt.NDArray[np.complexfloating[Any, Any]],
    ) -> npt.NDArray[np.complex128]:
        x: npt.NDArray[np.complex128] = np.asarray(psi, dtype=np.complex128)
        for numerator, lu in factors:
            rhs: npt.NDArray[np.complex128] = np.asarray(numerator @ x, dtype=np.complex128)
            x = lu.solve(rhs)
        return x

    return stepper
