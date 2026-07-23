"""Kronecker sum of operators over an arbitrary number of dimensions.

The construction behind every separable-kinetic Hamiltonian on a tensor
product of grids:

    kron_sum([A_0, ..., A_{D-1}]) = sum_d  I x ... x A_d x ... x I

Pure linear algebra -- this module knows nothing about grids, potentials or
physics, and accepts ANY square sparse matrices. That is deliberate: a future
angular-DVR, finite-difference or B-spline dimension composes with FEM-DVR-ECS
dimensions at no extra cost.

Index convention: numpy-native C order, i.e. the LAST axis is fastest, so the
result acts on `psi.ravel()` for `psi` of shape `(n_0, ..., n_{D-1})`. Note
eMoScat uses the OPPOSITE convention (first coordinate fastest,
`idx = i_r + i_R*N_r`, `FemDvrEcsGrid2d.cpp:169`); the two are physically
identical and differ only in basis ordering.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import scipy.sparse as sp

__all__ = ["kron_sum"]


def kron_sum(ops: Sequence[sp.spmatrix]) -> sp.csr_matrix:
    """Assemble `sum_d I x ... x ops[d] x ... x I` as a CSR matrix.

    Each `ops[d]` must be square. The result is square with dimension
    `prod(n_d)`. `D == 1` returns `ops[0]` unchanged (as CSR).
    """
    mats = list(ops)
    if not mats:
        raise ValueError("kron_sum requires at least one operator")
    for d, m in enumerate(mats):
        if m.shape[0] != m.shape[1]:
            raise ValueError(f"operator {d} is not square: shape {m.shape}")

    sizes = [int(m.shape[0]) for m in mats]

    # D == 1 falls out of the same loop below (left = right = 1, so the
    # single term is I_1 (x) m (x) I_1 == m): no special case needed, and it
    # sidesteps scipy-stubs' incomplete `spmatrix` mixin (no `.tocsr()` in
    # the stubs) that a literal `return ops[0]` conversion would hit.
    total: sp.csr_matrix | None = None
    for d, m in enumerate(mats):
        left = int(np.prod(sizes[:d])) if d else 1
        right = int(np.prod(sizes[d + 1 :])) if d < len(mats) - 1 else 1
        term = sp.kron(
            sp.identity(left, format="csr", dtype=complex),
            sp.kron(m, sp.identity(right, format="csr", dtype=complex), format="csr"),
            format="csr",
        )
        total = term if total is None else total + term

    assert total is not None  # len(mats) >= 1 guarantees at least one iteration
    return sp.csr_matrix(total)
