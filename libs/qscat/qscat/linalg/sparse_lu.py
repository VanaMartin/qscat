"""Cached sparse LU factorization: factor once, solve many right-hand sides.

A thin, typed wrapper over `scipy.sparse.linalg.splu` that adds the two things
the bare function lacks for our use: an explicit CSC conversion (splu warns
otherwise), and fill-in / memory diagnostics.

Those diagnostics are not decoration. At the sizes this library targets, the
factorization -- not the solve -- is the whole cost, and fill-in decides whether
a problem fits in RAM at all. A measured spike on the production N2 2-D deck
(N = 143,380, nnz = 3,276,450) gave x93 fill-in, 3.05e8 nonzeros in L+U, and
13.6 GB peak RSS with the default COLAMD ordering, against a 440 ms
back-substitution. Choosing an ordering is therefore a real decision, and
`ordering` + `fill_factor` + `memory_bytes` exist so it can be MEASURED rather
than guessed. It cannot affect correctness -- only speed and memory.

Reusing one factorization across right-hand sides is the point: in a scattering
calculation every final channel at a given energy shares the same matrix.
"""

from __future__ import annotations

from typing import Literal, cast

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
import scipy.sparse.linalg as spla

__all__ = ["SparseLU"]

_Ordering = Literal["NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD"]


class SparseLU:
    """LU factorization of a square sparse matrix, reusable across solves.

    `ordering` is scipy's `permc_spec`: one of `"NATURAL"`, `"MMD_ATA"`,
    `"MMD_AT_PLUS_A"`, `"COLAMD"` (the default). For a structurally symmetric
    pattern -- which a Kronecker-sum Hamiltonian has -- `"MMD_AT_PLUS_A"` is
    often the better choice; measure with `fill_factor` before assuming.
    """

    def __init__(self, A: sp.spmatrix, *, ordering: _Ordering = "COLAMD") -> None:
        if A.shape[0] != A.shape[1]:
            raise ValueError(f"matrix must be square, got shape {A.shape}")
        csc: sp.csc_matrix[np.complex128] = sp.csc_matrix(A, dtype=np.complex128)
        self._shape: tuple[int, int] = (int(csc.shape[0]), int(csc.shape[1]))
        self._nnz: int = int(csc.nnz)
        self._ordering = ordering
        self._lu: spla.SuperLU[np.complex128] = spla.splu(csc, permc_spec=ordering)

    @property
    def shape(self) -> tuple[int, int]:
        return self._shape

    @property
    def ordering(self) -> str:
        return self._ordering

    @property
    def fill_factor(self) -> float:
        """`(L.nnz + U.nnz) / A.nnz` -- how much denser the factors are."""
        if self._nnz == 0:
            return 1.0
        return float(self._lu.L.nnz + self._lu.U.nnz) / float(self._nnz)

    @property
    def memory_bytes(self) -> int:
        """Bytes actually held by the L and U factors (data + index arrays)."""
        total = 0
        for factor in (self._lu.L, self._lu.U):
            fcsc = factor.tocsc()
            total += fcsc.data.nbytes + fcsc.indices.nbytes + fcsc.indptr.nbytes
        return int(total)

    def solve(self, b: npt.NDArray[np.complex128]) -> npt.NDArray[np.complex128]:
        """Solve `A x = b` for one `(N,)` or several `(N, k)` right-hand sides."""
        rhs = np.asarray(b)
        if rhs.shape[0] != self._shape[0]:
            raise ValueError(
                f"right-hand side has leading dimension {rhs.shape[0]}, "
                f"expected {self._shape[0]}"
            )
        result = self._lu.solve(rhs.astype(np.complex128, copy=False))
        # mypy note: an inline `out: npt.NDArray[...] = self._lu.solve(...)` annotation here
        # pushes an expected-type context into SuperLU.solve's overload resolution that picks
        # the wrong (float64) overload despite a complex128 argument -- a scipy-stubs/mypy
        # interaction, not a real type error. `cast` sidesteps it; the dtype is guaranteed by
        # the explicit `.astype(np.complex128, ...)` just above.
        return cast(npt.NDArray[np.complex128], result)
