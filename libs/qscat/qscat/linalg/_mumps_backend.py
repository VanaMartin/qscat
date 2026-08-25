"""MUMPS complex-symmetric (``SYM=2``) sparse backend for `qscat.linalg.SparseLU`.

Import-guarded on purpose: this module imports fine WITHOUT ``python-mumps`` /
system MUMPS installed, so it lints and type-checks on a MUMPS-less dev box (the
Mac). `mumps_available()` reports whether the real backend can actually run;
`SparseLU` gates on it so ``backend="auto"`` cleanly falls back to SuperLU when
MUMPS is absent. All MUMPS execution happens in the Docker ``test`` image where
the system MUMPS libraries and the ``qscat[mumps]`` extra are installed.

Recipe (verified against SuperLU to machine precision, rel err 7.25e-16 on
an N=400 complex-symmetric system)::

    ctx = mumps.Context()
    ctx.set_matrix(sp.triu(A).tocsc(), symmetric=True)   # SYM=2: UPPER TRIANGLE
    ctx.analyze()                                         # ordering='auto'
    ctx.factor()                                          # factor once
    x = ctx.solve(b)                                      # solve many (reuse ctx)

THE CORRECTNESS TRAP: ``SYM=2`` reads ONLY the upper triangle. For a symmetric
``A`` we supply ``sp.triu(A)``; supplying the full matrix would double-count the
off-diagonals. A genuinely non-symmetric ``A`` must use the ``SYM=0`` (full
matrix, ``symmetric=False``) path instead.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.exceptions import BackendError

# Runtime import guard: typed as ``Any`` so this module type-checks with MUMPS
# absent (no stubs, no module). ``mumps_available()`` is the single gate every
# caller uses; the class body never runs unless it returned True.
try:
    import mumps as _mumps  # type: ignore[import-not-found]

    _MUMPS: Any = _mumps
except ImportError:  # pragma: no cover - the MUMPS-less path (e.g. the Mac)
    _MUMPS = None


def mumps_available() -> bool:
    """True iff ``python-mumps`` imported; NEVER raises (returns False instead).

    Callers use this to decide whether the MUMPS backend can run at all, so a
    missing MUMPS degrades to SuperLU rather than erroring.
    """
    return _MUMPS is not None


# A stored sparsity pattern: ``(indices, indptr)`` of a canonicalized CSC matrix.
_Pattern = tuple[npt.NDArray[Any], npt.NDArray[Any]]


def _pattern_of(csc: sp.csc_matrix[np.complex128]) -> _Pattern:
    """The canonical ``(indices, indptr)`` sparsity pattern of ``csc`` (copied).

    Canonicalizes first (sorted column indices, summed duplicates) so the
    stored pattern is well-defined and comparable to any later matrix put
    through the same canonicalization in `_check_pattern`.
    """
    canon = csc.sorted_indices()
    canon.sum_duplicates()
    return (canon.indices.copy(), canon.indptr.copy())


def _check_pattern(pattern: _Pattern, csc: sp.csc_matrix[np.complex128]) -> None:
    """Raise ``ValueError`` if ``csc``'s sparsity pattern differs from ``pattern``.

    ``pattern`` is the stored ``(indices, indptr)`` of the matrix that was
    symbolically analyzed. Reusing that analysis (MUMPS ``reuse_analysis=True``)
    is valid ONLY for an identical nonzero structure, so a mismatch must raise
    rather than silently reuse a wrong ordering and produce garbage. ``csc`` is
    canonicalized (sorted indices, summed duplicates) before comparison so the
    check is well-defined regardless of the caller's index ordering.
    """
    indices, indptr = pattern
    canon = csc.sorted_indices()
    canon.sum_duplicates()
    if canon.indptr.shape != indptr.shape or canon.indices.shape != indices.shape:
        raise ValueError(
            "refactor pattern mismatch: nonzero-structure shape "
            f"(indptr {canon.indptr.shape}, indices {canon.indices.shape}) "
            f"differs from the analyzed pattern "
            f"(indptr {indptr.shape}, indices {indices.shape})"
        )
    if not (np.array_equal(canon.indptr, indptr) and np.array_equal(canon.indices, indices)):
        raise ValueError(
            "refactor pattern mismatch: nonzero structure differs from the "
            "analyzed matrix (reuse_analysis requires an identical pattern)"
        )


# python-mumps exposes INFOG 1-based-aligned: infog[k] == INFOG(k), with index 0
# unused (its `mumps_int_array` mirrors the Fortran 1-based indexing directly, it
# is NOT the usual 0-based shift). Verified against the live array in the
# container: for a symmetric factorization the analysis/factorization pairs
# INFOG(3)==INFOG(9), INFOG(20)==INFOG(29) line up exactly under this convention.
_INFOG_ORDERING = 7  # ordering actually used during analysis (echo of ICNTL(7))
_INFOG_ENTRIES_IN_FACTORS = 29  # effective number of entries in the factors
_INFOG_MEMORY_MB = 22  # memory (MB) effectively used during factorization

# INFOG(7) ordering codes -> human names (MUMPS manual, ICNTL(7)/INFOG(7)).
# INFOG(7) reports the ordering MUMPS ACTUALLY used (0-6); it never reports 7
# ("auto" is an ICNTL(7) *request* value that MUMPS resolves to a concrete
# 0-6 code by analysis time), so there is deliberately no `7` key here -- an
# unexpected code falls through to `_ORDERING_NAMES.get(code, f"code-{code}")`.
_ORDERING_NAMES: dict[int, str] = {
    0: "amd",
    1: "user",
    2: "amf",
    3: "scotch",
    4: "pord",
    5: "metis",
    6: "qamd",
}


class _MumpsBackend:
    """MUMPS factorization of a square sparse matrix, reusable across solves.

    Mirrors `_ScipyBackend`'s shape: analyze+factor at construction,
    `solve(b)` on demand, and `fill_factor` / `memory_bytes()` /
    `ordering_used` sourced from MUMPS's INFOG diagnostics array. Uses
    ``symmetric`` (from `SparseLU`) to pick ``SYM=2`` (complex-symmetric, upper
    triangle only) versus ``SYM=0`` (general unsymmetric, full matrix).
    """

    name = "mumps"

    def __init__(self, csc: sp.csc_matrix[np.complex128], *, symmetric: bool) -> None:
        if _MUMPS is None:  # pragma: no cover - callers gate on mumps_available()
            raise BackendError("python-mumps is not installed")
        self._symmetric = symmetric
        self._ctx = _MUMPS.Context()
        # SYM=2 takes ONLY the upper triangle; SYM=0 takes the full matrix.
        a = sp.triu(csc).tocsc() if symmetric else csc
        # Store the SUPPLIED matrix's pattern (triu for SYM=2, full for SYM=0):
        # `refactor` must guard against a different structure before reusing
        # this analysis (see `_check_pattern`).
        self._pattern = _pattern_of(a)
        self._ctx.set_matrix(a, symmetric=symmetric)
        self._ctx.analyze()
        self._ctx.factor()

    def refactor(self, csc: sp.csc_matrix[np.complex128]) -> None:
        """Re-factorize ``csc`` reusing the existing symbolic analysis.

        ``csc`` MUST share the analyzed matrix's sparsity pattern (a diagonal
        shift ``E*I - H`` does). Re-supplies the same-pattern matrix and calls
        ``factor(reuse_analysis=True)``, so the persistent Context keeps its
        analysis and the (SCOTCH/METIS/...) ordering is NOT recomputed. Raises
        ``ValueError`` on a pattern mismatch, since reusing the analysis for a
        different structure would silently produce a wrong factorization.
        """
        a = sp.triu(csc).tocsc() if self._symmetric else csc
        _check_pattern(self._pattern, a)
        self._ctx.set_matrix(a, symmetric=self._symmetric)
        self._ctx.factor(reuse_analysis=True)

    def _infog(self, one_based: int) -> int:
        """Read INFOG(one_based), applying MUMPS's negative-million convention.

        The array is 1-based-aligned (``infog[k] == INFOG(k)``), so no index
        shift. MUMPS stores counts/memory that can overflow a 32-bit int as a
        negative value ``v`` meaning ``-v * 1e6`` (millions); fields that never
        overflow (e.g. the ordering code) are small and non-negative, so the
        convention is a no-op for them.
        """
        v = int(self._ctx.mumps_instance.infog[one_based])
        return -v * 1_000_000 if v < 0 else v

    @property
    def ordering_used(self) -> str:
        code = self._infog(_INFOG_ORDERING)
        return _ORDERING_NAMES.get(code, f"code-{code}")

    def fill_factor(self, nnz: int) -> float:
        return float(self._infog(_INFOG_ENTRIES_IN_FACTORS)) / float(nnz)

    def memory_bytes(self) -> int:
        return int(self._infog(_INFOG_MEMORY_MB)) * 1024 * 1024

    def solve(self, rhs: npt.NDArray[np.complex128]) -> npt.NDArray[np.complex128]:
        # python-mumps solves a single 1-D right-hand side per call; loop over
        # the columns of an (N, k) block to match `_ScipyBackend`'s contract.
        if rhs.ndim == 2:
            cols = [self._ctx.solve(rhs[:, j]) for j in range(rhs.shape[1])]
            out = np.column_stack(cols)
        else:
            out = self._ctx.solve(rhs)
        return np.asarray(out, dtype=np.complex128)
