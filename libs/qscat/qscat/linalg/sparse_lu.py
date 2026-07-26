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
`ordering` + `fill_factor` + `memory_bytes()` exist so it can be MEASURED
rather than guessed. Neither can affect correctness -- only speed and memory.

Reusing one factorization across right-hand sides is the point: in a scattering
calculation every final channel at a given energy shares the same matrix.

FILL-IN DIAGNOSTICS COST MEMORY DIFFERENTLY -- READ BEFORE CALLING
`memory_bytes()` ON A PRODUCTION MATRIX:

  `fill_factor` is free at any scale: SuperLU's own `nnz` attribute reports
  the total L+U nonzero count directly from its internal factorization
  structure, with NO conversion to sparse arrays and NO extra memory. Call it
  as often as you like.

  `memory_bytes()` is NOT free, which is why it is a method, not a property.
  Computing it forces scipy to materialize `self._lu.L` and `self._lu.U` as
  full CSC arrays (data + indices + indptr) -- and, measured directly on this
  class (see `docs/physics/nd-tensor-hamiltonian.md`), scipy's `SuperLU`
  object then CACHES those arrays internally for its own lifetime: a second
  access costs no extra memory (proof the first access is cached, not
  rebuilt), and the memory is NOT released by deleting your own references to
  the returned arrays -- only deleting this whole `SparseLU` object (and
  therefore the factorization itself) frees it. At production scale (L+U
  nnz = 3.05e8, complex128 data + int32 indices) that cache is on the order
  of **+6 GB on top of the already-documented 13.6 GB peak** -- enough to OOM
  a 32 GB laptop that would otherwise finish the factorization. Measure
  `memory_bytes()` on a reduced grid to characterize the scaling, NOT on the
  production matrix itself; use the cheap `fill_factor` (and the nnz-based
  estimate `fill_factor * A.nnz * 16` bytes for data alone) to reason about
  the production case without ever paying this cost there.
"""

from __future__ import annotations

from typing import Literal, cast

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
import scipy.sparse.linalg as spla

__all__ = ["SparseLU"]

_Ordering = Literal["NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD"]
_Backend = Literal["auto", "scipy", "mumps"]


class _ScipyBackend:
    """The original `scipy.sparse.linalg.splu` path, unchanged.

    Holds exactly the factorization object and semantics `SparseLU` used
    before backend dispatch existed: `fill_factor` reads SuperLU's own `nnz`
    with no materialization, `memory_bytes()` materializes (and permanently
    caches, on the `SuperLU` object) the `L`/`U` CSC factors. See the module
    docstring for the memory caveat.
    """

    name = "scipy"

    def __init__(self, csc: sp.csc_matrix[np.complex128], ordering: _Ordering) -> None:
        self._ordering = ordering
        self._lu: spla.SuperLU[np.complex128] = spla.splu(csc, permc_spec=ordering)

    @property
    def ordering_used(self) -> str:
        return self._ordering

    def fill_factor(self, nnz: int) -> float:
        return float(self._lu.nnz) / float(nnz)

    def memory_bytes(self) -> int:
        total = 0
        for factor in (self._lu.L, self._lu.U):
            fcsc = factor.tocsc()
            total += fcsc.data.nbytes + fcsc.indices.nbytes + fcsc.indptr.nbytes
        return int(total)

    def solve(self, rhs: npt.NDArray[np.complex128]) -> npt.NDArray[np.complex128]:
        result = self._lu.solve(rhs)
        # mypy note: an inline `out: npt.NDArray[...] = self._lu.solve(...)` annotation here
        # pushes an expected-type context into SuperLU.solve's overload resolution that picks
        # the wrong (float64) overload despite a complex128 argument -- a scipy-stubs/mypy
        # interaction, not a real type error. `cast` sidesteps it; the dtype is guaranteed by
        # the caller's explicit `.astype(np.complex128, ...)` before this is called.
        return cast(npt.NDArray[np.complex128], result)


class SparseLU:
    """LU factorization of a square sparse matrix, reusable across solves.

    `ordering` is scipy's `permc_spec`: one of `"NATURAL"`, `"MMD_ATA"`,
    `"MMD_AT_PLUS_A"`, `"COLAMD"` (the default). For a structurally symmetric
    pattern -- which a Kronecker-sum Hamiltonian has -- `"MMD_AT_PLUS_A"` is
    often the better choice; measure with `fill_factor` before assuming.

    A real-valued `A` is silently promoted: the internal CSC conversion always
    uses `dtype=np.complex128`, so values are preserved but memory doubles.

    `fill_factor` is cheap at any scale (reads SuperLU's own `nnz` count, no
    array materialization). `memory_bytes()` is NOT cheap -- it is a method,
    not a property, precisely so that its cost is opt-in rather than hidden
    behind attribute access -- and its cache is permanent for this object's
    lifetime. See the module docstring before calling it on a production-size
    matrix.

    `backend` selects the factorization engine: `"scipy"` is the SuperLU path
    above (the only one implemented so far); `"mumps"` is reserved for a
    complex-symmetric MUMPS factorization (not yet implemented -- provisioned
    in Docker but not wired in here); `"auto"` (the default) prefers MUMPS
    when available and falls back to scipy otherwise -- today that fallback
    is unconditional, since no MUMPS path exists yet, so `"auto"` and
    `"scipy"` are currently identical in every observable way. `backend_used`
    reports which one actually ran.

    `symmetric`, if left `None`, is auto-detected as `A == A.T` (an O(nnz)
    sparse comparison: `(abs(A - A.T)).max() == 0`, cheap relative to the
    factorization itself). It is informational only on the scipy path --
    SuperLU does not exploit symmetry -- but is stored for the future MUMPS
    path (Task 3), which will use it to select the complex-symmetric MUMPS
    matrix type instead of the general unsymmetric one.
    """

    def __init__(
        self,
        A: sp.spmatrix,
        *,
        ordering: _Ordering = "COLAMD",
        backend: _Backend = "auto",
        symmetric: bool | None = None,
    ) -> None:
        if A.shape[0] != A.shape[1]:
            raise ValueError(f"matrix must be square, got shape {A.shape}")
        csc: sp.csc_matrix[np.complex128] = sp.csc_matrix(A, dtype=np.complex128)
        self._shape: tuple[int, int] = (int(csc.shape[0]), int(csc.shape[1]))
        self._nnz: int = int(csc.nnz)
        self._ordering = ordering

        self._backend_used: str
        if backend == "mumps":
            # Task 3 seam: wire the real MUMPS complex-symmetric factorization
            # in here (using `self.symmetric` to pick the MUMPS matrix type).
            # Until then, an explicit request for MUMPS must fail loudly
            # rather than silently falling back to scipy -- and before doing
            # any (wasted) symmetry detection on this error path.
            raise RuntimeError(
                "MUMPS backend requested but not available "
                "(qscat[mumps] / system MUMPS missing)"
            )

        if symmetric is None:
            # O(nnz) sparse comparison -- cheap relative to the factorization
            # that follows, but not free, hence computed once and cached.
            symmetric = bool((abs(csc - csc.T)).max() == 0)
        self._symmetric = symmetric

        if backend == "scipy":
            self._impl: _ScipyBackend = _ScipyBackend(csc, ordering)
            self._backend_used = "scipy"
        else:  # backend == "auto"
            # Task 3 seam: try MUMPS first here (guarded by an availability
            # check), falling back to scipy on ImportError/unavailability.
            # For now there is no MUMPS path at all, so "auto" always falls
            # straight through to scipy.
            self._impl = _ScipyBackend(csc, ordering)
            self._backend_used = "scipy"

    @property
    def shape(self) -> tuple[int, int]:
        return self._shape

    @property
    def ordering(self) -> str:
        return self._ordering

    @property
    def symmetric(self) -> bool:
        """Whether `A` was treated as (complex-)symmetric `A == A.T`.

        Auto-detected from `A` when `symmetric=None` was passed (the default),
        or the explicit override. Informational on the scipy path; the MUMPS
        path (Task 3) uses it to select the complex-symmetric (`SYM=2`) matrix
        type instead of the general unsymmetric one.
        """
        return self._symmetric

    @property
    def backend_used(self) -> str:
        """Which backend actually factorized `A`: `"scipy"` (only option today)."""
        return self._backend_used

    @property
    def ordering_used(self) -> str:
        """The ordering the active backend actually used.

        On the scipy path this is scipy's `permc_spec` (identical to
        `ordering`); the MUMPS path (Task 3) will report MUMPS's own chosen
        ordering here instead.
        """
        return self._impl.ordering_used

    @property
    def fill_factor(self) -> float:
        """`(L.nnz + U.nnz) / A.nnz` -- how much denser the factors are.

        Cheap at any scale: on the scipy backend, `self._lu.nnz` is SuperLU's
        own reported L+U nonzero count, read directly off the internal
        factorization -- this NEVER materializes `L` or `U` as arrays and
        costs no extra memory (measured delta < 0.1 MB on an N=6000 matrix
        with a x300 fill-in). Contrast `memory_bytes()`, which does
        materialize them and is priced accordingly -- see the module
        docstring.
        """
        return self._impl.fill_factor(self._nnz)

    def memory_bytes(self) -> int:
        """Bytes actually held by the L and U factors (data + index arrays).

        NOT CHEAP -- a method, not a property, because computing this forces
        scipy to materialize `self._lu.L` and `self._lu.U` as full CSC arrays,
        which `SuperLU` then caches for this object's lifetime (measured: a
        second call allocates no further memory, and the cache is not
        released by dropping your own references to the result -- only
        deleting this `SparseLU` instance does). Read the module docstring's
        production-scale estimate (+6 GB on the N2 2-D deck) before calling
        this on anything but a reduced/test-scale matrix.
        """
        return self._impl.memory_bytes()

    def solve(self, b: npt.NDArray[np.complex128]) -> npt.NDArray[np.complex128]:
        """Solve `A x = b` for one `(N,)` or several `(N, k)` right-hand sides."""
        rhs = np.asarray(b)
        if rhs.ndim == 0:
            raise ValueError(
                f"right-hand side must be at least 1-D (an (N,) vector or an "
                f"(N, k) block of right-hand sides), got a 0-d scalar with "
                f"shape {rhs.shape}"
            )
        if rhs.shape[0] != self._shape[0]:
            raise ValueError(
                f"right-hand side has leading dimension {rhs.shape[0]}, "
                f"expected {self._shape[0]}"
            )
        return self._impl.solve(rhs.astype(np.complex128, copy=False))
