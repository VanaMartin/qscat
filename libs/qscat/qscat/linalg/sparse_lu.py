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

from collections.abc import Iterator
from contextlib import contextmanager
from contextvars import ContextVar
from typing import Literal, cast

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from qscat.exceptions import BackendError

from ._mumps_backend import _check_pattern, _MumpsBackend, _pattern_of, mumps_available

__all__ = ["Ordering", "SparseLU", "default_backend", "get_default_backend", "set_default_backend"]

# scipy splu's permc_spec -- the public name solver modules re-use.
Ordering = Literal["NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD"]
_Backend = Literal["auto", "scipy", "mumps"]

# Relative tolerance for the `symmetric=None` auto-detect: `A` is treated as
# (complex-)symmetric when `max|A - A.T| <= _SYM_RTOL * max|A|`. This is a
# SCALED tolerance, not exact equality, because the matrices this class exists
# to factor -- the N2 driven `(E_tot*I - H)` and Crank-Nicolson `(I + iH*dt/2)`
# decks -- are mathematically `A = A.T` but only symmetric to ROUND-OFF: they
# are assembled by Kronecker-sum reordering of float arrays, so `A - A.T` is
# not bit-zero but sits at the floating-point noise floor (measured on the real
# N2 working deck: `max|A - A.T| = 4.5e-13`, `max|A| = 1.3e4`, i.e. a RELATIVE
# asymmetry ~3.6e-17 -- essentially one ULP). Exact equality
# (`(abs(A - A.T)).max() == 0`) rejected every such matrix and silently forced
# the MUMPS path onto SYM=0 (general unsymmetric) instead of the SYM=2
# (complex-symmetric, single-triangle) mode that is the whole point of the
# backend. `1e-12` sits ~5 orders of magnitude above the real matrices' ~3.6e-17
# relative asymmetry (an enormous accept margin) yet ~12 orders below the O(1)
# relative asymmetry of a genuinely non-symmetric matrix (a decisive reject
# margin), so it cannot misclassify a truly asymmetric matrix as symmetric.
# That safety matters because SYM=2 treats the upper triangle as truth and
# reconstructs the lower from it; accepting a truly-asymmetric matrix would give
# a WRONG answer. The tight bound keeps that from happening. Callers can always
# override the auto-detect with an explicit `symmetric=True`/`False`.
_SYM_RTOL = 1e-12

# Context-local default that `backend="auto"` resolves against. Lets a caller
# force every internal `SparseLU(...)` -- e.g. the ones `ve_cross_section_2d`
# creates without exposing a `backend=` kwarg -- onto one engine, for
# differential backend-equivalence checks. An EXPLICIT `backend="scipy"` /
# `backend="mumps"` at a call site always wins over this default; only the
# `"auto"` sites (the default, and every call that does not name a backend)
# consult it. Itself defaults to `"auto"` (prefer MUMPS when available, else
# SuperLU), so absent any override the behaviour is exactly as before.
# A ContextVar rather than a module global so a scoped `default_backend(...)`
# block in one thread (or async task) cannot leak into another's `"auto"`
# resolution: a fresh `threading.Thread` copies the context at start, and
# `ContextVar.set`/`.reset` only ever mutate the CURRENT context -- concurrent
# threads flipping the default no longer race each other.
_DEFAULT_BACKEND: ContextVar[_Backend] = ContextVar(
    "qscat_sparse_lu_default_backend", default="auto"
)


def _validate_backend(name: _Backend) -> None:
    if name not in ("auto", "scipy", "mumps"):
        raise ValueError(f"unknown backend {name!r}; expected auto/scipy/mumps")


def set_default_backend(name: _Backend) -> None:
    """Set the default backend `SparseLU(backend="auto")` resolves to.

    HAZARD: this mutates the CURRENT context for the rest of the process
    (or thread/task) lifetime and is easy to leave flipped -- prefer the
    `default_backend` context manager, which restores the previous value
    on exit (including on exception). Only `"auto"` call sites consult
    this; an explicit `backend="scipy"`/`"mumps"` argument always wins.
    """
    _validate_backend(name)
    _DEFAULT_BACKEND.set(name)


def get_default_backend() -> _Backend:
    """The current default backend (see `set_default_backend`)."""
    return _DEFAULT_BACKEND.get()


@contextmanager
def default_backend(name: _Backend) -> Iterator[None]:
    """Temporarily force the `"auto"` backend to `name` within a `with` block.

    The recommended way to steer internal `SparseLU(...)` construction --
    e.g. forcing a whole computation that builds `SparseLU` internally
    (`projects.n2_2d_cross_section.ve_cross_section_2d`) through one specific
    factorization backend, so two backends can be compared for physics
    equivalence without threading a `backend=` kwarg through every call site.
    Scoped, exception-safe, and context-local, so concurrent threads cannot
    race each other's defaults.
    """
    _validate_backend(name)
    token = _DEFAULT_BACKEND.set(name)
    try:
        yield
    finally:
        _DEFAULT_BACKEND.reset(token)


class _ScipyBackend:
    """The original `scipy.sparse.linalg.splu` path, unchanged.

    Holds exactly the factorization object and semantics `SparseLU` used
    before backend dispatch existed: `fill_factor` reads SuperLU's own `nnz`
    with no materialization, `memory_bytes()` materializes (and permanently
    caches, on the `SuperLU` object) the `L`/`U` CSC factors. See the module
    docstring for the memory caveat.
    """

    name = "scipy"

    def __init__(self, csc: sp.csc_matrix[np.complex128], ordering: Ordering) -> None:
        self._ordering = ordering
        # Store the analyzed pattern for the `refactor` guard. scipy has no
        # symbolic-reuse hook, so `refactor` re-runs `splu`; the guard still
        # holds `refactor` to the same-pattern contract as the MUMPS path.
        self._pattern = _pattern_of(csc)
        self._lu: spla.SuperLU[np.complex128] = spla.splu(csc, permc_spec=ordering)

    def refactor(self, csc: sp.csc_matrix[np.complex128]) -> None:
        """Re-factorize ``csc`` (scipy: a fresh ``splu``, no symbolic reuse).

        scipy exposes no clean symbolic-reuse hook, so this simply re-runs
        ``splu`` with the original ordering -- correct, but with no speedup over
        constructing a new `SparseLU`. The pattern guard is kept so the scipy
        and MUMPS paths share one contract (same sparsity pattern required).
        """
        _check_pattern(self._pattern, csc)
        self._lu = spla.splu(csc, permc_spec=self._ordering)

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
    above; `"mumps"` is the complex-symmetric MUMPS factorization (available
    only where system MUMPS + the `qscat[mumps]` extra are installed -- the
    Docker image, not a bare Mac -- and raising `RuntimeError` if forced when
    absent); `"auto"` (the default) prefers MUMPS when available and falls back
    to scipy otherwise, so on a MUMPS-less box `"auto"` and `"scipy"` are
    identical in every observable way. `backend_used` reports which one
    actually ran. An `"auto"` call site also consults the context-local
    override set by `set_default_backend` / the `default_backend` context
    manager (an explicit `"scipy"`/`"mumps"` here overrides it) -- the seam
    used to force an entire computation that builds `SparseLU` internally
    onto one engine for a backend-equivalence check.

    `symmetric`, if left `None`, is auto-detected as `A == A.T` to a SCALED
    tolerance (an O(nnz) sparse comparison:
    `(abs(A - A.T)).max() <= _SYM_RTOL * abs(A).max()`, cheap relative to the
    factorization itself). The tolerance -- not exact equality -- is deliberate:
    the N2 decks this class factors are `A = A.T` mathematically but symmetric
    only to round-off (Kronecker-sum float reordering; ~3.6e-17 relative
    asymmetry), so exact equality would reject them and forfeit the whole point
    of the MUMPS backend (see `_SYM_RTOL`). The flag is informational only on
    the scipy path -- SuperLU does not exploit symmetry -- but on the MUMPS path
    it selects the complex-symmetric `SYM=2` matrix type (upper triangle only)
    instead of the general unsymmetric `SYM=0` one. Pass an explicit
    `symmetric=True`/`False` to override the auto-detect entirely.
    """

    def __init__(
        self,
        A: sp.spmatrix,
        *,
        ordering: Ordering = "COLAMD",
        backend: _Backend = "auto",
        symmetric: bool | None = None,
    ) -> None:
        if A.shape[0] != A.shape[1]:
            raise ValueError(f"matrix must be square, got shape {A.shape}")
        csc: sp.csc_matrix[np.complex128] = sp.csc_matrix(A, dtype=np.complex128)
        self._shape: tuple[int, int] = (int(csc.shape[0]), int(csc.shape[1]))
        self._nnz: int = int(csc.nnz)
        self._ordering = ordering

        # Resolve `"auto"` against the context-local default; an explicit
        # `"scipy"`/`"mumps"` at the call site is honored verbatim and never
        # consults the override.
        resolved: _Backend = get_default_backend() if backend == "auto" else backend

        self._backend_used: str
        if resolved == "mumps" and not mumps_available():
            # An explicit (or defaulted-to) request for MUMPS must fail loudly
            # rather than silently falling back to scipy -- and before doing
            # any (wasted) symmetry detection on this error path.
            raise BackendError(
                "MUMPS backend requested but not available (qscat[mumps] / system MUMPS missing)"
            )

        if symmetric is None:
            # O(nnz) sparse comparison -- cheap relative to the factorization
            # that follows, but not free, hence computed once and cached. A
            # SCALED tolerance (`_SYM_RTOL`), not exact equality: the real N2
            # decks are symmetric only to round-off (see `_SYM_RTOL`'s comment).
            scale = abs(csc).max() if csc.nnz else 0.0
            if scale == 0.0:
                symmetric = True  # a zero matrix is trivially symmetric
            else:
                symmetric = bool((abs(csc - csc.T)).max() <= _SYM_RTOL * scale)
        self._symmetric = symmetric

        self._impl: _ScipyBackend | _MumpsBackend
        if resolved == "scipy":
            self._impl = _ScipyBackend(csc, ordering)
            self._backend_used = "scipy"
        elif resolved == "mumps":
            # Availability already checked above; select SYM=2 vs SYM=0 from
            # the (auto-detected or overridden) symmetry flag.
            self._impl = _MumpsBackend(csc, symmetric=self._symmetric)
            self._backend_used = "mumps"
        else:  # resolved == "auto": prefer MUMPS when available, else scipy.
            if mumps_available():
                self._impl = _MumpsBackend(csc, symmetric=self._symmetric)
                self._backend_used = "mumps"
            else:
                self._impl = _ScipyBackend(csc, ordering)
                self._backend_used = "scipy"

    def refactor(self, A_new: sp.spmatrix) -> None:
        """Re-factorize `A_new` reusing this object's symbolic analysis.

        `A_new` MUST share the original matrix's sparsity pattern (e.g. a
        diagonal shift `E*I - H` across energies). On the MUMPS backend this
        reuses the analysis (skips re-ordering); on scipy it re-runs `splu`
        (correct, no reuse). Keeps the original backend and symmetry decision.
        Raises `ValueError` on a shape/pattern mismatch.
        """
        if A_new.shape != self._shape:
            raise ValueError(f"refactor shape {A_new.shape} != {self._shape}")
        csc: sp.csc_matrix[np.complex128] = sp.csc_matrix(A_new, dtype=np.complex128)
        self._impl.refactor(csc)
        self._nnz = int(csc.nnz)

    @property
    def shape(self) -> tuple[int, int]:
        return self._shape

    @property
    def ordering(self) -> str:
        return self._ordering

    @property
    def symmetric(self) -> bool:
        """Whether `A` was treated as (complex-)symmetric `A == A.T`.

        Auto-detected from `A` (to the scaled `_SYM_RTOL` tolerance) when
        `symmetric=None` was passed (the default), or the explicit override.
        Informational on the scipy path; the MUMPS path uses it to select
        the complex-symmetric (`SYM=2`) matrix type instead of the general
        unsymmetric one.
        """
        return self._symmetric

    @property
    def backend_used(self) -> str:
        """Which backend actually factorized `A`: `"scipy"` or `"mumps"`."""
        return self._backend_used

    @property
    def ordering_used(self) -> str:
        """The ordering the active backend actually used.

        On the scipy path this is scipy's `permc_spec` (identical to
        `ordering`); on the MUMPS path it is MUMPS's own chosen ordering
        (e.g. `"scotch"`/`"metis"`/`"amd"`), read from `INFOG(7)`.
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
                f"right-hand side has leading dimension {rhs.shape[0]}, expected {self._shape[0]}"
            )
        return self._impl.solve(rhs.astype(np.complex128, copy=False))
