"""Sparse shift-invert eigensolver: interior eigenvalues near a complex shift.

Krylov methods converge to EXTREMAL eigenvalues, but a resonance is an INTERIOR
one -- it sits in the lower half of the complex plane surrounded by the
discretized (rotated) continuum. The standard remedy is the shift-invert
spectral transform: run Arnoldi on `(A - sigma*I)^-1`, whose extremal
eigenvalues are `A`'s eigenvalues nearest `sigma`.

That transform costs one sparse solve per matrix-vector product, which is what
`SparseLU` provides -- and because `A - sigma*I` has the SAME sparsity pattern
for every `sigma` (a diagonal shift), `SparseLU.refactor` reuses the symbolic
analysis across shifts. A shift sweep therefore gets the same discount as the
time-independent energy sweep (see `docs/physics/ti-energy-sweep-reuse.md`).

See `docs/physics/shift-invert-eigensolver.md`.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from qscat.exceptions import ConvergenceError

from .sparse_lu import Ordering, SparseLU, _Backend

__all__ = ["ShiftInvertEigs"]


class ShiftInvertEigs:
    """The `k` eigenpairs of a sparse matrix nearest a complex shift.

    Holds the factorization: the first `near` call factors `A - sigma*I`, and
    every later call `refactor`s it, reusing the symbolic analysis (the shifted
    matrix's sparsity pattern does not depend on `sigma`). On the MUMPS backend
    that skips the SCOTCH ordering; on scipy it re-runs `splu` (correct, no
    reuse) -- the same trade `SparseLU.refactor` documents.

    `ordering`, `backend` and `symmetric` are forwarded verbatim to `SparseLU`.
    `k`, `ncv`, `tol` and `maxiter` are ARPACK controls; `k` may be overridden
    per call.

    Parameters
    ----------
    A : scipy.sparse.spmatrix
        Square matrix. Promoted to complex CSC internally.
    k : int, optional
        Number of eigenpairs to return per call (default 6).
    ordering : {"COLAMD", "NATURAL", "MMD_ATA", "MMD_AT_PLUS_A"}, optional
        SuperLU column ordering, forwarded to `SparseLU`.
    backend : {"auto", "scipy", "mumps"}, optional
        Factorization backend, forwarded to `SparseLU`.
    symmetric : bool or None, optional
        Complex-symmetry flag, forwarded to `SparseLU` (auto-detected if None).
    ncv : int or None, optional
        Krylov subspace size. ARPACK's default is used when None.
    tol : float, optional
        ARPACK relative tolerance; 0.0 means machine precision.
    maxiter : int or None, optional
        Maximum ARPACK restarts.

    Notes
    -----
    Eigenvectors are returned with numpy's Euclidean (`v^dagger v = 1`)
    normalization, exactly as `qscat.dvr.eigen` returns them. For
    exterior-complex-scaling observables re-normalize under the bilinear
    `c_product`.

    Examples
    --------
    >>> import numpy as np, scipy.sparse as sp
    >>> from qscat.linalg import ShiftInvertEigs
    >>> A = sp.diags(np.array([0.0, 1.0, 2.0, 10.0, 11.0], dtype=complex))
    >>> vals, vecs = ShiftInvertEigs(A, k=2).near(9.5 + 0.0j)
    >>> np.round(np.sort_complex(vals).real, 6)
    array([10., 11.])
    """

    def __init__(
        self,
        A: sp.spmatrix,
        *,
        k: int = 6,
        ordering: Ordering = "COLAMD",
        backend: _Backend = "auto",
        symmetric: bool | None = None,
        ncv: int | None = None,
        tol: float = 0.0,
        maxiter: int | None = None,
    ) -> None:
        if A.shape[0] != A.shape[1]:
            raise ValueError(f"matrix must be square, got shape {A.shape}")
        self._A: sp.csc_matrix[np.complex128] = sp.csc_matrix(A, dtype=np.complex128)
        self._n = int(self._A.shape[0])
        self._eye = sp.identity(self._n, format="csc", dtype=np.complex128)
        self._k = int(k)
        self._ordering: Ordering = ordering
        self._backend: _Backend = backend
        self._symmetric = symmetric
        self._ncv = ncv
        self._tol = float(tol)
        self._maxiter = maxiter
        self._lu: SparseLU | None = None
        self._n_factorizations = 0

    def near(
        self, sigma: complex, *, k: int | None = None
    ) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
        """The `k` eigenpairs of `A` nearest `sigma`, nearest first.

        Parameters
        ----------
        sigma : complex
            The shift: eigenvalues are returned in order of increasing
            `|E - sigma|`.
        k : int or None, optional
            Number of eigenpairs for this call; the constructor's `k` if None.

        Returns
        -------
        energies : ndarray of complex128, shape (k,)
            Eigenvalues of `A`, sorted by `|E - sigma|` ascending.
        vectors : ndarray of complex128, shape (n, k)
            `vectors[:, i]` is the eigenvector of `energies[i]`, Euclidean-
            normalized.

        Raises
        ------
        ValueError
            If `k` is not in `1 <= k < n - 1` (an ARPACK constraint).
        qscat.exceptions.ConvergenceError
            If ARPACK fails to converge.
        """
        k_eff = self._k if k is None else int(k)
        if k_eff < 1:
            raise ValueError(f"k must be >= 1, got {k_eff}")
        if k_eff >= self._n - 1:
            raise ValueError(
                f"k={k_eff} must be < n-1 = {self._n - 1} (an ARPACK constraint); "
                "for a matrix this small use the dense qscat.dvr.eigen instead"
            )

        shift = complex(sigma)
        # scipy's convention: OPinv solves (A - sigma*I) x = b. NOT the driven
        # solver's (E*I - H); passing that instead silently returns eigenvalues
        # reflected about sigma. Adding the identity explicitly also forces the
        # diagonal to be structurally present, so the pattern -- and hence the
        # symbolic analysis reused by `refactor` -- is identical for every shift.
        shifted = (self._A - shift * self._eye).tocsc()
        if self._lu is None:
            self._lu = SparseLU(
                shifted,
                ordering=self._ordering,
                backend=self._backend,
                symmetric=self._symmetric,
            )
        else:
            self._lu.refactor(shifted)
        self._n_factorizations += 1

        op_inv = spla.LinearOperator((self._n, self._n), matvec=self._lu.solve, dtype=np.complex128)
        try:
            vals, vecs = spla.eigs(
                self._A,
                k=k_eff,
                sigma=shift,
                OPinv=op_inv,
                ncv=self._ncv,
                tol=self._tol,
                maxiter=self._maxiter,
            )
        except spla.ArpackNoConvergence as exc:
            raise ConvergenceError(
                f"ARPACK did not converge at sigma={shift!r} with k={k_eff}: "
                f"{np.size(exc.eigenvalues)} of {k_eff} eigenvalues converged. "
                "Raise k, raise ncv, or move sigma closer to the target."
            ) from exc

        order = np.argsort(np.abs(vals - shift))
        return (
            np.asarray(vals[order], dtype=np.complex128),
            np.asarray(vecs[:, order], dtype=np.complex128),
        )

    @property
    def shape(self) -> tuple[int, int]:
        """Operator shape `(n, n)`."""
        return (self._n, self._n)

    @property
    def n_factorizations(self) -> int:
        """How many shifted matrices have been factored (analysis + refactors)."""
        return self._n_factorizations

    def _require_lu(self) -> SparseLU:
        if self._lu is None:
            raise RuntimeError(
                "no factorization yet -- call near(sigma) before reading diagnostics"
            )
        return self._lu

    @property
    def backend_used(self) -> str:
        """Which factorization engine ran (`"scipy"` or `"mumps"`)."""
        return self._require_lu().backend_used

    @property
    def ordering_used(self) -> str:
        """The ordering the factorization actually used."""
        return self._require_lu().ordering_used

    @property
    def fill_factor(self) -> float:
        """Factor nnz relative to the shifted matrix's nnz."""
        return self._require_lu().fill_factor

    def memory_bytes(self) -> int:
        """Factor memory. NOT cheap -- a method, not a property, so the cost is opt-in."""
        return self._require_lu().memory_bytes()
