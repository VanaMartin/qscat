"""Tests for `qscat.linalg.ShiftInvertEigs`.

The oracle throughout is dense `np.linalg.eig` on the same matrix: this class
must return exactly the eigenpairs the dense solver finds nearest the shift.
Matrices are COMPLEX SYMMETRIC (A == A.T, not Hermitian) -- what exterior
complex scaling produces and what every real use of this class will be.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from qscat.linalg import ShiftInvertEigs, c_product


def _complex_symmetric(n: int, seed: int) -> sp.csc_matrix:
    """A well-conditioned, sparse, complex-symmetric test matrix."""
    rng = np.random.default_rng(seed)
    nnz = 5 * n
    rows = rng.integers(0, n, size=nnz)
    cols = rng.integers(0, n, size=nnz)
    vals = rng.standard_normal(nnz) + 1j * rng.standard_normal(nnz)
    m = sp.coo_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex).tocsr()
    m = m + m.T  # symmetric, NOT conjugated
    m = m + sp.identity(n, format="csr", dtype=complex) * (10.0 + 3.0j)
    return sp.csc_matrix(m)


def _nearest_dense(A: sp.csc_matrix, sigma: complex, k: int) -> np.ndarray:
    """The k dense eigenvalues nearest `sigma`, ascending in |E - sigma|."""
    w = np.linalg.eigvals(A.toarray())
    return w[np.argsort(np.abs(w - sigma))][:k]


def test_eigenvalues_match_dense_nearest_the_shift() -> None:
    A = _complex_symmetric(200, seed=0)
    sigma = 9.0 + 2.0j
    k = 6
    vals, _ = ShiftInvertEigs(A, k=k).near(sigma)
    assert np.allclose(vals, _nearest_dense(A, sigma, k), rtol=1e-9, atol=1e-12)


def test_eigenvalues_are_sorted_by_distance_from_the_shift() -> None:
    A = _complex_symmetric(200, seed=1)
    sigma = 11.0 - 1.0j
    vals, _ = ShiftInvertEigs(A, k=5).near(sigma)
    d = np.abs(vals - sigma)
    assert np.all(np.diff(d) >= 0.0)


def test_eigenvectors_match_dense_up_to_scale() -> None:
    """Compared under the c-norm sqrt(v@v): the right notion of 'same vector up
    to scale' for a complex-symmetric operator (v@v != 0 away from an
    exceptional point)."""
    A = _complex_symmetric(150, seed=2)
    sigma = 8.5 + 1.5j
    vals, vecs = ShiftInvertEigs(A, k=3).near(sigma)
    w, V = np.linalg.eig(A.toarray())
    for i, val in enumerate(vals):
        j = int(np.argmin(np.abs(w - val)))
        u, v = V[:, j], vecs[:, i]
        u = u / np.sqrt(c_product(u, u))
        v = v / np.sqrt(c_product(v, v))
        assert abs(abs(c_product(u, v)) - 1.0) < 1e-6


def test_shift_sign_convention_is_A_minus_sigma_I() -> None:
    """A spectrum deliberately asymmetric about sigma: passing `sigma*I - A` as
    OPinv instead of `A - sigma*I` returns eigenvalues reflected about sigma,
    which this pins down."""
    diag = np.array([0.0, 1.0, 2.0, 3.0, 10.0, 11.0, 12.0], dtype=complex)
    A = sp.csc_matrix(sp.diags(diag, dtype=complex))
    sigma = 9.0 + 0.0j
    vals, _ = ShiftInvertEigs(A, k=2).near(sigma)
    assert np.allclose(np.sort_complex(vals), np.array([10.0, 11.0], dtype=complex))


def test_k_too_large_raises_and_points_at_the_dense_route() -> None:
    A = _complex_symmetric(20, seed=3)
    with pytest.raises(ValueError, match="qscat.dvr.eigen"):
        ShiftInvertEigs(A, k=19).near(1.0 + 0.0j)


def test_non_square_raises() -> None:
    A = sp.csc_matrix((5, 7), dtype=complex)
    with pytest.raises(ValueError, match="square"):
        ShiftInvertEigs(A)


def test_repeated_shifts_reuse_the_factorization_object() -> None:
    """Reuse must be invisible in the results and visible in the diagnostics."""
    A = _complex_symmetric(200, seed=4)
    s1, s2 = 9.0 + 2.0j, 12.0 - 1.0j
    solver = ShiftInvertEigs(A, k=4)
    v1, _ = solver.near(s1)
    v2, _ = solver.near(s2)
    assert solver.n_factorizations == 2
    fresh1, _ = ShiftInvertEigs(A, k=4).near(s1)
    fresh2, _ = ShiftInvertEigs(A, k=4).near(s2)
    assert np.allclose(v1, fresh1, rtol=1e-9, atol=1e-12)
    assert np.allclose(v2, fresh2, rtol=1e-9, atol=1e-12)


def test_diagnostics_delegate_to_the_factorization() -> None:
    A = _complex_symmetric(100, seed=5)
    solver = ShiftInvertEigs(A, k=3)
    with pytest.raises(RuntimeError, match="near"):
        _ = solver.backend_used  # nothing factored yet
    solver.near(9.0 + 1.0j)
    assert solver.backend_used in {"scipy", "mumps"}
    assert solver.ordering_used
    assert solver.fill_factor > 0.0
    assert solver.memory_bytes() > 0
    assert solver.shape == (100, 100)


def test_non_convergence_raises_convergence_error() -> None:
    """maxiter=1 is far too few restarts: ARPACK bails, and we translate."""
    from qscat.exceptions import ConvergenceError

    A = _complex_symmetric(300, seed=6)
    solver = ShiftInvertEigs(A, k=8, ncv=10, maxiter=1)
    with pytest.raises(ConvergenceError, match="ARPACK"):
        solver.near(0.0 + 0.0j)
