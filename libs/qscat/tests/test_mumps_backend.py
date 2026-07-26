"""MUMPS complex-symmetric (SYM=2) backend tests.

Every test is skipped unless `python-mumps` / system MUMPS is importable, so the
whole module SKIPS on a MUMPS-less box (the Mac) and RUNS in the Docker `test`
image where the `qscat[mumps]` extra and system MUMPS are installed.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from qscat.linalg import SparseLU
from qscat.linalg._mumps_backend import mumps_available

pytestmark = pytest.mark.skipif(
    not mumps_available(), reason="system MUMPS / qscat[mumps] not installed"
)


def _complex_symmetric(n: int, seed: int) -> sp.csc_matrix:
    rng = np.random.default_rng(seed)
    nnz = 5 * n
    r = rng.integers(0, n, nnz)
    c = rng.integers(0, n, nnz)
    v = rng.standard_normal(nnz) + 1j * rng.standard_normal(nnz)
    m = sp.coo_matrix((v, (r, c)), shape=(n, n), dtype=complex).tocsr()
    m = m + m.T  # complex SYMMETRIC (A == A.T, not Hermitian)
    m = m + sp.identity(n, format="csr", dtype=complex) * (10.0 + 3.0j)
    return sp.csc_matrix(m)


@pytest.mark.parametrize("n", [50, 400])
def test_mumps_solve_matches_scipy_to_roundoff(n: int) -> None:
    """THE gate: MUMPS SYM=2 (upper triangle) == SuperLU on the full matrix."""
    A = _complex_symmetric(n, seed=100 + n)
    rng = np.random.default_rng(7)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    x_mumps = SparseLU(A, backend="mumps").solve(b)
    x_scipy = SparseLU(A, backend="scipy").solve(b)
    assert np.linalg.norm(x_mumps - x_scipy) / np.linalg.norm(x_scipy) < 1e-10
    assert np.linalg.norm(A @ x_mumps - b) / np.linalg.norm(b) < 1e-10


def test_mumps_used_and_reports_sym2() -> None:
    # n large enough that MUMPS's factorization memory rounds to >= 1 MB (it is
    # reported in whole MB, so a tiny system reports 0 MB, not a failure).
    A = _complex_symmetric(500, seed=200)
    lu = SparseLU(A, backend="mumps")
    assert lu.backend_used == "mumps"
    assert lu.fill_factor >= 1.0
    assert lu.ordering_used  # non-empty; the ordering MUMPS chose
    assert lu.memory_bytes() > 0


def test_mumps_auto_prefers_mumps_when_available() -> None:
    """With MUMPS importable, backend='auto' (the default) selects MUMPS."""
    A = _complex_symmetric(100, seed=201)
    assert SparseLU(A).backend_used == "mumps"


def test_mumps_multi_rhs_matches_scipy() -> None:
    """An (N, k) block of right-hand sides solves column-wise, matching SuperLU."""
    n = 120
    A = _complex_symmetric(n, seed=202)
    rng = np.random.default_rng(9)
    B = rng.standard_normal((n, 4)) + 1j * rng.standard_normal((n, 4))
    x_mumps = SparseLU(A, backend="mumps").solve(B)
    x_scipy = SparseLU(A, backend="scipy").solve(B)
    assert x_mumps.shape == (n, 4)
    assert np.linalg.norm(x_mumps - x_scipy) / np.linalg.norm(x_scipy) < 1e-10


def _asymmetric(n: int, seed: int) -> sp.csc_matrix:
    """A well-conditioned complex matrix with A != A.T (no `m + m.T`)."""
    rng = np.random.default_rng(seed)
    nnz = 5 * n
    r = rng.integers(0, n, nnz)
    c = rng.integers(0, n, nnz)
    v = rng.standard_normal(nnz) + 1j * rng.standard_normal(nnz)
    m = sp.coo_matrix((v, (r, c)), shape=(n, n), dtype=complex).tocsr()
    m = m + sp.identity(n, format="csr", dtype=complex) * (10.0 + 3.0j)
    return sp.csc_matrix(m)


def test_mumps_upper_triangle_trap_would_be_caught() -> None:
    """A deliberately asymmetric matrix must NOT be silently treated as symmetric:
    forcing symmetric=True on a non-symmetric A must differ from the truth,
    proving the symmetric path really uses only the upper triangle."""
    n = 60
    A = _asymmetric(n, seed=3)
    assert abs(A - A.T).max() > 1e-6  # guard: the fixture really is asymmetric
    rng = np.random.default_rng(30)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    x_true = SparseLU(A, backend="scipy").solve(b)
    x_sym_wrong = SparseLU(A, backend="mumps", symmetric=True).solve(b)  # wrongly claims sym
    # forcing SYM=2 on a non-symmetric A (upper triangle only) gives a DIFFERENT,
    # wrong answer -- confirming the symmetric path genuinely drops the lower triangle
    assert np.linalg.norm(x_sym_wrong - x_true) / np.linalg.norm(x_true) > 1e-3
    # and the correct (unsymmetric SYM=0) MUMPS path matches scipy
    x_unsym = SparseLU(A, backend="mumps", symmetric=False).solve(b)
    assert np.linalg.norm(x_unsym - x_true) / np.linalg.norm(x_true) < 1e-10
