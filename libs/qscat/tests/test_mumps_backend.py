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


def _roundoff_symmetric(n: int, seed: int) -> sp.csc_matrix:
    """`A == A.T` mathematically but only to ROUND-OFF, like the real N2 decks.

    Built from a true symmetric matrix with its strictly-lower triangle
    perturbed by ~1e-13 (relative asymmetry ~1e-14) -- inside `_SYM_RTOL`=1e-12,
    so the auto-detect returns True, but NOT bit-exact, so the old
    exact-equality detect returned False and forced SYM=0.
    """
    A = _complex_symmetric(n, seed).tolil()
    rng = np.random.default_rng(seed + 999)
    for _ in range(3 * n):
        i = int(rng.integers(1, n))
        j = int(rng.integers(0, i))
        A[i, j] += (rng.standard_normal() + 1j * rng.standard_normal()) * 1e-13
    return sp.csc_matrix(A)


def test_mumps_engages_sym2_on_roundoff_symmetric_matrix() -> None:
    """The fix, proven on a realistic (round-off-symmetric) matrix, not just an
    exactly-symmetric synthetic one: the auto-detect now flags it symmetric, so
    the MUMPS path runs SYM=2, and it STILL matches SuperLU on the FULL matrix
    to tight rtol. rtol is 1e-9 (slightly looser than the exactly-symmetric
    1e-10) because SYM=2 takes the upper triangle as truth and reconstructs the
    lower from it, perturbing the matrix at the ~1e-13 round-off level."""
    n = 400
    A = _roundoff_symmetric(n, seed=300)
    assert abs(A - A.T).max() != 0.0  # NOT bit-exact (exact-equality would say False)
    assert abs(A - A.T).max() / abs(A).max() < 1e-12  # but inside the tolerance

    lu = SparseLU(A, backend="mumps")  # symmetric=None => auto-detect
    assert lu.symmetric is True  # SYM=2 genuinely engaged (was False before the fix)

    rng = np.random.default_rng(13)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    x_mumps = lu.solve(b)
    x_scipy = SparseLU(A, backend="scipy").solve(b)
    assert np.linalg.norm(x_mumps - x_scipy) / np.linalg.norm(x_scipy) < 1e-9
    assert np.linalg.norm(A @ x_mumps - b) / np.linalg.norm(b) < 1e-9


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


def test_mumps_refactor_reuses_analysis_matches_fresh() -> None:
    """analyze once, refactor(A_shift) per shift == fresh SuperLU each time."""
    n = 400
    A0 = _complex_symmetric(n, seed=300)
    lu = SparseLU(A0, backend="mumps")
    rng = np.random.default_rng(9)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    for shift in (2.0 + 1.0j, -3.0 + 0.5j, 5.0 - 2.0j):
        A = (A0 + shift * sp.identity(n, dtype=complex)).tocsc()
        lu.refactor(A)
        x = lu.solve(b)
        x_ref = SparseLU(A, backend="scipy").solve(b)
        assert np.linalg.norm(x - x_ref) / np.linalg.norm(x_ref) < 1e-9
    assert lu.backend_used == "mumps"
