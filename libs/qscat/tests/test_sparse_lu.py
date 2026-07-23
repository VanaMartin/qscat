"""Tests for `qscat.linalg.SparseLU` (V5).

Exercised on a COMPLEX SYMMETRIC matrix (H = H^T, not Hermitian), which is what
exterior complex scaling produces and what every real use of this class will be.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from qscat.linalg import SparseLU


def _complex_symmetric(n: int, seed: int) -> sp.csc_matrix:
    """A well-conditioned, sparse, complex-symmetric test matrix."""
    rng = np.random.default_rng(seed)
    nnz = 5 * n
    rows = rng.integers(0, n, size=nnz)
    cols = rng.integers(0, n, size=nnz)
    vals = rng.standard_normal(nnz) + 1j * rng.standard_normal(nnz)
    m = sp.coo_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex).tocsr()
    m = m + m.T                                             # complex SYMMETRIC, no conjugate
    m = m + sp.identity(n, format="csr", dtype=complex) * (10.0 + 3.0j)  # diagonally dominant
    return sp.csc_matrix(m)


def test_solve_residual_is_at_round_off() -> None:
    n = 200
    A = _complex_symmetric(n, seed=0)
    rng = np.random.default_rng(10)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    x = SparseLU(A).solve(b)
    assert np.linalg.norm(A @ x - b) / np.linalg.norm(b) < 1e-12


def test_matrix_is_complex_symmetric_not_hermitian() -> None:
    """Guard the fixture itself -- a Hermitian matrix would not exercise the point."""
    A = _complex_symmetric(50, seed=1)
    assert abs(A - A.T).max() < 1e-14
    assert abs(A - A.conj().T).max() > 1e-3


def test_multi_rhs_matches_looped_single_solves() -> None:
    n = 150
    A = _complex_symmetric(n, seed=2)
    rng = np.random.default_rng(11)
    B = rng.standard_normal((n, 4)) + 1j * rng.standard_normal((n, 4))
    lu = SparseLU(A)
    together = lu.solve(B)
    assert together.shape == (n, 4)
    for j in range(4):
        assert np.allclose(together[:, j], lu.solve(B[:, j]), rtol=0, atol=1e-12)


def test_factorization_is_reused_across_solves() -> None:
    """One factorization, many right-hand sides -- the whole reason this class exists."""
    n = 120
    A = _complex_symmetric(n, seed=3)
    lu = SparseLU(A)
    rng = np.random.default_rng(12)
    for _ in range(5):
        b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        assert np.linalg.norm(A @ lu.solve(b) - b) / np.linalg.norm(b) < 1e-12


def test_diagnostics_are_reported() -> None:
    A = _complex_symmetric(200, seed=4)
    lu = SparseLU(A)
    assert lu.shape == (200, 200)
    assert lu.fill_factor >= 1.0
    assert lu.memory_bytes() > 0


def test_fill_factor_is_a_sane_proxy_for_materialized_lu_nnz() -> None:
    """`fill_factor` reads SuperLU's own raw storage `nnz` -- cheap, no
    materialization -- rather than `(L.nnz + U.nnz) / A.nnz` from the
    materialized CSC factors. The two are NOT expected to match exactly
    (SuperLU's internal supernodal storage carries some explicit padding that
    compresses away when converted to CSC -- measured ~7% higher on this
    fixture), but they must stay within a sane factor of each other, as a
    guard against a wrong attribute or denominator. Computed independently
    via a bare `scipy.sparse.linalg.splu` call, not by reaching into
    `SparseLU`'s internals.
    """
    A = _complex_symmetric(150, seed=7)
    lu = SparseLU(A)
    bare = spla.splu(sp.csc_matrix(A, dtype=np.complex128), permc_spec="COLAMD")
    materialized = float(bare.L.nnz + bare.U.nnz) / float(A.nnz)
    assert lu.fill_factor / materialized == pytest.approx(1.0, rel=0.5)


def test_ordering_is_configurable_and_changes_fill() -> None:
    """Every ordering must solve correctly, and fill-in must actually differ.

    Measured `fill_factor` on this exact matrix (n=300, seed=5):
    NATURAL=19.5417, COLAMD=18.7879, MMD_AT_PLUS_A=9.6030. Only the large,
    structural gap (MMD_AT_PLUS_A roughly halving the fill relative to either
    other ordering) is asserted -- a strict 3-way ordering between NATURAL and
    COLAMD is SuperLU-version dependent and plausibly brittle, so it is not
    asserted here even though it held on this run.
    """
    n = 300
    A = _complex_symmetric(n, seed=5)
    rng = np.random.default_rng(13)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    fill_factors: dict[str, float] = {}
    for ordering in ("COLAMD", "MMD_AT_PLUS_A", "NATURAL"):
        lu = SparseLU(A, ordering=ordering)
        assert np.linalg.norm(A @ lu.solve(b) - b) / np.linalg.norm(b) < 1e-12
        fill_factors[ordering] = lu.fill_factor
    assert fill_factors["MMD_AT_PLUS_A"] < 0.7 * fill_factors["COLAMD"]
    assert fill_factors["MMD_AT_PLUS_A"] < 0.7 * fill_factors["NATURAL"]


def test_rejects_non_square() -> None:
    with pytest.raises(ValueError, match="square"):
        SparseLU(sp.csc_matrix(np.zeros((3, 4))))


def test_accepts_csr_input_without_warning() -> None:
    """CSR in must be converted internally, not warned about."""
    A = sp.csr_matrix(_complex_symmetric(80, seed=6))
    rng = np.random.default_rng(14)
    b = rng.standard_normal(80) + 1j * rng.standard_normal(80)
    with warnings.catch_warnings():
        warnings.simplefilter("error", sp.SparseEfficiencyWarning)
        x = SparseLU(A).solve(b)
    assert np.linalg.norm(A @ x - b) / np.linalg.norm(b) < 1e-12


def test_solve_rejects_0d_right_hand_side() -> None:
    """A 0-d `b` used to raise a bare `IndexError` (`rhs.shape[0]` on an empty
    shape tuple) -- an unhelpful, implementation-accident error. It must raise
    `ValueError` with a message that says what was wrong.
    """
    A = _complex_symmetric(10, seed=8)
    lu = SparseLU(A)
    with pytest.raises(ValueError, match="0-d"):
        lu.solve(np.array(1.0 + 0j))
