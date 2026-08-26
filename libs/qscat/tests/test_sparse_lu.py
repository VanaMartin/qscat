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
from qscat.linalg import SparseLU, default_backend, get_default_backend, set_default_backend
from qscat.linalg._mumps_backend import mumps_available


def _complex_symmetric(n: int, seed: int) -> sp.csc_matrix:
    """A well-conditioned, sparse, complex-symmetric test matrix."""
    rng = np.random.default_rng(seed)
    nnz = 5 * n
    rows = rng.integers(0, n, size=nnz)
    cols = rng.integers(0, n, size=nnz)
    vals = rng.standard_normal(nnz) + 1j * rng.standard_normal(nnz)
    m = sp.coo_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex).tocsr()
    m = m + m.T  # complex SYMMETRIC, no conjugate
    m = m + sp.identity(n, format="csr", dtype=complex) * (10.0 + 3.0j)  # diagonally dominant
    return sp.csc_matrix(m)


def test_solve_residual_is_at_round_off() -> None:
    n = 200
    A = _complex_symmetric(n, seed=0)
    rng = np.random.default_rng(10)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    x = SparseLU(A).solve(b)
    assert np.linalg.norm(A @ x - b) / np.linalg.norm(b) < 1e-12


def test_ordering_literal_is_public_and_single_sourced() -> None:
    """lib-M12: one public Ordering type; the solver modules import it."""
    from typing import get_args

    from qscat import linalg
    from qscat.core import dissociation, driven, lcp, problem

    assert "Ordering" in linalg.__all__
    assert get_args(linalg.Ordering) == ("NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD")
    # the four former private copies are gone
    for mod in (driven, dissociation, lcp, problem):
        assert not hasattr(mod, "_Ordering"), mod.__name__


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
    # Pinned to scipy: this asserts the SuperLU diagnostics contract (exact
    # materialized-factor bytes, always > 0). MUMPS reports memory in whole MB
    # -- 0 for a matrix this small -- and has its own diagnostics test.
    A = _complex_symmetric(200, seed=4)
    lu = SparseLU(A, backend="scipy")
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
    lu = SparseLU(A, backend="scipy")  # SuperLU-specific proxy, pin the backend
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
        # `ordering` is scipy's permc_spec; MUMPS chooses its own ordering, so
        # this fill-in comparison is meaningful only on the scipy backend.
        lu = SparseLU(A, ordering=ordering, backend="scipy")
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


def test_default_backend_prefers_mumps_when_available_else_scipy() -> None:
    """backend='auto' picks MUMPS when available, else SuperLU; either way it
    solves to roundoff. With MUMPS absent (e.g. the Mac) it IS the old SuperLU
    behaviour exactly."""
    A = _complex_symmetric(200, seed=20)
    rng = np.random.default_rng(21)
    b = rng.standard_normal(200) + 1j * rng.standard_normal(200)
    lu = SparseLU(A)  # unchanged call site
    assert lu.backend_used == ("mumps" if mumps_available() else "scipy")
    x = lu.solve(b)
    assert np.linalg.norm(A @ x - b) / np.linalg.norm(b) < 1e-12


def test_force_scipy_backend() -> None:
    A = _complex_symmetric(120, seed=22)
    lu = SparseLU(A, backend="scipy")
    assert lu.backend_used == "scipy"
    assert lu.ordering_used == "COLAMD"  # the default permc_spec on the scipy path


def test_ordering_still_applies_on_scipy_path() -> None:
    A = _complex_symmetric(300, seed=23)
    lu = SparseLU(A, ordering="MMD_AT_PLUS_A", backend="scipy")
    assert lu.ordering == "MMD_AT_PLUS_A"
    assert lu.backend_used == "scipy"


def test_symmetric_autodetect_flag_is_recorded() -> None:
    """A == A.T is genuinely detected (the MUMPS path uses it to pick SYM=2)."""
    A = _complex_symmetric(80, seed=24)  # symmetric fixture (A == A.T)
    lu = SparseLU(A)  # symmetric=None => auto-detect
    assert lu.symmetric is True  # detected symmetric
    assert lu.backend_used == ("mumps" if mumps_available() else "scipy")


def _asymmetric(n: int, seed: int) -> sp.csc_matrix:
    """A well-conditioned sparse complex matrix with A != A.T (built like
    `_complex_symmetric` but WITHOUT the `m + m.T` symmetrization)."""
    rng = np.random.default_rng(seed)
    nnz = 5 * n
    rows = rng.integers(0, n, size=nnz)
    cols = rng.integers(0, n, size=nnz)
    vals = rng.standard_normal(nnz) + 1j * rng.standard_normal(nnz)
    m = sp.coo_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex).tocsr()
    m = m + sp.identity(n, format="csr", dtype=complex) * (10.0 + 3.0j)
    return sp.csc_matrix(m)


def test_symmetric_autodetect_false_on_asymmetric_and_overridable() -> None:
    """A genuinely non-symmetric matrix is detected False; the flag is overridable."""
    A = _asymmetric(80, seed=25)
    assert abs(A - A.T).max() > 1e-6  # guard: the fixture really is asymmetric
    assert SparseLU(A).symmetric is False  # auto-detect
    assert SparseLU(A, symmetric=True).symmetric is True  # explicit override honored


def _roundoff_symmetric(n: int, seed: int) -> sp.csc_matrix:
    """A matrix that is `A == A.T` mathematically but only to round-off.

    Mirrors the real N2 decks: `(E_tot*I - H_2D)` is symmetric by construction
    yet its `max|A - A.T|` sits at the float noise floor (~1e-13 absolute on a
    matrix whose entries are O(1e4)) because it is assembled by Kronecker-sum
    float reordering, not bit-exact mirroring. Reproduced here by building a
    true symmetric matrix and perturbing its lower triangle by ~1e-13 -- a
    relative asymmetry ~1e-14, far below `_SYM_RTOL` = 1e-12.
    """
    A = _complex_symmetric(n, seed).tolil()
    rng = np.random.default_rng(seed + 999)
    for _ in range(3 * n):
        i = int(rng.integers(1, n))
        j = int(rng.integers(0, i))  # strictly-lower entry
        A[i, j] += (rng.standard_normal() + 1j * rng.standard_normal()) * 1e-13
    return sp.csc_matrix(A)


def test_roundoff_symmetric_is_detected_symmetric_under_tolerance() -> None:
    """The whole fix: an N2-like round-off-symmetric matrix (relative asymmetry
    ~1e-14, NOT bit-exact) is detected `symmetric is True` under `_SYM_RTOL`,
    where the old exact-equality check would have wrongly returned False and
    forfeited SYM=2. A genuinely asymmetric matrix (relative asymmetry O(1)) is
    still rejected -- pinning both sides of the tolerance WITHOUT needing MUMPS.
    """
    A = _roundoff_symmetric(120, seed=40)
    rel = abs(A - A.T).max() / abs(A).max()
    assert rel != 0.0  # guard: NOT bit-exactly symmetric (exact-equality would say False)
    assert rel < 1e-12  # but comfortably inside the tolerance
    assert SparseLU(A).symmetric is True  # detected symmetric under _SYM_RTOL

    B = _asymmetric(120, seed=41)
    assert abs(B - B.T).max() / abs(B).max() > 1e-2  # relative asymmetry O(1)
    assert SparseLU(B).symmetric is False  # still rejected -- tolerance is tight


def test_zero_matrix_scale_edge_does_not_crash_detection() -> None:
    """The `abs(A).max() == 0` edge must be handled by the auto-detect (a zero
    matrix is trivially symmetric) rather than crashing on a zero scale. A zero
    matrix is singular, so construction reaches -- and fails at -- factorization
    with SuperLU's `Factor is exactly singular` RuntimeError; the point is that
    the symmetry auto-detect BEFORE it does not raise (no divide/empty-reduce)."""
    Z = sp.csc_matrix((5, 5), dtype=complex)
    with pytest.raises(RuntimeError, match="singular"):
        SparseLU(Z, backend="scipy")


# --- V4: fallback / absence (Mac-runnable; these assert the MUMPS-absent path) ---


@pytest.mark.skipif(
    mumps_available(), reason="asserts the MUMPS-ABSENT RuntimeError path (Mac only)"
)
def test_forced_mumps_without_mumps_raises_clear_error() -> None:
    """`backend="mumps"` with MUMPS not installed must fail loudly, naming the
    missing extra, rather than silently falling back to SuperLU (Task 2's
    contract; confirmed here on the MUMPS-less box)."""
    A = _complex_symmetric(40, seed=30)
    with pytest.raises(RuntimeError, match=r"MUMPS backend requested but not available"):
        SparseLU(A, backend="mumps")


@pytest.mark.skipif(
    mumps_available(), reason="asserts the MUMPS-ABSENT auto==scipy path (Mac only)"
)
def test_auto_is_bit_identical_to_scipy_when_mumps_absent() -> None:
    """With MUMPS absent, `backend="auto"` IS the SuperLU path in every
    observable way -- same backend label and bit-identical solve."""
    A = _complex_symmetric(120, seed=31)
    rng = np.random.default_rng(32)
    b = rng.standard_normal(120) + 1j * rng.standard_normal(120)
    auto = SparseLU(A)  # backend defaults to "auto"
    scipy = SparseLU(A, backend="scipy")
    assert auto.backend_used == "scipy"
    assert np.array_equal(auto.solve(b), scipy.solve(b))  # bit-for-bit


# --- default-backend override (the seam V2 uses to force a whole computation) ---


def test_default_backend_override_is_scoped_and_restored() -> None:
    """`default_backend("scipy")` forces every `"auto"` site to SuperLU inside
    the block, and the previous default is restored on exit."""
    A = _complex_symmetric(80, seed=33)
    before = get_default_backend()
    with default_backend("scipy"):
        assert SparseLU(A).backend_used == "scipy"  # "auto" resolved to scipy
    assert get_default_backend() == before  # restored


def test_default_backend_is_context_local() -> None:
    """lib-m18: a default_backend(...) block in one thread must not leak
    into a concurrently running thread's "auto" resolution."""
    import threading

    from qscat.linalg import default_backend, get_default_backend

    seen: list[str] = []
    inside = threading.Event()
    release = threading.Event()

    def forcer() -> None:
        with default_backend("scipy"):
            inside.set()
            release.wait(timeout=10.0)

    t = threading.Thread(target=forcer)
    t.start()
    assert inside.wait(timeout=10.0)
    seen.append(get_default_backend())  # main thread, while forcer holds "scipy"
    release.set()
    t.join()
    assert seen == ["auto"]  # a process-global would have leaked "scipy"
    assert get_default_backend() == "auto"


def test_default_backend_override_restored_on_exception() -> None:
    before = get_default_backend()
    with pytest.raises(ValueError, match="boom"):
        with default_backend("scipy"):
            raise ValueError("boom")
    assert get_default_backend() == before


def test_explicit_backend_wins_over_default_override() -> None:
    """An explicit `backend="scipy"` ignores the override; only `"auto"` sites
    consult it."""
    A = _complex_symmetric(80, seed=34)
    with default_backend("scipy"):
        # explicit scipy honored (trivially), and the override does not force
        # anything an explicit arg already pins.
        assert SparseLU(A, backend="scipy").backend_used == "scipy"


def test_explicit_scipy_wins_over_mumps_default_override() -> None:
    """The DISTINGUISHING case: with the process default flipped to `"mumps"`,
    an explicit `backend="scipy"` still runs SuperLU -- proving the explicit arg
    overrides a NON-matching override, not just a coincidentally-equal one. Runs
    on any box: the explicit-scipy path never consults MUMPS availability, so
    this holds whether or not MUMPS is installed (the `"mumps"` override would
    only bite an `"auto"` site)."""
    A = _complex_symmetric(80, seed=36)
    with default_backend("mumps"):
        assert SparseLU(A, backend="scipy").backend_used == "scipy"


def test_set_default_backend_rejects_unknown_name() -> None:
    with pytest.raises(ValueError, match="unknown backend"):
        set_default_backend("nope")  # type: ignore[arg-type]


@pytest.mark.skipif(
    mumps_available(), reason="asserts MUMPS-ABSENT behaviour of the override (Mac only)"
)
def test_default_backend_mumps_override_raises_when_absent() -> None:
    """Forcing the default to `"mumps"` with MUMPS absent makes even an
    `"auto"` call site raise the clear error (the override resolves `"auto"`
    to `"mumps"`, which then fails loudly)."""
    A = _complex_symmetric(40, seed=35)
    with default_backend("mumps"):
        with pytest.raises(RuntimeError, match="MUMPS backend requested but not available"):
            SparseLU(A)


def test_refactor_scipy_matches_fresh_factorization() -> None:
    """refactor(A1) then solve == a fresh SparseLU(A1).solve, on the scipy path."""
    n = 200
    A0 = _complex_symmetric(n, seed=40)
    # A1 = A0 with a different diagonal shift -> SAME sparsity pattern
    A1 = (A0 + (2.0 + 1.0j) * sp.identity(n, dtype=complex)).tocsc()
    rng = np.random.default_rng(41)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    lu = SparseLU(A0, backend="scipy")
    lu.refactor(A1)
    x = lu.solve(b)
    x_fresh = SparseLU(A1, backend="scipy").solve(b)
    assert np.linalg.norm(x - x_fresh) / np.linalg.norm(x_fresh) < 1e-10
    assert np.linalg.norm(A1 @ x - b) / np.linalg.norm(b) < 1e-10


def test_refactor_rejects_pattern_mismatch() -> None:
    """reuse_analysis is only valid for an identical pattern -> guard raises."""
    A0 = _complex_symmetric(80, seed=42)
    B = _complex_symmetric(80, seed=43)  # different random pattern, same shape
    lu = SparseLU(A0, backend="scipy")
    with pytest.raises(ValueError, match="pattern"):
        lu.refactor(B)
    with pytest.raises(ValueError, match=r"shape|pattern"):
        lu.refactor(_complex_symmetric(70, seed=44))  # different shape too


def test_refactor_reuses_backend_and_symmetry() -> None:
    A0 = _complex_symmetric(100, seed=45)
    A1 = (A0 + (1.0 + 1.0j) * sp.identity(100, dtype=complex)).tocsc()
    lu = SparseLU(A0, backend="scipy")
    assert lu.symmetric is True
    lu.refactor(A1)
    assert lu.backend_used == "scipy" and lu.symmetric is True
