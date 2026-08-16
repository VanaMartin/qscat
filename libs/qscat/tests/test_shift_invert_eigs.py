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


def _n2_electronic(angle_deg: float, R: float = 2.02):
    """Dense and sparse builds of the SAME fixed-R electronic ECS Hamiltonian."""
    from qscat.core.grids import electronic_grid
    from qscat.dvr import kinetic, kinetic_sparse
    from qscat.model import N2

    g = electronic_grid(r_max=16.0, order=7, n_complex=6, angle_deg=angle_deg)
    V = np.asarray(N2.surface(g.points, R), dtype=np.complex128)
    H_dense = kinetic(g, 1.0) + np.diag(V)
    H_sparse = sp.csc_matrix(kinetic_sparse(g, 1.0) + sp.diags(V))
    return g, H_dense, H_sparse


# The N2 pole sits near the neutral curve: an ECS eigenvalue is ABSOLUTE, so it
# carries `v0(R)` (= -0.751 Ha at R = 2.02), while the literature quotes E_res
# RELATIVE to it. Hence a seed around -0.66 Ha, not around +0.09 Ha.
_N2_R = 2.02
_N2_SEED = -0.60 - 0.02j  # deliberately offset from the pole at -0.6613 - 0.0083j
_N2_WINDOW = (-1.0, 0.0, -0.1, 0.0)


def _sparse_two_angle_pole(k: int = 8) -> tuple[complex, float]:
    """The N2 electronic pole via the sparse path at both ECS angles.

    Note what this does NOT do: pick the pole out of one spectrum. A single
    shift-invert window contains the pole AND its neighbouring rotated-continuum
    eigenvalues, several of which are narrower in |Im E| than the pole itself.
    Only the two-angle criterion separates them.
    """
    from qscat.ecs import find_resonance_pole

    _, _, Hs_a = _n2_electronic(35.0)
    _, _, Hs_b = _n2_electronic(44.0)
    va, _ = ShiftInvertEigs(Hs_a, k=k).near(_N2_SEED)
    vb, _ = ShiftInvertEigs(Hs_b, k=k).near(_N2_SEED)
    return find_resonance_pole(va, vb, _N2_WINDOW)


def test_sparse_eigenvalues_match_dense_at_the_same_ecs_angle() -> None:
    """The differential oracle, per angle: every eigenvalue the sparse solver
    returns is one the dense solver returns, on the same matrix."""
    from qscat.dvr import eigen

    _, Hd, Hs = _n2_electronic(35.0)
    vals, _ = ShiftInvertEigs(Hs, k=8).near(_N2_SEED)
    w = eigen(Hd)[0]
    for val in vals:
        j = int(np.argmin(np.abs(w - val)))
        assert abs(w[j] - val) <= 1e-9 * abs(val) + 1e-12


def test_sparse_reproduces_the_dense_electronic_resonance_pole() -> None:
    """End to end: the two-angle pole built from sparse spectra equals the one
    built from dense spectra. Both are midpoints of the same two Hamiltonians,
    so they must agree to solver precision -- unlike a single-angle eigenvalue,
    which differs from the midpoint by half the (physical) angle residual."""
    from qscat.dvr import eigen
    from qscat.ecs import find_resonance_pole

    _, Hd_a, _ = _n2_electronic(35.0)
    _, Hd_b, _ = _n2_electronic(44.0)
    E_dense, res_dense = find_resonance_pole(eigen(Hd_a)[0], eigen(Hd_b)[0], _N2_WINDOW)
    assert res_dense < 1e-3  # a genuine angle-stable pole

    E_sparse, res_sparse = _sparse_two_angle_pole()
    assert abs(E_sparse - E_dense) <= 1e-9 * abs(E_dense) + 1e-12
    assert res_sparse == pytest.approx(res_dense, rel=1e-6)


def test_sparse_pole_carries_the_published_n2_resonance_parameters() -> None:
    """The pole the sparse path lands on is the physical one: E_res and Gamma
    fall inside the experimental band quoted in docs/physics/n2-resonance.md
    (E_res 2.3-2.5 eV, Gamma 0.35-0.55 eV; Schulz, Berman/Domcke).

    An ECS eigenvalue is ABSOLUTE, so `v0(R)` must be subtracted before
    comparing with a literature E_res, which is measured from the neutral curve.
    """
    from qscat.model import N2
    from qscat.units import HARTREE_TO_EV

    E_pole, _ = _sparse_two_angle_pole()
    e_res_ev = (E_pole.real - complex(N2.v0(_N2_R)).real) * HARTREE_TO_EV
    gamma_ev = -2.0 * E_pole.imag * HARTREE_TO_EV
    assert 2.3 < e_res_ev < 2.5
    assert 0.35 < gamma_ev < 0.55


@pytest.mark.parametrize("offset", [0.001, 0.01, 0.05, 0.2])
def test_pole_is_found_from_an_offset_shift(offset: float) -> None:
    """How far may the seed shift sit from the pole and still find it at k=8?
    The measured boundary is quoted in docs/physics/shift-invert-eigensolver.md."""
    from qscat.dvr import eigen
    from qscat.ecs import find_resonance_pole

    _, Hd_a, Hs_a = _n2_electronic(35.0)
    _, Hd_b, _ = _n2_electronic(44.0)
    E_pole, _ = find_resonance_pole(eigen(Hd_a)[0], eigen(Hd_b)[0], _N2_WINDOW)
    vals, _ = ShiftInvertEigs(Hs_a, k=8).near(E_pole + offset * (1.0 + 1.0j))
    assert np.min(np.abs(vals - E_pole)) <= 1e-6 * abs(E_pole)


def test_continuum_adjacent_shift_returns_angle_unstable_eigenvalues() -> None:
    """A shift parked on the rotated continuum returns eigenvalues that MOVE
    when the ECS angle changes, while a pole does not.

    This is what makes a shift-invert window on its own useless as a resonance
    finder: the window is a local slice of the spectrum, and being in it says
    nothing about being a pole. Only comparing two ECS angles does."""
    from qscat.ecs import match_angle_stable

    _, _, Hs_a = _n2_electronic(35.0)
    _, _, Hs_b = _n2_electronic(44.0)
    sigma = -0.30 - 0.40j  # deep in the rotated continuum, far from the pole
    va, _ = ShiftInvertEigs(Hs_a, k=10).near(sigma)
    vb, _ = ShiftInvertEigs(Hs_b, k=10).near(sigma)
    stable, _, _ = match_angle_stable(va, vb, (-1.5, 0.5, -1.0, 0.0), rel_tol=1e-4, atol=1e-8)
    assert stable.size == 0  # continuum: nothing is angle-stable


def test_sparse_eigenvectors_match_dense_on_the_electronic_hamiltonian() -> None:
    from qscat.dvr import eigen

    _, Hd, Hs = _n2_electronic(35.0)
    vals, vecs = ShiftInvertEigs(Hs, k=4).near(_N2_SEED)
    w, V = eigen(Hd)
    for i, val in enumerate(vals):
        j = int(np.argmin(np.abs(w - val)))
        assert abs(w[j] - val) <= 1e-9 * abs(val) + 1e-12
        u = V[:, j] / np.sqrt(c_product(V[:, j], V[:, j]))
        v = vecs[:, i] / np.sqrt(c_product(vecs[:, i], vecs[:, i]))
        assert abs(abs(c_product(u, v)) - 1.0) < 1e-6
