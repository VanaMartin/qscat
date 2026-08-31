"""Spec gate 1, the one that replaces Houfek's certification: at s = 0 the
coupled Hamiltonian is block diagonal and its l = 1 block IS the shipped
model's Hamiltonian."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from qscat.core.grids import segmented_grid
from qscat.dvr import FemDvrEcsGrid, TensorGrid
from qscat.model import NO

from projects.no_coupled_channels.anisotropy import TwoCentreWell
from projects.no_coupled_channels.model import CoupledModel, DiagonalChannelModel

N_CH = 4
# Round-off bound, RELATIVE to the reference matrix magnitude. Absolute
# bounds are wrong here: the centrifugal term l(l+1)/(2 r^2) makes the
# largest matrix entries grow without limit as the first radial node
# approaches the origin, so an absolute tolerance silently tightens or
# loosens with the grid. Measured on the grid below, both differences sit
# at ~2e-17 relative; 1e-14 leaves ample headroom while staying far tighter
# than any real coupling (~1e-1 relative).
RTOL = 1e-14


def _electronic() -> FemDvrEcsGrid:
    """A deliberately small electronic grid. `electronic_grid` cannot be
    used: it hardcodes inner segments out to 10 bohr and rejects any
    `r_max` below that, so the smallest grid it can build is far larger
    than an identity test needs."""
    return segmented_grid(((4, 8.0),), ((2, 12.0),), angle_deg=35.0, quadrature=6)


def _tensor_grid() -> TensorGrid:
    """A deliberately small 2-D grid: this test is about identity, not
    physics. 29 x 24 = 696 points."""
    nu = segmented_grid(((3, 4.0),), ((2, 6.0),), angle_deg=30.0, quadrature=6, x_min=1.0)
    return TensorGrid([_electronic(), nu])


def _magnitude(H: sp.csr_matrix) -> float:
    """Largest entry of `H` -- the scale a round-off claim is relative to."""
    return float(np.max(np.abs(H.data)))


def _block(H: sp.csr_matrix, i: int, j: int, n: int) -> sp.csr_matrix:
    return sp.csr_matrix(H[i * n : (i + 1) * n, j * n : (j + 1) * n])


def test_s0_hamiltonian_is_block_diagonal() -> None:
    model = CoupledModel(well=TwoCentreWell(base=NO, s=0.0, kappa=0.3), n_channels=N_CH)
    tg = _tensor_grid()
    H = model.hamiltonian(tg)
    n = tg.size
    bound = RTOL * _magnitude(sp.csr_matrix(NO.hamiltonian(tg)))
    for i in range(N_CH):
        for j in range(N_CH):
            if i != j:
                blk = _block(H, i, j, n)
                assert blk.nnz == 0 or float(np.max(np.abs(blk.data))) < bound


def test_s0_first_block_is_the_shipped_model_exactly() -> None:
    model = CoupledModel(well=TwoCentreWell(base=NO, s=0.0, kappa=0.3), n_channels=N_CH)
    tg = _tensor_grid()
    n = tg.size
    got = _block(model.hamiltonian(tg), 0, 0, n)
    ref = sp.csr_matrix(NO.hamiltonian(tg))
    diff = (got - ref).tocoo()
    assert diff.nnz == 0 or float(np.max(np.abs(diff.data))) < RTOL * _magnitude(ref)


def test_one_channel_is_the_fixed_l_model_through_the_same_code() -> None:
    """n_channels = 1 must be the SAME assembly path, not a special case."""
    well = TwoCentreWell(base=NO, s=0.0, kappa=0.3)
    tg = _tensor_grid()
    one = CoupledModel(well=well, n_channels=1).hamiltonian(tg)
    ref = sp.csr_matrix(NO.hamiltonian(tg))
    diff = (one - ref).tocoo()
    assert diff.nnz == 0 or float(np.max(np.abs(diff.data))) < RTOL * _magnitude(ref)


def test_electronic_hamiltonian_is_complex_symmetric() -> None:
    model = CoupledModel(well=TwoCentreWell(base=NO, s=1.0, kappa=0.3), n_channels=3)
    g = _electronic()
    H = model.electronic_hamiltonian(g, 2.3 + 0.0j).toarray()
    np.testing.assert_allclose(H, H.T, atol=1e-13)


def test_coupling_is_present_once_the_anisotropy_is_on() -> None:
    model = CoupledModel(well=TwoCentreWell(base=NO, s=1.0, kappa=0.3), n_channels=2)
    g = _electronic()
    n = g.n
    H = sp.csr_matrix(model.electronic_hamiltonian(g, 2.3 + 0.0j))
    off = _block(H, 0, 1, n)
    assert off.nnz > 0
    assert float(np.max(np.abs(off.data))) > 1e-6


def test_coupling_table_matches_potential_nd_c_order() -> None:
    """`_coupling_table`'s off-diagonal flattening must use the SAME C-order
    (last-axis-fastest) convention `potential_nd` uses for the diagonal
    blocks -- checked on `_tensor_grid`'s 29 x 24 grid, where the two axes
    have UNEQUAL size.

    That is deliberate, not incidental: `potential_nd`'s own docstring warns
    that a transposed potential is invisible whenever the two axes happen to
    have the SAME size (`np.meshgrid`'s default `"xy"` indexing swaps the
    first two axes, and a same-size swap produces a same-shape array that
    looks fine). `_coupling_table` builds its off-diagonal blocks from
    `well.v_block(l, lp, r, R)` on the same broadcastable `r, R` that
    `TensorGrid.points()` hands to `potential_nd`, then ravels with the same
    default (C) order -- this test is the one place that equivalence is
    checked directly, on a grid shaped so a transposition could not hide.
    """
    from qscat.dvr import potential_nd

    well = TwoCentreWell(base=NO, s=0.3, kappa=0.3)
    model = CoupledModel(well=well, n_channels=2)
    tg = _tensor_grid()
    assert len(set(tg.shape)) == tg.ndim, f"grid axes must be unequal, got {tg.shape}"
    r, R = tg.points()
    table = model._coupling_table(r, R)
    l0, l1 = model.channel_ells()
    expected = potential_nd(tg, lambda rr, RR: well.v_block(l0, l1, rr, RR))
    np.testing.assert_allclose(table[0][1], expected, rtol=0, atol=1e-14)


def test_diagonal_channel_model_matches_the_shipped_surface_at_s0() -> None:
    """DiagonalChannelModel(l=1) at s = 0 must be NO itself -- that identity is
    what makes it a valid cross-check against local_complex_potential."""
    dm = DiagonalChannelModel(well=TwoCentreWell(base=NO, s=0.0, kappa=0.3), l=1)
    r = np.linspace(0.4, 10.0, 37)[:, None]
    R = np.array([1.9, 2.4, 3.1])[None, :]
    np.testing.assert_allclose(dm.surface(r, R), NO.surface(r, R), rtol=0, atol=1e-14)
    assert dm.ell == NO.ell
    assert dm.mu == NO.mu
    assert dm.charge == NO.charge


def test_interaction_matrix_is_the_perturbation_alone() -> None:
    """At s = 0 the coupled interaction must be block-diagonal, and its l = 1
    block must be exactly the shipped model's interaction -- v0 and the
    centrifugal term belong to the free Hamiltonian, not to the perturbation."""
    well = TwoCentreWell(base=NO, s=0.0, kappa=0.3)
    tg = _tensor_grid()
    n = tg.size
    V = sp.csr_matrix(CoupledModel(well=well, n_channels=3).interaction_matrix(tg))
    assert V.shape == (3 * n, 3 * n)

    # RELATIVE to the interaction's own magnitude, not absolute. The
    # off-diagonal blocks vanish by quadrature exactness, so what survives is
    # round-off on the potential -- measured, 2.9e-15 against a magnitude of
    # 6.0, i.e. 4.8e-16 relative. An absolute bound would silently tighten or
    # loosen with lambda(R), and on a weaker interaction it would start
    # rejecting a correct implementation.
    scale = float(np.max(np.abs(NO.interaction_diag(tg))))
    for i in range(3):
        for j in range(3):
            if i != j:
                blk = _block(V, i, j, n)
                assert blk.nnz == 0 or float(np.max(np.abs(blk.data))) < 1e-13 * scale

    got = np.asarray(_block(V, 0, 0, n).diagonal())
    np.testing.assert_allclose(got, NO.interaction_diag(tg), rtol=0, atol=1e-14)


def test_interaction_matrix_is_complex_symmetric() -> None:
    """Every operator here is complex symmetric, never Hermitian -- the ECS
    contour makes it so, and a Hermitian-only routine downstream would be
    silently wrong rather than loudly."""
    well = TwoCentreWell(base=NO, s=0.6, kappa=0.5)
    tg = _tensor_grid()
    V = sp.csr_matrix(CoupledModel(well=well, n_channels=3).interaction_matrix(tg))
    diff = (V - V.T).tocoo()
    assert diff.nnz == 0 or float(np.max(np.abs(diff.data))) < 1e-14


def test_interaction_matrix_couples_once_the_anisotropy_is_on() -> None:
    well = TwoCentreWell(base=NO, s=1.0, kappa=0.3)
    tg = _tensor_grid()
    n = tg.size
    V = sp.csr_matrix(CoupledModel(well=well, n_channels=2).interaction_matrix(tg))
    off = _block(V, 0, 1, n)
    assert off.nnz > 0
    assert float(np.max(np.abs(off.data))) > 1e-6
