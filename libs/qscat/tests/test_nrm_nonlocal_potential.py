"""Tests for the nonlocal potential F(E,R,R') (PRA 77 Eq. 60-61)."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import segmented_grid
from qscat.core.nrm.ingredients import NrmIngredients
from qscat.core.nrm.nonlocal_potential import continue_to_tail, nonlocal_operator
from qscat.dvr import kinetic
from qscat.model import F2


@pytest.fixture(scope="module")
def nuc():
    return segmented_grid(
        ((6, 2.0), (8, 6.0)), ((1, 7.0), (2, 16.0)), angle_deg=45.0, quadrature=10
    )


def _ingredients(nuc, n_states=3, coupling=0.05):
    """A hand-built ingredient set: flat couplings decaying to ~0 by R0
    (Eq. 67 -- a properly decaying discrete state has decoupled by the ECS
    boundary), well-separated E_n. The decay is what keeps
    `nonlocal_operator`'s tail-coupling guard satisfied: a genuinely flat
    coupling extended into the ECS tail is exactly the ill-posed input that
    guard exists to reject (see `test_rejects_a_discrete_state_that_leaks_
    into_the_ecs_tail` below).
    """
    real = nuc.points.imag == 0.0
    R = nuc.points[real].real[::-1]  # descending
    e_n = np.tile(np.array([0.2 - 0.01j, 0.5 - 0.02j, 0.9 - 0.03j])[:n_states], (R.size, 1))
    r0 = float(R[0])
    decay_start = r0 - 1.0
    t = np.clip((R - decay_start) / (r0 - decay_start), 0.0, 1.0)
    envelope = (1.0 - t) ** 6  # 1 in the interior, exactly 0 at R0
    v_dn = coupling * envelope[:, None] * np.ones((1, n_states), dtype=np.complex128)
    return NrmIngredients(
        R=R,
        v_d_discrete=F2.v0(R).astype(np.complex128),
        E_n=e_n,
        V_dn=v_dn,
    )


def test_shape_and_symmetry(nuc):
    """F is complex SYMMETRIC (not Hermitian) -- H = H^T under ECS."""
    f = nonlocal_operator(_ingredients(nuc), nuc, F2, e_total=0.05)
    assert f.shape == (nuc.n, nuc.n)
    assert np.allclose(f, f.T, rtol=1e-10, atol=1e-14)


def test_scales_quadratically_with_the_coupling(nuc):
    """F is bilinear in V_dn: doubling the coupling quadruples F."""
    f1 = nonlocal_operator(_ingredients(nuc, coupling=0.05), nuc, F2, e_total=0.05)
    f2 = nonlocal_operator(_ingredients(nuc, coupling=0.10), nuc, F2, e_total=0.05)
    assert np.allclose(f2, 4.0 * f1, rtol=1e-10, atol=1e-14)


def test_zero_coupling_gives_zero(nuc):
    f = nonlocal_operator(_ingredients(nuc, coupling=0.0), nuc, F2, e_total=0.05)
    assert np.all(f == 0.0)


def test_states_add_independently(nuc):
    """Eq. (60) is a plain sum over n -- the 3-state F is the sum of three."""
    total = nonlocal_operator(_ingredients(nuc, n_states=3), nuc, F2, e_total=0.05)
    parts = sum(
        nonlocal_operator(_ingredients(nuc, n_states=3), nuc, F2, e_total=0.05, n_states=k)
        - nonlocal_operator(_ingredients(nuc, n_states=3), nuc, F2, e_total=0.05, n_states=k - 1)
        for k in (1, 2, 3)
    )
    assert np.allclose(total, parts, rtol=1e-10, atol=1e-14)


def test_single_state_matches_an_explicit_green_solve(nuc):
    """One state, built term by term against Eq. (60)-(61) directly.

    `v_dn`/`e_n` are mapped onto the full nuclear grid via `continue_to_tail`
    (the module's own, separately-exercised, ECS-tail continuation rule) so
    this test's expected value tracks whatever coupling profile `_ingredients`
    uses; the Green's-function sandwich itself (`diag @ inv @ diag`) is built
    independently of `nonlocal_operator`.
    """
    ing = _ingredients(nuc, n_states=1, coupling=0.05)
    e_tot = 0.05
    f = nonlocal_operator(ing, nuc, F2, e_total=e_tot, n_states=1)

    v_dn = continue_to_tail(ing.V_dn[:, 0], ing.R, nuc)
    e_n = continue_to_tail(ing.E_n[:, 0], ing.R, nuc)
    m = (
        e_tot * np.eye(nuc.n, dtype=np.complex128)
        - kinetic(nuc, F2.mu)
        - np.diag(F2.v0(nuc.points) + e_n)
    )
    expected = np.diag(v_dn) @ np.linalg.inv(m) @ np.diag(v_dn)
    assert np.allclose(f, expected, rtol=1e-9, atol=1e-14)


def test_rejects_more_states_than_available(nuc):
    with pytest.raises(ValueError, match="n_states"):
        nonlocal_operator(_ingredients(nuc, n_states=3), nuc, F2, e_total=0.05, n_states=99)


def test_rejects_a_discrete_state_that_leaks_into_the_ecs_tail(nuc):
    """A discrete state that has NOT decoupled by R0 (Eq. 67) must be
    rejected, not silently smeared across the ECS tail as a spurious
    long-range coupling. This is the guard `_ingredients`'s decay envelope
    is built to satisfy; this test pins the guard actually firing when it
    doesn't -- undetected before this test, since the brief's own reference
    fixture (flat coupling, no decay) tripped it on every other test here.
    """
    real = nuc.points.imag == 0.0
    R = nuc.points[real].real[::-1]
    e_n = np.full((R.size, 1), 0.2 - 0.01j, dtype=np.complex128)
    v_dn = np.full((R.size, 1), 0.05, dtype=np.complex128)  # flat, no decay
    ing = NrmIngredients(R=R, v_d_discrete=F2.v0(R).astype(np.complex128), E_n=e_n, V_dn=v_dn)
    with pytest.raises(ValueError, match="tail"):
        nonlocal_operator(ing, nuc, F2, e_total=0.05)


def test_continue_to_tail_maps_by_position_not_by_array_order(nuc):
    """Pins `continue_to_tail`'s OUTPUT against an independently built array
    -- not one constructed by calling `continue_to_tail` itself -- using an
    R-ASYMMETRIC profile.

    A flat or symmetric-under-reversal profile cannot distinguish a correct
    nearest-R mapping from one that copies values by array position (i.e.
    end-for-end mirrored: R0's value lands at R~0 and vice versa). Every
    other test in this file builds its "expected" arrays either from a flat
    profile or by calling `continue_to_tail` itself, so none of them would
    catch that mistake -- this test is the one that does.
    """
    real_mask = nuc.points.imag == 0.0
    pts_real = nuc.points[real_mask].real
    assert np.array_equal(pts_real, np.sort(pts_real))  # grid convention: ascending

    R = pts_real[::-1]  # descending, R[0] = R0 (largest)
    values = np.arange(R.size, dtype=np.complex128)  # strictly monotonic in R -> asymmetric

    out = continue_to_tail(values, R, nuc)

    # Built independently of continue_to_tail: pts_real is ascending, so its
    # i-th point is R's (size-1-i)-th (descending) point -- plain reversal.
    expected_real = values[::-1]
    np.testing.assert_array_equal(out[real_mask], expected_real)
    # Tail extends at R0's value, i.e. R's index 0 (the largest R).
    np.testing.assert_array_equal(out[~real_mask], values[0])


def test_rejects_ingredient_nodes_that_do_not_match_the_nuclear_grid(nuc):
    """`ing.R` must be exactly `nuclear_grid`'s real nodes -- a subsampled or
    otherwise mismatched set is rejected rather than silently
    piecewise-constant-continued (E_n in particular has no self-limiting
    guard the way V_dn's tail check provides).
    """
    real = nuc.points.imag == 0.0
    R_full = nuc.points[real].real[::-1]
    R_subsampled = R_full[::4]  # every 4th node -- mismatched with nuc's real nodes
    e_n = np.full((R_subsampled.size, 1), 0.2 - 0.01j, dtype=np.complex128)
    v_dn = np.zeros((R_subsampled.size, 1), dtype=np.complex128)
    ing = NrmIngredients(
        R=R_subsampled,
        v_d_discrete=F2.v0(R_subsampled).astype(np.complex128),
        E_n=e_n,
        V_dn=v_dn,
    )
    with pytest.raises(ValueError, match="coincide"):
        nonlocal_operator(ing, nuc, F2, e_total=0.05)


def test_single_state_matches_an_independent_loop_reference(nuc):
    """A SECOND, independently-coded check that the weights are not
    reapplied -- deleting `test_single_state_matches_an_explicit_green_solve`
    must not silently restore the sqrt(W) bug. Builds F entry-by-entry with
    plain Python loops (no vectorized broadcasting, and an asymmetric
    coupling profile distinct from `_ingredients`'s), so it does not share
    an implementation mistake with either `nonlocal_operator` or the other
    comparator test.
    """
    real = nuc.points.imag == 0.0
    R = nuc.points[real].real[::-1]
    r0 = float(R[0])
    # Asymmetric, non-flat coupling: linear ramp from 0 at R0 to 0.08 at the
    # innermost node, distinct in shape from _ingredients' envelope.
    v_dn_col = 0.08 * np.clip((r0 - R) / r0, 0.0, 1.0)
    e_n_col = np.full(R.size, 0.3 - 0.015j, dtype=np.complex128)
    v_dn_col = v_dn_col.astype(np.complex128)

    ing = NrmIngredients(
        R=R,
        v_d_discrete=F2.v0(R).astype(np.complex128),
        E_n=e_n_col[:, None],
        V_dn=v_dn_col[:, None],
    )
    e_tot = 0.07
    f = nonlocal_operator(ing, nuc, F2, e_total=e_tot, n_states=1)

    v_dn = continue_to_tail(v_dn_col, R, nuc)
    e_n = continue_to_tail(e_n_col, R, nuc)
    t = kinetic(nuc, F2.mu)
    v0 = F2.v0(nuc.points)
    m = np.empty((nuc.n, nuc.n), dtype=np.complex128)
    for i in range(nuc.n):
        for j in range(nuc.n):
            m[i, j] = (e_tot if i == j else 0.0) - t[i, j] - ((v0[i] + e_n[i]) if i == j else 0.0)
    g = np.linalg.inv(m)
    expected = np.empty((nuc.n, nuc.n), dtype=np.complex128)
    for i in range(nuc.n):
        for j in range(nuc.n):
            expected[i, j] = v_dn[i] * g[i, j] * v_dn[j]  # no sqrt(W) anywhere

    assert np.allclose(f, expected, rtol=1e-9, atol=1e-14)
