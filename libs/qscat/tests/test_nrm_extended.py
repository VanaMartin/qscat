"""Tests for the extended-space form of the nonlocal model (PRA 47 Eq. 2.1)."""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from qscat.core.grids import electronic_grid, segmented_grid
from qscat.core.nrm.discrete_state import AsymptoticDiscreteState
from qscat.core.nrm.extended import extended_hamiltonian
from qscat.core.nrm.ingredients import NrmIngredients, nrm_ingredients
from qscat.core.nrm.nonlocal_potential import continue_to_tail, nonlocal_operator
from qscat.dvr import kinetic
from qscat.model import F2
from scipy.sparse.linalg import spsolve

N_STATES = 3


@pytest.fixture(scope="module")
def nuc():
    return segmented_grid(
        ((9, 1.8), (1, 2.0), (5, 2.5), (4, 2.7), (20, 10.7)),
        ((1, 11.0), (1, 12.5), (1, 14.0), (3, 30.0)),
        angle_deg=45.0,
        quadrature=14,
    )


@pytest.fixture(scope="module")
def elec():
    return electronic_grid(r_max=16.0, order=8, n_complex=6)


@pytest.fixture(scope="module")
def ing(nuc, elec):
    phi_d = AsymptoticDiscreteState(elec, F2, R_inf=nuc.R0)
    real = nuc.points.imag == 0.0
    R = nuc.points[real].real[::-1]  # strictly descending, as nrm_ingredients requires
    return nrm_ingredients(elec, F2, phi_d, R)


def test_eliminating_the_arms_reproduces_the_nonlocal_operator(nuc, ing):
    """(E - H_ext)^-1 restricted to the d-block IS (E - T - V_d - F(E))^-1.

    This is the whole design: the extended space is a RESUMMATION of Eq.
    (2.1), not an approximation of it, so the two resolvents must agree to
    solver precision -- no physics tolerance is involved.
    """
    e_total = 0.05
    n_r = nuc.n
    h_ext = extended_hamiltonian(ing, nuc, F2, n_states=N_STATES)
    a_ext = (e_total * sp.identity(h_ext.shape[0], dtype=np.complex128) - h_ext).tocsc()

    rhs = np.zeros((h_ext.shape[0], n_r), dtype=np.complex128)
    rhs[:n_r, :] = np.eye(n_r, dtype=np.complex128)
    d_block = np.asarray(spsolve(a_ext, rhs))[:n_r, :]

    f = nonlocal_operator(ing, nuc, F2, e_total, n_states=N_STATES)
    v_d = continue_to_tail(ing.v_d_discrete, ing.R, nuc)
    a_ti = (
        e_total * np.eye(n_r, dtype=np.complex128)
        - kinetic(nuc, F2.mu)
        - np.diag(v_d)
        - f
    )
    want = np.linalg.inv(a_ti)
    # Scale-consistent: an EXACT identity, so gate on the achieved relative
    # error (measured 4.45e-14, six orders of headroom) rather than an
    # absolute atol that would pass an implementation ~1e5x worse than this
    # one on the fixture's large (~1.65e2) entries.
    assert np.linalg.norm(d_block - want) / np.linalg.norm(want) < 1e-10
    assert np.allclose(d_block, want, rtol=1e-8, atol=1e-12)


def test_the_matrix_is_complex_symmetric(nuc, ing):
    """ECS makes it symmetric, NOT Hermitian -- a `.conj().T` here is a bug."""
    h = extended_hamiltonian(ing, nuc, F2, n_states=N_STATES)
    assert abs(h - h.T).max() < 1e-12


def test_n_states_zero_drops_the_nonlocality(nuc, ing):
    """With no arms the block Hamiltonian is the bare T + V_d."""
    h = extended_hamiltonian(ing, nuc, F2, n_states=0)
    v_d = continue_to_tail(ing.v_d_discrete, ing.R, nuc)
    expected = kinetic(nuc, F2.mu) + np.diag(v_d)
    assert h.shape == (nuc.n, nuc.n)
    assert np.allclose(h.toarray(), expected, rtol=1e-10, atol=1e-12)


def test_rejects_more_states_than_available(nuc, ing):
    n_avail = ing.E_n.shape[1]
    with pytest.raises(ValueError, match="n_states"):
        extended_hamiltonian(ing, nuc, F2, n_states=n_avail + 1)


def test_rejects_ingredient_nodes_that_do_not_match_the_nuclear_grid(nuc):
    """`ing.R` must be exactly `nuclear_grid`'s real nodes -- mirrors
    `test_nrm_nonlocal_potential.py`'s equivalent guard, which
    `check_nodes_coincide` is now shared with (see `nonlocal_potential.py`).
    """
    real = nuc.points.imag == 0.0
    R_full = nuc.points[real].real[::-1]
    R_subsampled = R_full[::4]  # every 4th node -- mismatched with nuc's real nodes
    e_n = np.full((R_subsampled.size, 1), 0.2 - 0.01j, dtype=np.complex128)
    v_dn = np.zeros((R_subsampled.size, 1), dtype=np.complex128)
    mismatched = NrmIngredients(
        R=R_subsampled,
        v_d_discrete=F2.v0(R_subsampled).astype(np.complex128),
        E_n=e_n,
        V_dn=v_dn,
    )
    with pytest.raises(ValueError, match="coincide"):
        extended_hamiltonian(mismatched, nuc, F2)


def test_rejects_a_discrete_state_that_leaks_into_the_ecs_tail(nuc):
    """A discrete state that has NOT decoupled by R0 (Eq. 67) must be
    rejected -- mirrors `test_nrm_nonlocal_potential.py`'s equivalent guard,
    which `check_tail_coupling` is now shared with (see
    `nonlocal_potential.py`).
    """
    real = nuc.points.imag == 0.0
    R = nuc.points[real].real[::-1]
    e_n = np.full((R.size, 1), 0.2 - 0.01j, dtype=np.complex128)
    v_dn = np.full((R.size, 1), 0.05, dtype=np.complex128)  # flat, no decay
    leaky = NrmIngredients(R=R, v_d_discrete=F2.v0(R).astype(np.complex128), E_n=e_n, V_dn=v_dn)
    with pytest.raises(ValueError, match="tail"):
        extended_hamiltonian(leaky, nuc, F2)
