"""Tests for the extended-space form of the nonlocal model (PRA 47 Eq. 2.1)."""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from qscat.core.grids import electronic_grid, segmented_grid
from qscat.core.nrm.coupling import v_dk_plus
from qscat.core.nrm.discrete_state import AsymptoticDiscreteState
from qscat.core.nrm.dissociation import _boundary_node
from qscat.core.nrm.extended import extended_hamiltonian
from qscat.core.nrm.ingredients import NrmIngredients, nrm_ingredients
from qscat.core.nrm.nonlocal_potential import continue_to_tail, nonlocal_operator
from qscat.core.vibrational import vibrational_states
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
    a_ti = e_total * np.eye(n_r, dtype=np.complex128) - kinetic(nuc, F2.mu) - np.diag(v_d) - f
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


def test_launch_state_drives_the_time_independent_solve(nuc, elec, ing):
    """The launch state reconstructs Eq. (52)'s right-hand side.

    Feed the reconstruction to the shipped `solve_nuclear` and the shipped
    sigma_DA must come back -- a tautology-free check that the launch basis
    matches the TI route, not just that it round-trips through itself.
    """
    from qscat.core.nrm.dissociation import (
        da_sigma_from_psi,
        nrm_da_cross_section,
        solve_nuclear,
    )
    from qscat.core.nrm.extended import initial_packet

    phi_d = AsymptoticDiscreteState(elec, F2, R_inf=nuc.R0)
    eps, chi = vibrational_states(nuc, F2.mu, 4, F2.v0)
    e_kin, v_init = 0.03, 0

    launch = initial_packet(
        nuc,
        elec,
        F2,
        phi_d,
        ing,
        eps,
        chi,
        v_init,
        np.array([e_kin]),
        n_states=N_STATES,
        rank_tol=1e-10,
    )
    assert launch.vectors.shape[0] == (1 + N_STATES) * nuc.n
    assert np.allclose(launch.vectors[nuc.n :, :], 0.0)  # the arms start empty

    xi = (launch.vectors @ launch.coeffs)[: nuc.n, 0]
    e_total = e_kin + float(eps[v_init])
    f = nonlocal_operator(ing, nuc, F2, e_total, n_states=N_STATES)
    v_d = continue_to_tail(ing.v_d_discrete, ing.R, nuc)
    psi_d = solve_nuclear(nuc, F2.mu, v_d, f, xi, e_total)

    # Same expression `nrm_da_cross_section` uses for the anion electronic
    # threshold -- not re-derived here.
    v_d_full = continue_to_tail(ing.v_d_discrete, ing.R, nuc)
    eps_e = float(v_d_full[_boundary_node(nuc)].real)
    got = da_sigma_from_psi(nuc, F2.mu, psi_d, e_total, eps_e, e_kin)
    want = float(
        nrm_da_cross_section(
            nuc, elec, F2, phi_d, eps, chi, v_init, e_kin, ingredients=ing, n_states=N_STATES
        )
    )
    assert got == pytest.approx(want, rel=1e-10)


def test_the_launch_matrix_is_numerically_low_rank(nuc, elec, ing):
    """The economy of Task 3 rests on this; measure it, do not assume it.

    `truncation_error` (`sigma_{r+1}/sigma_1`) bounds error relative to the
    LARGEST column, not to each energy's own -- see
    `LaunchBasis.truncation_error`'s docstring. Task 3's accuracy depends on
    the actual PER-COLUMN reconstruction error, so gate on that directly
    rather than trusting `truncation_error` as a per-energy bound (it can't
    be one: `keep` is chosen so `sv[keep] <= rank_tol * sv[0]` by
    construction, so `truncation_error < rank_tol` holds regardless of
    whether the truncation is actually accurate).
    """
    from qscat.core.nrm.extended import initial_packet

    phi_d = AsymptoticDiscreteState(elec, F2, R_inf=nuc.R0)
    eps, chi = vibrational_states(nuc, F2.mu, 4, F2.v0)
    energies = np.linspace(0.010, 0.050, 9)
    launch = initial_packet(
        nuc, elec, F2, phi_d, ing, eps, chi, 0, energies, n_states=N_STATES, rank_tol=1e-6
    )
    assert launch.rank <= 4
    assert launch.truncation_error < 1e-6

    # Exact (rank_tol=0.0) reconstruction as the per-column reference --
    # test_full_rank_reconstructs_every_column_exactly already gates that it
    # matches the direct v_dk_plus construction to round-off, so comparing
    # the truncated reconstruction against it (rather than against a second
    # direct build) isolates the truncation's own error.
    exact = initial_packet(
        nuc, elec, F2, phi_d, ing, eps, chi, 0, energies, n_states=N_STATES, rank_tol=0.0
    )
    m_exact = (exact.vectors @ exact.coeffs)[: nuc.n, :]
    m_trunc = (launch.vectors @ launch.coeffs)[: nuc.n, :]
    col_err = np.linalg.norm(m_trunc - m_exact, axis=0) / np.linalg.norm(m_exact, axis=0)
    # Measured worst-column error 1.15e-6, 2.2x above `truncation_error`
    # (5.26e-7) -- gated with headroom above the measurement, not at the
    # `truncation_error` value (which would pass trivially, per this test's
    # own docstring above).
    assert col_err.max() < 5e-6


def test_full_rank_reconstructs_every_column_exactly(nuc, elec, ing):
    """`rank_tol=0.0` must reproduce the per-energy launch vectors to
    round-off, so the truncation is the ONLY approximation the factorization
    introduces."""
    from qscat.core.nrm.extended import initial_packet

    phi_d = AsymptoticDiscreteState(elec, F2, R_inf=nuc.R0)
    eps, chi = vibrational_states(nuc, F2.mu, 4, F2.v0)
    energies = np.array([0.015, 0.025, 0.035])
    launch = initial_packet(
        nuc, elec, F2, phi_d, ing, eps, chi, 0, energies, n_states=N_STATES, rank_tol=0.0
    )
    assert launch.rank == energies.size
    assert launch.truncation_error == 0.0

    recon = launch.vectors @ launch.coeffs
    for j, e_kin in enumerate(energies):
        direct = (
            continue_to_tail(v_dk_plus(elec, F2, phi_d, ing.R, float(e_kin)), ing.R, nuc) * chi[0]
        )
        got = recon[: nuc.n, j]
        assert np.linalg.norm(got - direct) / np.linalg.norm(direct) < 1e-10
