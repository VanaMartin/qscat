"""`sigma_DA(E)` by propagation, against the per-energy solve (PRA 77 Eq. 54).

The gate is `test_td_da_matches_the_time_independent_cross_section`; the rest
cover the conventions the gate cannot separate (closed channels, scalar
shape, the truncated-transform warning) at a cost that does not need
`@slow`.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import electronic_grid, segmented_grid
from qscat.core.lcp import lcp_da_cross_section, local_complex_potential
from qscat.core.nrm.discrete_state import AsymptoticDiscreteState
from qscat.core.nrm.dissociation import nrm_da_cross_section, solve_nuclear
from qscat.core.nrm.ingredients import nrm_ingredients
from qscat.core.nrm.nonlocal_potential import continue_to_tail
from qscat.core.nrm.td_cross_section import (
    td_nrm_da_cross_section,
    td_nrm_ve_cross_section,
)
from qscat.core.nrm.vibrational_excitation import nrm_ve_cross_section, t_resonant
from qscat.core.vibrational import vibrational_states
from qscat.model import F2, N2

# The eMoScat F2 NUCLEAR deck (`reference/eMoScat/input/F2/grids.txt`, 2nd
# declaration), transcribed as literals for the same reason
# `test_nrm_propagation.py` transcribes its N2 one: `libs/qscat` must not
# import `validation` or `apps`. `validation/diatomic/config.py`'s
# `MoleculeConfig` and `qscat_run.presets` carry the same numbers, locked to
# each other by `validation/diatomic/test_da_grid.py`.
#
# It has to be THIS deck and not a hand-rolled coarse one. F2's dissociating
# wave carries `K_R = 55.6-64.2` over `E = 0.02-0.05 Ha` (wavelength
# 0.098-0.113 bohr); on a coarse N2-style nuclear grid the packet cannot
# propagate out at all -- measured 2026-08-19, the centroid creeps
# 2.66 -> 2.77 bohr and then oscillates while `S(t)` plateaus at 0.673, i.e.
# grid-trapped flux masquerading as a bound component. This deck carries 65
# points/bohr over `R = 3-10.7` (largest node spacing 0.023 bohr) and on it
# the centroid climbs monotonically while `<P>_t` rises to `K_R`.
_F2_REAL = ((9, 1.8), (1, 2.0), (5, 2.5), (4, 2.596908), (4, 2.7), (40, 10.7))
_F2_COMPLEX = (
    (1, 10.8),
    (1, 11.0),
    (1, 11.5),
    (1, 12.5),
    (1, 14.0),
    (1, 18.0),
    (4, 30.0),
    (2, 101.0),
)


@pytest.fixture(scope="module")
def f2_deck():
    """F2 on the production NUCLEAR deck (974 points) with a REDUCED
    electronic grid (55 points against the production 132).

    The nuclear deck cannot be reduced -- it is what makes the packet leave
    (comment above). The electronic one can: the comparison this file makes
    is DIFFERENTIAL, both routes running on the same `ing` on the same grids,
    so it does not require a physically converged electronic basis, only one
    on which `Psi_d^+` is well defined. What it costs to keep the production
    electronic deck is the whole budget: `H_ext` is `(1 + n_states) * N_R`,
    so 132 electronic points give a 128568-square matrix at 0.83 s/step
    against 53570 and 0.30 s/step here (measured 2026-08-19, 12-core dev
    machine), and the gate needs thousands of steps either way.

    What the reduction does NOT move is the physics the extraction depends
    on: `eps_e` reads -0.126942 here against -0.126931 on the production
    electronic deck (1.1e-5 Ha), and `min V_d` -0.149272 against -0.149263.
    What it DOES move is the absolute cross section -- `sigma_TI` reads
    [3.443, 1.559, 0.296] here against [3.352, 1.682, 0.279] at
    `order=6, n_complex=3, r_max=16` -- so these numbers are NOT the
    converged F2 cross section and must not be quoted as one;
    `validation/diatomic` owns that. The TD/TI RATIO is what this fixture
    measures, and it is insensitive to the reduction: the two decks agree on
    it to 3-4 significant figures at every checkpoint of the convergence
    scan (e.g. 0.6609/0.7648/0.4351 against 0.6506/0.7607/0.4413 at
    T=3000).

    `n_states` is left at `None` throughout: the COMPLETE projected-state sum
    is a correctness requirement of the time-dependent route, not a
    performance knob (`propagation._warn_if_diverging`).
    """
    nuc = segmented_grid(_F2_REAL, _F2_COMPLEX, angle_deg=35.0, quadrature=14)
    elec = electronic_grid(r_max=13.0, order=5, n_complex=2)
    phi_d = AsymptoticDiscreteState(elec, F2, R_inf=nuc.R0)
    # nrm_ingredients requires strictly DESCENDING R.
    r_values = nuc.points[nuc.points.imag == 0.0].real[::-1]
    ing = nrm_ingredients(elec, F2, phi_d, r_values)
    eps, chi = vibrational_states(nuc, F2.mu, 4, F2.v0)
    return nuc, elec, phi_d, ing, eps, chi


@pytest.fixture(scope="module")
def small_deck():
    """A deliberately small F2 deck for the convention tests.

    Nothing here is physically converged and nothing here is compared
    against the time-independent route -- these tests check argument
    handling, shape, and the warning, all of which are grid-independent.
    """
    nuc = segmented_grid(((3, 2.5), (4, 6.0)), ((2, 20.0),), angle_deg=35.0, quadrature=8)
    elec = electronic_grid(r_max=11.0, order=6, n_complex=2)
    phi_d = AsymptoticDiscreteState(elec, F2, R_inf=nuc.R0)
    r_values = nuc.points[nuc.points.imag == 0.0].real[::-1]
    ing = nrm_ingredients(elec, F2, phi_d, r_values)
    eps, chi = vibrational_states(nuc, F2.mu, 3, F2.v0)
    return nuc, elec, phi_d, ing, eps, chi


def test_non_positive_energies_are_zero_like_the_time_independent_route(small_deck):
    """`E <= 0` gives `0.0` on both routes -- and costs no propagation.

    `initial_packet` RAISES on a non-positive energy while
    `nrm_da_cross_section` silently returns zero, so a batch mixing open and
    closed entries is exactly where the two routes could stop being a
    one-identifier substitution.
    """
    nuc, elec, phi_d, ing, eps, chi = small_deck
    kw = {"ingredients": ing, "dt": 1.0, "n_steps": 4}
    got = td_nrm_da_cross_section(nuc, elec, F2, phi_d, eps, chi, 0, np.array([-0.01, 0.0]), **kw)
    assert got.shape == (2,)
    assert np.array_equal(got, np.zeros(2))


def test_a_closed_entry_does_not_disturb_the_open_ones(small_deck):
    """The open columns of a mixed batch are what they are on their own.

    `open_idx` has to map propagation columns back onto the OUTPUT slots; an
    off-by-one there would shift every cross section by an energy and still
    return finite, positive, plausible numbers.
    """
    nuc, elec, phi_d, ing, eps, chi = small_deck
    kw = {"ingredients": ing, "dt": 1.0, "n_steps": 6}
    with pytest.warns(UserWarning):
        mixed = td_nrm_da_cross_section(
            nuc, elec, F2, phi_d, eps, chi, 0, np.array([-0.01, 0.03, 0.05]), **kw
        )
        alone = td_nrm_da_cross_section(
            nuc, elec, F2, phi_d, eps, chi, 0, np.array([0.03, 0.05]), **kw
        )
    assert mixed[0] == 0.0
    assert np.allclose(mixed[1:], alone, rtol=1e-10)


def test_a_scalar_energy_returns_a_scalar(small_deck):
    nuc, elec, phi_d, ing, eps, chi = small_deck
    with pytest.warns(UserWarning):
        got = td_nrm_da_cross_section(
            nuc, elec, F2, phi_d, eps, chi, 0, 0.03, ingredients=ing, dt=1.0, n_steps=4
        )
    assert np.ndim(got) == 0


def test_a_truncated_transform_is_warned_about(small_deck):
    """Four steps is not a converged transform, and the caller is told so.

    The failure mode this guards is not a blow-up: a packet that has not
    finished leaving returns a finite, positive, plausible `sigma_DA` that
    is simply wrong. Nothing else in the returned array says so.
    """
    nuc, elec, phi_d, ing, eps, chi = small_deck
    with pytest.warns(UserWarning, match="initial real-region norm"):
        td_nrm_da_cross_section(
            nuc, elec, F2, phi_d, eps, chi, 0, 0.03, ingredients=ing, dt=1.0, n_steps=4
        )


@pytest.mark.slow
def test_td_da_matches_the_time_independent_cross_section(f2_deck):
    """GATE: same model, same grids, two routes -- `sigma_DA` must agree.

    EXPENSIVE: ~31 min on the 12-core dev machine (6000 steps at 0.30 s,
    plus ~10 s of ingredients, ~17 s of factorization, and the three-energy
    TI oracle). There is no cheaper honest version -- see the convergence
    table below for what shortening `T` actually costs.

    Measured 2026-08-19 at these settings: `sigma_TD/sigma_TI` =
    1.0097 / 1.0138 / 1.0102 at E = 0.02 / 0.03 / 0.05 Ha, i.e. agreement to
    1.4%, with `unabsorbed/S(0)` = 0.007 / 0.008 / 0.009. The gate sits at
    5e-2, ~3.6x the achieved error; the extra headroom is not idle, because
    the residual is an OSCILLATION and not a bias (see below), so another
    platform's round-off can land at a different phase of it.

    WHY `T = 12000` AND NOT LESS. The transform's truncation error decays
    with the packet, and on this deck the packet needs ~2500 a.u. just to
    REACH the ECS absorber at R = 10.7. Measured worst-of-three-energies
    `max |sigma_TD/sigma_TI - 1|` against `T` (dt=2, this fixture):

        T =  4000   0.29        T = 10000   0.024
        T =  6000   0.13        T = 12000   0.014
        T =  8000   0.065       T = 14000   0.022

    -- so a `T = 8000` gate would have to sit near 1.5e-1 and a `T = 6000`
    one near 3e-1. Beyond T = 12000 the curve stops improving: from there
    through T = 15500 every ratio lies inside [0.976, 1.025] with no trend,
    which is residual oscillation, not remaining truncation. That is also
    why the gate must not be read as "1.4% is the method's accuracy" -- 1.4%
    is where this run's oscillation happened to sit.

    WHY `dt = 2`. The `dt`-error and the `T`-error are separable and add in
    quadrature; at order 3 the propagation error falls as `dt^6`, so dt=2 is
    expected to sit far below the truncation floor here. Measured rather
    than assumed: a dt=1 rerun of this fixture over T <= 6000 reproduces the
    dt=2 ratios to 3-4 significant figures at every checkpoint (at T=6000,
    0.9872/1.0407/0.8729 against 0.9872/1.0406/0.8733), for twice the cost.

    WHY THE PACKET IS KNOWN TO HAVE LEFT, rather than assumed: `<R>_t`
    climbs monotonically 2.66 -> 9.33 bohr through T = 4000 while `<P>_t`
    rises toward `K_R`, and `S(t)` then falls 0.94 -> 0.007 as the wave
    crosses into the absorber. It does not reach zero, and cannot: F2's
    `V_d(R)` well holds >= 24 near-real modes carrying ~5e-3 of the launch
    norm (module docstring of `td_cross_section.py`). `S(t)` plateaus at
    0.006-0.009 from T ~ 12000 on -- which is why this gate is keyed to
    `sigma_DA` being stationary in `T`, not to an absolute survival floor.
    """
    nuc, elec, phi_d, ing, eps, chi = f2_deck
    energies = np.array([0.02, 0.03, 0.05])
    want = nrm_da_cross_section(nuc, elec, F2, phi_d, eps, chi, 0, energies, ingredients=ing)
    got = td_nrm_da_cross_section(
        nuc, elec, F2, phi_d, eps, chi, 0, energies, ingredients=ing, dt=2.0, n_steps=6000
    )
    assert np.all(want > 0.0), "the TI oracle is zero -- pick open energies"
    rel = np.abs(got - want) / want
    assert np.all(rel < 5e-2), f"sigma_TD/sigma_TI = {got / want}"


@pytest.fixture(scope="module")
def f2_lcp(f2_deck):
    """`(V_d(R), Gamma(R))` for the same F2 deck the nonlocal route uses.

    Built on the fixture's own REDUCED electronic grid, plus a second copy of
    it at a different ECS angle (`local_complex_potential` needs two to match
    the pole). That is not a compromise here: every comparison below is
    DIFFERENTIAL -- the shipped `lcp_da_cross_section` and the Markovian
    propagation consume this same curve -- so the curve only has to be a
    curve, not a converged one. It costs ~3 s.
    """
    nuc, elec, _, _, _, _ = f2_deck
    elec_b = electronic_grid(r_max=13.0, order=5, n_complex=2, angle_deg=40.0)
    return local_complex_potential(F2, nuc, elec, elec_b)


def test_markovian_requires_the_local_potential(small_deck):
    nuc, elec, phi_d, ing, eps, chi = small_deck
    with pytest.raises(ValueError, match="needs both Vd and Gamma"):
        td_nrm_da_cross_section(
            nuc, elec, F2, phi_d, eps, chi, 0, 0.03, markovian=True, dt=1.0, n_steps=2
        )


def test_the_local_potential_is_rejected_without_markovian(small_deck):
    """Silently ignoring `Vd`/`Gamma` would return a nonlocal answer to a local call."""
    nuc, elec, phi_d, ing, eps, chi = small_deck
    with pytest.raises(ValueError, match="markovian=True arguments"):
        td_nrm_da_cross_section(
            nuc,
            elec,
            F2,
            phi_d,
            eps,
            chi,
            0,
            0.03,
            Vd=np.zeros(nuc.n, dtype=np.complex128),
            Gamma=np.zeros(nuc.n),
            dt=1.0,
            n_steps=2,
        )


def test_markovian_rejects_an_arm_count(small_deck):
    nuc, elec, phi_d, ing, eps, chi = small_deck
    with pytest.raises(ValueError, match="no projected-state arms"):
        td_nrm_da_cross_section(
            nuc,
            elec,
            F2,
            phi_d,
            eps,
            chi,
            0,
            0.03,
            markovian=True,
            n_states=3,
            Vd=np.zeros(nuc.n, dtype=np.complex128),
            Gamma=np.ones(nuc.n),
            dt=1.0,
            n_steps=2,
        )


@pytest.mark.slow
def test_markovian_limit_reproduces_the_lcp_cross_section(f2_deck, f2_lcp):
    """GATE: PRA 47 Eq. (2.11) -- the LCP *is* the Markovian limit.

    Eq. (2.11) collapses the memory kernel to `i[Delta_L - (i/2)Gamma_L]
    delta(R-R') delta(t)`, leaving Eq. (2.15)'s one-curve propagation. That is
    the model `qscat.core.lcp.lcp_da_cross_section` already solves in the
    frequency domain, so the two must return the same number -- and the
    agreement is a DIFFERENTIAL one (same deck, same `(V_d, Gamma)`, same
    Eq. (54) extraction), which is why it can be gated tightly.

    Measured 2026-08-24 at `dt=2, T=12000`: `sigma_TD/sigma_TI` =
    1.000215 / 1.000198 / 0.999892 at E = 0.02 / 0.03 / 0.05 Ha. The residual
    is transform truncation and nothing else -- dt = 1, 2 and 4 reproduce
    those three ratios to all six digits printed, so the `dt^6` propagation
    error is far below it, and extending `T` shrinks it (max |ratio-1| =
    2.4e-2 / 1.3e-3 / 2.2e-4 / 3.6e-4 / 7.7e-5 at T = 4000 / 8000 / 12000 /
    16000 / 20000, the last two being oscillation about the converged value).
    The gate sits at 5e-3, ~14x the worst measured residual in that
    stationary range -- far tighter than the nonlocal route's 5e-2, which it
    should be: this comparison has no model difference in it at all.

    It is CHEAP, unlike the nonlocal gate: no arms means `H_ext` is `N_R`
    square (974) rather than `(1 + n_states) * N_R` (53570), so 6000 steps
    cost ~4 s against ~30 min. It carries `@slow` for the fixtures it shares
    (the electronic pole walk and the deck build), not for the propagation.
    """
    nuc, elec, phi_d, _, eps, chi = f2_deck
    Vd, Gamma = f2_lcp
    energies = np.array([0.02, 0.03, 0.05])

    want = lcp_da_cross_section(nuc, F2.mu, Vd, Gamma, eps, chi, 0, energies)
    got = td_nrm_da_cross_section(
        nuc,
        elec,
        F2,
        phi_d,
        eps,
        chi,
        0,
        energies,
        markovian=True,
        Vd=Vd,
        Gamma=Gamma,
        dt=2.0,
        n_steps=6000,
    )
    assert np.all(want > 0.0), "the LCP reference is zero -- pick open energies"
    rel = np.abs(got - want) / want
    assert np.all(rel < 5e-3), f"sigma_TD/sigma_LCP = {got / want}"


@pytest.mark.slow
def test_the_discrete_state_potential_does_not_reproduce_the_lcp(f2_deck, f2_lcp):
    """The OTHER `V_d` -- and why Eq. (2.15) cannot be read as taking it.

    PRA 77 Eq. (20)'s `V_d = V_0 + <phi_d|H_el|phi_d>` (`NrmIngredients.
    v_d_discrete`) and `qscat.core.lcp`'s `Vd` (`E_res(R) + V_0`) are the two
    candidates for Eq. (2.15)'s bracket, and the paper's own Eq. (2.14),
    `V_d + Delta_L = E_res + V_0`, says the second one is what belongs there:
    the first is missing the level shift `Delta_L`. Measured rather than
    argued -- swapping it in gives `sigma/sigma_LCP` = 0.346 / 0.419 / 7.14 at
    E = 0.02 / 0.03 / 0.05 Ha (2026-08-24, this fixture, dt=2, T=12000), a
    disagreement that is large AND energy-dependent, i.e. not a normalization
    anyone could absorb.

    The two potentials agree exactly where `Gamma = 0` and separate only where
    it does not, which is `Delta_L`'s own support: measured on this deck,
    `V_d(Eq.20) - V_d(LCP)` = 2e-6 / 4e-5 / 9.5e-4 / 0.042 / 0.27 / 1.17 Ha at
    R = 3.99 / 3.50 / 3.01 / 2.49 / 2.20 / 1.51 bohr, against `Gamma` = 0 for
    the first three and 0.0095 for the last three.

    This test exists so that reading is not re-litigated by a later
    "simplification" that reaches for the `ing` already in scope.
    """
    nuc, elec, phi_d, ing, eps, chi = f2_deck
    Vd, Gamma = f2_lcp
    energies = np.array([0.02, 0.03, 0.05])

    want = lcp_da_cross_section(nuc, F2.mu, Vd, Gamma, eps, chi, 0, energies)
    vd_eq20 = continue_to_tail(ing.v_d_discrete, ing.R, nuc)
    got = td_nrm_da_cross_section(
        nuc,
        elec,
        F2,
        phi_d,
        eps,
        chi,
        0,
        energies,
        markovian=True,
        Vd=vd_eq20,
        Gamma=Gamma,
        dt=2.0,
        n_steps=6000,
    )
    ratio = got / want
    assert np.all(np.abs(ratio - 1.0) > 0.3), f"sigma/sigma_LCP = {ratio}"


# --- vibrational excitation -------------------------------------------------
#
# The N2 deck below is the one `test_nrm_propagation.py` gates the
# vector-to-vector identity on (179 nuclear x 74 electronic, all 73 projected
# states), transcribed rather than imported because these test modules are not
# an importable package. Using the SAME deck is the point: the VE cross
# section is a contraction of the very vector that gate checks, so the two
# residuals are directly comparable.
#
# VE and DA differ in WHERE the observable lives, and it shows in what `T`
# each needs. `t_resonant` contracts `Psi_d` against `chi_f`, which sits in
# the interaction region, so the transform converges as the amplitude THERE
# decays; `da_sigma_from_psi` reads the wavefunction value at the outermost
# real node, so its packet has to cross the whole box first. On F2's deck
# below that is worth two orders at a sixth of the propagation: VE converges
# by T = 2000 while DA needs T = 12000 -- and at T = 2000 that packet still
# holds 94% of its initial real-region norm, so nothing has LEFT.

_N2_REAL = ((3, 1.5), (8, 3.0), (2, 4.0), (4, 8.0))
_N2_COMPLEX = ((3, 20.0),)


@pytest.fixture(scope="module")
def n2_deck():
    """N2 on a small-but-COMPLETE deck: 179 nuclear x 74 electronic, 73 arms.

    Deliberately coarser than the production `N2:emoscat` deck. Every
    comparison below is DIFFERENTIAL -- both routes run on the same `ing` on
    the same grids -- so it does not need a physically converged
    discretisation, only one where the packet decays inside an affordable
    propagation. The absolute cross sections here are NOT converged N2
    numbers and must not be quoted as such; `validation/diatomic/ve_nrm.py`
    owns those.

    `n_states` is left at `None`: the COMPLETE projected-state sum is a
    correctness requirement of the time-dependent route, not a performance
    knob (`propagation._warn_if_diverging`).
    """
    nuc = segmented_grid(_N2_REAL, _N2_COMPLEX, angle_deg=35.0, quadrature=10)
    elec = electronic_grid(r_max=11.0, order=6, n_complex=3)
    phi_d = AsymptoticDiscreteState(elec, N2, R_inf=nuc.R0)
    r_values = nuc.points[nuc.points.imag == 0.0].real[::-1]
    ing = nrm_ingredients(elec, N2, phi_d, r_values)
    eps, chi = vibrational_states(nuc, N2.mu, 4, N2.v0)
    return nuc, elec, phi_d, ing, eps, chi


def test_ve_non_positive_energies_are_zero_like_the_time_independent_route(small_deck):
    """`E <= 0` gives a row of zeros on both routes, and costs no propagation."""
    nuc, elec, phi_d, ing, eps, chi = small_deck
    got = td_nrm_ve_cross_section(
        nuc,
        elec,
        F2,
        phi_d,
        eps,
        chi,
        0,
        [0, 1],
        np.array([-0.01, 0.0]),
        ingredients=ing,
        dt=1.0,
        n_steps=4,
    )
    assert got.shape == (2, 2)
    assert np.array_equal(got, np.zeros((2, 2)))


def test_ve_a_scalar_energy_returns_one_entry_per_final_channel(small_deck):
    """`driven.ve_cross_section`'s shape convention, which `nrm_ve_cross_section`
    follows and this route must too: scalar `E` drops the energy axis."""
    nuc, elec, phi_d, ing, eps, chi = small_deck
    with pytest.warns(UserWarning):
        got = td_nrm_ve_cross_section(
            nuc, elec, F2, phi_d, eps, chi, 0, [0, 1], 0.03, ingredients=ing, dt=1.0, n_steps=4
        )
    assert got.shape == (2,)
    assert np.all(got > 0.0)


def test_ve_a_closed_final_channel_is_zero(small_deck):
    """A `v'` above the total energy contributes `0.0`, not a negative-`k` value.

    The open channels of the same call are unaffected, which is what a
    mis-slotted `vprimes` loop would break while still returning finite,
    positive, plausible numbers.
    """
    nuc, elec, phi_d, ing, eps, chi = small_deck
    e_kin = 0.5 * float(eps[2] - eps[0])  # opens v'=1, leaves v'=2 shut
    with pytest.warns(UserWarning):
        got = td_nrm_ve_cross_section(
            nuc, elec, F2, phi_d, eps, chi, 0, [0, 1, 2], e_kin, ingredients=ing, dt=1.0, n_steps=4
        )
    assert got[2] == 0.0
    assert np.all(got[:2] > 0.0)


def test_markovian_ve_rejects_the_background(small_deck):
    """`T^bg` is built from `phi_d`, which the local model does not have.

    Allowing it would silently assemble a third model -- a local resonant
    term plus a nonlocal background -- and return it under a call that reads
    as the LCP.
    """
    nuc, elec, phi_d, ing, eps, chi = small_deck
    with pytest.raises(ValueError, match="include_background=False"):
        td_nrm_ve_cross_section(
            nuc,
            elec,
            F2,
            phi_d,
            eps,
            chi,
            0,
            [0, 1],
            0.03,
            markovian=True,
            Vd=np.zeros(nuc.n, dtype=np.complex128),
            Gamma=np.ones(nuc.n),
            dt=1.0,
            n_steps=2,
        )


def test_markovian_ve_requires_the_local_potential(small_deck):
    nuc, elec, phi_d, ing, eps, chi = small_deck
    with pytest.raises(ValueError, match="needs both Vd and Gamma"):
        td_nrm_ve_cross_section(
            nuc,
            elec,
            F2,
            phi_d,
            eps,
            chi,
            0,
            [0, 1],
            0.03,
            markovian=True,
            include_background=False,
            dt=1.0,
            n_steps=2,
        )


@pytest.mark.slow
def test_td_ve_matches_the_time_independent_cross_section(n2_deck):
    """GATE: same model, same grids, two routes -- `sigma_VE` must agree.

    EXPENSIVE: ~8-11 min, almost all of it TWO 4000-step propagations (one per
    `include_background` setting; the propagation cannot be shared between
    them through the public API). The time-independent oracle costs ~4 s and
    the deck ~1 s.

    BOTH `include_background` settings are checked, and that pair is a
    diagnostic rather than a duplicate. `T^bg` (Eq. 37) contains no `Psi_d`:
    it is energy-domain, static, and produced by the SAME `t_background` call
    on both routes, so it is bit-identical between them. A discrepancy that
    appeared only with `include_background=True` could therefore not be a
    background error -- it would mean the resonant term is being combined
    with it wrongly (a sign, a conjugation, a channel energy).

    Measured 2026-08-24 at these settings, `sigma_TD/sigma_TI` over
    E = 0.06 / 0.10 / 0.15 Ha and v' = 0 / 1:

        include_background=True   0.999765  0.999793
                                  1.000073  0.999828
                                  0.999932  1.000271
        include_background=False  0.999786  0.999794
                                  1.000098  0.999830
                                  0.999911  1.000268

    i.e. agreement to 2.71e-4 (background) and 2.68e-4 (no
    background), against a gate at 1e-3 -- 3.7x the achieved error.

    WHY `T = 4000`. Worst-of-six `max |ratio - 1|` against `T` on this deck
    (dt = 1, `include_background=True`):

        T = 1000   4.2e-1        T = 4000   2.7e-4
        T = 2000   5.0e-2

    The residual at T = 4000 sits at the level of the vector-to-vector
    identity gate (`test_nrm_propagation.py`, 1.73e-4 on this same deck at
    the same `T`/`dt`), which is what it should be: the cross section is a
    contraction of that vector, and the two worst channels here are the
    E = 0.06 ones -- the near-threshold final channels converge slowest.

    WHY `dt = 1` AND NOT 2, unlike the DA gate. The error budget measured in
    `test_nrm_propagation.py` -- `truncation(T) = 0.40*sqrt(S(T)/S(0))` and
    `propagation(dt=1) = 1.43e-4`, adding in quadrature -- makes the two
    terms COMPARABLE here, so the `dt^6` term is no longer negligible:
    dt = 2 costs a factor 64 in it. Measured on this fixture at T = 4000,
    `include_background=True`: max |ratio - 1| = 1.52e-2 at dt = 2 against
    2.71e-4 at dt = 1 -- a factor of 56 where `dt^6` predicts 64. The F2 DA gate can afford dt = 2
    because its truncation floor is 1.4e-2, two orders above its own
    propagation error; this one cannot.
    """
    nuc, elec, phi_d, ing, eps, chi = n2_deck
    energies = np.array([0.06, 0.10, 0.15])
    vprimes = [0, 1]

    for include_background in (True, False):
        want = nrm_ve_cross_section(
            nuc,
            elec,
            N2,
            phi_d,
            eps,
            chi,
            0,
            vprimes,
            energies,
            ingredients=ing,
            include_background=include_background,
        )
        got = td_nrm_ve_cross_section(
            nuc,
            elec,
            N2,
            phi_d,
            eps,
            chi,
            0,
            vprimes,
            energies,
            ingredients=ing,
            include_background=include_background,
            dt=1.0,
            n_steps=4000,
        )
        assert np.all(want > 0.0), "the TI oracle is zero -- pick open channels"
        rel = np.abs(got - want) / want
        assert np.all(rel < 1e-3), (
            f"include_background={include_background}: sigma_TD/sigma_TI = {got / want}"
        )


@pytest.mark.slow
def test_markovian_ve_reproduces_the_local_cross_section(f2_deck, f2_lcp):
    """GATE: PRA 47 Eq. (2.11) for VE -- the LCP *is* the Markovian limit.

    The DA sibling of this test
    (`test_markovian_limit_reproduces_the_lcp_cross_section`) compares against
    the shipped `lcp_da_cross_section`. There is no `lcp_ve_cross_section` in
    `qscat` -- the repository's LCP VE route lives in
    `projects/n2_ti_cross_section/cross_section.py`, which `libs/qscat` must
    not import -- so the reference is assembled here from shipped `qscat`
    pieces: `solve_nuclear` with `F -> diag(-(i/2)Gamma)` (the substitution
    `dissociation`'s module docstring documents `solve_nuclear` as being
    public FOR), the doorway `sqrt(Gamma/2pi) chi_v` on both ends, and
    `sigma = 4 pi^3 |S|^2 / 2E`. That is the same formula the projects module
    implements, verified equal to every printed digit (2026-08-24, this deck,
    ratio 1.000000000000 at all six entries).

    THE DOORWAY IS SUBSTITUTED AT BOTH ENDS, which is the whole content of
    `markovian=True` for VE. Eq. (2.11) localizes the kernel; the LCP this
    repository ships also localizes the DOORWAY (`lcp_initial_packet`), and
    for VE that doorway appears twice -- as the launch state AND as the
    coupling the exit channel is contracted against. Keeping the nonlocal
    `V^+_dk_f` at the exit while localizing the kernel would be a third
    model, and it would NOT reproduce this reference.

    Measured 2026-08-24 at `dt = 2, T = 4000`: every one of the six
    `sigma_TD/sigma_LCP` entries lies within 3.4e-6 of 1 (E = 0.02/0.03/0.05
    Ha, v' = 0/1). It is converged well before that and stays there --
    max |ratio-1| = 5.6e-6 / 3.3e-6 / 4.3e-6 at T = 2000 / 4000 / 12000 --
    against the DA Markovian gate on this SAME deck needing T = 12000 to
    reach 2.2e-4. That gap is the point of `td_nrm_ve_cross_section`'s
    docstring: `T^res` weights `Psi_d` by `chi_f`, which lives in the
    interaction region, so the transform converges as the packet decays out
    of that region rather than waiting for it to cross the whole box.

    The gate sits at 5e-5, ~15x the worst residual in the stationary range,
    matching the DA Markovian gate's own headroom philosophy: this comparison
    contains no model difference at all, so it can be gated far tighter than
    the nonlocal route's.

    CHEAP: no arms means `H_ext` is `N_R` square (974) rather than 53570, so
    2000 steps cost a few seconds. It carries `@slow` for the fixtures it
    shares (the deck build and the electronic pole walk), not for the
    propagation.
    """
    nuc, elec, phi_d, _, eps, chi = f2_deck
    Vd, Gamma = f2_lcp
    energies = np.array([0.02, 0.03, 0.05])
    vprimes = [0, 1]

    doorway = np.sqrt(Gamma / (2.0 * np.pi)).astype(np.complex128)
    no_kernel = np.zeros((nuc.n, nuc.n), dtype=np.complex128)
    want = np.zeros((energies.size, len(vprimes)), dtype=np.float64)
    for ie, e_kin in enumerate(energies):
        e_total = float(e_kin) + float(eps[0])
        psi = solve_nuclear(nuc, F2.mu, Vd - 0.5j * Gamma, no_kernel, doorway * chi[0], e_total)
        for jv, vp in enumerate(vprimes):
            if e_total - float(eps[vp]) <= 0.0:
                continue
            t = t_resonant(chi[vp], doorway, psi)
            want[ie, jv] = 4.0 * np.pi**3 * abs(t) ** 2 / (2.0 * float(e_kin))

    got = td_nrm_ve_cross_section(
        nuc,
        elec,
        F2,
        phi_d,
        eps,
        chi,
        0,
        vprimes,
        energies,
        markovian=True,
        include_background=False,
        Vd=Vd,
        Gamma=Gamma,
        dt=2.0,
        n_steps=2000,
    )
    assert np.all(want > 0.0), "the LCP reference is zero -- pick open channels"
    rel = np.abs(got - want) / want
    assert np.all(rel < 5e-5), f"sigma_TD/sigma_LCP = {got / want}"


@pytest.mark.slow
def test_td_ve_matches_the_time_independent_cross_section_on_f2(f2_deck):
    """GATE: the same VE comparison on the molecule whose DA channel is OPEN.

    N2's VE packet decays by autodetachment alone. F2's also DISSOCIATES -- at
    these energies the DA channel is open, and this is the same `Psi_d` the
    31-minute DA gate above propagates, on the same fixture. So this is not a
    second copy of the N2 test: it asks whether the VE contraction converges
    while the DA one, from the very same packet, is still nowhere near
    converged.

    It does. Measured 2026-08-24 at `dt = 2, T = 2000`, E = 0.03/0.05 Ha,
    v' = 0/1:

        include_background=True   1.000034  1.000027  1.000039  1.000059
        include_background=False  1.000059  1.000022  1.000098  1.000061

    against this deck's DA gate at 0.29 (T = 4000) and 1.4e-2 (T = 12000).
    And it is NOT because the packet has left: `unabsorbed/S(0)` here is
    0.938, so the truncation warning fires (expected -- `_UNABSORBED_TOL` is
    calibrated on DA, where the observable IS the far-field amplitude). What
    has decayed is the amplitude under `chi_f`, which is all `t_resonant`
    integrates.

    WHY THE GATE IS 1e-2 AND NOT 1e-4. The residual does NOT keep falling: at
    T = 4000 it reads 2.5e-3 against T = 2000's 5.9e-5. That is not a
    regression, it is the OSCILLATION F2's near-real modes leave in the
    transform -- the >=24 modes with `|Im E| = 1.5e-7 ... 7.7e-6` living in
    the `V_d` well at R ~ 3.36 (module docstring, sec. 4.2 of
    `docs/physics/nrm-time-dependent.md`), which sit UNDER `chi_f` and
    contribute a term that oscillates in `T` rather than decaying. No
    affordable `T` removes them. So the defensible statement about F2 is an
    AMPLITUDE, <= 2.5e-3 over the measured range, and the gate is 4x that
    rather than a multiple of whichever phase T = 2000 happened to land on.
    N2, which has no such well, converges monotonically and is gated at 1e-3
    accordingly.

    COST: ~10-20 min, two 1000-step propagations of the 53570-square `H_ext`
    (plus ~22 s per `include_background` setting for the time-independent
    oracle, whose `t_background` runs an electronic scattering solve at each
    of 819 real nuclear nodes).
    """
    nuc, elec, phi_d, ing, eps, chi = f2_deck
    energies = np.array([0.03, 0.05])
    vprimes = [0, 1]

    for include_background in (True, False):
        want = nrm_ve_cross_section(
            nuc,
            elec,
            F2,
            phi_d,
            eps,
            chi,
            0,
            vprimes,
            energies,
            ingredients=ing,
            include_background=include_background,
        )
        got = td_nrm_ve_cross_section(
            nuc,
            elec,
            F2,
            phi_d,
            eps,
            chi,
            0,
            vprimes,
            energies,
            ingredients=ing,
            include_background=include_background,
            dt=2.0,
            n_steps=1000,
        )
        assert np.all(want > 0.0), "the TI oracle is zero -- pick open channels"
        rel = np.abs(got - want) / want
        assert np.all(rel < 1e-2), (
            f"include_background={include_background}: sigma_TD/sigma_TI = {got / want}"
        )
