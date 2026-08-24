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
from qscat.core.nrm.dissociation import nrm_da_cross_section
from qscat.core.nrm.ingredients import nrm_ingredients
from qscat.core.nrm.nonlocal_potential import continue_to_tail
from qscat.core.nrm.td_cross_section import td_nrm_da_cross_section
from qscat.core.vibrational import vibrational_states
from qscat.model import F2

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
