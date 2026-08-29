"""`local_width`, and the identity that says what the auxiliary packets mean.

The second half of this file is a MEASUREMENT that had not been made before.
The time-dependent nonlocal model resums PRA 47 Eq. (2.1)'s memory integral
into auxiliary nuclear packets, and it is tempting to read `‖phi_n‖^2` as the
population of electronic channel `n`. Under exterior complex scaling `H_ext` is
complex SYMMETRIC, so nothing conserves that norm, and the one exact statement
available -- what the coupling takes out of `Psi_d` it must put into the arms --
holds only to the extent `V_dn` is real. These tests record how far that is
from true on the N2 gate deck.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import pytest
import scipy.sparse as sp
from qscat.core.grids import electronic_grid, segmented_grid
from qscat.core.nrm.coupling import gamma_from_coupling, v_dk_plus
from qscat.core.nrm.discrete_state import AsymptoticDiscreteState
from qscat.core.nrm.extended import LaunchBasis, extended_hamiltonian, initial_packet
from qscat.core.nrm.ingredients import nrm_ingredients
from qscat.core.nrm.memory import MemoryRecorder, MemorySpec, local_width
from qscat.core.nrm.nonlocal_potential import continue_to_tail, nonlocal_operator
from qscat.core.nrm.propagation import propagate_nrm
from qscat.core.vibrational import vibrational_states
from qscat.dvr import FemDvrEcsGrid
from qscat.dvr.spec import ElementSpec, GridSpec
from qscat.evolution import make_pade_stepper
from qscat.model import N2

# The nuclear half of `test_nrm_propagation.py`'s N2 gate deck, verbatim: 179
# nuclear points, 153 of them real. Both fixtures below share it, so the
# `local_width` identity and the imbalance measurement differ only in the
# ELECTRONIC box -- which is the one thing the identity turns out to be
# sensitive to (see `width_deck` and the module docstring of `memory.py`).
_N2_REAL_SEGMENTS = ((3, 1.5), (8, 3.0), (2, 4.0), (4, 8.0))
_N2_COMPLEX_SEGMENTS = ((3, 20.0),)

_E_KIN = 0.10
_V_INIT = 0


def _nuclear_grid():
    return segmented_grid(_N2_REAL_SEGMENTS, _N2_COMPLEX_SEGMENTS, angle_deg=35.0, quadrature=10)


@pytest.fixture(scope="module")
def width_deck():
    """The gate deck's nuclear grid with an `r_max = 16` ELECTRONIC box.

    The gate deck's own `r_max = 11` box is too small for the local-limit
    identity below -- measured median ratio 1.218 there against 0.996 here,
    converged (0.997 at `r_max` = 24 and 30). That is not a defect of the gate
    deck: it was sized for `test_nrm_propagation.py`'s algebraic transform
    identity, which runs both routes on the same ingredients and so cannot see
    an under-converged electronic box at all. `n_complex=4` rather than 3 keeps
    `nrm_ingredients`' adiabatic tracking clean (at `n_complex=3` it warns
    about a `_sign_align` overlap).
    """
    nuc = _nuclear_grid()
    elec = electronic_grid(r_max=16.0, order=6, n_complex=4)
    phi_d = AsymptoticDiscreteState(elec, N2, R_inf=nuc.R0)
    # nrm_ingredients requires strictly DESCENDING R.
    r_values = nuc.points[nuc.points.imag == 0.0].real[::-1]
    ing = nrm_ingredients(elec, N2, phi_d, r_values)
    eps, _chi = vibrational_states(nuc, N2.mu, 4, N2.v0)
    e_total = _E_KIN + float(eps[_V_INIT])
    f_matrix = nonlocal_operator(ing, nuc, N2, e_total, n_states=None)
    return nuc, elec, phi_d, f_matrix, e_total


def test_local_width_reproduces_the_eq_68_width_at_the_local_energy(width_deck):
    """`F`'s local limit IS a width, and it matches the independently computed one.

    `nonlocal-resonance-model.md` Sec. 9 records this ratio as median 0.977 on
    NO and 1.011 on F2. It is the reason the Markovian reference can be taken
    from `F` rather than from a pole walk that freezes.

    THE COMPARISON IS AGAINST Eq. (68), NOT AGAINST `local_complex_potential`.
    `F(E)` carries one total energy, so its local limit is the width at the
    LOCAL electron energy `eps_loc(R) = E - V_0(R)`; the LCP's `Gamma` is the
    width at the RESONANCE position `eps_res(R) = E_res(R) - V_0(R)`. On this
    deck those are nowhere near each other: `V_0` is a well, so `eps_loc` is
    positive only BETWEEN R = 1.7426 and R = 2.4284, and at the outer crossing
    `eps_res` is still ~0.07 Ha. The ratio of `local_width` to the pole walk's
    `Gamma` therefore runs 0.12 to 8.9 over the R-range where `Gamma` is
    nonzero, monotonically, with the walk's frozen region playing no part. That
    is a difference of energy argument, not an error in either quantity, and no
    frozen-region cut can repair it.
    `gamma_from_coupling(v_dk_plus(...))` evaluated per node at `eps_loc(R)` is
    an INDEPENDENT construction (P-projected scattering states, not a resolvent
    sum) at the energy `F` was actually built at. Eq. (68) against the pole
    width AT `eps_res` is gated separately, in `test_nrm_coupling.py`.

    MEASURED (2026-08-27, this deck): median 0.9958 over the 28 nodes with
    `eps_loc > 0.02` Ha, and 1.00 +- 0.009 across the 19-node plateau
    `1.88 <= R <= 2.26`. The worst node is 0.880 at R = 1.8261 -- the SECOND
    innermost of the window, not the innermost (R = 1.7967 reads 0.897): the
    degradation tracks proximity to the inner `eps_loc` crossing at R = 1.7426
    but is not monotone node by node. The `eps_loc > 0.02` cut is the same
    near-threshold exclusion `test_nrm_coupling.py` applies for the same
    reason: at `eps_loc` = 0.016 and 0.006 the ratio reads 1.82 and 8.49, where
    the width varies by orders across the ~10 nodes the kernel spans and no
    local limit exists to compare to.
    """
    nuc, elec, phi_d, f_matrix, e_total = width_deck
    gamma_loc = local_width(f_matrix, nuc)

    real = np.flatnonzero(nuc.points.imag == 0.0)
    r_real = nuc.points[real].real
    eps_loc = e_total - np.real(N2.v0(r_real))
    open_ = np.flatnonzero(eps_loc > 0.02)
    assert open_.size > 20, "no usable comparison window -- check the deck"

    gamma_68 = np.array(
        [
            gamma_from_coupling(
                v_dk_plus(elec, N2, phi_d, np.array([r_real[j]]), float(eps_loc[j]))
            )[0]
            for j in open_
        ]
    )
    ratio = gamma_loc[real][open_] / gamma_68
    assert np.all(np.isfinite(ratio))
    median = float(np.median(ratio))
    worst = float(np.max(np.abs(ratio - 1.0)))
    assert abs(median - 1.0) < 0.03, f"median Gamma_loc/Gamma_68 = {median:.4f}"
    # 12.0% measured at R = 1.8261, the second-innermost node of the window;
    # 0.20 is headroom for that one point, not a claim that 20% would be
    # acceptable at the plateau.
    assert worst < 0.20, f"worst node |ratio - 1| = {worst:.4f} (median {median:.4f})"


def test_the_diagonal_of_f_is_not_the_local_limit(width_deck):
    """`diag F` is 0.16x the row sum here -- the claim `local_width` rests on.

    The kernel spans ~10 nuclear nodes, so most of each row sits off the
    diagonal. Reading `diag F` as the local potential (the obvious mistake, and
    the one the `sqrt(w)` row sum exists to avoid) would under-report the width
    by a factor of six on this deck, and by seven on the time-independent decks
    of `nonlocal-resonance-model.md` Sec. 9 (0.14x).
    """
    nuc, _elec, _phi_d, f_matrix, e_total = width_deck
    gamma_loc = local_width(f_matrix, nuc)
    gamma_diag = -2.0 * np.diag(f_matrix).imag

    real = np.flatnonzero(nuc.points.imag == 0.0)
    r_real = nuc.points[real].real
    open_ = np.flatnonzero(e_total - np.real(N2.v0(r_real)) > 0.02)
    share = float(np.median(gamma_diag[real][open_] / gamma_loc[real][open_]))
    assert 0.05 < share < 0.35, f"median diag/row-sum = {share:.4f}"


def test_local_width_rejects_a_matrix_that_is_not_the_grid_s(width_deck):
    """A `F` from another grid is caught, not broadcast against these weights."""
    nuc, _elec, _phi_d, f_matrix, _e_total = width_deck
    with pytest.raises(ValueError, match="one row and column per nuclear DVR node"):
        local_width(f_matrix[:-1, :-1], nuc)


# --- the identity that decides what `‖phi_n‖^2` may be called ----------------


@pytest.fixture(scope="module")
def propagation_deck():
    """`test_nrm_propagation.py`'s N2 gate deck, unchanged, plus `V_dn`.

    Unchanged deliberately: the numbers below are meant to be comparable with
    the prototype's and with the transform identity's, so the electronic box
    stays at `r_max = 11` even though `width_deck` above needs a larger one.
    """
    nuc = _nuclear_grid()
    elec = electronic_grid(r_max=11.0, order=6, n_complex=3)
    phi_d = AsymptoticDiscreteState(elec, N2, R_inf=nuc.R0)
    r_values = nuc.points[nuc.points.imag == 0.0].real[::-1]
    ing = nrm_ingredients(elec, N2, phi_d, r_values)
    eps, chi = vibrational_states(nuc, N2.mu, 4, N2.v0)
    v_dn = np.array([continue_to_tail(ing.V_dn[:, n], ing.R, nuc) for n in range(ing.E_n.shape[1])])
    e_total = _E_KIN + float(eps[_V_INIT])
    gamma_loc = local_width(nonlocal_operator(ing, nuc, N2, e_total, n_states=None), nuc)
    return nuc, elec, phi_d, ing, eps, chi, v_dn, gamma_loc


# The finite-difference probe: a 3-point central difference of step `_FD_H`
# about `t = _FD_TIME`, on the packet the loop above has already propagated
# there. `t = 20` is past the launch transient (the arms fill within ~10 a.u.)
# and long before the packet reaches the ECS tail, so both block norms are
# smooth and no flux has left the real region -- the FULL-grid and real-region
# variants of every number below agree to all digits shown. FULL grid is what
# is asserted, because only there is `d‖psi‖^2/dt = 2 Im<psi|H psi>` exact;
# restricting to the real region turns genuine outgoing flux into an apparent
# residual, which is the difference between this scheme and a real-region one.
_FD_TIME = 20
_FD_H = 0.02
# Measured 2026-08-27 at `_FD_H`: 4.92e-5 (d) and 1.71e-5 (arm). The residual
# is FINITE-DIFFERENCE TRUNCATION, not a defect -- halving `_FD_H` divides it
# by ~4.4 through 1.15e-3 / 3.70e-4 (h=0.08), 2.31e-4 / 7.77e-5, 4.92e-5 /
# 1.71e-5, 1.11e-5 / 3.94e-6, 2.63e-6 / 9.42e-7 (h=0.005), i.e. clean h^2.
# 5e-4 is ~10x the measured d value and ~29x the arm one, and still fails a
# 0.1% error in the coupling block (9.8e-4, mutation-tested).
_FD_TOL = 5e-4


def _assert_the_rates_are_what_the_propagation_does(
    h_ext: sp.csr_matrix,
    psi_fd: npt.NDArray[np.complex128],
    n_r: int,
    n_arm: int,
    v_dn: npt.NDArray[np.complex128],
) -> None:
    """`d/dt‖block‖^2` minus the block's own term IS the coupling term.

    The one statement that ties `exch_d` / `exch_arm` to the Hamiltonian
    actually being propagated rather than to the expressions the test writes
    down. The block-diagonal terms come from `h_ext` (they are not what is
    under suspicion); the coupling terms come from `ing.V_dn` via
    `continue_to_tail`, independently of `h_ext`'s assembly -- which is what
    makes a changed coupling block show up here.

    DO NOT ROUTE THIS THROUGH `MemoryRecorder`, however tempting it looks now
    that the caller above does. The recorder reads its coupling OUT OF `h_ext`
    (`h_ext[:n_r, n_r:]` and the transpose block), on purpose, so that the
    observables measure the operator actually being propagated. Feeding that
    same coupling into the check that `h_ext` propagates what the coupling
    describes closes the loop: `fd - own` and `coupling` would both derive from
    `h_ext`, the residual would vanish identically, and a wrong coupling block
    would pass silently. That is precisely the circularity the Task 1 review
    found in the algebraic assertion and asked to close, so the two paths stay
    two -- `ing.V_dn` here, `h_ext` in the recorder. Consolidating them is not
    a simplification; it deletes the test.
    """
    fine = make_pade_stepper(h_ext, _FD_H, 3)
    s0 = psi_fd
    s1 = fine(s0)
    s2 = fine(s1)

    def block_norms(x: npt.NDArray[np.complex128]) -> npt.NDArray[np.float64]:
        # CONJUGATING product again, over the FULL grid this time -- see the
        # comment above `_FD_TIME` for why the real-region variant is not used.
        return np.array([np.vdot(x[:n_r], x[:n_r]).real, np.vdot(x[n_r:], x[n_r:]).real])

    fd = (block_norms(s2) - block_norms(s0)) / (2.0 * _FD_H)
    own = np.array(
        [
            2.0 * float(np.vdot(s1[:n_r], h_ext[:n_r, :n_r] @ s1[:n_r]).imag),
            2.0 * float(np.vdot(s1[n_r:], h_ext[n_r:, n_r:] @ s1[n_r:]).imag),
        ]
    )
    d = s1[:n_r]
    arms = s1[n_r:].reshape(n_arm, n_r)
    coupling = np.array(
        [
            2.0 * float(np.vdot(d, np.einsum("nr,nr->r", v_dn, arms)).imag),
            2.0 * float(np.sum((arms.conj() * v_dn * d[None, :]).imag)),
        ]
    )
    rel = np.abs((fd - own) - coupling) / np.abs(coupling)
    assert rel[0] < _FD_TOL, f"d-block rate is not its coupling term: rel {rel[0]:.3e}"
    assert rel[1] < _FD_TOL, f"arm-block rate is not its coupling term: rel {rel[1]:.3e}"
    # And the d block's OWN term is ~0 (measured 4.7e-6 of the coupling): all of
    # `Psi_d`'s norm rate at this time is exchange with the arms, which is why
    # the exchange observable can be read as the discrete state's whole budget
    # while the packet is still inside the real region. The arms are NOT like
    # that -- their own term is -21x their coupling term here, because `E_n` is
    # complex and the arm blocks are strongly dissipative. That asymmetry is
    # the mechanism that makes `H_ext` dissipative at all, so it is asserted in
    # the direction it actually holds rather than for both blocks.
    assert abs(own[0] / coupling[0]) < 1e-3, f"d-block own term {own[0] / coupling[0]:.3e}"


@pytest.mark.slow
def test_the_coupling_exchange_balances_to_the_extent_V_dn_is_real(propagation_deck):
    """What the coupling removes from `Psi_d` it adds to the arms -- exactly only
    where `V_dn` is real, which under ECS it is not.

    For `i d/dt Psi = H Psi`, `d‖psi‖^2/dt = 2 Im<psi|H psi>`. The coupling
    block contributes `2 Im<Psi_d|sum_n V_dn phi_n>` to the discrete state and
    `2 Im<phi_n|V_dn Psi_d>` to the arms, and those two cancel iff `V_dn` is
    real; in general their sum is `4 sum_n Re[conj(Psi_d) phi_n] Im[V_dn]`.
    Both rates use the CONJUGATING product restricted to the real nuclear
    region -- a population is not a c-product, and that restriction is the
    whole point of the measurement.

    MEASURED (2026-08-27, N2 gate deck, `n_states=None`, E_kin = 0.10 Ha,
    dt = 1, 200 steps): the residual is not small. Against the LARGER of the
    two one-sided rates it is median 0.82, max 1.06 (min 0.005); against the
    discrete state's own rate `|exchange_d|` it is median 4.6 and reaches 1.1e4
    at the step where that rate passes through zero. Over the full 4000-step
    run the same numbers are median 0.90 / max 1.06 and median 9.5 / max 1.1e4.
    The arms GAIN conjugating norm from the coupling several times faster than
    `Psi_d` loses it, so the two rates do not describe one transfer.

    WHERE IT LIVES: entirely in the real nuclear region -- the ECS tail carries
    6e-71 of the run-summed `|density|` over this 200-step window (2e-55 over
    4000 steps), since `V_dn` there is ~1e-13 by Eq. (67) -- and within it, in
    the autodetachment region where `Gamma_loc` peaks. Split by R, and the
    split MOVES OUTWARD with time as the packet does, so the window has to be
    named with the number: over these 200 steps it is [2.0, 2.2] 68.6%,
    [1.8, 2.0] 17.2%, [2.2, 2.5] 14.2%, and nothing at all beyond 2.5; over
    4000 steps it is [2.0, 2.2] 57%, [2.2, 2.5] 32%, [1.8, 2.0] 11%, with
    under 0.5% (0.43%) outside [1.8, 2.5]. This is the electronic rotation
    leaking into nuclear-space bookkeeping exactly where the physics is, not an
    artifact of the nuclear absorber.

    So `‖phi_n‖^2` is NOT a population, and the partition observable must be
    reported as a relative channel decomposition. The exchange RATE is
    unaffected: it is a rate, not an amplitude, and needs no such reading.

    WHAT COMPUTES THE RATES. `MemoryRecorder` -- the shipped one. This test kept
    its own per-step copy of the arithmetic until the recorder existed; a gate
    on a private copy measures the copy, not the code that ships. Switching it
    over left every number above unchanged (checked 2026-08-27: the medians and
    maxima agree to every digit quoted, and the two implementations' `exch_d`
    and `exch_arm` agree to 2.1e-16 and 2.3e-15 relative).

    TWO THINGS DELIBERATELY DO NOT GO THROUGH THE RECORDER, and both look like
    duplication until you try removing them. `residual` is computed here
    because it is what the recorder's `imbalance` is tested AGAINST -- taking
    both from the recorder would be numpy checking one expression against
    itself. It is worth more than that now: `residual` is built from
    `ing.V_dn` while `imbalance` comes from `h_ext`'s assembled coupling block,
    so `imbalance == residual` stopped being an algebraic tautology when this
    test was re-pointed and became a CROSS-CHECK of the assembled operator
    against the ingredients it was assembled from. Mutation-tested 2026-08-27:
    scaling `h_ext[:n_r, n_r:]` by 1.001 takes `max|imbalance - residual|` from
    5.9e-15 to 2.6e-3 of `scale` -- twelve orders, and nine past the 1e-12
    gate. Do not "simplify" it by sourcing both sides alike.
    And the finite-difference block below keeps its own `ing.V_dn`,
    because the recorder's coupling comes out of `h_ext` and that check exists
    to catch a wrong `h_ext`; see its docstring for what merging the two paths
    would delete.

    WHY THE FINITE-DIFFERENCE BLOCK BELOW IS PART OF THIS TEST. The algebraic
    assertion (`imbalance == residual`) is numpy recomputing one expression two
    ways: it holds for ANY `psi` and ANY `v_real` and touches no qscat code, so
    on its own it would keep passing if `extended_hamiltonian` stopped
    propagating what these rates describe -- the natural "make it Hermitian"
    edit `off_diag = coup.T` -> `coup.conj().T`, say. The finite-difference
    check closes that: it measures `d/dt` of each block's norm along the ACTUAL
    propagation and requires the remainder, after subtracting that block's own
    diagonal-block term, to equal the coupling term built from `ing.V_dn`.
    Mutation-tested 2026-08-27: conjugating the lower coupling block takes the
    arm residual from 1.7e-5 to 3.3e-1, and merely rescaling it by 1.001 takes
    it to 9.8e-4 -- both far outside the 5e-4 gate.
    """
    nuc, elec, model_phi_d, ing, eps, chi, v_dn, gamma_loc = propagation_deck
    n_r = nuc.n
    n_arm = ing.E_n.shape[1]
    h_ext = extended_hamiltonian(ing, nuc, N2)
    launch = initial_packet(
        nuc, elec, N2, model_phi_d, ing, eps, chi, _V_INIT, np.array([_E_KIN]), rank_tol=1e-10
    )

    dt, n_steps = 1.0, 200
    step = make_pade_stepper(h_ext, dt, 3)
    psi = (launch.vectors @ launch.coeffs)[:, 0].astype(np.complex128)
    real = nuc.points.imag == 0.0
    v_real = v_dn[:, real]

    # THE SHIPPED RECORDER, not a copy of it. Everything below that is not the
    # independent `residual` algebra comes from `MemoryRecorder`, so what this
    # gate measures is the code `propagate_nrm(memory=...)` actually runs. It is
    # fed the propagation this test steps itself rather than being reached
    # through `propagate_nrm`, because the finite-difference block needs a
    # mid-run snapshot that no returned `TdNrmResult` carries.
    rec = MemoryRecorder(
        h_ext, nuc, MemorySpec(gamma_local=gamma_loc, n_channels=3), n_steps=n_steps, n_energies=1
    )
    residual = np.empty(n_steps + 1)
    psi_fd = psi.copy()
    for m in range(n_steps + 1):
        col = psi[:, None]
        rec.record(m, col[:n_r], col[n_r:])
        d = psi[:n_r][real]
        arms = psi[n_r:].reshape(n_arm, n_r)[:, real]
        # The one expression the recorder does NOT provide, kept independent on
        # purpose: `4 sum_n Re[conj(Psi_d) phi_n] Im[V_dn]` is what the two
        # rates are claimed to sum to, so computing it from the recorder would
        # make the identity below vacuous. CONJUGATING product, restricted to
        # the real nuclear region: a probability rate, not a c-product.
        residual[m] = 4.0 * float(np.sum((d.conj()[None, :] * arms).real * v_real.imag))
        if m == _FD_TIME:
            psi_fd = psi.copy()
        if m < n_steps:
            psi = step(psi)

    exch_d = rec.exchange[:, 0]
    imbalance = rec.imbalance[:, 0]
    exch_arm = imbalance - exch_d  # the recorder reports the sum, not the term
    scale = float(np.abs(imbalance).max())

    # The Markovian rate on a REAL deck, with `local_width`'s unclamped
    # round-off negatives in it: it stays strictly negative anyway (measured
    # max -3.85e-6 over this window), which is what makes a POSITIVE nonlocal
    # `exchange` a result and not a sign convention.
    # `test_the_markovian_exchange_is_never_positive` gates the arithmetic and
    # the absence of a clamp; this gates it against the actual `Gamma_loc`.
    # NOTE the negatives are bigger here than the 32-node / -3.8e-9 set
    # `local_width`'s docstring records: that was measured on `width_deck`'s
    # converged r_max = 16 electronic box, and this deck's r_max = 11 one puts
    # `Gamma_loc` 22% high (13 nodes, worst -1.8e-6). Same conclusion, and the
    # reason the campaign decks must use a converged box.
    assert np.all(rec.exchange_local < 0.0), f"max {rec.exchange_local.max():.3e}"
    # The algebraic half of the identity, and the only part that is exact: the
    # sum of the two rates IS the Im(V_dn) term. Round-off, nothing else.
    assert np.abs(imbalance - residual).max() < 1e-12 * scale

    live = np.maximum(np.abs(exch_d), np.abs(exch_arm)) > 0.0
    frac = np.abs(imbalance[live]) / np.maximum(np.abs(exch_d[live]), np.abs(exch_arm[live]))
    median, worst = float(np.median(frac)), float(np.max(frac))
    # Bounds with headroom over the measured 0.82 / 1.06. The LOWER bound is
    # the load-bearing one: it is what forbids calling `‖phi_n‖^2` a
    # population, and a change that made this residual small would be a change
    # in what the observables of `memory.py` are allowed to claim -- it must
    # not land silently.
    assert median > 0.5, f"median imbalance fraction {median:.4f} -- re-read memory.py"
    assert worst < 1.5, f"max imbalance fraction {worst:.4f}"

    _assert_the_rates_are_what_the_propagation_does(h_ext, psi_fd, n_r, n_arm, v_dn)


# --- the recorder ------------------------------------------------------------


def _absorbing_symmetric(rng, n):
    """A complex-SYMMETRIC, strictly absorbing `H_ext` of size `n`.

    Symmetric because that is the structure ECS gives (`.T`, never
    `.conj().T`). The `-i` shift is computed from the spectrum rather than
    guessed: a fixed shift that is ample at one size leaves a genuinely growing
    mode at another (`a + a.T`'s largest imaginary eigenvalue grows like
    `sqrt(n)`), and `propagate_nrm` would then warn -- correctly -- that the
    fixture, not the recorder, is broken.
    """
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    h = a + a.T
    shift = float(np.max(np.linalg.eigvals(h).imag)) + 1.0
    return sp.csr_matrix(h - 1j * shift * np.eye(n))


def _synthetic_deck(seed=5, n_arm=3):
    """A small grid + a random absorbing arrow-free `H_ext` on it.

    The Hamiltonian is deliberately NOT `extended_hamiltonian`'s arrow shape:
    the recorder reads the coupling blocks straight out of whatever `h_ext` it
    is given (`h_ext[:n_r, n_r:]` and `h_ext[n_r:, :n_r]`), so a dense random
    coupling exercises that path rather than the diagonal special case, and a
    recorder that quietly assumed `sp.diags(V_dn)` would fail here.
    """
    rng = np.random.default_rng(seed)
    nuc = segmented_grid(((2, 2.0), (2, 4.0)), ((2, 8.0),), angle_deg=35.0, quadrature=5)
    n_ext = (1 + n_arm) * nuc.n
    h = _absorbing_symmetric(rng, n_ext)
    psi0 = np.zeros((n_ext, 2), dtype=np.complex128)
    psi0[: nuc.n, :] = rng.normal(size=(nuc.n, 2)) + 1j * rng.normal(size=(nuc.n, 2))
    launch = LaunchBasis(
        vectors=psi0,
        coeffs=np.eye(2, dtype=np.complex128),
        energies=np.array([0.05, 0.09]),
        e_total=np.array([0.05, 0.09]),
        truncation_error=0.0,
    )
    return nuc, h, launch, n_arm


def test_memory_off_changes_nothing():
    """The default path is untouched, and switching memory ON does not perturb it.

    Two claims in one, because "before" is not a thing a test can hold: the
    `memory=None` run must report every new field as `None` while still
    producing the diagnostics it always did, and the `memory=spec` run must
    reproduce `psi_d`/`survival`/`centroid`/`momentum` BIT FOR BIT. Bitwise,
    not `allclose` -- the recorder is meant to observe the propagation, and an
    observation that moves the last digit of `psi_d` is a recorder that has
    got into the arithmetic. Every gate in this suite propagates, so a
    regression here is a regression everywhere.
    """
    nuc, h, launch, n_arm = _synthetic_deck()
    off = propagate_nrm(h, launch, nuc, dt=0.05, n_steps=60)
    assert off.arm_norm is None
    assert off.arm_norm_by_channel is None
    assert off.arm_peak is None
    assert off.exchange is None
    assert off.exchange_local is None
    assert off.imbalance is None

    gamma = np.linspace(0.0, 0.02, nuc.n)
    on = propagate_nrm(h, launch, nuc, dt=0.05, n_steps=60, memory=MemorySpec(gamma))
    assert np.array_equal(on.psi_d, off.psi_d)
    assert np.array_equal(on.survival, off.survival)
    assert np.array_equal(on.centroid, off.centroid)
    assert np.array_equal(on.momentum, off.momentum)
    assert np.array_equal(on.time, off.time)

    assert on.exchange is not None
    assert on.exchange.shape == (61, 2)
    assert on.arm_norm_by_channel is not None
    # `n_channels=None` (the default) resolves EVERY channel.
    assert on.arm_norm_by_channel.shape == (61, n_arm, 2)
    assert on.arm_peak is not None
    assert on.arm_peak.shape == (n_arm, 2)

    # Truncating the time series must not truncate `arm_peak`, which is the
    # field that tells a caller whether the kept blocks were the right ones.
    few = propagate_nrm(h, launch, nuc, dt=0.05, n_steps=60, memory=MemorySpec(gamma, 2))
    assert few.arm_norm_by_channel is not None and few.arm_peak is not None
    assert few.arm_norm_by_channel.shape == (61, 2, 2)
    assert few.arm_peak.shape == (n_arm, 2)
    assert np.array_equal(few.arm_peak, on.arm_peak)
    assert np.array_equal(few.arm_norm_by_channel, on.arm_norm_by_channel[:, :2])
    # And it IS the running maximum of the per-channel series.
    assert np.allclose(on.arm_peak, on.arm_norm_by_channel.max(axis=0))

    # WHICH BLOCK IS WHICH. The only assertion in the file that ties a
    # per-channel series to a SPECIFIC arm block of `h_ext`, and it is
    # load-bearing for `arm_peak`, whose whole job is to say which blocks
    # carried the flux -- Task 3 reads it to decide whether a positional
    # truncation kept the right ones.
    #
    # Everything else here is permutation-invariant and cannot see a mis-set
    # channel index: `arm_norm` is a total over all arms; `exchange` and
    # `imbalance` use the same index set on both sides of the coupling; and the
    # `arm_peak` check above compares it against `arm_norm_by_channel`, i.e. to
    # itself. Mutation-tested 2026-08-27: reversing the recorder's channel
    # order (`np.arange(n_arm)` -> `np.arange(n_arm)[::-1]`, memory.py's
    # `arm_real_idx`) passed all 8 tests before this line and fails here after
    # it. n = 2 is the block a reversal moves; step 1 rather than 0, since
    # every arm is identically zero at launch.
    arms_1 = make_pade_stepper(h, 0.05, 3)(launch.vectors)[nuc.n :, :] @ launch.coeffs
    real_idx = np.flatnonzero(nuc.real_points <= nuc.R0)
    block_2 = arms_1[2 * nuc.n + real_idx, :]
    assert np.allclose(on.arm_norm_by_channel[1, 2], (np.abs(block_2) ** 2).sum(axis=0), rtol=1e-12)


def _two_state_deck(v_coup=0.3):
    """`H_ext = [[0, V], [V, 0]]` on a ONE-node nuclear grid: n_r = 1, one arm.

    A single GLL element at `quadrature=3` retains exactly one basis function
    after the Dirichlet drop, and its node sits at R = 0.5 < R0 = 1.0, so it is
    a real-region node and the recorder's real-region restriction keeps it.
    That is what makes the closed form below the WHOLE observable rather than
    one term of it.
    """
    grid = FemDvrEcsGrid(GridSpec(quadrature=3, elements=[ElementSpec(1.0, 0.0)]))
    assert grid.n == 1 and grid.real_points[0] <= grid.R0
    h = sp.csr_matrix(np.array([[0.0, v_coup], [v_coup, 0.0]], dtype=np.complex128))
    psi0 = np.array([[1.0 + 0.0j], [0.0 + 0.0j]])
    launch = LaunchBasis(
        vectors=psi0,
        coeffs=np.eye(1, dtype=np.complex128),
        energies=np.array([0.0]),
        e_total=np.array([0.0]),
        truncation_error=0.0,
    )
    return grid, h, launch


def test_exchange_matches_a_direct_two_state_calculation():
    """A 2x2 `H_ext` with one arm has a closed form, and the recorder must hit it.

    For `i d/dt Psi = H Psi` with `H = [[0, V], [V, 0]]` real and
    `Psi(0) = [1, 0]`, `Psi(t) = [cos(Vt), -i sin(Vt)]`. Then

        exchange(t)   = 2 Im<Psi_d|V phi> = -V sin(2Vt),
        arm_norm(t)   = sin^2(Vt),
        exchange_local(t) = -Gamma cos^2(Vt),
        imbalance(t)  = 0 exactly, since V is REAL.

    and `-V sin(2Vt)` is `d/dt cos^2(Vt)`, i.e. the exchange really is the rate
    at which the discrete block's conjugating norm changes. This is arithmetic:
    it fixes the factor of 2, the sign, and the conjugation side, and it does
    not depend on any physics being right. `t` stops at `V t = 1.5 < pi/2` so
    `survival` falls monotonically and `_warn_if_diverging` -- which is looking
    at an undamped, genuinely oscillating norm here -- has nothing to say.
    """
    v_coup = 0.3
    gamma = 0.017
    grid, h, launch = _two_state_deck(v_coup)
    spec = MemorySpec(gamma_local=np.array([gamma]), n_channels=1)
    res = propagate_nrm(h, launch, grid, dt=0.05, n_steps=100, memory=spec)

    t = res.time
    assert res.exchange is not None
    assert res.arm_norm is not None
    assert res.exchange_local is not None
    assert res.imbalance is not None
    assert np.allclose(res.exchange[:, 0], -v_coup * np.sin(2 * v_coup * t), atol=1e-9)
    assert np.allclose(res.arm_norm[:, 0], np.sin(v_coup * t) ** 2, atol=1e-9)
    assert np.allclose(res.exchange_local[:, 0], -gamma * np.cos(v_coup * t) ** 2, atol=1e-9)
    assert np.max(np.abs(res.imbalance[:, 0])) < 1e-14
    assert res.arm_peak is not None
    assert res.arm_peak[0, 0] == pytest.approx(np.sin(v_coup * t[-1]) ** 2, abs=1e-9)
    # And the exchange IS d/dt of the discrete block's norm here (nothing else
    # touches it): a central difference of `survival` reproduces it.
    fd = (res.survival[2:, 0] - res.survival[:-2, 0]) / (2 * 0.05)
    assert np.allclose(fd, res.exchange[1:-1, 0], atol=1e-4)


def test_the_markovian_exchange_is_never_positive():
    """`-<Psi_d|Gamma_loc|Psi_d>` with `Gamma_loc >= 0` is a loss, never a gain.

    This is what makes a POSITIVE nonlocal `exchange` a result rather than a
    sign convention, so it is asserted rather than assumed: the two curves are
    read against each other, and if the Markovian one could be positive the
    comparison would say nothing.

    The second half is a NO-CLAMP gate. `local_width` deliberately does not
    clamp its round-off negatives (see its docstring: 32 of 179 nodes on the N2
    deck, worst -3.8e-9, every one where the true width is zero), because a
    negative entry is the one diagnostic worth keeping. A recorder that clamped
    -- "to be safe" -- would hide exactly that, so a deliberately negative
    `gamma_local` is required to come back POSITIVE here.
    """
    nuc, h, launch, _n_arm = _synthetic_deck(seed=9)
    gamma = np.abs(np.linspace(-0.01, 0.03, nuc.n))
    res = propagate_nrm(
        h, launch, nuc, dt=0.05, n_steps=60, memory=MemorySpec(gamma_local=gamma, n_channels=2)
    )
    assert res.exchange_local is not None
    assert np.all(res.exchange_local <= 0.0)
    assert np.all(res.exchange_local[0] < 0.0), "a live packet against a nonzero width"

    flipped = propagate_nrm(
        h, launch, nuc, dt=0.05, n_steps=60, memory=MemorySpec(gamma_local=-gamma, n_channels=2)
    )
    assert flipped.exchange_local is not None
    assert np.all(flipped.exchange_local >= 0.0)
    assert np.allclose(flipped.exchange_local, -res.exchange_local, rtol=1e-12)


def test_the_recorder_rejects_what_it_cannot_measure():
    """Every guard, since each one stands for a way of getting silently wrong
    numbers rather than an error: a `gamma_local` from another grid would be
    broadcast against the wrong nodes; a complex one would make the Markovian
    "rate" complex; and `nuclear_grid=None` has no real region to restrict to,
    so the observables would quietly include the ECS absorber."""
    nuc, h, launch, _n_arm = _synthetic_deck()
    with pytest.raises(ValueError, match="one per nuclear DVR node"):
        propagate_nrm(h, launch, nuc, dt=0.05, n_steps=2, memory=MemorySpec(np.zeros(nuc.n - 1)))
    with pytest.raises(ValueError, match="must be real"):
        MemorySpec(np.zeros(nuc.n, dtype=np.complex128))
    with pytest.raises(ValueError, match="must be 1-D"):
        MemorySpec(np.zeros((nuc.n, 2)))
    with pytest.raises(ValueError, match="n_channels must be >= 0 or None"):
        MemorySpec(np.zeros(nuc.n), n_channels=-1)
    with pytest.raises(ValueError, match="requires a nuclear_grid"):
        propagate_nrm(h, launch, None, dt=0.05, n_steps=2, memory=MemorySpec(np.zeros(nuc.n)))
