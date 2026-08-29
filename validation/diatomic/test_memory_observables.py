"""What the memory-observable campaign asserts, and what it deliberately does not.

The campaign itself (`memory_observables.py`) is a ~1.5 hour run on three
production decks; nothing here reruns it. What is gated is the part that is
cheap and DURABLE:

  * the SIGN STRUCTURE of the exchange rate -- the Markovian limit can only
    lose, by construction, and it does so on real decks with `local_width`'s
    unclamped round-off negatives in it;
  * that `Gamma_loc` is a real, positive-peaked width on all three campaign
    decks, so `exchange_local` and `Gamma_eff` are defined for each;
  * that `summarize` computes the normalizations its FIELD NAMES claim, since
    "the number was quoted without its normalization" is the specific error
    this sub-project keeps making;
  * that a command-line override cannot overwrite the campaign `.npz` it is
    supposed to be checked against.

NOT gated: exact peak times, the number of returning steps, or the ordering
itself. Those are measurements of a particular deck at a particular energy over
a particular window, they take an hour and a half to reproduce, and freezing
them here would turn a finding into a fixture.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import pytest
from qscat.core.nrm.extended import extended_hamiltonian, initial_packet
from qscat.core.nrm.memory import MemorySpec, local_width
from qscat.core.nrm.nonlocal_potential import nonlocal_operator
from qscat.core.nrm.propagation import propagate_nrm

from validation.diatomic.memory_observables import (
    COARSE_GRAINED_RETURN,
    DECK_COST,
    DECKS,
    E_BOX_LADDER,
    ENERGY_LADDER,
    ENERGY_WINDOWS,
    N2_RESOLUTION,
    N2_RESOLUTION_FLOOR,
    RESOLVED_RETURN,
    _build,
    _stem,
    replace_energy,
    replace_window,
    summarize,
)

# A short N2 propagation: 40 steps is past the launch transient (the arms fill
# within ~10 a.u.) and nowhere near the late-time regime where the campaign's
# returning flux appears. That is on purpose -- the claim gated below is the
# one that holds at EVERY time, not the one that needs 4000 steps.
_SHORT_STEPS = 40


@pytest.fixture(scope="module")
def n2_short() -> dict[str, npt.NDArray[np.float64]]:
    """40 steps of the N2 campaign deck, in `summarize`'s input shape."""
    deck = DECKS["N2"]
    nuc, elec, phi_d, ing, eps, chi = _build(deck)
    e_total = deck.e_kin + float(eps[0])
    gamma = local_width(nonlocal_operator(ing, nuc, deck.model, e_total, n_states=None), nuc)
    h_ext = extended_hamiltonian(ing, nuc, deck.model)
    launch = initial_packet(
        nuc, elec, deck.model, phi_d, ing, eps, chi, 0, np.array([deck.e_kin]), rank_tol=1e-10
    )
    res = propagate_nrm(
        h_ext,
        launch,
        nuc,
        dt=deck.dt,
        n_steps=_SHORT_STEPS,
        order=3,
        # `n_states=None` above and every channel here: a truncated arm set
        # makes `H_ext` non-dissipative, and the recorder's cost is measured
        # identical at `None` and at 4.
        memory=MemorySpec(gamma_local=gamma, n_channels=None),
    )
    assert res.arm_norm is not None
    assert res.arm_norm_by_channel is not None
    assert res.arm_peak is not None
    assert res.exchange is not None
    assert res.exchange_local is not None
    assert res.imbalance is not None
    n_arm = res.arm_peak.shape[0]
    return {
        "molecule": np.array("N2"),
        "e_kin": np.array(deck.e_kin),
        "dt": np.array(deck.dt),
        "min_overlap": np.array(ing.min_overlap),
        "n_arm": np.array(n_arm),
        "n_ext": np.array(h_ext.shape[0]),
        "R_real": nuc.points[nuc.points.imag == 0.0].real,
        "gamma_local": gamma,
        "time": res.time,
        "survival": res.survival[:, 0],
        "arm_norm": res.arm_norm[:, 0],
        "arm_channel_index": np.arange(n_arm),
        "arm_norm_by_channel": res.arm_norm_by_channel[:, :, 0],
        "arm_peak": res.arm_peak[:, 0],
        "exchange": res.exchange[:, 0],
        "exchange_local": res.exchange_local[:, 0],
        "imbalance": res.imbalance[:, 0],
    }


@pytest.mark.slow
def test_the_markovian_limit_can_only_lose(n2_short) -> None:
    """`exchange_local <= 0` at every step -- the claim the whole comparison rests on.

    `-<Psi_d|Gamma_loc|Psi_d>` is non-positive wherever `Gamma_loc` is, and
    `local_width` is deliberately NOT clamped: on this deck 32 of 179 nodes
    come back negative from round-off, the worst at -3.8e-9 against a peak of
    2.5e-2. This gates that they stay round-off -- that the Markovian rate is
    negative on a REAL deck rather than only in the algebra -- which is what
    makes a POSITIVE nonlocal `exchange` a result and not a sign convention.
    Step 0 is included: the arms are empty there, so `exchange` is exactly 0
    while `exchange_local` already is not.
    """
    exl = n2_short["exchange_local"]
    assert np.all(exl <= 0.0), f"Markovian rate went positive: max {exl.max():.3e}"
    assert exl.min() < 0.0, "Markovian rate is identically zero -- Gamma_loc is not being read"


@pytest.mark.slow
def test_the_exchange_starts_at_exactly_zero_and_then_loses(n2_short) -> None:
    """`exchange[0] == 0` (empty arms) and the launch transient is a LOSS.

    The first is exact, not approximate: at `t = 0` every `phi_n` is zero, so
    `2 Im<Psi_d|sum_n V_dn phi_n>` has nothing to sum. It matters because a
    "positive step" count must not include it -- the campaign counts strict
    positives for exactly this reason.
    """
    ex = n2_short["exchange"]
    assert ex[0] == 0.0
    assert ex[1:].max() <= 0.0, (
        "the nonlocal exchange is already positive within the launch transient, "
        "which no deck has done -- the returning flux is a LATE-time feature"
    )


@pytest.mark.slow
def test_the_imbalance_stays_O_1_so_the_arms_are_not_a_population(n2_short) -> None:
    """Task 1's residual is the size of the exchange, not a correction to it.

    `4 sum_n Re[conj(Psi_d) phi_n] Im[V_dn]` vanishes only where `V_dn` is
    real, and under ECS it is not. Measured median 0.822 of the larger
    one-sided rate on the gate deck over 200 steps; this gates that it is
    O(1) rather than the exact value, because the value moves with the window
    and with the electronic box. If this ever came back small, `arm_norm`
    could be read as a population and every label in the campaign would be
    wrong in the other direction.
    """
    sm = summarize(n2_short)
    assert sm.imbalance_median_rel > 0.1, (
        f"imbalance is only {sm.imbalance_median_rel:.3f} of the larger rate -- "
        "re-measure before calling arm_norm a relative decomposition"
    )


@pytest.mark.slow
def test_every_campaign_deck_has_a_usable_local_width() -> None:
    """`Gamma_loc` is real, peaked and positive on all three decks, and its
    negatives are round-off.

    No propagation: this is the ingredient build plus one `F(E)`, which is what
    makes it affordable for all three. It gates the precondition every
    observable in the campaign shares -- `exchange_local` and `Gamma_eff` both
    divide the run by this array -- and it gates it per molecule rather than
    carrying N2's answer across.
    """
    for name, deck in DECKS.items():
        nuc, _elec, _phi_d, ing, eps, _chi = _build(deck)
        e_total = deck.e_kin + float(eps[0])
        gamma = local_width(nonlocal_operator(ing, nuc, deck.model, e_total, n_states=None), nuc)
        assert not np.iscomplexobj(gamma)
        peak = float(gamma.max())
        assert peak > 0.0, f"{name}: Gamma_loc has no positive peak"
        worst = float(gamma.min())
        assert worst > -1e-4 * peak, (
            f"{name}: most negative Gamma_loc is {worst:.3e} against a peak of "
            f"{peak:.3e} -- that is too large to be round-off"
        )
        assert ing.min_overlap >= 0.5, (
            f"{name}: nrm_ingredients tracking overlap {ing.min_overlap:.3f} -- the "
            "adiabatic walk paired the wrong P-space state, so E_n/V_dn are mixed"
        )


def test_each_deck_uses_a_box_the_ladder_measured() -> None:
    """Every campaign deck's electronic box is a rung of its own ladder, and a
    rung that tracks cleanly.

    The failure this catches is the one the sub-project already made once:
    inheriting a gate deck's electronic box because it was there. The gate
    decks' boxes are IN the ladder (N2 at `r_max = 11` reads 1.218 x the Eq.
    (68) width, F2's TD-figure box 1.545 x), so this also fails if a deck is
    quietly moved back to one.
    """
    for name, deck in DECKS.items():
        rung = (deck.e_r_max, deck.e_order, deck.e_n_complex)
        ladder = E_BOX_LADDER[name]
        assert rung in ladder, f"{name}: box {rung} was never measured -- run --converge"
        peak, ratio, overlap = ladder[rung]
        assert overlap >= 0.5, f"{name}: chosen rung tracks at {overlap}"
        assert 0.9 < ratio < 1.1, f"{name}: chosen rung reads {ratio} x the Eq. (68) width"
        # And it agrees with the largest rung measured -- the convergence
        # claim itself, read off the recorded table rather than re-run.
        largest = max(ladder, key=lambda k: (k[0], k[1], k[2]))
        assert abs(peak - ladder[largest][0]) < 0.02 * ladder[largest][0], (
            f"{name}: peak Gamma_loc {peak:.4e} at the chosen box against "
            f"{ladder[largest][0]:.4e} at {largest}"
        )


def _synthetic() -> dict[str, npt.NDArray[np.float64]]:
    """A hand-built campaign record with known answers.

    `survival` halves every 10 a.u.; `exchange` is negative except for one
    positive step; the arm norm peaks at step 5 with a known channel split.
    """
    t = np.arange(0.0, 41.0)
    s = 2.0 * 0.5 ** (t / 10.0)
    ex = -np.ones_like(t) * 1e-3
    ex[0] = 0.0
    ex[20] = +4e-4
    exl = -np.ones_like(t) * 2e-3
    arm_c = np.zeros((t.size, 4))
    arm_c[:, 0] = 0.1
    arm_c[:, 1] = 0.2
    arm_c[:, 2] = 0.3
    arm_c[:, 3] = 0.4
    arm_c[5] = [1.0, 2.0, 3.0, 4.0]
    return {
        "molecule": np.array("SYN"),
        "e_kin": np.array(0.1),
        "dt": np.array(1.0),
        "min_overlap": np.array(1.0),
        "n_arm": np.array(4),
        "n_ext": np.array(40),
        "R_real": np.zeros(7),
        "gamma_local": np.ones(7),
        "time": t,
        "survival": s,
        "arm_norm": arm_c.sum(axis=1),
        # Block order deliberately NOT 0,1,2,3: the largest are 3 and 2, and
        # `summarize` has to find blocks 0-3 through the index map rather than
        # by taking the first four columns.
        "arm_channel_index": np.array([3, 2, 1, 0]),
        "arm_norm_by_channel": arm_c[:, ::-1],
        "arm_peak": np.array([1.0, 2.0, 3.0, 4.0]),
        "exchange": ex,
        "exchange_local": exl,
        "imbalance": ex * 1.5,
    }


def test_summarize_computes_the_normalizations_its_field_names_claim() -> None:
    """`_raw`, `_over_s0` and `_over_s` are three different numbers, and each
    field is the one it says.

    The specific error this gates is the one the design document made for
    several revisions: quoting `+2.420e-4` (a `/S_d(0)` figure) as if it were
    the raw `+8.776e-7`. Two agents then reproduced the normalized number
    while checking against a figure that never said it was normalized.
    """
    sm = summarize(_synthetic())
    assert sm.s_d0 == 2.0
    assert sm.exchange_max_raw == pytest.approx(4e-4)
    assert sm.exchange_max_over_s0 == pytest.approx(4e-4 / 2.0)
    # /S_d(t) at t = 20, where S has fallen to a quarter of S(0).
    assert sm.exchange_max_over_s == pytest.approx(4e-4 / (2.0 * 0.25))
    assert sm.exchange_min_raw == pytest.approx(-1e-3)
    assert sm.exchange_min_over_s0 == pytest.approx(-5e-4)
    assert sm.n_positive == 1
    assert sm.t_first_positive == 20.0
    # `Gamma_eff` is the golden-rule rate from the SAME `Gamma_loc` the
    # Markovian curve uses: -exchange_local(0)/S_d(0).
    assert sm.gamma_eff == pytest.approx(2e-3 / 2.0)


def test_summarize_finds_blocks_0_to_3_through_the_stored_index_map() -> None:
    """`arm_first_four_share` is about `h_ext`'s block ORDER, not column order.

    `MemorySpec.n_channels=4` keeps the first four arm BLOCKS. The `.npz`
    stores a union of the largest and the first four, in ascending block
    order, so reading its first four columns would silently answer a different
    question on any deck where the two sets differ.
    """
    sm = summarize(_synthetic())
    # Every block 0-3 is present, so the share is the whole arm norm.
    assert sm.arm_first_four_share == pytest.approx(1.0)
    assert sm.arm_top_blocks == (3, 2, 1, 0)
    assert sm.arm_peak_time == 5.0
    assert sm.arm_peak_over_s0 == pytest.approx(10.0 / 2.0)


def test_the_decay_law_is_read_against_one_rate_constant() -> None:
    """`decay_times` are first crossings and `decay_ratio_at_levels` compares
    them to `exp(-Gamma_eff t)`.

    The synthetic `S` halves every 10 a.u. (rate ln2/10 = 0.0693) against a
    `Gamma_eff` of 1e-3, so the packet decays far FASTER than the golden rule
    here and every ratio is below 1. The campaign's molecules do the opposite
    at late times, and the ratio is the number that says by how much.
    """
    sm = summarize(_synthetic())
    assert sm.decay_times[0] == pytest.approx(10.0)
    assert sm.decay_times[1] == pytest.approx(34.0)  # first t with S/S0 <= 0.1
    assert np.isnan(sm.decay_times[2])  # 0.01 is not reached inside 40 a.u.
    assert sm.decay_ratio_at_levels[0] == pytest.approx(0.5 / np.exp(-1e-3 * 10.0))
    assert np.isnan(sm.decay_ratio_at_levels[2])


def test_an_override_cannot_overwrite_the_campaign_npz() -> None:
    """An energy rung or a `dt` check writes its own file.

    Both are things this campaign actually runs -- the F2 energy ladder either
    side of the LCP's unity crossing, and the `dt`-halving check on N2 -- and
    both would otherwise land on the campaign deck's own stem and be
    indistinguishable from it afterwards.
    """
    deck = DECKS["F2"]
    assert _stem(deck) == "f2-nrm-memory-observables"
    assert _stem(replace_energy(deck, 0.010)) == "f2-nrm-memory-observables-e0.01"
    assert _stem(replace_window(deck, 8000, 0.5)) == "f2-nrm-memory-observables-dt0.5"
    assert _stem(deck, smoke=True) == "f2-nrm-memory-observables-n20"


def test_the_deck_energies_are_where_the_lcp_comparison_lives() -> None:
    """Each campaign energy sits inside the window its molecule's LCP failure
    was measured over.

    N2 at 0.10 Ha is inside the 0.06-0.16 Ha VE window (its DA channel opens
    only at +0.5016 Ha); F2 at 0.030 Ha is inside 0.010-0.050, next to the
    ratio's unity crossing at ~0.032; NO at 0.200 Ha is inside 0.150-0.300 and
    above NO's +0.1719 Ha DA threshold. An energy outside these is not
    comparable to anything the energy domain measured.
    """
    windows = {"N2": (0.06, 0.16), "F2": (0.010, 0.050), "NO": (0.1719, 0.300)}
    for name, (lo, hi) in windows.items():
        assert lo <= DECKS[name].e_kin <= hi, (
            f"{name} at {DECKS[name].e_kin} Ha is outside {lo}-{hi}"
        )


def test_replace_energy_keeps_everything_else() -> None:
    """`replace_energy` changes one field. It exists because a multi-energy
    `LaunchBasis` would share one `gamma_local` across columns at different
    energies, so an energy ladder must be separate propagations."""
    deck = DECKS["F2"]
    other = replace_energy(deck, 0.050)
    assert other.e_kin == 0.050
    assert other.nuc_real == deck.nuc_real
    assert other.e_r_max == deck.e_r_max
    assert other.n_steps == deck.n_steps


def test_the_recorded_deck_cost_is_self_consistent() -> None:
    """`DECK_COST`'s `H_ext` order is `(1 + n_arm) * N_R` for every molecule.

    The table exists so "run these on a big machine" is a number someone can
    plan against, and an inconsistent row is worse than no row -- it would be
    quoted. This checks the one relation the table cannot get wrong by
    accident, and that every campaign deck has an entry.
    """
    assert set(DECK_COST) == set(DECKS)
    for name, (n_r, n_arm, n_ext, nnz, peak_gb, _minutes) in DECK_COST.items():
        assert n_ext == (1 + n_arm) * n_r, f"{name}: {n_ext} != (1 + {n_arm}) * {n_r}"
        assert nnz > n_ext, f"{name}: fewer nonzeros than rows"
        assert 0.0 < peak_gb < 68.7, f"{name}: peak RSS {peak_gb} GB is not a laptop number"


def test_the_n2_returning_flux_is_below_its_own_resolution_floor() -> None:
    """The retraction, asserted from the recorded numbers rather than restated.

    This is the one place a reader can check the claim in `N2_RESOLUTION`'s
    comment without rerunning four propagations: the largest positive
    excursion any of the four N2 runs produced, as a fraction of that run's
    own exchange peak, against the floor between the two FINEST runs. If a
    future change to the propagator or the deck lifts the signal over that
    floor, this test fails and the retraction is due for revisiting -- which
    is the outcome it exists to notice.
    """
    finest = ("order 3, dt=0.5", "order 3, dt=0.25")
    floor, sign_agreement = N2_RESOLUTION_FLOOR[finest]
    worst_positive = max(N2_RESOLUTION[k][0] for k in N2_RESOLUTION)
    assert worst_positive < floor / 10.0, (
        f"N2's largest positive excursion is {worst_positive:.3e} of peak against a "
        f"{floor:.3e} floor -- it now clears it, so re-measure before quoting the "
        "retraction"
    )
    # And the sign pattern does not converge: every pair, including the
    # finest, disagrees on better than a quarter of the steps.
    assert all(agree < 0.8 for _floor, agree in N2_RESOLUTION_FLOOR.values())
    assert sign_agreement < 0.8


def test_the_integral_observables_do_converge_on_n2() -> None:
    """The other half of the verdict: what the retraction does NOT touch.

    `nonlocality` and the arm-norm peak agree across all four propagations,
    so the campaign keeps them. Asserting this beside the retraction is what
    stops the retraction being read as "the N2 run is unusable".
    """
    nonlocality = [v[2] for v in N2_RESOLUTION.values()]
    arm_peak = [v[3] for v in N2_RESOLUTION.values()]
    assert max(nonlocality) - min(nonlocality) < 0.01 * min(nonlocality)
    assert max(arm_peak) - min(arm_peak) < 0.01 * min(arm_peak)
    # The onset time, by contrast, moves by a factor of six -- the contrast is
    # the point, so it is asserted rather than left to the comment.
    onsets = [v[1] for v in N2_RESOLUTION.values()]
    assert max(onsets) / min(onsets) > 5.0


def test_the_return_clears_its_own_floor_only_on_no() -> None:
    """The pointwise discriminator, on all three molecules, measured alike.

    Each molecule's largest return is read against the floor for its OWN
    returning window, so the three ratios are comparable even though the raw
    rates differ by four orders. They separate cleanly:

      NO 2.2-2.6x ABOVE its floor -- resolved
      F2 0.65x and 1.02x -- AT it, and the 1.02 is not a pass, it is a tie
      N2 0.17x and 0.07x -- an order of magnitude below

    F2's marginality is the honest reading and is asserted as marginality
    rather than rounded either way: refining `dt` twice does not move it off
    the floor, because the floor and the feature shrink together.

    THE CONCORDANCES ARE READ AGAINST THEIR NULLS, which this test once got
    wrong by asserting a coin-toss null of 0.5 on a CONDITIONAL fraction. With
    the right null -- the comparison run's own positive-step rate -- N2's 0.427
    is a lift of +0.171, a real signal rather than a failure, and F2's
    larger-looking 0.622 is +0.048. On this metric N2 outranks F2. The verdict
    never rested on that pair, so it does not move: NO's lift is +0.513 to
    +0.568, three to twelve times either of them.

    If this flips on N2 the retraction is due for revisiting; if it flips on
    NO, the one molecule whose returning flux is fully resolved has gone.
    """
    ratios: dict[str, list[float]] = {}
    lifts: dict[str, list[float]] = {}
    disagreements: dict[str, list[float]] = {}
    for (molecule, _), (peak_share, floor, conc, null, rel) in RESOLVED_RETURN.items():
        ratios.setdefault(molecule, []).append(peak_share / floor)
        lifts.setdefault(molecule, []).append(conc - null)
        disagreements.setdefault(molecule, []).append(rel)

    assert min(ratios["NO"]) > 2.0, f"NO no longer clears its floor: {ratios['NO']}"
    assert max(ratios["N2"]) < 0.5, f"N2's return rose toward its floor: {ratios['N2']}"
    # Marginal on both sides of unity -- neither "clears" nor "fails".
    assert 0.5 < max(ratios["F2"]) < 1.5, f"F2 stopped being marginal: {ratios['F2']}"

    # NO is the only molecule whose returning steps reproduce well above its
    # own null. N2 and F2 both show a small positive lift; the test does NOT
    # rank them against each other, because on this metric they rank the
    # opposite way to the binned one and neither gap is meaningful.
    assert min(lifts["NO"]) > 0.4, f"NO's lift fell to {lifts['NO']}"
    for molecule in ("N2", "F2"):
        assert max(lifts[molecule]) < 0.25, f"{molecule}'s lift rose to {lifts[molecule]}"
        assert min(lifts[molecule]) > 0.0, (
            f"{molecule}'s lift went non-positive ({lifts[molecule]}) -- it is small, "
            "but it has never been zero and 'at chance' would be the wrong word for it"
        )
    assert min(lifts["NO"]) > 2.0 * max(max(lifts["N2"]), max(lifts["F2"]))

    assert max(disagreements["NO"]) < min(min(disagreements["N2"]), min(disagreements["F2"])), (
        "NO's disagreement at its returning steps is no longer the smallest"
    )


def test_the_return_columns_are_frozen_in_energy() -> None:
    """Everything read as a RETURN describes the molecule, not the energy.

    Thirteen rungs across three molecules. Over a 4-6x change in the decay rate
    the onset does not move and `max positive / peak` moves in the fourth
    figure. Without this the per-molecule return numbers would be quotable only
    at their campaign energies.
    """
    by_molecule: dict[str, list[tuple[float, ...]]] = {}
    for (molecule, _), row in ENERGY_LADDER.items():
        by_molecule.setdefault(molecule, []).append(row)
    assert set(by_molecule) == {"N2", "NO", "F2"}
    for molecule, rows in by_molecule.items():
        peak_share = [r[3] for r in rows]
        assert max(peak_share) - min(peak_share) < 0.01 * min(peak_share), molecule


def test_nonlocality_inflates_near_threshold_and_the_post_peak_window_fixes_it() -> None:
    """The observable's near-threshold failure, and that the remedy removes it.

    With empty arms `|X - X_loc| = |X_loc|` identically, so the ratio is pinned
    near 1 regardless of the kernel. Near a threshold the Markovian REFERENCE
    collapses -- `int|X_loc|` falls by more than an order of magnitude across
    N2's and F2's ladders -- while that floor does not, so the FULL-RUN ratio
    inflates for reasons that are not about nonlocality.

    Integrating from the arm-norm peak removes the contaminated WINDOW rather
    than the contaminated rungs. The test of that is that the lowest rung of a
    ladder stops being an automatic outlier: on the full-run column it is the
    maximum for its molecule, and post-peak N2's is not.
    """
    for molecule in ("N2", "F2"):
        ref = [row[2] for (mol, _), row in ENERGY_LADDER.items() if mol == molecule]
        assert max(ref) / min(ref) > 10.0, f"{molecule}'s Markovian reference stopped collapsing"

    # Read over the IN-WINDOW rungs, since the out-of-window ones are not part
    # of any claim. On the full-run column each ladder's LOWEST energy is its
    # maximum -- that is the inflation, visible with no criterion at all.
    in_window: dict[str, list[tuple[float, tuple[float, ...]]]] = {}
    for (molecule, e_kin), row in ENERGY_LADDER.items():
        lo, hi = ENERGY_WINDOWS[molecule]
        if lo <= e_kin <= hi:
            in_window.setdefault(molecule, []).append((e_kin, row))
    for molecule in ("N2", "F2"):
        rows = sorted(in_window[molecule])
        assert rows[0][1][0] == max(r[1][0] for r in rows), (
            f"{molecule}'s lowest in-window rung is no longer its full-run maximum"
        )

    # Post-peak the picture inverts on N2: its maximum moves to its HIGHEST
    # energy. F2's stays at its lowest -- that rung sits at the DA threshold
    # and is the most nonlocal even post-peak, which is physics rather than
    # arithmetic, and it is why F2's floor rather than its ceiling faces NO.
    n2 = sorted(in_window["N2"])
    assert n2[-1][1][1] == max(r[1][1] for r in n2), (
        "N2's post-peak maximum left its highest energy -- the launch-transient "
        "inflation may be back"
    )


def test_nonlocality_orders_the_three_molecules_in_window() -> None:
    """N2 < NO < F2 on every in-window rung, read post-peak.

    The only exclusions are the two rungs this ladder added outside the
    molecules' own declared energy windows -- a criterion this module has
    carried since before the ladder existed, applied to BOTH offenders rather
    than to whichever was inconvenient. An earlier construction excluded four
    rungs by a cut invented after the fact; it is gone.

    `NO < F2` is asserted twice: here, and on the RAW full-run column over all
    seventeen rungs where it holds with no window argument at all. It is the
    inequality to quote if only one can be.
    """
    valid: dict[str, list[float]] = {}
    for (molecule, e_kin), row in ENERGY_LADDER.items():
        lo, hi = ENERGY_WINDOWS[molecule]
        if lo <= e_kin <= hi:
            valid.setdefault(molecule, []).append(row[1])
    assert {k: len(v) for k, v in sorted(valid.items())} == {"F2": 5, "N2": 6, "NO": 4}

    assert max(valid["N2"]) < min(valid["NO"]), f"N2 {valid['N2']} no longer below NO {valid['NO']}"
    assert max(valid["NO"]) < min(valid["F2"]), f"NO {valid['NO']} no longer below F2 {valid['F2']}"
    # Margins, not just signs. N2/NO is the narrow one and the first to watch.
    assert min(valid["NO"]) / max(valid["N2"]) - 1.0 > 0.05, "the N2/NO margin closed"
    assert min(valid["F2"]) / max(valid["NO"]) - 1.0 > 0.10, "the NO/F2 margin closed"

    raw_no = [row[0] for (mol, _), row in ENERGY_LADDER.items() if mol == "NO"]
    raw_f2 = [row[0] for (mol, _), row in ENERGY_LADDER.items() if mol == "F2"]
    assert max(raw_no) < min(raw_f2), (
        "NO < F2 no longer holds on the raw full-run column, which is the one "
        "statement here needing no argument about windows"
    )


def test_nos_flatness_is_not_evidence_because_nothing_moved() -> None:
    """NO's ladder never applied the perturbation, and that is recorded.

    NO's `nonlocality` is flat to 0.3% over its whole declared window -- but its
    `Gamma_eff` and Markovian reference are flat too, where N2's and F2's move
    by factors of 4-35. Flatness of an output under an input that did not vary
    is not a measurement, and this gate exists so that "NO's memory is
    energy-independent" cannot be reintroduced from the output column alone.
    """
    ref = [row[2] for (mol, _), row in ENERGY_LADDER.items() if mol == "NO"]
    assert max(ref) / min(ref) < 1.1, (
        f"NO's Markovian reference now varies by {max(ref) / min(ref):.2f}x -- if this is "
        "real, its flat nonlocality becomes a measurement and the note should say so"
    )


def test_the_coarse_grained_return_separates_no_from_the_other_two() -> None:
    """Two bands, at every bin width, against every rerun.

    The pointwise sign of `exchange` is not converged on any of the three
    molecules -- refining `dt` keeps the sign-flip period at roughly two steps
    rather than fixing it in atomic units -- so the returning flux is compared
    on time-averaged bins instead (`coarse_grained_return`). Two things this
    gate exists to hold fixed, both of which were once got wrong here:

    THE NULL IS NOT ONE HALF. The concordance is conditional, so its null is
    the comparison run's own positive-bin rate. F2's is ~0.59; reading its raw
    0.638-0.750 against 0.5 put it in a middle band it does not occupy.

    SIGN IS NOT ENOUGH. Two runs can agree about which bins return and disagree
    by a factor of five about how much. Both columns are asserted.
    """
    for width in (5.0, 10.0, 20.0, 50.0):
        n2_lift = COARSE_GRAINED_RETURN[("N2", width)][1]
        f2_lift = COARSE_GRAINED_RETURN[("F2", width)][1]
        no_lift = COARSE_GRAINED_RETURN[("NO", width)][0]
        # NO is resolved; the other two are at chance and are not ranked
        # against each other, because the measurement does not separate them.
        assert no_lift > 0.5, f"bin {width}: NO's worst lift fell to {no_lift}"
        assert n2_lift < 0.2, f"bin {width}: N2's best lift rose to {n2_lift}"
        assert f2_lift < 0.2, f"bin {width}: F2's best lift rose to {f2_lift}"
        assert no_lift > 3.0 * max(n2_lift, f2_lift), (
            f"bin {width}: the gap closed -- NO {no_lift} against N2 {n2_lift} / F2 {f2_lift}"
        )

        # Magnitude, the half that sign agreement cannot see.
        no_mag = COARSE_GRAINED_RETURN[("NO", width)][3]
        n2_mag = COARSE_GRAINED_RETURN[("N2", width)][2]
        f2_mag = COARSE_GRAINED_RETURN[("F2", width)][2]
        assert no_mag < 0.25, f"bin {width}: NO's binned returns now differ by {no_mag}"
        assert n2_mag > 0.5, f"bin {width}: N2's binned returns now agree to {n2_mag}"
        assert f2_mag > 0.5, f"bin {width}: F2's binned returns now agree to {f2_mag}"

    # A lift of zero can mean "agrees at chance" or "the other run has no
    # positive bins at all", so the bin counts travel with the table.
    for key, row in COARSE_GRAINED_RETURN.items():
        assert row[4] > 0, f"{key} records no returning bins"
