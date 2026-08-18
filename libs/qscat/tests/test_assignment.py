"""Tests for `qscat.core.assignment` -- pairing exact poles to BO states.

The end-to-end test here is the one that matters, and the separable limit makes
it exact. With `V(r, R) = v_el(r) + v_nuc(R)`:

- `exact_resonance_states` finds poles at `e_pole_el + eps_vib_v`, whose
  eigenvectors are exactly `phi_res (x) chi_v`;
- `resonance_curve` + `bo_basis` build those same products independently, from
  a pole walk rather than from a 2-D solve.

So `overlap` must return 1, `pair_by_overlap` must return `ok`, and the shift
must be zero -- not approximately, but to solver precision on the identical
grid. A pipeline that mislabels states, loses phase alignment, or pairs by the
wrong index fails this outright.

The verdict tests around it are unit tests with hand-built states, because a
spurious state is by definition one no physical model produces on demand.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core import exact_resonance_states
from qscat.core.assignment import (
    OverlapPair,
    overlap,
    pair_by_overlap,
    pair_one_to_one,
    peak_alignment,
    peak_positions,
)
from qscat.core.bo import BoBasis, BoState, bo_basis_from_levels, resonance_curve
from qscat.core.grids import ecs_angle_family, electronic_grid, nuclear_grid
from qscat.core.lcp import lcp_resonance_levels
from qscat.dvr import hamiltonian_nd
from qscat.linalg import c_product
from qscat.model import N2

_R_FIXED = 2.02
_MU = 12528.0
_D_MORSE, _A_MORSE, _RE_MORSE = 0.30, 1.0, 2.02


def _v_el(r):
    return np.asarray(N2.surface(r, _R_FIXED), dtype=np.complex128)


def _v_nuc(R):
    x = 1.0 - np.exp(-_A_MORSE * (np.asarray(R) - _RE_MORSE))
    return np.asarray(_D_MORSE * x * x - _D_MORSE, dtype=np.complex128)


class _SeparableModel:
    """`V(r, R) = v_el(r) + v_nuc(R)`, with a genuine electronic RESONANCE.

    Unlike `test_bo.py`'s binding well, this one carries N2's shape resonance at
    fixed `R`, which is what `resonance_curve` is built to follow.
    """

    charge = 0
    mu = _MU
    ell = 0

    def surface(self, r, R):
        return _v_el(r) + _v_nuc(R)

    def v0(self, R):
        return _v_nuc(R)

    def hamiltonian(self, tgrid):
        return hamiltonian_nd(tgrid, [1.0, self.mu], self.surface)


def _elec(theta: float):
    return electronic_grid(r_max=14.0, order=6, n_complex=5, angle_deg=theta)


def _nuc(theta: float):
    return nuclear_grid(r_max=16.0, quadrature=6, n_complex=3, angle_deg=theta)


_WINDOW = (-1.5, 0.5, -0.5, 0.0)
_SEED_WINDOW = (-1.0, 0.0, -0.1, 0.0)


@pytest.fixture(scope="module")
def separable_basis():
    """BO products on the separable model's RESONANCE curve.

    The nuclear side comes from `lcp_resonance_levels`, not `vibrational_states`:
    a resonance curve carries a width, so its levels are quasi-bound and complex.
    `bo_basis_from_levels` is the seam that lets the same comparator serve both
    the bound (ion / Rydberg) and the resonance (neutral / anion) case.
    """
    g_R_a, g_R_b = _nuc(25.0), _nuc(30.0)
    model = _SeparableModel()
    el_a, el_b = _elec(35.0), _elec(44.0)
    cur = resonance_curve(model, el_a, el_b, g_R_a, _SEED_WINDOW, with_states=True)
    cur_b = resonance_curve(model, el_a, el_b, g_R_b, _SEED_WINDOW, with_states=False)
    levels = lcp_resonance_levels(
        g_R_a,
        g_R_b,
        _MU,
        cur.energies[0].real.astype(np.complex128),
        cur_b.energies[0].real.astype(np.complex128),
        np.maximum(0.0, -2.0 * cur.energies[0].imag),
        n_levels=3,
    )
    return cur, bo_basis_from_levels(cur, levels.energies, levels.states)


@pytest.fixture(scope="module")
def separable_poles(separable_basis):
    """`exact_resonance_states` seeded at the BO levels of the same model."""
    _, basis = separable_basis
    base, moved_el, moved_nu = ecs_angle_family(
        _elec, _nuc, electronic_angles=(35.0, 44.0), nuclear_angles=(25.0, 30.0)
    )
    e_levels, _ = basis.flat()
    shifts = [complex(e, -0.001) for e in e_levels]
    return exact_resonance_states(
        _SeparableModel(), base, moved_el, moved_nu, shifts=shifts, k=10, window=_WINDOW
    )


# --- the end-to-end oracle ---------------------------------------------------


def test_every_exact_pole_pairs_to_its_bo_state_with_unit_overlap(separable_basis, separable_poles):
    """In the separable limit the exact state IS the BO product, so overlap = 1.

    This is the whole pipeline at once: the pole walk that builds the electronic
    curve, the phase alignment across `R`, the vibrational ladder, the c-product
    normalization, and the pairing. Any one of them wrong and the overlap falls
    away from 1 -- a lost phase alignment in particular drives it toward 0.
    """
    _, basis = separable_basis
    res = separable_poles
    assert res.energies.size >= 3, "the search should find the seeded levels"

    matched = 0
    for i in range(res.energies.size):
        pair = pair_by_overlap(res.energies[i], res.states[:, i], basis)
        if pair.verdict == "basis-limited":
            continue  # a pole above the three levels built here
        assert pair.verdict == "ok", f"{pair.pole_energy}: {pair.verdict} @ {pair.overlap}"
        assert pair.overlap == pytest.approx(1.0, abs=1e-5)
        # The separable limit has NO Born-Oppenheimer shift, so the only
        # residual is the pole walk's own: its window recentring and residual
        # tolerance leave ~1e-7 Ha, the same scale as the Gamma noise the
        # `local_complex_potential` comparison sees. Anything larger would be a
        # real disagreement between the two routes.
        assert abs(pair.shift_mev) < 1e-2, "the separable limit has no BO shift"
        matched += 1
    assert matched >= 3, f"only {matched} poles paired cleanly"


def test_the_pairing_is_a_bijection_on_the_separable_model(separable_basis, separable_poles):
    """Distinct poles must claim distinct levels -- no level assigned twice."""
    _, basis = separable_basis
    res = separable_poles
    claimed = [
        p.level
        for p in (
            pair_by_overlap(res.energies[i], res.states[:, i], basis)
            for i in range(res.energies.size)
        )
        if p.level is not None
    ]
    assert len(claimed) == len(set(claimed)), f"a level was claimed twice: {claimed}"


def test_overlap_agrees_with_the_energy_pairing_where_both_are_defined(
    separable_basis, separable_poles
):
    """`pair_one_to_one` is a cross-check, and here it must agree.

    With no BO error the two criteria cannot disagree; where they do on a real
    model, that disagreement is the finding.
    """
    _, basis = separable_basis
    res = separable_poles
    level_e, keys = basis.flat()
    by_energy = pair_one_to_one(res.energies.real, level_e)
    for pole_i, level_i in by_energy.items():
        pair = pair_by_overlap(res.energies[pole_i], res.states[:, pole_i], basis)
        if pair.level is not None:
            assert pair.level == keys[level_i]


# --- overlap -----------------------------------------------------------------


def test_overlap_of_a_state_with_itself_is_one(separable_basis):
    _, basis = separable_basis
    for _, s in basis.items():
        assert overlap(s.psi, s.psi) == pytest.approx(1.0, abs=1e-12)


def test_overlap_is_scale_invariant(separable_basis):
    """Both sides are normalized inside, so a rescaled input changes nothing."""
    _, basis = separable_basis
    a = basis[(0, 0)].psi
    b = basis[(0, 1)].psi
    assert overlap(a, b) == pytest.approx(overlap(3.7 * a, -0.2j * b), abs=1e-12)


def test_distinct_vibrational_levels_are_nearly_orthogonal(separable_basis):
    """Different `v` in one curve: orthogonal by construction, so overlap ~ 0."""
    _, basis = separable_basis
    assert overlap(basis[(0, 0)].psi, basis[(0, 1)].psi) < 1e-6


def test_overlap_of_a_self_orthogonal_vector_is_zero():
    """A c-product norm can genuinely vanish; returning 0 beats dividing by it."""
    v = np.array([1.0, 1.0j], dtype=np.complex128)  # c_product(v, v) == 0
    assert complex(c_product(v, v)) == pytest.approx(0.0, abs=1e-15)
    assert overlap(v, np.array([1.0, 0.0], dtype=np.complex128)) == 0.0


# --- verdicts ----------------------------------------------------------------


def _toy_basis(n: int = 6) -> BoBasis:
    """An orthonormal toy basis, so overlaps are exactly computable by hand."""
    states = {}
    for j in range(2):
        for v in range(2):
            psi = np.zeros(n, dtype=np.complex128)
            psi[2 * j + v] = 1.0
            states[(j, v)] = BoState(psi=psi, energy=-0.10 + 0.001 * (2 * j + v), curve=j, vib=v)
    return BoBasis(
        energies=np.array([[-0.100, -0.099], [-0.098, -0.097]], dtype=np.float64), states=states
    )


def test_a_clean_match_is_ok():
    basis = _toy_basis()
    psi = basis[(1, 0)].psi
    pair = pair_by_overlap(-0.098 + 0j, psi, basis)
    assert pair.verdict == "ok"
    assert pair.level == (1, 0)
    assert pair.overlap == pytest.approx(1.0)
    assert pair.shift_mev == pytest.approx(0.0, abs=1e-9)
    assert pair.is_resonance


def test_an_equal_blend_is_mixed_and_quotes_no_shift_worth_trusting():
    """Past a certain coupling neither BO label describes the state."""
    basis = _toy_basis()
    psi = (basis[(0, 0)].psi + basis[(1, 1)].psi) / np.sqrt(2.0)
    pair = pair_by_overlap(-0.099 + 0j, psi, basis)
    assert pair.verdict == "mixed"
    assert pair.second_overlap / pair.overlap == pytest.approx(1.0, abs=1e-9)


def test_a_state_outside_the_basis_is_basis_limited_without_thresholds():
    """Without the energy argument the honest verdict withholds judgement."""
    basis = _toy_basis()
    orphan = np.zeros(6, dtype=np.complex128)
    orphan[5] = 1.0
    pair = pair_by_overlap(-0.098 + 0j, orphan, basis)
    assert pair.verdict == "basis-limited"
    assert pair.level is None
    assert pair.is_resonance, "basis-limited withholds a label, it does not reject"
    assert np.isnan(pair.shift_ev)


def test_a_state_outside_a_PROVABLY_COMPLETE_basis_is_spurious():
    """With the closed-channel check satisfied, a vanishing overlap convicts.

    This is the H2+ finding in miniature: angle stability found the state, and
    only the admissibility argument turns "no partner here" into "not a
    resonance".
    """
    basis = _toy_basis()
    orphan = np.zeros(6, dtype=np.complex128)
    orphan[5] = 1.0
    # No threshold lies above this energy, so nothing is admissible and the
    # basis trivially covers it -- the strongest form of "the basis is not the
    # problem".
    thresholds = np.array([-0.20, -0.15])
    pair = pair_by_overlap(-0.098 + 0j, orphan, basis, thresholds)
    assert pair.verdict == "spurious"
    assert pair.level is None
    assert not pair.is_resonance


def test_a_poor_but_nonzero_match_is_weak_not_ok():
    """A moderate score against a wrong partner must not pass as a label.

    Two H2+ poles reached 0.21 and 0.11 against levels they are not, simply
    because nothing competed with them. Without the `weak` band those would have
    been reported as clean identifications with a quotable shift.
    """
    basis = _toy_basis()
    psi = 0.3 * basis[(0, 0)].psi
    psi[5] = np.sqrt(1.0 - 0.09)
    thresholds = np.array([-0.20])  # nothing above: the basis provably covers
    pair = pair_by_overlap(-0.100 + 0j, psi, basis, thresholds)
    assert pair.verdict == "weak"
    assert pair.level == (0, 0)
    assert 0.10 < pair.overlap < 0.50


def test_a_clear_match_far_away_in_energy_is_distant():
    """A big overlap across a big gap is mixing, not a shift worth quoting."""
    basis = _toy_basis()
    psi = basis[(0, 0)].psi
    pair = pair_by_overlap(-0.090 + 0j, psi, basis)  # 10 mHa = 272 meV away
    assert pair.verdict == "distant"
    assert pair.overlap == pytest.approx(1.0)


def test_thresholds_can_be_tightened_by_keyword():
    basis = _toy_basis()
    psi = basis[(0, 0)].psi
    off = -0.1005 + 0j  # 0.5 mHa = 13.6 meV below the level: inside the default
    assert pair_by_overlap(off, psi, basis).verdict == "ok"
    assert pair_by_overlap(off, psi, basis, max_shift_ev=0.010).verdict == "distant"


def test_pair_by_overlap_refuses_a_basis_without_states():
    basis = BoBasis(energies=np.array([[0.0]]), states={})
    with pytest.raises(ValueError, match="with_states=True"):
        pair_by_overlap(0.0 + 0j, np.zeros(2, dtype=np.complex128), basis)


def test_shift_mev_is_a_thousand_times_shift_ev():
    """The unit that was once got wrong by a factor of 27.2.

    An earlier revision computed `1000 * dE(Ha)` and called it meV -- that is
    milli-HARTREE. The conversion now goes through `qscat.units.HARTREE_TO_EV`,
    and this pins the two accessors to each other.
    """
    p = OverlapPair(0.0, (0, 0), 1.0, None, 0.0, 0.0123, "ok")
    assert p.shift_mev == pytest.approx(12.3)


# --- pair_one_to_one ---------------------------------------------------------


def test_bijection_resolves_a_crossing_that_nearest_neighbour_gets_wrong():
    """The measured H2+ geometry: two close levels, two farther-apart poles.

    Nearest-neighbour hands BOTH poles to the lower level and leaves the upper
    unclaimed. Requiring a bijection uses the constraint the physics supplies --
    distinct states occupy distinct levels -- and recovers the right answer.
    """
    levels = np.array([0.0, 20e-6])  # 20 uHa apart
    poles = np.array([-4e-6, 6e-6])  # both nearer the LOWER level
    nearest = [int(np.argmin(np.abs(levels - p))) for p in poles]
    assert nearest == [0, 0], "the failure mode this test exists to rule out"

    assignment = pair_one_to_one(poles, levels)
    assert assignment == {0: 0, 1: 1}


def test_pairs_beyond_max_distance_are_dropped_not_forced():
    levels = np.array([0.0, 1.0])
    poles = np.array([1e-6, 0.5])
    assert pair_one_to_one(poles, levels, max_distance=1e-3) == {0: 0}


def test_empty_inputs_pair_to_nothing():
    assert pair_one_to_one(np.array([]), np.array([1.0])) == {}
    assert pair_one_to_one(np.array([1.0]), np.array([])) == {}


# --- peaks -------------------------------------------------------------------


def test_peak_positions_finds_prominent_maxima_only():
    e = np.linspace(0.0, 1.0, 101)
    s = np.full_like(e, 0.01)
    s[30] = 10.0  # towering
    s[70] = 0.011  # a ripple, not a peak
    found = peak_positions(e, s)
    assert found.size == 1
    assert found[0] == pytest.approx(e[30])


def test_peak_positions_needs_a_long_enough_sweep():
    assert peak_positions(np.arange(4.0), np.arange(4.0)).size == 0
    with pytest.raises(ValueError, match="equal length"):
        peak_positions(np.arange(5.0), np.arange(4.0))


def test_peak_alignment_measures_in_widths_not_in_hartree():
    """The unit is the point: the same energy gap is near or far by width.

    A 2e-5 Ha displacement is a third of a width for a 6e-5 Ha resonance and ten
    widths for a 2e-6 Ha one. Reporting Hartree would make the two look alike.
    """
    peaks = np.array([0.0, 1e-3])
    marks = np.array([2e-5, 1e-3])
    wide = peak_alignment(marks, peaks, width=6e-5)
    narrow = peak_alignment(marks, peaks, width=2e-6)
    assert wide.median_widths < 1.0 < narrow.median_widths
    assert wide.within_one_width == 2
    assert narrow.within_one_width == 1
    assert wide.n_marks == 2 and wide.n_peaks == 2
    assert "2 marks vs 2 peaks" in str(wide)


def test_peak_alignment_rejects_degenerate_inputs():
    with pytest.raises(ValueError, match="width must be positive"):
        peak_alignment([0.0], [0.0], width=0.0)
    with pytest.raises(ValueError, match="needs both"):
        peak_alignment([], [0.0], width=1e-5)


def test_basis_complete_lets_a_neutral_caller_assert_coverage():
    """The neutral escape hatch: no Rydberg series to count, so no thresholds.

    An anion resonance curve carries ONE electronic state, so its basis is
    complete once every vibrational level is built -- but `admissible_levels`
    has no way to know that, and would report every orphan `basis-limited`.
    Asserting completeness is a claim the caller has to be able to defend,
    because it is what turns `basis-limited` into `spurious`.
    """
    basis = _toy_basis()
    orphan = np.zeros(6, dtype=np.complex128)
    orphan[5] = 1.0
    assert pair_by_overlap(-0.098 + 0j, orphan, basis).verdict == "basis-limited"
    assert (
        pair_by_overlap(-0.098 + 0j, orphan, basis, basis_complete=True).verdict == "spurious"
    )
    # And it overrides thresholds rather than being combined with them.
    covering = np.array([-0.20])
    assert (
        pair_by_overlap(
            -0.098 + 0j, orphan, basis, covering, basis_complete=False
        ).verdict
        == "basis-limited"
    )
