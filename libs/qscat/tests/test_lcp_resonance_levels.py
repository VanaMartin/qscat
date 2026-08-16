"""Analytic oracles for `qscat.core.lcp.lcp_resonance_levels`.

Two exact benchmarks, no convergence hand-waving:

1. Gamma = 0 with a bare Morse curve -> the analytic Morse spectrum,
   E_n = -D (1 - alpha (n + 1/2) / sqrt(2 mu D))^2, with Im E ~ 0.
2. Gamma = Gamma_0 CONSTANT -> a constant imaginary term commutes with
   everything, so the spectrum must shift rigidly by exactly -i Gamma_0/2.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import segmented_grid
from qscat.core.lcp import lcp_resonance_levels

# A deep, heavy Morse well (N2-like) so the analytic infinite-domain formula
# holds to high accuracy: V(0) ~ 64 Ha is an effectively infinite inner wall.
MU, D0, ALPHA0, RE = 12766.36, 0.75102, 1.1535, 2.01943

# Lowest five analytic levels; the sixth (~ -0.6838) lies outside WINDOW.
WINDOW = (-0.76, -0.69, -1e-6, 1e-6)


def morse_levels(n: np.ndarray) -> np.ndarray:
    x = ALPHA0 * (n + 0.5) / np.sqrt(2.0 * MU * D0)
    return -D0 * (1.0 - x) ** 2


def morse(R: np.ndarray) -> np.ndarray:
    z = np.asarray(R, dtype=np.complex128)
    return D0 * (np.exp(-2 * ALPHA0 * (z - RE)) - 2 * np.exp(-ALPHA0 * (z - RE)))


def grid_pair(angle_a: float = 35.0, angle_b: float = 25.0, n_real: int = 30):
    """Two nuclear grids sharing every real node, differing only in tail angle."""
    real, complex_ = [(n_real, 6.0)], [(6, 16.0)]
    return (
        segmented_grid(real, complex_, angle_deg=angle_a, quadrature=10),
        segmented_grid(real, complex_, angle_deg=angle_b, quadrature=10),
    )


def test_zero_width_reproduces_the_analytic_morse_spectrum():
    ga, gb = grid_pair()
    Vd_a, Vd_b = morse(ga.points), morse(gb.points)
    Gamma = np.zeros(ga.n, dtype=np.float64)

    out = lcp_resonance_levels(ga, gb, MU, Vd_a, Vd_b, Gamma, window=WINDOW)

    assert out.energies.shape == (5,)
    np.testing.assert_allclose(out.energies.real, morse_levels(np.arange(5)), rtol=1e-5)
    assert np.all(np.abs(out.energies.imag) < 1e-8)
    assert np.all(out.widths < 1e-8)
    # Bound states live entirely in the real region.
    assert np.all(out.real_weight > 0.999)


def test_constant_width_shifts_the_spectrum_rigidly():
    ga, gb = grid_pair()
    Vd_a, Vd_b = morse(ga.points), morse(gb.points)
    g0 = 0.01
    zero = lcp_resonance_levels(ga, gb, MU, Vd_a, Vd_b, np.zeros(ga.n), window=WINDOW)
    shifted = lcp_resonance_levels(
        ga,
        gb,
        MU,
        Vd_a,
        Vd_b,
        np.full(ga.n, g0),
        window=(-0.76, -0.69, -0.5 * g0 - 1e-6, -0.5 * g0 + 1e-6),
    )
    # H(Gamma_0) = H(0) - i (Gamma_0/2) I exactly -- round-off, not tolerance.
    np.testing.assert_allclose(shifted.energies, zero.energies - 0.5j * g0, atol=1e-10)
    np.testing.assert_allclose(shifted.widths, g0, atol=1e-10)


def test_states_are_c_product_normalized():
    from qscat.linalg import c_product

    ga, gb = grid_pair()
    out = lcp_resonance_levels(
        ga, gb, MU, morse(ga.points), morse(gb.points), np.zeros(ga.n), window=WINDOW
    )
    for state in out.states:
        assert abs(c_product(state, state) - 1.0) < 1e-10

    # The bound levels above have negligible ECS-tail amplitude, so the raw
    # LAPACK v^dagger v = 1 eigenvector already satisfies |sum c^2 - 1| ~
    # 1e-15 there too -- c_product-normalizing is a no-op and this alone
    # does not discriminate the bilinear norm from the conjugated one (nor
    # would it catch a `vdot`-based implementation). A level with real
    # amplitude in the ECS-rotated tail does discriminate: near/above the
    # anion dissociation limit (Re E -> 0) a state's eigenvector has
    # sizeable weight on the complex-rotated nodes, where the
    # (non-conjugated) bilinear c-product and the conjugated v^dagger v
    # genuinely disagree.
    near_threshold_window = (-0.0165, -0.005, -1e-4, 1e-4)
    out_tail = lcp_resonance_levels(
        ga,
        gb,
        MU,
        morse(ga.points),
        morse(gb.points),
        np.zeros(ga.n),
        window=near_threshold_window,
        rel_tol=1e-3,
        golden_rule=False,
    )
    assert out_tail.states.shape[0] == 1
    assert out_tail.real_weight[0] < 0.96  # meaningful (~5%) tail weight
    state = out_tail.states[0]
    assert abs(c_product(state, state) - 1.0) < 1e-10
    assert abs(np.vdot(state, state) - 1.0) > 1e-6


def test_n_levels_truncates_to_the_lowest():
    ga, gb = grid_pair()
    out = lcp_resonance_levels(
        ga,
        gb,
        MU,
        morse(ga.points),
        morse(gb.points),
        np.zeros(ga.n),
        window=WINDOW,
        n_levels=2,
    )
    assert out.energies.shape == (2,)
    np.testing.assert_allclose(out.energies.real, morse_levels(np.arange(2)), rtol=1e-5)


def test_mismatched_real_regions_raise():
    ga, _ = grid_pair()
    _, gb_wrong = grid_pair(n_real=20)  # different real discretization
    with pytest.raises(ValueError, match="real nodes"):
        lcp_resonance_levels(
            ga,
            gb_wrong,
            MU,
            morse(ga.points),
            morse(gb_wrong.points),
            np.zeros(ga.n),
            window=WINDOW,
        )


def test_non_positive_mass_raises():
    ga, gb = grid_pair()
    with pytest.raises(ValueError, match="mu must be positive"):
        lcp_resonance_levels(
            ga,
            gb,
            0.0,
            morse(ga.points),
            morse(gb.points),
            np.zeros(ga.n),
            window=WINDOW,
        )


def test_golden_rule_returns_nan_when_comparator_window_is_empty():
    """A diagnostic-only failure must not take down the primary result.

    Near/above the anion dissociation limit, a level's Gamma=0 comparator
    already has a nonzero Im E purely from V_d's own complex ECS-tail
    continuation (no explicit width needed) -- large enough that the
    golden-rule branch's tight `[-atol, atol]` Im band catches NO angle-
    stable state there, even though the primary (wider-window) solve above
    correctly keeps the level as physical. `match_angle_stable` raises
    `ValueError` in that situation (window catches nothing); the primary
    energies must still come back, with `golden_rule` reporting `nan`.
    """
    ga, gb = grid_pair()
    Vd_a, Vd_b = morse(ga.points), morse(gb.points)
    Gamma = np.zeros(ga.n)
    window = (-0.0165, -0.005, -1e-4, 1e-4)

    out = lcp_resonance_levels(ga, gb, MU, Vd_a, Vd_b, Gamma, window=window, rel_tol=1e-3)

    assert out.energies.shape == (1,)
    assert np.all(np.isnan(out.golden_rule))


def test_golden_rule_nans_a_level_with_no_plausible_comparator():
    """The distance guard, not just the whole-call try/except.

    A wider window mixes well-separated bound levels (their comparator IS
    findable -- `golden_rule` should closely match `energies`, same as
    `test_zero_width_reproduces_the_analytic_morse_spectrum`'s Gamma=0
    case) with the same near-dissociation-limit levels as the test above
    (comparator excluded by the tight Im band). `E0` is non-empty here (the
    bound levels' comparators exist), so unconditional nearest-neighbor
    `argmin` pairing would silently glue the excluded levels to a distant
    bound-level comparator instead of reporting `nan` -- this is the
    per-level guard, exercised independently of the all-nan crash path.
    """
    ga, gb = grid_pair()
    Vd_a, Vd_b = morse(ga.points), morse(gb.points)
    Gamma = np.zeros(ga.n)
    window = (-0.03, -0.005, -1e-4, 1e-4)

    with pytest.warns(UserWarning, match="dropped 1 level"):
        out = lcp_resonance_levels(ga, gb, MU, Vd_a, Vd_b, Gamma, window=window, rel_tol=1e-3)

    assert out.energies.shape == (6,)
    # The four well-separated bound levels: a real, findable comparator.
    np.testing.assert_allclose(out.golden_rule[:4], out.energies[:4], atol=1e-8)
    # The two near-dissociation-limit levels: no plausible comparator.
    assert np.all(np.isnan(out.golden_rule[4:]))


def _lcp_grids():
    """F2's LCP decks: fine nuclear grid at two tail angles + two electronic angles."""
    from qscat.core.grids import electronic_grid

    real = [(9, 1.8), (1, 2.0), (5, 2.5), (4, 2.596908), (4, 2.7), (4, 10.7)]
    cx = [(15, 30.0)]
    nuc_a = segmented_grid(real, cx, angle_deg=25.0, quadrature=12)
    nuc_b = segmented_grid(real, cx, angle_deg=15.0, quadrature=12)
    ea = electronic_grid(r_max=16.0, order=8, n_complex=6, angle_deg=35.0)
    eb = electronic_grid(r_max=16.0, order=8, n_complex=6, angle_deg=44.0)
    return nuc_a, nuc_b, ea, eb


@pytest.mark.slow
def test_f2_levels_are_bound_and_narrow():
    from qscat.core.lcp import resonance_levels
    from qscat.model import F2

    nuc_a, nuc_b, ea, eb = _lcp_grids()
    # Same reason as test_gamma_support_condition_holds_for_f2: the library
    # default re_half_width=im_half_width=0.05 under-resolves the walk
    # through the R~2.6 crossing, leaking a spurious ~2e-5 Gamma INSIDE the
    # well (where the vibrational wavefunctions have real support) -- that
    # would contaminate the widths/golden-rule assertions below. 0.01
    # resolves the walk cleanly (see the sibling test's docstring).
    out = resonance_levels(
        F2, nuc_a, nuc_b, ea, eb, n_levels=6, re_half_width=0.01, im_half_width=0.01
    )

    assert out.energies.size > 0
    assert np.all(np.diff(out.energies.real) > 0)  # ascending, non-degenerate
    assert np.all(out.widths >= 0.0)  # clamped
    assert np.all(out.residuals < 1e-3)  # genuinely angle-stable
    assert np.all(out.real_weight > 0.5)  # localized, not continuum
    # The comparator must agree with the complex result on the narrow levels:
    # for Gamma_v << level spacing the shift is first-order.
    narrow = out.widths < 1e-4
    assert narrow.any()
    np.testing.assert_allclose(out.energies.real[narrow], out.golden_rule.real[narrow], atol=1e-4)


@pytest.mark.slow
def test_gamma_support_condition_holds_for_f2():
    """Vana 2017 Sec. 1.5: Im V_res is nonzero ONLY where v0(R) < E_res(R).

    With the DEFAULT `re_half_width=im_half_width=0.05`, the walk's
    window-recentering step is too coarse to track the pole precisely
    through the crossing: the worst in-bound-region Gamma is 2.17e-5,
    localized to ~5 points at R~2.5976-2.608 (just above the 2.596908
    segment boundary), where Im(E_pole) turns on steeply over ~0.01 bohr in
    R. That is a numerical-resolution artifact of the walk's step size, not
    the physical answer -- the support condition says Gamma is EXACTLY zero
    there. Tightening to 0.01 resolves the walk through the crossing and
    drops the worst value to 3e-14, confirming the artifact (not the
    physics) was the culprit. The half-widths below are there to resolve
    the walk, not to make this assertion pass.
    """
    from qscat.core.lcp import local_complex_potential
    from qscat.model import F2

    nuc_a, _nuc_b, ea, eb = _lcp_grids()
    Vd, Gamma = local_complex_potential(F2, nuc_a, ea, eb, re_half_width=0.01, im_half_width=0.01)
    real = nuc_a.points.imag == 0.0
    R = nuc_a.points[real].real
    bound_region = Vd[real].real < F2.v0(R).real  # anion below neutral: no autodetachment
    assert np.all(Gamma[real][bound_region] < 1e-6)


def test_grid_angle_bound_is_enforced():
    from qscat.core.grids import electronic_grid
    from qscat.core.lcp import resonance_levels
    from qscat.model import H2P

    real, cx = [(5, 1.0), (90, 14.05)], [(3, 30.0)]
    too_steep = segmented_grid(real, cx, angle_deg=30.0, quadrature=8)  # > 22.5
    ok = segmented_grid(real, cx, angle_deg=20.0, quadrature=8)
    e = electronic_grid(r_max=16.0, order=8, n_complex=6, angle_deg=35.0)
    with pytest.raises(ValueError, match="max_nuclear_ecs_angle_deg"):
        resonance_levels(H2P, too_steep, ok, e, e)


def test_levels_are_independent_of_the_ecs_angle_pair():
    """A physical level does not move when the rotation angles change."""
    ga1, gb1 = grid_pair(angle_a=35.0, angle_b=25.0)
    ga2, gb2 = grid_pair(angle_a=45.0, angle_b=20.0)
    kw = dict(window=WINDOW, n_levels=5)
    one = lcp_resonance_levels(
        ga1, gb1, MU, morse(ga1.points), morse(gb1.points), np.zeros(ga1.n), **kw
    )
    two = lcp_resonance_levels(
        ga2, gb2, MU, morse(ga2.points), morse(gb2.points), np.zeros(ga2.n), **kw
    )
    np.testing.assert_allclose(one.energies, two.energies, rtol=1e-6, atol=1e-9)


def test_levels_converge_under_h_refinement():
    """Refining the real region must not move the levels.

    NOTE on the fixture: the brief's original (n_real=20, n_real=40) pair
    fails this assertion -- not from a defect in `lcp_resonance_levels`, but
    because a single-segment FEM-DVR converges EXPONENTIALLY in n_real for
    this smooth Morse well (verified directly against the analytic
    spectrum: n_real=20 is already accurate to ~7.6e-6, n_real=40 to
    ~1.3e-10), so two levels differ by ~1.1e-5 relative -- past rtol=1e-6.
    (n_real=30, n_real=60) keeps the same coarse/fine h-refinement intent
    (30 is this file's baseline resolution used throughout) while landing
    in a regime where rtol=1e-6 is actually meaningful (max relative
    difference ~1.0e-7 here, verified against the analytic oracle too) --
    see the Task 5 report for the full convergence table.

    This test does NOT check `residuals` shrinking under refinement:
    `residuals` measures ECS-TAIL angle stability, and `nuclear_grid_a`/
    `nuclear_grid_b` share every real node by construction, so real-region
    discretization error is common to both spectra and cancels out of the
    residual -- it stays pinned near machine precision (~1e-15) regardless
    of `n_real` (see `ResonanceLevels.residuals`'s docstring). `energies`
    actually moving under refinement, checked above, is the real
    convergence signal.
    """
    coarse_a, coarse_b = grid_pair(n_real=30)
    fine_a, fine_b = grid_pair(n_real=60)
    kw = dict(window=WINDOW, n_levels=5)
    coarse = lcp_resonance_levels(
        coarse_a,
        coarse_b,
        MU,
        morse(coarse_a.points),
        morse(coarse_b.points),
        np.zeros(coarse_a.n),
        **kw,
    )
    fine = lcp_resonance_levels(
        fine_a, fine_b, MU, morse(fine_a.points), morse(fine_b.points), np.zeros(fine_a.n), **kw
    )
    np.testing.assert_allclose(coarse.energies.real, fine.energies.real, rtol=1e-6)


def test_golden_rule_matches_the_complex_result_for_a_weak_constant_width():
    """First-order perturbation theory is EXACT for a constant Gamma, so the
    comparator must reproduce the complex levels to round-off."""
    ga, gb = grid_pair()
    g0 = 1e-5
    out = lcp_resonance_levels(
        ga,
        gb,
        MU,
        morse(ga.points),
        morse(gb.points),
        np.full(ga.n, g0),
        window=(-0.76, -0.69, -0.5 * g0 - 1e-6, -0.5 * g0 + 1e-6),
    )
    np.testing.assert_allclose(out.golden_rule, out.energies, atol=1e-9)


def test_golden_rule_can_be_switched_off():
    ga, gb = grid_pair()
    out = lcp_resonance_levels(
        ga,
        gb,
        MU,
        morse(ga.points),
        morse(gb.points),
        np.zeros(ga.n),
        window=WINDOW,
        golden_rule=False,
    )
    assert np.all(np.isnan(out.golden_rule.real))


def test_energies_diverge_from_golden_rule_for_a_broad_r_dependent_level():
    """A large, R-DEPENDENT Gamma is the one case that can actually break
    first-order perturbation theory -- unlike constant Gamma (see
    `test_golden_rule_matches_the_complex_result_for_a_weak_constant_width`),
    which commutes with H(0) and shifts the whole spectrum rigidly by
    construction (the comparator is then exact to round-off, not a real
    test of divergence). `golden_rule` freezes the eigenstate at its
    Gamma=0 shape and adds only the linear-in-Gamma width
    `<chi_v|Gamma|chi_v>`; it cannot see the wavefunction REARRANGING away
    from a large, spatially localized Gamma, which is a second-order (and
    higher) effect that only the full complex diagonalization captures.

    Gamma here is a Gaussian bump straddling the outer part of the well
    (centered at RE + 0.4, far from the compact ground state but well
    inside the reach of the more extended upper levels), so within ONE run
    the ground level stays perturbative while the topmost level in the
    window goes genuinely non-perturbative.
    """
    ga, gb = grid_pair()
    Vd_a, Vd_b = morse(ga.points), morse(gb.points)
    R = ga.points.real
    center, sigma, g0 = RE + 0.4, 0.1, 1.5
    Gamma = g0 * np.exp(-((R - center) ** 2) / sigma**2)

    out = lcp_resonance_levels(
        ga, gb, MU, Vd_a, Vd_b, Gamma, window=(-0.76, -0.69, -0.3, 1e-6), n_levels=5
    )

    # Pin the count exactly (not >=): the assertions below index [0] and
    # [-1] to mean "narrowest" and "broadest", which only holds for this
    # specific 4-level result -- a 5th in-window level would silently
    # retarget them onto the wrong level.
    assert out.energies.size == 4
    assert not np.any(np.isnan(out.golden_rule))  # every level found a comparator
    # Angle-stable throughout -- this is a real level, not continuum leakage.
    assert np.all(out.residuals < 1e-10)

    diff = np.abs(out.energies - out.golden_rule)
    # Ground level: compact, barely samples the bump -- small, near-
    # perturbative deviation (nowhere near the round-off of the constant-
    # Gamma case, but two-plus orders of magnitude below the broad level).
    assert diff[0] < 5e-5
    # Topmost level: the deviation itself EXCEEDS the golden-rule width
    # prediction -- a qualitative breakdown, not just a quantitative one.
    assert diff[-1] > out.widths[-1]
    assert diff[-1] > 100.0 * diff[0]


def test_two_angle_selection_isolates_the_resonance_from_the_rotated_continuum():
    """The angle-stability selection is only trivially exercised by the
    other oracles above: their levels are bound-state limits with
    `real_weight` > 0.999, so any reasonable criterion would keep them.
    Pick a window near the anion dissociation limit instead, where the
    physical level has genuine amplitude on the ECS-rotated tail
    (`real_weight` meaningfully below 1, as in
    `test_lcp.py::test_states_are_c_product_normalized`'s near-threshold
    state) and sits among several raw discretized-continuum eigenvalues in
    the same window. `match_angle_stable` must keep exactly the one
    angle-invariant level and reject the rest -- which DO move with the
    ECS angle, so their count in a fixed window is not itself invariant.
    """
    from qscat.dvr import eigen, kinetic

    window = (-0.0165, -0.003, -5e-3, 1e-4)
    angle_pairs = [(35.0, 25.0), (45.0, 20.0), (40.0, 30.0)]
    selected_counts = []
    raw_counts = []
    for angle_a, angle_b in angle_pairs:
        ga, gb = grid_pair(angle_a=angle_a, angle_b=angle_b)
        Vd_a, Vd_b = morse(ga.points), morse(gb.points)
        # rel_tol=1e-3 (10x the 1e-4 default): this near-threshold level's
        # residual runs 2e-6..1e-5, which the default `max(rel_tol*|E|,
        # atol)` threshold (~1e-6 at |E|~0.01) would reject outright -- not
        # tolerance-fudging, this level genuinely needs the looser cut to
        # be selected at all.
        out = lcp_resonance_levels(
            ga,
            gb,
            MU,
            Vd_a,
            Vd_b,
            np.zeros(ga.n),
            window=window,
            rel_tol=1e-3,
            golden_rule=False,
        )
        assert out.energies.shape == (1,)  # exactly the one physical level
        assert out.residuals[0] < 1e-4  # tightly angle-stable
        assert out.real_weight[0] < 0.96  # real amplitude in the rotated tail
        selected_counts.append(out.energies.size)

        E_a, _ = eigen(kinetic(ga, MU) + np.diag(Vd_a))
        re_lo, re_hi, im_lo, im_hi = window
        raw_counts.append(
            int(
                np.count_nonzero(
                    (E_a.real >= re_lo)
                    & (E_a.real <= re_hi)
                    & (E_a.imag >= im_lo)
                    & (E_a.imag <= im_hi)
                )
            )
        )

    # The selected level count is stable across angle pairs...
    assert len(set(selected_counts)) == 1
    # ...while the raw eigenvalue count in the same window is NOT: the
    # discretized continuum rotates with theta, so how much of it lands in
    # a fixed box changes with the angle pair. This is exactly what the
    # selection criterion is discarding.
    assert len(set(raw_counts)) > 1
    assert all(raw > sel for raw, sel in zip(raw_counts, selected_counts, strict=True))
