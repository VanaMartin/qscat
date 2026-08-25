"""GATE: NRM vibrational excitation against the exact 2-D oracle.

Every band below is RECORDED from the first converged run, not a preset
target -- the design spec requires the success criterion here to be a
measurement. The run that produced each one is written into the test that
asserts it; their physical reading is in
`docs/physics/nonlocal-resonance-model.md`.

EVERY NUMBER IN THIS FILE IS MEASURED ON A GATE ENERGY GRID, not on every
energy in the plotted window. For N2 that grid is the one the gate still runs
(11 energies). For F2 the recorded numbers come from the ORIGINAL SIX-ENERGY
sweep; the gate now runs only the two binding anchors of it (0.02, 0.04 -- see
`_ENERGIES` for why those two hold the bands). Every F2 figure below is
labelled with the grid it was recorded on. A denser sweep of the same window
lands on
N2's ~0.01 Ha boomerang peaks instead of straddling them, and the
approximations are worst exactly there, so the dense bands are wider. The
comments on `_BANDS`, `_BAND_N2_01` and `_N2_ABS_DEV_CEILING_01` carry the
101-energy figure-grid values alongside each recorded one, and explain why the
gate keeps its 11 energies rather than re-recording. Do not read a failure on
a densified sweep as a regression.

THE HEADLINE, `sigma_NRM(B)+bg / sigma_exact` (`n_states=100`, 0->0 and 0->1
together):

    N2   [0.997062, 1.000647]     11 energies, 0.06-0.16 Ha (the gate's grid)
    F2   [0.996225, 1.006923]      6 energies, 0.02-0.09 Ha (recorded sweep;
                                   the gate runs its two binding anchors,
                                   0.02 and 0.04, which define this band)

i.e. the nonlocal model with the R-independent discrete state and the Eq. (37)
background terms reproduces the exact 2-D solver to better than 0.7% on both
molecules, in both the elastic and the first inelastic channel. That is the
result the spec exists for, and unlike the DA comparison it is not a
one-molecule claim.

The other three routes over the same runs, for contrast:

              N2 (0.06-0.16 Ha)        F2 (0.02-0.09 Ha)
    A + bg    [0.85398, 1.05868]       [0.56528, 1.14013]
    B, no bg  [0.51715, 1.33786]       [0.04227, 0.90109]
    LCP       [0.10629, 4.56775]       [0.000177, 0.35574]

READ THE LCP ROW AS ABSOLUTE DEVIATION, NOT AS A RATIO. PRA 77 draws these on
LINEAR axes (p. 012710-8, Fig. 4), so a ratio far from 1 where sigma is small
is invisible in print and is NOT a disagreement with the paper. On N2 the two
channels behave completely differently once measured that way:

    N2 0->1   max |sigma_LCP - sigma_exact| = 0.53 bohr^2  (axis runs to 14)
    N2 0->0   max |sigma_LCP - sigma_exact| = 8.71 bohr^2  (peak sigma 35)

The 0->1 row is ~4% of its panel's axis -- the LCP is within 5.2% of exact
wherever sigma exceeds half its peak, and the wide ratios come only from the
wings. The 0->0 row is a genuine, visible failure at LARGE sigma, and it is
this paper's own missing-background claim. See
`test_n2_0to1_agrees_on_the_scale_fig_8_asserts`.

On the 101-energy figure grid those two absolute deviations are 0.71 and
11.65 bohr^2 -- the factor of 16 survives -- and the "within 5.2%" figure
becomes 9.9% (23 energies clear the half-peak cut there), against the elastic
channel's 61% under the identical cut. Quote the dense numbers when describing
the printed figure; see docs/physics/nonlocal-resonance-model.md Sec. 8.6.

VIBRATIONAL EXCITATION IS THE CHANNEL WITH PUBLISHED CURVES BEHIND IT. From
PRA 77's own panel inventory (`reference/literature/houfek-2008-pra77-012710.md`):
N2 0->0 appears in Fig. 4 (physical `phi_d`, choice A) and Fig. 8
(R-independent `phi_d`, choice B); F2 0->1 in Figs. 6 and 8; **N2 0->1 only in
Fig. 4** -- Fig. 8's N2 panels are 0->0 and 0->8; and **F2 0->0 is plotted
nowhere**. All four (molecule, v') pairs are gated here against OUR exact 2-D
solver, which exists for every transition -- what F2 0->0 lacks is external
corroboration, not an oracle.

ENERGY WINDOWS ARE THE PAPER'S: N2 VE is plotted over 0.05-0.17 Ha (Fig. 4
top) and F2 VE 0->1 over 0-0.10 Ha (Fig. 6 top); `_ENERGIES` stays inside
both. They are per molecule because the published windows are.

KNOWN, QUANTIFIED WEAKNESS IN THE N2 INGREDIENTS. `nrm_ingredients` warns on
N2's deck (`min_overlap = 0.0148 < 0.5`, both discrete-state choices) that
`_sign_align` paired different P-space states between adjacent nuclear nodes.
It is real but demonstrably not load-bearing here. Measured: the mispairing is
**three R-steps involving five states, swapped in PAIRS** -- (77, 78) at
R = 1.7608 -> 1.7374, (62, 63) at 1.3695 -> 1.2913 and (63, 64) at
0.9799 -> 0.9347 bohr, identical for both discrete-state choices. Three
independent reasons it cannot move the numbers gated below:

1. Those states are essentially DECOUPLED. Their largest `|V_dn|` anywhere on
   the walk is 1.04e-4 (choice A) / 7.08e-4 (choice B), against deck maxima of
   0.341 / 0.683 -- a relative coupling of 3e-4 / 1e-3. Eq. (60) is bilinear
   in `V_dn`, so they enter `F` at the ~1e-7 / ~1e-6 relative level, three to
   four orders below the agreement asserted below.
2. Each swapped pair is NEAR-DEGENERATE -- `E_n` = 212.6833/212.7180 and
   41.2234/41.2333 and 41.5195/41.5493 hartree. Relabelling the two members of
   a degenerate pair is very nearly a no-op in `F(E)`, which sums over the
   pair; only the split between them can matter, and it is 1.6e-4 relative.
3. Those `E_n` sit at 41-213 hartree against collision energies of 0.06-0.16,
   i.e. 260x to 3500x, so `F(E)`'s energy denominators suppress them further.

`phi_d` itself is continuous across N2's whole walk (minimum adjacent overlap
0.9895 for choice A, exactly 1.0 for the R-independent choice B) -- this is NOT
NO's choice-A pathology, where `phi_d` was genuinely discontinuous at a branch
switch.
"""

from __future__ import annotations

import functools

import numpy as np
import numpy.typing as npt
import pytest

from validation.diatomic.ve_nrm import VeComparison, compare

MOLECULES = ["N2", "F2"]

VPRIMES = [0, 1]

# MEASURED state-sum truncation (Eq. 60), one value for both molecules.
# Laddered over n_states = 10/25/40/55/70/85/100/(120)/all with the
# ingredients built once and reused, at two energies per molecule and for
# BOTH discrete-state choices and BOTH background settings. The number quoted
# per rung is the WORST |sigma(n) - sigma(n_prev)| / sigma(n) over the two
# energies and both channels, for the with-background curve:
#   N2 (E = 0.08/0.14 Ha; 106 states available, elec.n - 1)
#     A: 2.1e-2 at n=40, 1.0e-4 at 55, 3.0e-12 at 85, 2.6e-15 at 100,
#        bit-identical to the untruncated sum.
#     B: converges LATER and non-monotonically -- 1.9 (i.e. 190%) at n=40,
#        1.1e-1 at 55, 1.2e-5 at 70, 1.8e-11 at 85, 5.3e-15 at 100.
#   F2 (E = 0.03/0.08 Ha; 131 states available)
#     A: 2.0e-1 at n=40, 6.2e-4 at 55, 5.5e-6 at 70, 1.0e-10 at 100.
#     B: 1.1e-1 at n=40, 7.2e-2 at 55, 2.5e-2 at 70, 1.4e-3 at 85,
#        5.3e-7 at 100, 1.2e-14 at 120.
#
# The two molecules' shapes differ (F2/B is still moving at 1e-3 where N2/B
# has already reached 1e-11), which is why both were laddered rather than
# inheriting one number -- the same lesson spec 1 learned from NO. 100 is the
# smallest round value at which every one of the eight combinations sits
# within 1e-6 of its own untruncated sum, and it fits inside N2's 106-state
# ceiling.
_N_STATES = 100

# Energy grids, inside the paper's own plotted windows (module docstring).
#
# N2 keeps the full 11-point 0.06-0.16 Ha sweep. It is not the expensive one
# (95 s against F2's 746 s), and `test_n2_0to1_agrees_on_the_scale_fig_8_asserts`
# genuinely needs the CURVE: its claim is about where sigma sits relative to its
# own peak, which two points cannot locate.
#
# F2 is GATED AT TWO ANCHORS, reduced from the six that were originally run
# (0.02, 0.03, 0.04, 0.05, 0.07, 0.09 -- all six recorded in the docstrings
# below and in docs/physics/nonlocal-resonance-model.md). Nothing here needs
# F2's curve: sigma_VE varies smoothly and monotonically over the window, with
# no boomerang structure, and every F2 assertion in this module is either a
# per-energy band or a worst-over-energies extreme.
#
# The two are chosen so that every BINDING extreme survives -- the gate is as
# tight at two points as it was at six, not merely cheaper:
#
#   0.02  choice B's band MAXIMUM      (0->1 ratio 1.00692, the 0.006923 that
#                                       sets _CHOICE_A_ERROR_FLOOR's B side)
#         choice A's WORST error       (0->1 ratio 0.56528, i.e. 0.43472 --
#                                       the floor `test_choice_a_is_worse_
#                                       than_choice_b` pins)
#   0.04  choice B's band MINIMUM      (0->1 ratio 0.99623)
#
# 0.03/0.05/0.07/0.09 each sat strictly inside the bands those two define, so
# dropping them removes cost and no constraint. Re-measuring the full sweep is
# a deliberate act (widen this array), not something CI does on every run.
#
# MEASURED, both grids run back to back on the same machine with nothing else
# on it: six energies 1818 s, two energies 705 s -- 61% off, a factor of 2.58.
# That decomposes as ~278 s per energy against ~149 s of fixed setup (grid,
# LCP pole walk, ingredients), i.e. this really is per-energy work and the
# saving scales with the count. Note the 746 s in `_comparison`'s docstring is
# an older figure that did NOT reproduce here; treat 1818 s as the current
# six-energy cost on this hardware.
_ENERGIES: dict[str, npt.NDArray[np.float64]] = {
    "N2": np.array([0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16]),
    "F2": np.array([0.02, 0.04]),
}

# RECORDED `sigma_NRM(B)+bg / sigma_exact` bands (0->0 and 0->1 pooled).
# Measured N2 [0.997062, 1.000647] (max deviation from unity 0.294%, at
# E=0.09, 0->1) and F2 [0.996225, 1.006923] (max deviation 0.692%, at E=0.02,
# 0->1). Each band keeps roughly 50-70% headroom on that largest deviation,
# matching the convention `test_nrm.py` set for the DA gate. They are NOT
# equal because the two decks are not: F2's band was recorded on the
# six-energy sweep (the gate now runs its two binding anchors) against a
# larger exact-2D grid, and its own deviations run about twice N2's.
#
# THESE ARE ANCHOR-GRID BANDS, NOT WINDOW-WIDE BOUNDS. On the 101-energy grid
# `validation/diatomic/ve_nrm_figure.py` renders (0.06-0.16 Ha, step 0.001),
# N2's choice-B ratio spans [0.994539, 1.000647] -- below this band's 0.995
# floor. That is a grid-SAMPLING difference, not a regression: an 11-point
# sweep straddles N2's ~0.01 Ha boomerang peaks rather than landing on them,
# and the approximations are worst exactly at the peaks. The band is kept at
# the 11 energies the gate actually runs, where it is measured and sharp;
# re-recording it on 101 energies would put a ~14 min sweep inside a test and
# slacken the band on its own points for no gain in detection. If you densify
# the sweep, widen the band with it -- do not read the failure as a defect.
# See docs/physics/nonlocal-resonance-model.md Sec. 8.9.
_BANDS: dict[str, tuple[float, float]] = {
    "N2": (0.995, 1.005),
    "F2": (0.990, 1.012),
}

# RECORDED band for N2 0->1 across the NRM routes that INCLUDE the background
# (choices A and B). Measured A [0.853977, 1.058678], B [0.997062, 1.000647];
# the band keeps ~35% headroom on A's largest deviation (14.6%, at E=0.09).
# The LCP's RATIO leaves this band in the wings -- which is a statement about
# small sigma, not a disagreement with Fig. 8; see
# `test_n2_0to1_agrees_on_the_scale_fig_8_asserts`.
#
# ANCHOR-GRID BAND, same caveat as `_BANDS` above: on the 101-energy figure
# grid choice A's 0->1 ratio spans [0.828515, 1.132600] and so runs past this
# band's 1.10 ceiling. Kept at the gate's 11 energies for the same reasons.
_BAND_N2_01: tuple[float, float] = (0.80, 1.10)

# RECORDED N2 0->1 LCP ratios over the same 11 energies: [0.379276, 1.201764],
# worst at E=0.06 where sigma_exact is 0.287 bohr^2 (2.8% of the 10.23 peak
# ON THIS GRID; the dense grid resolves the peak at 11.48, making it 2.5%).
_N2_01_LCP_WORST_RATIO = 0.379276

# RECORDED max |sigma_LCP - sigma_exact| in bohr^2 on N2, the quantity Fig. 4's
# LINEAR axis actually shows (p. 012710-8): 0.531 for 0->1 against 8.71 for
# 0->0, a factor of 16. The 0->1 panel's axis runs to 14 bohr^2, so 0.531 is
# ~4% of it -- invisible in print, which is what lets the Fig. 8 caption say
# the calculations are "practically the same" there. The elastic channel's
# 8.71 bohr^2, against a 35 bohr^2 peak, is not invisible at all.
#
# GRID-SAMPLED, like the bands above. On the 101-energy figure grid the same
# two quantities are 0.707 (0->1) and 11.65 (0->0) bohr^2 -- the factor of 16
# survives, and both stay inside the ceiling/floor asserted here, but the
# ABSOLUTE numbers in this comment are the coarse-sweep ones. Likewise the
# "within 5.2% wherever sigma exceeds half its peak" figure in
# `test_n2_0to1_agrees_on_the_scale_fig_8_asserts`: densely that criterion
# covers 23 energies and the worst deviation is 9.9% (at E=0.073). Quote 9.9%
# when describing the printed figure, 5.2% only for this grid.
_N2_ABS_DEV_CEILING_01 = 1.0
_N2_ABS_DEV_FLOOR_00 = 5.0

# RECORDED floors under choice A's own error, so the A-vs-B ORDERING test
# cannot pass vacuously if choice A ever silently improved. Measured worst
# |sigma_A/sigma_exact - 1| = 0.146023 (N2) and 0.434720 (F2); the floors keep
# ~30% headroom below those. They are a claim in their own right -- that the
# Born-Oppenheimer breakdown of PRA 77 Sec. VI A is PRESENT at this size --
# not just a guard.
_CHOICE_A_ERROR_FLOOR: dict[str, float] = {"N2": 0.10, "F2": 0.30}


@functools.cache
def _comparison(molecule: str) -> VeComparison:
    """All six curves for `molecule`, computed once per session.

    One `compare` call is expensive -- MEASURED end to end on the 12-core dev
    machine: **N2 94.9 s** (11 energies; exact-2D sweep ~63 s, LCP pole walk
    4 s, four NRM curves ~2.1-2.6 s per energy each) and, for F2, **705 s at
    the two anchors the gate runs** -- against 1818 s for the full six-energy
    sweep (see `_ENERGIES`, where both were re-measured back to back; an
    earlier 745.9 s six-energy figure did NOT reproduce and is superseded).
    F2's cost is dominated by the exact-2D sweep at 128568 unknowns and ~10 GB
    peak RSS, on top of a 22 s LCP pole walk, 33 s of grid setup, 2 x 14 s of
    ingredients and ~10-12 s per energy per NRM curve.
    Every test in this module therefore shares one result per molecule rather
    than recomputing it.
    """
    return compare(molecule, _ENERGIES[molecule], VPRIMES, n_states=_N_STATES)


@pytest.mark.slow
@pytest.mark.parametrize("molecule", MOLECULES)
def test_every_route_is_finite_and_open(molecule: str) -> None:
    """No route may be zero or non-finite -- otherwise every ratio is 0/0.

    NOT a formality: `nrm_ve_cross_section` returns exactly `0.0` for a closed
    channel (`E_tot - eps_vf <= 0`), and `driven.ve_cross_section` does the
    same, so an energy grid edited below the 0->1 threshold would make every
    comparison below vacuously true. This is what stops that silently.
    """
    c = _comparison(molecule)
    for name, sigma in (
        ("exact", c.sigma_exact),
        ("lcp", c.sigma_lcp),
        ("nrm_a", c.sigma_nrm_a),
        ("nrm_a_nobg", c.sigma_nrm_a_nobg),
        ("nrm_b", c.sigma_nrm_b),
        ("nrm_b_nobg", c.sigma_nrm_b_nobg),
    ):
        assert np.all(np.isfinite(sigma)), f"{molecule}/{name} is not finite:\n{sigma}"
        assert np.all(sigma > 0.0), (
            f"{molecule}/{name} has a non-positive sigma at "
            f"{c.energies[np.any(sigma <= 0.0, axis=1)]} Ha -- a channel is "
            "closed there, so the comparison is vacuous"
        )


@pytest.mark.slow
@pytest.mark.parametrize("molecule", MOLECULES)
def test_nrm_with_background_reproduces_the_exact_oracle(molecule: str) -> None:
    """The headline: PRA 77 claims nonlocal+bg is essentially exact for VE.

    Choice B (`AsymptoticDiscreteState`, the R-independent bound state) plus
    the Eq. (37) background terms, against `qscat.core.driven.ve_cross_section`.

    RECORDED (`n_states=100`), `sigma_NRM(B)+bg / sigma_exact`. N2's rows are
    the 11 energies the gate runs; F2's rows are the ORIGINAL SIX-ENERGY
    measurement (0.02, 0.03, 0.04, 0.05, 0.07, 0.09), of which the gate now
    runs the first and third -- the two that define the band asserted below:

        N2  0->0  1.00019 0.99978 0.99953 0.99916 0.99886 0.99882
                  0.99902 0.99928 0.99950 0.99966 0.99977
            0->1  1.00065 1.00014 0.99843 0.99706 0.99920 0.99897
                  0.99871 0.99858 0.99857 0.99862 0.99868
        F2  0->0  0.99805 0.99901 0.99983 1.00032 1.00044 1.00014
            0->1  1.00692 0.99822 0.99623 0.99640 0.99766 0.99816

    Unlike the DA result this rests on the channel PRA 77 plots for every
    molecule in its study, and it holds for the elastic channel -- the one the
    paper says a bare LCP curve cannot get right -- as well as the inelastic.
    """
    lo, hi = _BANDS[molecule]
    c = _comparison(molecule)
    ratio = c.sigma_nrm_b / c.sigma_exact
    assert np.all((ratio > lo) & (ratio < hi)), (
        f"{molecule}: NRM(B)+bg / exact left the recorded band [{lo}, {hi}]:\n{ratio}"
    )


@pytest.mark.slow
@pytest.mark.parametrize("molecule", MOLECULES)
@pytest.mark.parametrize("choice", ["a", "b"])
def test_background_matters_more_for_elastic_than_inelastic(molecule: str, choice: str) -> None:
    """A falsifiable prediction: the background matters more for the elastic
    channel than for the first inelastic one.

    WHAT THE TRACKED REFERENCE NOTE SUPPORTS. `reference/literature/
    houfek-2008-pra77-012710.md` records "background contributions to VE
    T-matrix matter for BOTH elastic and inelastic channels, largest for F2
    (broadest resonance)" at p. 012710-1 (abstract) and p. 012710-10
    (Sec. VI B/C). The stronger, ordered form tested here -- that the
    importance DECREASES with increasing inelasticity, p. 012710-8 -- is not
    carried by that note as a quoted sentence with its own locator, so it is
    asserted here as a MEASUREMENT of this repo's own solvers rather than as
    a reproduction of a quoted claim. Adding it to the reference note (with
    a page-checked locator) is the way to promote it.

    Measured as the median relative change dropping `T^bg` makes,
    `|sigma_bg - sigma_nobg| / sigma_bg`:

                    0->0      0->1     ratio
        N2  A      0.7102    0.10655    6.7x
        N2  B      0.2221    0.007233  30.7x
        F2  A      0.99844   0.93392    1.07x
        F2  B      0.86964   0.23393    3.7x

    The direction holds in all four combinations. The MAGNITUDES also carry
    the paper's second claim -- that the background is largest for the
    broadest resonance: F2's background is the whole cross section (dropping
    it costs 87% of the elastic answer) while N2's costs 22%.
    """
    c = _comparison(molecule)
    with_bg = c.sigma_nrm_a if choice == "a" else c.sigma_nrm_b
    without_bg = c.sigma_nrm_a_nobg if choice == "a" else c.sigma_nrm_b_nobg
    rel_00 = np.abs(with_bg[:, 0] - without_bg[:, 0]) / with_bg[:, 0]
    rel_01 = np.abs(with_bg[:, 1] - without_bg[:, 1]) / with_bg[:, 1]
    assert np.median(rel_00) > np.median(rel_01), (
        f"{molecule}/{choice}: background not more important for elastic: "
        f"median 0->0 {np.median(rel_00):.4g} vs 0->1 {np.median(rel_01):.4g}"
    )


@pytest.mark.slow
def test_n2_0to1_agrees_on_the_scale_fig_8_asserts() -> None:
    """PRA 77's Fig. 8 caption drops the N2 0->1 panel "because results of all
    calculations are practically the same in this particular case"
    (p. 012710-10, quoted in `reference/literature/houfek-2008-pra77-012710.md`).
    Fig. 8's own curves are the four Fig. 4 defines -- exact, LCP, nonlocal,
    nonlocal+bg -- so "all calculations" includes the LOCAL one.

    That claim is about curves on a LINEAR axis running to 14 bohr^2
    (p. 012710-8, Fig. 4 middle panel), NOT about a bounded ratio. Both are
    measured here, and only the first is what the caption asserts.

    RECORDED, `sigma / sigma_exact` on N2 0->1 AT THE 11 GATE ENERGIES over
    0.06-0.16 Ha (the dense 101-energy figure grid in brackets -- see the
    comments on `_BAND_N2_01` and `_N2_ABS_DEV_CEILING_01`):

        A + bg    [0.85398, 1.05868]   (dense [0.82852, 1.13260])
        B + bg    [0.99706, 1.00065]   (dense [0.99454, 1.00065])
        LCP       [0.37928, 1.20176]   (dense [0.37928, 1.32884])

    The LCP ratio looks alarming and is not. Per energy it runs 0.379, 0.948,
    0.955, 1.202, 1.009, 0.963, 0.924, 0.869, 0.801, 0.729, 0.659 -- and
    wherever sigma exceeds half its 10.23 bohr^2 peak the LCP is within
    **5.2%** of exact. (5.2% is an ANCHOR-GRID figure. On the 101-energy
    figure grid the peak resolves to 11.48 bohr^2, the same criterion covers
    23 energies, and the worst deviation is 9.9% -- still small next to the
    ELASTIC channel's 61% under the identical cut. See the comment on
    `_N2_ABS_DEV_CEILING_01`.) The 0.379 is at E=0.06, where sigma_exact is
    0.287 bohr^2, 2.8% of peak; the absolute miss there is 0.18 bohr^2, about
    1% of the panel's axis. Across these 11 energies the LCP's worst ABSOLUTE
    deviation is 0.53 bohr^2 (densely 0.71). Fig. 8's caption survives intact.

    The mechanism is the one the paper names at p. 012710-9: the bare LCP's
    resonance width `Gamma(R)` is energy INDEPENDENT, so it degrades where the
    channel energy is far from `E_res` -- i.e. in the wings, which is exactly
    where sigma is small. It also sits inside the cross-model band this repo
    already gates the LCP at (`validation/n2/reference.py`,
    `ANCHOR_FACTOR = 3.0`; 1/0.379 = 2.64).

    CONTRAST THE ELASTIC CHANNEL, asserted below as the counter-example: there
    the LCP drifts 4.57 -> 0.106 monotonically and its worst absolute deviation
    is 8.71 bohr^2 against a 35 bohr^2 peak -- a failure at LARGE sigma, plainly
    visible in print, and it is this paper's own missing-background claim
    (p. 012710-1, 012710-10). A single ratio band cannot tell those two apart;
    the absolute-deviation assertions below can.

    None of this is an artifact of this driver. Rebuilding the LCP through the
    project's own independent N2 route (`projects.n2_ti_cross_section.vres.
    vres_on_grid` on its own 428-point nuclear grid, instead of
    `qscat.core.lcp.local_complex_potential` on the deck used here) reproduces
    it to ~0.1%: 0.10882 vs 0.10890 bohr^2 at E=0.06 and 5.8686 vs 5.8735 at
    E=0.08.
    """
    lo, hi = _BAND_N2_01
    c = _comparison("N2")
    for label, sigma in (("nrm-a", c.sigma_nrm_a), ("nrm-b", c.sigma_nrm_b)):
        ratio = sigma[:, 1] / c.sigma_exact[:, 1]
        assert np.all((ratio > lo) & (ratio < hi)), (
            f"{label}: N2 0->1 left the recorded band [{lo}, {hi}]:\n{ratio}"
        )

    # What Fig. 8's caption actually claims: on the paper's linear axis the
    # 0->1 curves lie on top of each other, INCLUDING the LCP's.
    abs_dev_01 = np.abs(c.sigma_lcp[:, 1] - c.sigma_exact[:, 1]).max()
    assert abs_dev_01 < _N2_ABS_DEV_CEILING_01, (
        f"N2 0->1: the LCP's worst absolute deviation is {abs_dev_01:.4g} bohr^2, "
        f"over the recorded {_N2_ABS_DEV_CEILING_01} -- it would now be visible on "
        "Fig. 4's linear axis, so the Fig. 8 caption's reading no longer holds"
    )

    # And the elastic channel is the counter-example that keeps the assertion
    # above from being vacuous: same molecule, same run, a visible failure.
    abs_dev_00 = np.abs(c.sigma_lcp[:, 0] - c.sigma_exact[:, 0]).max()
    assert abs_dev_00 > _N2_ABS_DEV_FLOOR_00, (
        f"N2 0->0: the LCP's worst absolute deviation is only {abs_dev_00:.4g} "
        f"bohr^2, under the recorded floor {_N2_ABS_DEV_FLOOR_00} -- the elastic "
        "channel was the RECORDED visible failure (8.71); if it has improved, the "
        "missing-background story in this module's docstring needs re-measuring"
    )

    # The ratio departure in the wings, recorded so a change is noticed.
    lcp_ratio = c.sigma_lcp[:, 1] / c.sigma_exact[:, 1]
    assert lcp_ratio.min() < lo, (
        "the LCP's N2 0->1 ratio no longer departs in the wings -- the RECORDED "
        f"observation is that it does (worst {_N2_01_LCP_WORST_RATIO:.6g}); "
        f"measured now: {lcp_ratio}"
    )


@pytest.mark.slow
@pytest.mark.parametrize("molecule", MOLECULES)
def test_nrm_b_beats_the_lcp_everywhere(molecule: str) -> None:
    """Choice B is closer to exact than the LCP, at every energy and channel.

    A comparison rather than a threshold: whatever the absolute agreement
    turns out to be, keeping the energy dependence, the nonlocality and the
    background terms must not leave us further from the exact answer than
    throwing them away did.

    RECORDED worst `|sigma/sigma_exact - 1|` over the window:

        N2   LCP 3.5678   B 0.002938    factor 1214
        F2   LCP 0.99982  B 0.006923    factor  144

    On F2 the gap is not close: the LCP elastic cross section is three to four
    orders of magnitude below exact at every energy measured -- ratio 1.8e-4
    to 2.1e-2 over the recorded six-energy sweep, of which the gate now checks
    the two anchors -- which is PRA 77's "largest for the broadest resonance" statement
    about the background terms, seen from the side of the model that omits
    them.
    """
    c = _comparison(molecule)
    err_lcp = np.abs(c.sigma_lcp / c.sigma_exact - 1.0)
    err_nrm = np.abs(c.sigma_nrm_b / c.sigma_exact - 1.0)
    assert np.all(err_nrm < err_lcp), (
        f"{molecule}: NRM(B)+bg is not closer to exact than the LCP everywhere.\n"
        f"  lcp errors:\n{err_lcp}\n  nrm errors:\n{err_nrm}"
    )


@pytest.mark.slow
@pytest.mark.parametrize("molecule", MOLECULES)
def test_choice_a_is_worse_than_choice_b(molecule: str) -> None:
    """The Born-Oppenheimer breakdown of PRA 77 Sec. VI A, in VE.

    Both choices get the SAME background terms here, so this isolates the
    discrete state itself. RECORDED worst `|sigma/sigma_exact - 1|`:

        N2   A 0.14602   B 0.002938
        F2   A 0.43472   B 0.006923

    Choice A's error is worst in the inelastic channel (N2 0->1 [0.85398,
    1.05868] against 0->0 [0.93356, 1.01208]; F2 0->1 [0.56528, 1.14013]
    against 0->0 [0.82841, 1.00192]) and, on F2, worst at the low-energy end
    -- the same shape the DA gate recorded for choice A on F2.
    """
    c = _comparison(molecule)
    err_a = np.abs(c.sigma_nrm_a / c.sigma_exact - 1.0)
    err_b = np.abs(c.sigma_nrm_b / c.sigma_exact - 1.0)
    assert err_b.max() < err_a.max(), (
        f"{molecule}: choice A is no longer worse than choice B "
        f"(worst A {err_a.max():.5g}, worst B {err_b.max():.5g})"
    )
    # The ordering alone would still hold if choice A improved to 0.9%, quietly
    # falsifying the 14.6%/43.5% recorded above. Pin the SIZE too.
    floor = _CHOICE_A_ERROR_FLOOR[molecule]
    assert err_a.max() > floor, (
        f"{molecule}: choice A's worst error is now {err_a.max():.5g}, under the "
        f"recorded floor {floor} -- the Born-Oppenheimer breakdown this test "
        "documents (RECORDED 0.146023 for N2, 0.434720 for F2) has shrunk, so the "
        "recorded numbers need re-measuring rather than this floor lowering"
    )
