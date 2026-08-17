"""GATE (validation check 3): NRM sigma_DA against the exact-2D oracle.

Every band below is RECORDED from the first converged run, not a preset
target -- the design spec requires the success criterion here to be a
measurement. The run that produced each one is written into the test that
asserts it; their physical reading is in
`docs/physics/nonlocal-resonance-model.md`.

The headline measurement, sigma/sigma_exact at the anchors:

    F2   LCP  0.263  0.465  0.889  1.425  1.733
         A    0.292  0.671  0.844  0.854  0.901
         B    1.019  1.003  0.9994 0.9987 0.9990
    NO   LCP  8.1e-3 8.5e-2 5.97   1.90   75.4
         A    6.4e-2 0.571  32.9   8.58   244
         B    4.7e-8 1.8e-7 4.7e-6 5.7e-7 3.4e-6

On F2 -- the only molecule for which PRA 77 publishes a DA cross section at
all -- the paper's prediction holds exactly: choice B reproduces the exact 2-D
oracle and choice A is degraded by the Born-Oppenheimer breakdown Sec. VI A
describes. `test_nrm_b_beats_the_lcp_on_f2` is that claim, and it is the
gate's one load-bearing physical assertion.

NO IS OUTSIDE THE PAPER'S TESTED RANGE. Its DA channel opens at +0.1719 Ha,
above every energy window PRA 77 plots for NO (0.01-0.08 Ha), so NO DA is
energetically shut throughout the published data -- Fig. 6's DA panel and
Fig. 8's DA panel are F2 only. Nothing measured here for NO contradicts the
literature; the literature is silent. What NO shows is recorded as an
observation (`test_no_approximations_are_all_flat_below_a_structured_exact`),
not asserted as a law, and choice B's absolute scale there is unexplained.

ENERGIES ARE PER MOLECULE because the DA thresholds are: measured
`eps_e - eps[0]` is **-0.0691 Ha for F2** (open at every positive `E`) and
**+0.1719 Ha for NO**. A shared grid would compare zeros for NO, and
`test_the_da_channel_is_open` exists to keep that from ever passing silently.

COVERAGE LIMIT worth knowing when reading a NO number here: the underlying
discrete-continuum coupling `v_dk_plus` is gated against the LCP's independent
ECS-pole width on **F2 only**
(`libs/qscat/tests/test_nrm_coupling.py::test_gamma_matches_the_lcp_width`),
and that gate's tight window -- `Gamma/E < 0.35` and `E_res > 0.02` Ha -- is
EMPTY on NO. Measured over NO's 47 genuinely-tracked open points, `Gamma/E`
rises monotonically inward from 0.124 and crosses 0.35 at R ~ 2.196 where
`E_res` is still 0.0158 Ha, so the last narrow point (0.338 at E_res=0.0150)
sits below the energy floor and the first point above the floor is already
broad (0.409 at E_res=0.0202). The two cuts do not overlap: NO's resonance is
narrow only where it is also near threshold. So on NO the coupling itself
carries no tight cross-check, only the order-of-magnitude one.
"""

from __future__ import annotations

import functools

import numpy as np
import numpy.typing as npt
import pytest

from validation.diatomic.nrm import NrmComparison, compare

MOLECULES = ["F2", "NO"]

ENERGIES: dict[str, npt.NDArray[np.float64]] = {
    "F2": np.array([0.010, 0.020, 0.030, 0.040, 0.050]),
    "NO": np.array([0.175, 0.180, 0.185, 0.190, 0.200]),
}

# MEASURED state-sum truncations (Eq. 60). Both molecules and both discrete-
# state choices were laddered over n_states = 10/25/40/55/70/80/90/100/110/120/
# all, with the ingredients built once and reused:
#
#   F2 / B  (task-8-report.md, E=0.03): converged to <1% at 75, bit-identical
#           to the untruncated 131-state sum from 95 up.
#   F2 / A  (E=0.02/0.03/0.04): 3e-4 at 40, plateau (<2e-12) from 90 up.
#   NO / B  (E=0.18/0.24/0.30): 4.7e-2 at 55, 5.7e-5 at 70, exact from 100 up.
#           Low n OVERSHOOTS here (rel. 312 at n=10, E=0.18) where F2
#           undershoots -- the two convergence shapes are qualitatively
#           different, which is why NO was re-laddered rather than inheriting
#           F2's number.
#   NO / A  (E=0.175/0.185/0.200): 3e-4 at 40, plateau (<4e-14) from 90 up.
#
# 100 is the smallest round value that reproduces the untruncated sum to
# numerical identity for all four combinations.
N_STATES: dict[str, int] = {"F2": 100, "NO": 100}

# RECORDED sigma_NRM/sigma_exact bands, choice B (`AsymptoticDiscreteState`).
#
# F2: measured [0.998669, 1.018735] -- deviations from unity of 1.9% (E=0.010,
# nearest threshold), 0.33% (0.020), then 0.059% / 0.13% / 0.097%. Those are
# quotable as stated, because the oracle's own floor is known and far below
# them: an INDEPENDENT 1000-point resonance-aware grid gives sigma_DA(F2, 0.03)
# = 1.6562 against this 974-point deck's 1.65611 (docs/physics/
# discretisation-tuning.md, validation/tuning/test_resonance_aware.py), i.e.
# the exact reference is converged to ~5e-5 relative. The deviations above are
# 11x to 340x that floor, so they are physics, not grid noise -- including
# their rise toward threshold, which is the physically interesting part. The
# band keeps ~50% headroom on the largest.
#
# NO: measured [4.672e-8, 4.719e-6]; the band is the enclosing decade pair,
# since what it records is an unexplained order-of-magnitude offset in a
# regime no published data covers, and pinning its mantissa would be false
# precision.
_BANDS_B: dict[str, tuple[float, float]] = {
    "F2": (0.99, 1.03),
    "NO": (1e-8, 1e-5),
}

# RECORDED sigma_NRM/sigma_exact band, choice A (`PhysicalDiscreteState`).
# F2 ONLY: measured [0.292258, 0.900562] -- a systematic under-prediction that
# worsens toward threshold, the Born-Oppenheimer breakdown of PRA 77 Sec. VI A.
#
# NO IS DELIBERATELY ABSENT. Its choice-A ingredients are built on a `phi_d`
# that is genuinely DISCONTINUOUS at R = 2.2657 bohr, the node where the state
# switches from the scattering branch to the bound one -- and that node carries
# the largest |chi_0| on the whole grid (0.1996). Eq. (60) is bilinear in
# `V_dn(R_i) V_dn(R_j)`, so the tracked P-space basis is unusable across that
# step (measured `|c_product(prev, cur)|` down to 3.3e-15, with 62 affected
# states inside the n_states=100 truncation). Gating those numbers would assert
# values the report says must not be quoted. F2 escapes this only because its
# deck happens to place the two crossing nodes 0.0005 bohr apart, so its
# `phi_d` barely changes across the switch (min overlap 0.891, nothing < 0.5).
# The separate small-R pole freeze below R = 1.5187 is NOT the reason: |chi_0|
# is 6.9e-16 there, so it cannot move sigma_DA.
_BANDS_A: dict[str, tuple[float, float]] = {
    "F2": (0.25, 0.95),
}


@functools.cache
def _comparison(molecule: str) -> NrmComparison:
    """All four routes for `molecule`, computed once per session.

    One `compare` call is ~10 minutes (measured: F2 454 s for the exact-2D
    sweep alone, plus 24 s for the LCP curve and 2 x ~95 s for the two NRM
    routes; NO 214 s + 17 s + 2 x ~60 s), so every test in this module shares
    one result rather than recomputing it.
    """
    return compare(molecule, ENERGIES[molecule], n_states=N_STATES[molecule])


@pytest.mark.slow
@pytest.mark.parametrize("molecule", MOLECULES)
def test_the_da_channel_is_open(molecule: str) -> None:
    """No route may be zero at any anchor -- otherwise every ratio is 0/0.

    NOT a formality. The DA channel opens at `E = eps_e - eps[0]`, which is
    +0.1719 Ha for NO: on the shared 0.010-0.050 Ha grid the task brief
    originally proposed, all four routes return exactly 0.0 for NO and every
    comparison below would be vacuously satisfied by `nan`s. This test is what
    stops a future energy edit from re-introducing that silently.
    """
    c = _comparison(molecule)
    for name, sigma in (
        ("exact", c.sigma_exact),
        ("lcp", c.sigma_lcp),
        ("nrm_a", c.sigma_nrm_a),
        ("nrm_b", c.sigma_nrm_b),
    ):
        assert np.all(np.isfinite(sigma)), f"{molecule}/{name} is not finite: {sigma}"
        assert np.all(sigma > 0.0), (
            f"{molecule}/{name} has a non-positive sigma_DA at "
            f"{c.energies[sigma <= 0.0]} Ha -- the DA channel is closed there, "
            "so the comparison is vacuous"
        )


@pytest.mark.slow
def test_nrm_b_beats_the_lcp_on_f2() -> None:
    """Choice B is closer to exact than the LCP is, at every F2 anchor.

    The spec's primary claim, and a comparison rather than a threshold:
    whatever the absolute agreement turns out to be, keeping the energy
    dependence and the nonlocality must not leave us further from the exact
    answer than throwing them away did.

    RECORDED (n_states=100, energies 0.010-0.050 Ha):

        |sigma/sigma_exact - 1|
          LCP  0.73718  0.53469  0.11092  0.42454  0.73328
          B    0.01873  0.00329  0.00059  0.00133  0.00097

    i.e. choice B recovers essentially the whole of the LCP's documented
    0.263 -> 1.736 error sweep, beating it at every anchor by a factor of 39
    (E=0.010) to 758 (E=0.050).
    """
    c = _comparison("F2")
    err_lcp = np.abs(c.sigma_lcp / c.sigma_exact - 1.0)
    err_nrm = np.abs(c.sigma_nrm_b / c.sigma_exact - 1.0)
    assert np.all(err_nrm <= err_lcp), (
        f"F2: NRM(B) is not closer to exact than the LCP at every anchor.\n"
        f"  lcp errors: {err_lcp}\n  nrm errors: {err_nrm}"
    )


@pytest.mark.slow
def test_no_approximations_are_all_flat_below_a_structured_exact() -> None:
    """RECORDED OBSERVATION on NO, in a regime no published data covers.

    This is a measurement, NOT a requirement. Nothing here asserts that an
    approximation *must* fail; if a later change makes one of them track the
    exact answer, this test should be updated to record that, not treated as
    a regression.

    NO's DA channel opens at E = +0.1719 Ha. PRA 77 publishes no NO (or N2)
    DA cross section at all -- its Fig. 6 DA panel and Fig. 8's DA panel are
    F2 only, and the NO windows it plots (0.01-0.08 Ha) lie entirely below
    that threshold, so NO DA is energetically shut throughout everything the
    paper tested. The Sec. VI B claim that choice B "gives exact results" for
    DA rests on the single F2 panel. Whatever happens here is therefore
    outside the model's tested range and contradicts nothing published.

    What is measured, sigma_DA (bohr^2) at E = 0.175 / 0.180 / 0.185 / 0.190
    / 0.200 Ha:

        exact  1.614e-2  1.566e-3  2.266e-5  7.080e-5  1.718e-6
        LCP    1.305e-4  1.330e-4  1.352e-4  1.347e-4  1.296e-4
        A      1.034e-3  8.945e-4  7.446e-4  6.075e-4  4.195e-4
        B      7.540e-10 2.833e-10 1.069e-10 4.055e-11 5.921e-12

    ONE phenomenon, not three failures: the converged exact cross section
    swings four orders of magnitude with real structure across a 0.025 Ha
    span, while ALL THREE approximations are nearly flat over it -- the LCP
    within a factor of 1.04, choice A within 2.5, choice B within 127 (and
    127 over four orders of exact variation is still flat by comparison).
    They differ in absolute scale, not in their inability to follow the
    structure. That is one observation about local, nearly-energy-independent
    resonance models in a steeply-varying near-threshold regime.

    Choice B's absolute scale (5-8 orders below exact) is UNEXPLAINED. It is
    not a discretisation artifact -- sigma_B is converged in the electronic
    grid to five significant figures (r_max 16 -> 24 bohr, DVR order 8 -> 10,
    n_complex 6 -> 8: 7.540e-10 / 7.540e-10 / 7.543e-10 / 7.541e-10 at
    E=0.175) and in the state sum to numerical identity, on the same nuclear
    deck the exact and LCP routes use, and the LCP local-limit bridge
    reproduces `lcp_da_cross_section` to 2.3e-14 on that deck. An
    equation-by-equation audit against Eq. (55)-(61) found no defect. A
    doorway-position explanation was proposed and REFUTED: raising `v_init`
    moves F2's doorway inside its own resonant region (overlap fraction
    0.37 -> 0.71) and choice B stays exact there (0.9994-1.0004).

    The assertion below records the flatness, which is the robust part of the
    observation, rather than the ordering of the three errors.
    """
    c = _comparison("NO")
    exact_swing = c.sigma_exact.max() / c.sigma_exact.min()
    assert exact_swing > 1e3, (
        f"NO: the exact sigma_DA no longer varies by orders across these "
        f"anchors (swing {exact_swing:.3g}); the observation this test "
        "records is about a steeply-varying exact reference"
    )
    for name, sigma, recorded in (
        ("lcp", c.sigma_lcp, 1.04),
        ("nrm_a", c.sigma_nrm_a, 2.5),
        ("nrm_b", c.sigma_nrm_b, 127.0),
    ):
        swing = sigma.max() / sigma.min()
        assert swing < 10.0 * recorded, (
            f"NO/{name}: swing {swing:.4g} over the anchors, far above the "
            f"recorded {recorded:.4g} -- this approximation is no longer flat "
            "and the recorded observation needs revisiting"
        )
        assert swing < exact_swing / 10.0, (
            f"NO/{name}: swing {swing:.4g} is no longer flat relative to the "
            f"exact reference's {exact_swing:.4g}"
        )


@pytest.mark.slow
@pytest.mark.parametrize("molecule", MOLECULES)
def test_nrm_b_ratio_band(molecule: str) -> None:
    """RECORDED band for choice B's sigma_NRM/sigma_exact.

    F2 measured: 1.018735 1.003286 0.999413 0.998669 0.999032
    NO measured: 4.672e-8 1.809e-7 4.719e-6 5.727e-7 3.447e-6
    """
    lo, hi = _BANDS_B[molecule]
    c = _comparison(molecule)
    ratio = c.sigma_nrm_b / c.sigma_exact
    assert np.all((ratio > lo) & (ratio < hi)), (
        f"{molecule}: choice-B ratio {ratio} left the recorded band [{lo}, {hi}]"
    )


@pytest.mark.slow
@pytest.mark.parametrize("molecule", sorted(_BANDS_A))
def test_nrm_a_ratio_band(molecule: str) -> None:
    """RECORDED band for choice A's sigma_NRM/sigma_exact. F2 only.

    F2 measured: 0.292258 0.671111 0.843850 0.853668 0.900562

    NO is excluded on purpose -- see `_BANDS_A`. Its choice-A numbers
    (0.064085 0.571017 32.8576 8.58031 244.229) are recorded in the Task 9
    report but are NOT gated, because the P-space basis they are built from is
    untrackable across NO's `phi_d` branch discontinuity.
    """
    lo, hi = _BANDS_A[molecule]
    c = _comparison(molecule)
    ratio = c.sigma_nrm_a / c.sigma_exact
    assert np.all((ratio > lo) & (ratio < hi)), (
        f"{molecule}: choice-A ratio {ratio} left the recorded band [{lo}, {hi}]"
    )


@pytest.mark.slow
def test_choice_a_degrades_against_choice_b_on_f2() -> None:
    """PRA 77's predicted Born-Oppenheimer breakdown of choice A, on F2.

    The paper (Sec. VI A) predicts the physical discrete state degrades for DA
    where the R-independent one does not. Measured on F2, choice A is farther
    from the exact answer than choice B at every anchor, by factors of 38
    (E=0.010) to 266 (E=0.030) in |sigma/sigma_exact - 1|:

        A  0.70774  0.32889  0.15615  0.14633  0.09944
        B  0.01873  0.00329  0.00059  0.00133  0.00097

    Asserted for F2 only, which is the only molecule PRA 77 publishes a DA
    cross section for. On NO both choices are flat against a four-decade exact
    swing and neither tracks it, so ordering them by |sigma/sigma_exact - 1|
    there would rank two approximations that fail the same way; see
    `test_no_approximations_are_all_flat_below_a_structured_exact`.
    """
    c = _comparison("F2")
    err_a = np.abs(c.sigma_nrm_a / c.sigma_exact - 1.0)
    err_b = np.abs(c.sigma_nrm_b / c.sigma_exact - 1.0)
    assert np.all(err_b < err_a), (
        f"F2: choice A is no longer the degraded one.\n  A errors: {err_a}\n  B errors: {err_b}"
    )
