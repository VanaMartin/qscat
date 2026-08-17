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

so the paper's prediction (B near-exact for DA, A degraded by a Born-
Oppenheimer breakdown) holds on F2 and INVERTS on NO. `test_nrm_b_beats_the_
lcp_on_f2` and `test_nrm_b_does_not_beat_the_lcp_on_no` are the two halves of
that finding, and both are load-bearing assertions.

ENERGIES ARE PER MOLECULE because the DA thresholds are: measured
`eps_e - eps[0]` is **-0.0691 Ha for F2** (open at every positive `E`) and
**+0.1719 Ha for NO**. A shared grid would compare zeros for NO, and
`test_the_da_channel_is_open` exists to keep that from ever passing silently.
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
# F2: measured [0.998669, 1.018735]; the band keeps ~50% headroom on the 1.9%
# maximum deviation. NO: measured [4.672e-8, 4.719e-6]; the band is the
# enclosing decade pair, since the finding there is an order-of-magnitude
# failure and pinning its mantissa would be false precision.
_BANDS_B: dict[str, tuple[float, float]] = {
    "F2": (0.99, 1.03),
    "NO": (1e-8, 1e-5),
}

# RECORDED sigma_NRM/sigma_exact bands, choice A (`PhysicalDiscreteState`).
# F2: measured [0.292258, 0.900562] -- a systematic under-prediction that
# worsens toward threshold. NO: measured [0.064085, 244.229] -- choice A is
# nearly energy-independent (1.03e-3 -> 4.19e-4) while the exact sigma_DA falls
# four decades across the anchors, so the ratio sweeps almost four decades.
_BANDS_A: dict[str, tuple[float, float]] = {
    "F2": (0.25, 0.95),
    "NO": (0.03, 500.0),
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
def test_nrm_b_does_not_beat_the_lcp_on_no() -> None:
    """RECORDED NEGATIVE RESULT: on NO, choice B is worse than the LCP.

    The mirror image of `test_nrm_b_beats_the_lcp_on_f2`, and the reason the
    spec's primary claim is asserted for F2 alone. Measured with the same code
    and the same `n_states`:

        sigma_DA (bohr^2), E = 0.175 / 0.180 / 0.185 / 0.190 / 0.200 Ha
          exact  1.614e-2  1.566e-3  2.266e-5  7.080e-5  1.718e-6
          LCP    1.305e-4  1.330e-4  1.352e-4  1.347e-4  1.296e-4
          B      7.540e-10 2.833e-10 1.069e-10 4.055e-11 5.921e-12

    Choice B under-predicts by five to eight orders of magnitude. The
    comparison is made on |log10(sigma/sigma_exact)| rather than the ratio
    error above: over-prediction is unbounded there while under-prediction
    saturates at 1, so the plain ratio error would score B (~1.0 everywhere)
    as "better" than the LCP wherever the LCP over-predicts, which is an
    artifact of the metric and not a result.

    NOT a discretisation artifact. sigma_B is converged in the electronic grid
    to five significant figures (r_max 16 -> 24 bohr, DVR order 8 -> 10,
    n_complex 6 -> 8: 7.540e-10 / 7.540e-10 / 7.543e-10 / 7.541e-10 at
    E=0.175) and in the state sum to numerical identity, on the same nuclear
    deck the exact and LCP routes use. The mechanism is measured: NO's
    doorway `chi_0` peaks at R=2.134 bohr, INSIDE the resonant region
    (crossing at R_c=2.288, Gamma_LCP=1.32e-2 there), where the R-independent
    `phi_b` is a poor discrete state -- |V_d(B) - V_d(LCP)| = 0.023 Ha at the
    peak. F2's `chi_0` peaks at R=2.745, on the BOUND side of its own crossing
    (R_c=2.597, Gamma_LCP=0), where `phi_b` is nearly the true discrete state
    and the same difference is 0.0053 Ha.
    """
    c = _comparison("NO")
    err_lcp = np.abs(np.log10(c.sigma_lcp / c.sigma_exact))
    err_nrm = np.abs(np.log10(c.sigma_nrm_b / c.sigma_exact))
    assert np.all(err_nrm > err_lcp), (
        f"NO: NRM(B) is no longer worse than the LCP at every anchor -- the "
        f"recorded negative result has changed and the physics note needs "
        f"revisiting.\n  lcp log-errors: {err_lcp}\n  nrm log-errors: {err_nrm}"
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
@pytest.mark.parametrize("molecule", MOLECULES)
def test_nrm_a_ratio_band(molecule: str) -> None:
    """RECORDED band for choice A's sigma_NRM/sigma_exact.

    F2 measured: 0.292258 0.671111 0.843850 0.853668 0.900562
    NO measured: 0.064085 0.571017 32.8576  8.58031  244.229
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

    Asserted for F2 only: on NO the ordering inverts (choice A lands within
    1-2 orders of the exact answer while choice B is 5-8 orders below it), and
    that inversion is recorded by the two ratio-band tests above.
    """
    c = _comparison("F2")
    err_a = np.abs(c.sigma_nrm_a / c.sigma_exact - 1.0)
    err_b = np.abs(c.sigma_nrm_b / c.sigma_exact - 1.0)
    assert np.all(err_b < err_a), (
        f"F2: choice A is no longer the degraded one.\n  A errors: {err_a}\n  B errors: {err_b}"
    )
