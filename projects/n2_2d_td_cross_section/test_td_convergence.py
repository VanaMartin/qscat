"""The sigma(E) "boomerang" curve + usable window (sub-project #7, Task 5).

ONE ~215s propagation at `TD_WORKING_GRID` (module-scope fixture, via
`convergence.sigma_curve`) produces the correlation functions `c_{v'}(t_n)`;
every energy in `E_GRID` is then a free Tannor-Weeks transform of that SAME
stored trajectory -- exactly the "one propagation, whole curve" claim this
task exists to prove.

**A caveat this file used to document, now withdrawn as measurement error.**
It previously described a "finite-T energy-resolution limit": below E~0.13 Ha
the exact `sigma_TI(E)` has boomerang sub-features spanning only ~0.004 Ha,
comparable to the propagation's `2*pi/T ~ 0.0042` Ha, and the transform was
said to average over them, so pointwise TD-vs-TI agreement there was claimed
to be unattainable (citing ratios 5.7 at E=0.09 and 0.37 at E=0.11) and
fixable only by a longer `T`.

That was wrong. Those ratios came from the order-1 Crank-Nicolson propagator
in use at the time, which under-converges badly over a long propagation. With
order-3 Pade (`qscat.evolution.make_pade_stepper`, the default since), the
same single T=1500 propagation tracks the exact oracle everywhere in the
tested range. Measured 2026-08-17 at `TD_WORKING_GRID`:

    E      sigma_TD    sigma_TI    ratio
    0.10    5.95946     6.12278    0.9733
    0.14    1.06012     1.05869    1.0013
    0.15    0.61850     0.62575    0.9884
    0.16    0.39993     0.40265    0.9932
    0.17    0.28170     0.27927    1.0087
    0.18    0.20843     0.20570    1.0132

Two lessons worth keeping. First, the fine structure is genuinely there in
`sigma_TI` -- what was wrong was the claim that the transform could not
resolve it. Second, the old per-energy tolerances were justified by a
tolerance that "widens with distance from the incident spectral peak
(p0**2/2 = 0.125 Ha)"; the measurements above run the other way -- E=0.18,
the furthest point, is among the best at 1.3%, while E=0.10 is the worst at
2.7%. The trend did not exist either. A single uniform tolerance replaces the
per-energy dict.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.dvr import TensorGrid

from projects.n2_2d_cross_section.cross_section_2d import ve_cross_section_2d
from projects.n2_2d_cross_section.hamiltonian2d import MU
from projects.n2_2d_td_cross_section.convergence import (
    TD_WORKING_GRID,
    sigma_curve,
    td_working_tgrid,
    usable_window,
)
from projects.n2_ti_cross_section.vibrational import vibrational_states

TG: TensorGrid = td_working_tgrid()
EPS, CHI = vibrational_states(TG.grids[1], MU, 4)

V_INIT = 0
VPRIMES = [1]

# One uniform tolerance, not a per-energy dict. Every energy here agrees with
# the exact TI oracle to <= 2.7% at the order-3 Pade default (measurements in
# the module docstring), so the old 0.10-0.35 per-energy tolerances -- fitted
# to order-1 Crank-Nicolson error and to a peak-distance trend that the
# measurements contradict -- were 4-26x looser than the code needs and could
# no longer detect a regression.
#
# 0.05 sits ~2x above the worst measured deviation (2.7% at E=0.10), leaving
# headroom for cross-architecture BLAS differences without going slack.
RTOL = 0.05
E_GRID = np.array([0.10, 0.14, 0.15, 0.16, 0.17, 0.18])


@pytest.fixture(scope="module")
def curve() -> np.ndarray:
    """The ONE propagation this whole file needs, transformed at every `E_GRID`
    energy for free (see module docstring)."""
    return sigma_curve(TG, EPS, CHI, V_INIT, VPRIMES, E_GRID)


@pytest.mark.slow
def test_sigma_curve_matches_ti_on_the_smooth_branch(curve: np.ndarray) -> None:
    assert curve.shape == (len(E_GRID), len(VPRIMES))

    sigma_ti = ve_cross_section_2d(TG, EPS, CHI, V_INIT, VPRIMES, E_GRID)
    for i, e in enumerate(E_GRID):
        sigma_td = float(curve[i, 0])
        assert sigma_td >= 0.0
        assert sigma_td == pytest.approx(float(sigma_ti[i, 0]), rel=RTOL), (
            f"E={e}: sigma_TD={sigma_td!r} vs sigma_TI={float(sigma_ti[i, 0])!r} "
            f"outside rtol={RTOL}"
        )


def test_usable_window_is_nonempty_and_sensible() -> None:
    """Cheap (no propagation): `|eta_incident(E)|` alone, evaluated on the
    electronic grid, brackets the spectral peak and both TI anchors.
    """
    e_scan = np.linspace(0.04, 0.22, 19)  # step 0.01
    (e_lo, e_hi), eta_abs = usable_window(
        TG.grids[0], e_scan, wp_in=TD_WORKING_GRID["wp_in"], frac=0.5
    )
    assert eta_abs.shape == e_scan.shape
    assert e_lo < e_hi

    # Peak should sit near p0**2/2 = 0.125 Ha (exact value up to the 0.01 Ha
    # discreteness of e_scan).
    peak_e = float(e_scan[int(np.argmax(eta_abs))])
    assert peak_e == pytest.approx(0.125, abs=0.015)

    # The window must bracket both of Task 4's validated anchors.
    assert e_lo <= 0.10
    assert e_hi >= 0.15
