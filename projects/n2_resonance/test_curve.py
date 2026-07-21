"""V3 -- R-scan smoothness test for `pole.resonance_curve` (Task 3,
sub-project #2).

Traces E_res(R)/Gamma(R)/V_d(R) across R in [1.6, 3.0] Bohr with the
continuation-window walk in `resonance_curve` (seeded near R0, walking
outward). Physically this range spans the resonance-to-bound-state crossing
(the shape resonance at short R narrows and, past R ~ 2.3-2.4, becomes a
real, angle-independent bound state -- the usual dissociative-attachment
picture), so the assertions below check for a mode-hop-free, physically
smooth curve rather than any single closed-form oracle:

  - Gamma(R) >= 0 everywhere (Gamma is defined as max(0, -2*Im(E_pole)), so
    this also guards against that clip being silently bypassed).
  - E_res(R) is monotonically non-increasing across the scan (matches the
    dissociative curve dropping as the bond stretches; a mode-hop onto a
    neighboring branch -- which sit ~0.02-0.05 Ha apart at these grids, per
    `test_pole.py`'s eigenvalue-window choice -- would break monotonicity).
  - The second difference of E_res (curvature) stays small relative to the
    first differences (no discontinuity a few x the local step): a genuine
    mode-hop shows up as a second-difference spike comparable to a full
    branch spacing, far above the smooth curve's actual curvature.
"""

from __future__ import annotations

import numpy as np

from projects.n2_resonance import pole
from projects.n2_resonance.grid_n2 import n2_electronic_grid

R_GRID = np.linspace(1.6, 3.0, 15)


def test_V3_curve_is_smooth_and_gamma_nonnegative():
    ga, gb = n2_electronic_grid(35.0), n2_electronic_grid(44.0)
    E_res, Gamma, V_d = pole.resonance_curve(R_GRID, ga, gb)

    assert E_res.shape == R_GRID.shape
    assert Gamma.shape == R_GRID.shape
    assert V_d.shape == R_GRID.shape

    # Gamma >= 0 everywhere (allow only roundoff below zero).
    assert np.all(Gamma >= -1e-12), Gamma

    # E_res(R) monotonically non-increasing across the scan (dissociative
    # curve; a mode-hop onto a neighboring branch would break this).
    dE = np.diff(E_res)
    assert np.all(dE <= 1e-6), (R_GRID, E_res, dE)

    # No jump much bigger than the local step: curvature (second difference)
    # stays a small multiple of the typical step size, not a branch-spacing
    # sized spike.
    ddE = np.diff(dE)
    step_scale = np.median(np.abs(dE))
    assert np.max(np.abs(ddE)) < 5 * step_scale, (ddE, step_scale)

    # Gamma should also vary smoothly: over this R range the resonance
    # narrows monotonically to zero as it closes into a bound state (no
    # discontinuous jump back up, which a mode-hop onto a neighboring branch
    # would produce). Tolerance accommodates matching-residual-scale noise
    # right at the resonance-to-bound-state crossing, where Gamma is already
    # numerically indistinguishable from zero.
    assert np.all(np.diff(Gamma) <= 1e-6), Gamma

    # Physically expected regimes: a genuine resonance (Gamma > 0) at the
    # short-R end, and a closed (bound, Gamma == 0) state at the long-R end.
    assert Gamma[0] > 0.05, Gamma[0]
    assert Gamma[-1] < 1e-6, Gamma[-1]
