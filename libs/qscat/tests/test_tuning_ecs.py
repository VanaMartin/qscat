from __future__ import annotations

import numpy as np
from qscat.tuning import max_stable_angle, tune_ecs_tail


def test_gaussian_angle_capped_below_45():
    # a Gaussian interaction exp(-alpha r^2) diverges on the rotated contour for theta>45
    V = lambda z: -np.exp(-0.4 * np.asarray(z) ** 2)
    ang = max_stable_angle(V, R0=12.0, tail_extent=40.0)
    assert ang <= 35.0 + 1e-9                    # never above the double-ECS cap
    # a bare -1/r (no Gaussian growth) is limited only by the cap
    ang2 = max_stable_angle(lambda z: -1.0 / np.asarray(z), R0=12.0, tail_extent=40.0)
    assert 30.0 <= ang2 <= 35.0 + 1e-9


def test_divergence_rejection_caps_gaussian_when_cap_is_high():
    # With a high angle_cap the scan reaches the Gaussian's ~45deg critical angle,
    # so the divergence-rejection branch MUST fire and cap below the cap -- this
    # actually exercises the logic (unlike the 35deg default, which returns the cap).
    gauss = lambda z: -np.exp(-0.4 * np.asarray(z) ** 2)
    ang = max_stable_angle(gauss, R0=12.0, tail_extent=40.0, angle_cap=60.0)
    assert ang < 60.0                      # divergence rejection fired (not just the cap)
    assert 40.0 <= ang <= 55.0             # near the ~45deg Gaussian critical angle
    # control: a non-diverging -1/r is limited only by the (high) cap
    ang2 = max_stable_angle(lambda z: -1.0 / np.asarray(z), R0=12.0, tail_extent=40.0, angle_cap=60.0)
    assert ang2 > 55.0                     # rises to ~the cap (no divergence)


def test_tail_absorbs_fast_wave():
    # K=58 (the F2 DA wave); at 35 deg the tail must reach ~1e-12 decay
    els = tune_ecs_tail(58.0, R0=10.7, angle=35.0, order=14, decay_target=1e-12)
    L = sum(els)
    assert np.exp(-58.0 * L * np.sin(np.deg2rad(35.0))) <= 1e-11   # decayed
    assert all(els[i] <= els[i + 1] + 1e-12 for i in range(len(els) - 1))  # exp-growth
