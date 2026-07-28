from __future__ import annotations

import numpy as np
from qscat.model import F2, N2
from qscat.tuning import interaction_region


def test_f2_interaction_region_brackets_the_coupling() -> None:
    R_lo, R_hi = interaction_region(F2)
    assert 0.5 < R_lo < 2.5 and 3.0 < R_hi < 8.0  # ~[1.5,4]-ish, where lambda(R) is significant
    assert R_lo < R_hi


def test_region_is_where_vint_is_still_transitioning() -> None:
    # N2's lambda(R) saturates to a substantial (non-zero) plateau on BOTH
    # sides -- verified numerically, see resonance.py's module docstring --
    # so "outside [R_lo, R_hi] is small" (the F2-style check) does not hold
    # here. What DOES hold: past R_hi, s(R) has essentially stopped
    # changing (it's within a small fraction of its peak of its own
    # asymptote) -- the property the downstream sampler actually relies on
    # to freeze at a single far value.
    R_lo, R_hi = interaction_region(N2, frac=0.05)
    r = np.linspace(0.1, 15, 200)
    s = lambda R: np.max(np.abs(np.real(N2.v_int(r[:, None], np.array([[R]])))))
    peak = max(s(R) for R in np.linspace(0.5, 6, 60))
    assert abs(s(R_hi + 1.0) - s(R_hi)) < 0.05 * peak * 1.5  # already saturated
    assert R_lo < R_hi
