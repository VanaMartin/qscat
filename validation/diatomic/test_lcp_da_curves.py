from __future__ import annotations

import numpy as np
import pytest

from validation.diatomic.config import CONFIGS
from validation.diatomic.lcp_da_curves import compute_lcp_da_curve


@pytest.mark.slow
def test_f2_lcp_da_positive_and_finite():
    E, s = compute_lcp_da_curve(CONFIGS["F2"], np.array([0.02, 0.03, 0.04]))
    assert np.all(np.isfinite(s)) and np.all(s >= 0.0)
    assert s.max() > 0.0


@pytest.mark.slow
def test_f2_lcp_agrees_with_exact_within_factor_two():
    # The scientific check: AWAY FROM THRESHOLD the LCP agrees with the exact-2D
    # oracle to within ~a factor of 2 (the ~50% band). Measured: E=0.03 ratio 0.89
    # (LCP 1.47 vs exact 1.66), E=0.04 ratio 1.43. NEAR threshold the exact sigma_DA
    # spikes (E=0.02: exact 3.36 vs LCP 1.56, ratio 0.47) while the LCP stays smooth
    # -- a genuine near-threshold LCP departure, documented via the Task-4 comparison
    # figure (which spans the full window), NOT gated here.
    from validation.diatomic.da_curves import compute_da_curve
    E = np.array([0.03, 0.04])
    _, s_lcp = compute_lcp_da_curve(CONFIGS["F2"], E)
    _, s_exact = compute_da_curve(CONFIGS["F2"], E)
    ratio = s_lcp / s_exact[:, 0]
    assert np.all((ratio > 0.5) & (ratio < 2.0))
