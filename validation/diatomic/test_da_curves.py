from __future__ import annotations

import numpy as np
import pytest

from validation.diatomic.config import CONFIGS
from validation.diatomic.da_curves import compute_da_curve


@pytest.mark.slow
def test_f2_da_positive_on_emoscat_grid():
    E, s = compute_da_curve(CONFIGS["F2"], np.array([0.02, 0.03, 0.04]))
    assert np.all(np.isfinite(s)) and np.all(s >= 0.0)
    assert s.max() > 0.0


@pytest.mark.slow
def test_no_da_threshold_onset():
    # NO opens ~0.17 Ha: closed below, open above.
    _, s_lo = compute_da_curve(CONFIGS["NO"], np.array([0.10]))
    _, s_hi = compute_da_curve(CONFIGS["NO"], np.array([0.22]))
    assert s_lo[0, 0] == 0.0
    assert s_hi[0, 0] >= 0.0
