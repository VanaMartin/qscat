"""Internal checks for the LCP VE cross section, re-homed onto validation.

the model-independent internal checks that shipped with the projects toy
model, re-homed onto the graduated solver: sigma real & >=0; a closed
channel gives exactly 0; the v=0->1 cross section is resonance-enhanced in
the ~2-3 eV region relative to near threshold.

the Houfek anchor comparison lives in `validation/n2/test_anchor_gate.py`.
"""

from __future__ import annotations

import numpy as np
from qscat.core.lcp import lcp_ve_cross_section
from qscat.model import N2

from validation.n2.cross_section import build_system


def test_sigma_real_and_nonnegative():
    grid, eps, chi, Vd, Gamma = build_system()
    for E in (0.02, 0.05, 0.1, 0.15, 0.2):
        sigma = lcp_ve_cross_section(grid, N2.mu, Vd, Gamma, eps, chi, 0, [0, 1, 2, 3], E)
        assert sigma.shape == (4,)
        assert np.all(np.isfinite(sigma)) and np.all(sigma >= 0.0)


def test_closed_channel_is_exactly_zero():
    grid, eps, chi, Vd, Gamma = build_system()
    E = 0.001
    assert E + eps[0] - eps[3] < 0  # sanity: v'=3 is closed at this E
    assert lcp_ve_cross_section(grid, N2.mu, Vd, Gamma, eps, chi, 0, [3], E)[0] == 0.0


def test_v0_to_v1_resonance_enhancement():
    grid, eps, chi, Vd, Gamma = build_system()
    near = lcp_ve_cross_section(grid, N2.mu, Vd, Gamma, eps, chi, 0, [1], 0.02)[0]
    res = lcp_ve_cross_section(grid, N2.mu, Vd, Gamma, eps, chi, 0, [1], 0.1)[0]
    assert res > near  # the ~2-3 eV Pi_g resonance enhances 0->1
