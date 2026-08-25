from __future__ import annotations

import numpy as np
import pytest
from qscat.core.dissociation import anion_electronic_states
from qscat.model import N2

from projects.potential_factory.extract import extract_target
from projects.potential_factory.target import CouplingTarget, Curve
from projects.potential_factory.tracker import ElectronicPair


def test_curve_from_table_interpolates_and_is_nan_outside():
    c = Curve.from_table(np.array([1.0, 2.0, 3.0]), np.array([1.0, 4.0, 9.0]))
    assert abs(c(2.5) - 6.25) < 0.05
    assert np.isnan(c(5.0))


def test_alt_houfek_coupling_has_the_threshold_exponent():
    ct = CouplingTarget.from_alt_houfek(
        a0=13.83669,
        a1=0.892095,
        a2=-0.935987,
        b0=3.015014,
        b1=0.718160,
        alpha=2.5,
        R_range=(1.8, 2.8),
    )
    e1, e2 = 1e-4, 2e-4
    slope = np.log(ct.gamma_tilde(e2, 2.0) / ct.gamma_tilde(e1, 2.0)) / np.log(e2 / e1)
    assert abs(slope - 2.5) < 1e-3


@pytest.mark.slow
def test_extract_target_from_n2_is_self_consistent():
    pair = ElectronicPair()
    R_desc = np.linspace(3.0, 1.6, 8)
    t = extract_target(N2, pair=pair, R_desc=R_desc, n_eps=4)
    assert t.ell == 2 and t.coordinates == ("R",)
    assert t.neutral is not None and t.resonance is not None and t.coupling is not None
    np.testing.assert_allclose(t.neutral.curve(R_desc), N2.v0(R_desc).real, atol=1e-12)
    eps_e, _ = anion_electronic_states(pair.grid_a, N2, 10.0, 1)
    assert abs(-t.resonance.ea - (eps_e[0] - N2.v0(10.0).real)) < 1e-10
    # A node with no gated pole is dropped, not frozen -- at least
    # 6 of the 8 R_desc nodes must survive, and Gamma(R) must show the real
    # near-equilibrium N2 shape-resonance width (~0.1 Ha at R=1.6 bohr), not
    # a frozen near-zero artifact.
    assert t.resonance.v_ion.x is not None and t.resonance.v_ion.x.size >= 6
    assert t.resonance.gamma(1.6) > 0.01
    # ResonanceTarget.R_range is the SURVIVING nodes' range, not
    # the full requested R_desc range.
    assert t.resonance.R_range == (
        float(t.resonance.v_ion.x.min()),
        float(t.resonance.v_ion.x.max()),
    )
    # threshold law of the extracted coupling: slope ~ l + 1/2 at small eps
    e1, e2 = 0.002, 0.004
    slope = np.log(t.coupling.gamma_tilde(e2, 1.9) / t.coupling.gamma_tilde(e1, 1.9)) / np.log(
        e2 / e1
    )
    assert abs(slope - 2.5) < 0.3
