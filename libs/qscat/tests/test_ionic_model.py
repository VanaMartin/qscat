"""Differential + protocol tests for the H2+ ionic model (`qscat.model.H2P`)."""

from __future__ import annotations

import numpy as np
from qscat.model import H2P, ResonanceModel


def test_h2p_is_a_resonance_model_with_charge():
    assert isinstance(H2P, ResonanceModel)
    assert H2P.charge == -1 and H2P.ell == 1
    # m_p/2, per Vana 2017 Table 1.2 and Hvizdos et al. PRA 97, 022704 (2018).
    assert H2P.mu == 918.076


def test_h2p_v0_is_the_ion_morse():
    R = np.array([2.0, 3.0, 8.0])
    V0, R0, a = 0.1027, 2.0, 0.69
    expect = V0 * (np.exp(-2 * a * (R - R0)) - 2 * np.exp(-a * (R - R0)))
    assert np.allclose(H2P.v0(R).real, expect, atol=1e-12)
    assert abs(H2P.v0(np.array([2.0]))[0].real + 0.1027) < 1e-12  # min -V0 at R0


def test_h2p_v_int_matches_sigma_capture():
    r, R = np.array([1.5]), np.array([2.5])
    a1, a2, a3, a4 = 1.6435, 6.2, 0.0125, 1.15
    Q = (a2 - R - a3 * R**4) / 7.0
    S = np.tanh(R / a4) ** 4
    E = np.exp(-r**2 / 3.0) / r
    expect = -a1 * (1 - np.tanh(Q)) * S * E
    assert np.allclose(H2P.v_int(r, R).real, expect, atol=1e-12)


def test_h2p_surface_has_coulomb_tail():
    # surface - (v0 + v_int + centrifugal) == charge/r == -1/r
    r, R = np.array([2.0]), np.array([3.0])
    cent = H2P.ell * (H2P.ell + 1) / (2.0 * r**2)
    tail = H2P.surface(r, R) - (H2P.v0(R) + H2P.v_int(r, R) + cent)
    assert np.allclose(tail, -1.0 / r, atol=1e-12)


def test_diatomic_models_are_neutral():
    from qscat.model import F2, N2

    assert N2.charge == 0 and F2.charge == 0
    assert isinstance(N2, ResonanceModel)
