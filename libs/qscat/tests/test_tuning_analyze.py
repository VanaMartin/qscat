from __future__ import annotations

import numpy as np
from qscat.tuning import analyze_potential


def test_harmonic_turning_points_and_k():
    # V = 1/2 m w^2 x^2, m=1, w=1, E_max=2 -> turning points at x=+-2, k(0)=sqrt(2*2)=2
    m, w, E = 1.0, 1.0, 2.0
    V = lambda x: 0.5 * m * w**2 * np.asarray(x) ** 2
    p = analyze_potential(V, -5.0, 5.0, m, E)
    assert np.isclose(abs(p.turning_points).max(), 2.0, atol=1e-2)
    assert np.isclose(p.k[np.argmin(np.abs(p.x))], np.sqrt(2 * m * E), atol=1e-2)  # k at x=0
    # forbidden beyond the turning points: k=0, kappa>0
    assert p.k[np.argmax(p.x)] == 0.0 and p.kappa[np.argmax(p.x)] > 0.0


def test_detects_coulomb_singularity():
    V = lambda x: -1.0 / np.asarray(x)
    p = analyze_potential(V, 1e-3, 50.0, 1.0, 0.1)
    assert p.singularities.size >= 1 and p.singularities.min() < 0.05


def test_k_scales_with_mass_and_energy():
    V = lambda x: np.zeros_like(np.asarray(x, dtype=float))
    p = analyze_potential(V, 0.1, 10.0, 918.25, 0.05)  # heavy: large k
    assert np.allclose(p.k, np.sqrt(2 * 918.25 * 0.05), atol=1e-6)
