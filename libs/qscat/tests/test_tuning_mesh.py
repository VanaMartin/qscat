from __future__ import annotations

import numpy as np
from qscat.tuning import analyze_potential, equidistribution_elements, optimal_real_mesh


def _flat_profile(k_const, x_max, m):
    E = k_const**2 / (2 * m)
    return analyze_potential(lambda x: np.zeros_like(np.asarray(x, float)), 0.0, x_max, m, E)


def test_uniform_k_gives_uniform_elements():
    # constant k -> equidistribution => (nearly) uniform element lengths
    p = _flat_profile(2.0, 20.0, 1.0)
    els = equidistribution_elements(p, 8, phase_per_element=5.0, min_len=1e-3, max_len=10.0)
    assert np.std(els) / np.mean(els) < 0.05                 # ~uniform
    assert abs(2.0 * np.mean(els) - 5.0) < 0.3               # k*h ~ phase_per_element


def test_finer_where_k_larger():
    # step: k large for x<5, small for x>5 -> smaller elements in the fast region
    def V(x):
        x = np.asarray(x, float); return np.where(x < 5.0, -2.0, 0.0)
    p = analyze_potential(V, 0.0, 20.0, 1.0, 0.5)
    els = equidistribution_elements(p, 8, phase_per_element=5.0, min_len=1e-3, max_len=10.0)
    # reconstruct boundaries; mean element length in [0,5] < mean in [5,20]
    b = np.concatenate([[0.0], np.cumsum(els)])
    lo = np.diff(b)[b[:-1] < 5.0].mean(); hi = np.diff(b)[b[:-1] >= 5.0].mean()
    assert lo < hi


def test_optimal_mesh_minimizes_points():
    p = _flat_profile(2.0, 20.0, 1.0)
    els, order = optimal_real_mesh(p, phase_coeff=1.0, min_len=1e-3, max_len=10.0)
    assert order in (6, 8, 10, 14)
    assert len(els) * (order - 1) <= 400                     # a sane point budget for this easy case
