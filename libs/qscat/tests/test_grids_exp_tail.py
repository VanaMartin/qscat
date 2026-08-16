from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import _ecs_tail, fem_grid_exp_tail


def test_real_endpoint_is_ecs_pivot():
    g = fem_grid_exp_tail([(2, 1.0), (3, 4.0)], angle_deg=10.0, quadrature=8, tail_n=5)
    assert g.R0 == pytest.approx(4.0)


def test_tail_element_lengths_grow_as_ecs_tail():
    # last real segment is (3, 4.0) starting at 1.0 -> element length 1.0
    base = 1.0
    alpha = 0.2
    skip = 2
    tail_n = 6
    g = fem_grid_exp_tail(
        [(2, 1.0), (3, 4.0)],
        angle_deg=10.0,
        quadrature=8,
        tail_n=tail_n,
        tail_alpha=alpha,
        tail_skip=skip,
    )
    expected = _ecs_tail(base, tail_n, skip=skip, alpha=alpha)
    assert expected[0] == pytest.approx(base)
    assert expected[skip - 1] == pytest.approx(base)
    assert expected[skip] > expected[skip - 1]
    # real region has 2 + 3 = 5 elements; tail elements follow, so the total
    # (unscaled) real coordinate span should reach the sum of real + tail
    # element lengths (within Dirichlet-endpoint-drop rounding).
    total_tail_span = sum(expected)
    assert g.real_points.max() == pytest.approx(4.0 + total_tail_span, rel=0.02)
    # monotonic growth beyond skip
    assert expected == sorted(expected)


def test_grid_has_complex_points_on_tail():
    g = fem_grid_exp_tail([(2, 1.0), (3, 4.0)], angle_deg=10.0, quadrature=8, tail_n=5)
    assert np.iscomplexobj(g.points)
    assert np.any(np.abs(g.points.imag) > 0)


def test_rejects_n_less_than_1():
    with pytest.raises(ValueError):
        fem_grid_exp_tail([(0, 1.0)], angle_deg=10.0, quadrature=8, tail_n=5)


def test_rejects_non_increasing_endpoints():
    with pytest.raises(ValueError):
        fem_grid_exp_tail([(2, 1.0), (1, 0.5)], angle_deg=10.0, quadrature=8, tail_n=5)


def test_rejects_low_quadrature():
    with pytest.raises(ValueError):
        fem_grid_exp_tail([(2, 1.0)], angle_deg=10.0, quadrature=1, tail_n=5)
