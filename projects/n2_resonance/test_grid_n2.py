"""Failing-first tests for the N2 electronic ECS grid factory (Task 1,
sub-project #2). Builds on the validated `qscat.dvr.FemDvrEcsGrid`.
"""

import numpy as np

from projects.n2_resonance.grid_n2 import n2_electronic_grid


def test_pivot_matches_r_pivot():
    g = n2_electronic_grid(35.0)
    assert np.isclose(g.R0, 10.0)


def test_real_region_points_are_real_and_below_pivot():
    g = n2_electronic_grid(35.0)
    real_mask = g.real_points < g.R0 - 1e-9
    assert real_mask.any()
    assert np.allclose(g.points[real_mask].imag, 0.0)


def test_tail_points_are_rotated_by_angle():
    angle = 35.0
    g = n2_electronic_grid(angle)
    tail_mask = g.real_points > g.R0 + 1e-9
    assert tail_mask.any()
    tail = g.points[tail_mask]
    ang = np.degrees(np.angle(tail - g.R0))
    assert np.allclose(ang, angle, atol=1e-6)


def test_n_positive():
    g = n2_electronic_grid(35.0)
    assert g.n > 0


def test_two_angles_share_identical_real_points():
    g35 = n2_electronic_grid(35.0)
    g44 = n2_electronic_grid(44.0)
    assert g35.real_points.shape == g44.real_points.shape
    assert np.array_equal(g35.real_points, g44.real_points)
    # but the tail complex points genuinely differ (angle actually applied)
    tail35 = g35.points[g35.real_points > g35.R0 + 1e-9]
    tail44 = g44.points[g44.real_points > g44.R0 + 1e-9]
    assert not np.allclose(tail35, tail44)


def test_default_element_counts():
    g = n2_electronic_grid(35.0)
    assert len(g.spec.elements) == 8 + 8
    real_elements = g.spec.elements[:8]
    complex_elements = g.spec.elements[8:]
    assert all(el.angle_deg == 0.0 for el in real_elements)
    assert all(el.angle_deg == 35.0 for el in complex_elements)
    assert np.allclose([el.length for el in real_elements], 10.0 / 8)
    assert np.allclose([el.length for el in complex_elements], (30.0 - 10.0) / 8)
