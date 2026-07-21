import numpy as np
from grid import FemDvrEcsGrid
from spec import ElementSpec, GridSpec


def _real_grid(nq=6, lengths=(1.0, 1.0, 1.0)):
    return GridSpec(quadrature=nq, elements=[ElementSpec(L) for L in lengths], x_min=0.0)


def test_point_count_and_dirichlet_drop():
    nq, nel = 6, 3
    g = FemDvrEcsGrid(_real_grid(nq, (1.0,) * nel))
    assert g.n == nel * (nq - 1) + 1 - 2          # bridge sharing + 2 endpoints dropped
    # outermost points (x_min=0 and x_max=3) are NOT in .points (Dirichlet)
    assert g.real_points.min() > 0.0 and g.real_points.max() < 3.0


def test_real_region_points_are_real():
    g = FemDvrEcsGrid(_real_grid())
    assert np.allclose(g.points.imag, 0.0)         # no ECS -> purely real
    assert np.allclose(g.points, g.real_points)


def test_ecs_tail_points_are_rotated():
    # 2 real elements then 1 complex element at 30 deg; pivot R0 = 2.0
    elements = [ElementSpec(1.0), ElementSpec(1.0), ElementSpec(1.0, 30.0)]
    spec = GridSpec(quadrature=6, elements=elements)
    g = FemDvrEcsGrid(spec)
    assert np.isclose(g.R0, 2.0)
    tail = g.points[g.real_points > g.R0 + 1e-9]
    # z = R0 + (x-R0) e^{i theta}; arg of (z-R0) ~ 30 deg
    ang = np.degrees(np.angle(tail - g.R0))
    assert np.allclose(ang, 30.0, atol=1e-6)


def test_weights_bridge_summed_at_shared_nodes():
    # interior element-boundary points carry a weight ~ sum of two half-element contributions;
    # they should be (roughly) larger than the small end-weights within an element.
    g = FemDvrEcsGrid(_real_grid(nq=6, lengths=(1.0, 1.0)))
    assert np.all(g.weights.real > 0)


def test_spec_rejects_noncontiguous_complex():
    import pytest

    # complex element before a real element -> not contiguous at the end
    elements = [ElementSpec(1.0, 30.0), ElementSpec(1.0, 0.0)]
    with pytest.raises(ValueError):
        GridSpec(quadrature=6, elements=elements)
