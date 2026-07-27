from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import segmented_grid


def test_reproduces_emoscat_n2_nuclear_deck():
    # N2 nuclear (input/N2/grids.txt, 2nd declaration): real to 12.0, tail to 55.
    g = segmented_grid(
        [(2, 1.0), (1, 1.5), (10, 3.0), (2, 4.0), (2, 6.0), (6, 12.0)],
        [(1, 13.0), (2, 16.0), (1, 18.0), (4, 55.0)],
        angle_deg=35.0,
        quadrature=14,
    )
    assert g.R0 == pytest.approx(12.0)            # ECS pivot = last real endpoint
    # real region ends at 12, tail runs onto the complex plane past it
    # (Dirichlet drop removes the endpoint, so real_points.max() < 55.0)
    assert float(g.real_points.max()) == pytest.approx(55.0, rel=0.01)
    assert np.iscomplexobj(g.points) and np.any(np.abs(g.points.imag) > 0)


def test_element_lengths_are_uniform_per_segment():
    # the 10 elements over [1.5, 3.0] are each 0.15 bohr
    g = segmented_grid([(1, 1.5), (10, 3.0)], [], angle_deg=35.0, quadrature=8)
    assert g.R0 == pytest.approx(3.0)             # no complex tail -> pivot at real end
    assert np.max(np.abs(g.points.imag)) == pytest.approx(0.0)  # pure real


@pytest.mark.parametrize(
    "real_seg",
    [[(0, 1.0)], [(2, 1.0), (1, 0.5)]],  # n<1 ; non-increasing endpoint
)
def test_rejects_bad_segments(real_seg):
    with pytest.raises(ValueError):
        segmented_grid(real_seg, [], angle_deg=35.0, quadrature=8)
