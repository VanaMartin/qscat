"""`FemDvrEcsGrid.real_index_near` — the shared nearest-real-region-DVR-index
primitive (was copy-pasted as `_index_near`/`_real_index_near`/`_nuclear_index_near`
across the runner and the validation TD drivers)."""

from __future__ import annotations

import numpy as np
from qscat.core.grids import nuclear_grid


def test_real_index_near_hits_the_nearest_real_node() -> None:
    g = nuclear_grid(r_max=20.0, quadrature=8, n_complex=4)
    real = g.real_points
    in_region = real <= g.R0
    # pick a genuine real-region node and ask for a value just beside it
    target_idx = int(np.flatnonzero(in_region)[len(np.flatnonzero(in_region)) // 2])
    target_val = float(real[target_idx])
    assert g.real_index_near(target_val + 1e-6) == target_idx


def test_real_index_near_never_returns_a_complex_tail_index() -> None:
    g = nuclear_grid(r_max=20.0, quadrature=8, n_complex=5)
    # ask for a value far out in the ECS tail; the answer must still be a
    # real-region index (tail points are masked out), i.e. the largest real node.
    idx = g.real_index_near(1e6)
    assert g.real_points[idx] <= g.R0
    real_region = g.real_points[g.real_points <= g.R0]
    assert np.isclose(g.real_points[idx], real_region.max())


def test_real_index_near_matches_the_legacy_masked_argmin() -> None:
    g = nuclear_grid(r_max=24.0, quadrature=10, n_complex=5)
    for r_value in (0.5, 2.0, 5.5, 9.9, 50.0):
        real = g.real_points
        masked = np.where(real <= g.R0, real, np.inf)
        expected = int(np.argmin(np.abs(masked - r_value)))
        assert g.real_index_near(r_value) == expected
