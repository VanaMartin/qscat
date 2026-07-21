"""Drift guard for the duplicated `n2_electronic_grid` factory.

`n2_electronic_grid` is hand-duplicated in two places: the toy-model
`projects/n2_resonance/grid_n2.py` (validated in
`projects/n2_resonance/test_grid_n2.py`) and this package's own
`resonance.py` (re-implemented in-line per `resonance.py`'s module
docstring, since `validation/` deliberately does not import from
`projects/`). Nothing stops the two copies from silently drifting apart on
a future edit to one but not the other -- this test cross-checks the grid
objects they build are identical, mirroring how
`projects/n2_resonance/test_potential.py` cross-checks `potential.py`
against `validation/n2/model.py`.
"""

from __future__ import annotations

import numpy as np

from projects.n2_resonance import grid_n2 as ref_grid_n2
from validation.n2 import resonance


def _assert_grids_identical(angle_deg: float) -> None:
    grid = resonance.n2_electronic_grid(angle_deg)
    ref_grid = ref_grid_n2.n2_electronic_grid(angle_deg)

    assert grid.n == ref_grid.n
    assert grid.R0 == ref_grid.R0
    np.testing.assert_array_equal(grid.points, ref_grid.points)
    np.testing.assert_array_equal(grid.weights, ref_grid.weights)


def test_grid_matches_reference_toy_model_at_35_deg():
    _assert_grids_identical(35.0)


def test_grid_matches_reference_toy_model_at_44_deg():
    _assert_grids_identical(44.0)
