"""Guard: qscat_run's N2 preset decks and the projects convergence grids must
stay identical. The two exist separately because layering forbids qscat_run
importing projects (see qscat_run/presets.py's docstring); this test is what
keeps the transcription from drifting -- the N2 sibling of
validation/diatomic/test_da_grid.py::test_diatomic_decks_match_presets.
Byte-identical (np.array_equal), not allclose: both sides build through the
same qscat.core.grids factories from the same literals."""

from __future__ import annotations

import numpy as np
from qscat.dvr import TensorGrid

from projects.n2_2d_cross_section.convergence import working_tgrid
from projects.n2_2d_td_cross_section.convergence import td_working_tgrid


def _assert_identical(a: TensorGrid, b: TensorGrid) -> None:
    for ga, gb in zip(a.grids, b.grids, strict=True):
        assert ga.n == gb.n
        assert ga.R0 == gb.R0
        assert np.array_equal(ga.points, gb.points)
        assert np.array_equal(ga.weights, gb.weights)
        assert np.array_equal(ga.real_points, gb.real_points)


def test_n2_ti_preset_matches_projects_working_grid():
    import qscat_run.presets as presets

    _assert_identical(presets._n2_ti_grid(), working_tgrid())


def test_n2_td_preset_matches_projects_td_working_grid():
    import qscat_run.presets as presets

    _assert_identical(presets._n2_td_grid(), td_working_tgrid())
