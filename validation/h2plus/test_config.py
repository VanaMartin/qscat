from __future__ import annotations

import numpy as np
import pytest
from qscat.dvr import TensorGrid

from validation.h2plus.config import N_CHANNELS, full_grid, proxy_grid


def test_full_grid_electronic_and_nuclear_pivots():
    g = full_grid()
    assert isinstance(g, TensorGrid)
    assert g.ndim == 2
    assert g.grids[0].R0 == pytest.approx(1300.0)  # electronic real region -> 1300
    assert g.grids[1].R0 == pytest.approx(14.0)  # nuclear real region -> 14
    assert np.iscomplexobj(g.grids[0].points)
    assert np.iscomplexobj(g.grids[1].points)


def test_proxy_grid_electronic_and_nuclear_pivots():
    g = proxy_grid()
    assert isinstance(g, TensorGrid)
    assert g.ndim == 2
    assert g.grids[0].R0 == pytest.approx(60.0)  # reduced electronic real region
    assert g.grids[1].R0 == pytest.approx(14.0)  # nuclear real region unchanged
    assert np.iscomplexobj(g.grids[0].points)
    assert np.iscomplexobj(g.grids[1].points)


def test_n_channels_is_three():
    assert N_CHANNELS == 3
