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


@pytest.mark.parametrize("which", ["full", "proxy"])
def test_h2p_decks_match_presets(which: str) -> None:
    """Guard: the eMoScat H2+ decks here and `qscat_run.presets`' copies must
    stay byte-identical (the exact-2D DR curves now run through qscat-run, so the
    two copies -- kept separate because `qscat_run` must not import `validation`
    -- must not drift). Same invariant as the diatomic deck guard."""
    import qscat_run.presets as presets

    here = {"full": full_grid, "proxy": proxy_grid}[which]()
    there = {"full": presets._h2p_full_grid, "proxy": presets._h2p_proxy_grid}[which]()
    for axis in (0, 1):
        g1, g2 = here.grids[axis], there.grids[axis]
        assert g1.n == g2.n
        assert g1.R0 == g2.R0
        assert np.allclose(g1.points, g2.points)
        assert np.allclose(g1.weights, g2.weights)
