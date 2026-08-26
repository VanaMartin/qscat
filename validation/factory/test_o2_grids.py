"""The O2 preset in `apps/qscat-run` carries the tuner's decks verbatim (to the
6 decimals they are written with): `o2_grids.o2_decks()` is the source, the
preset the copy, and this test is the lock -- the same layering as
`validation/diatomic/test_da_grid.py` for the eMoScat decks."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.model import O2
from qscat_run import presets

from validation.factory.o2_grids import NUCLEAR_R_CUT, deck, o2_decks, o2_grids, truncate_real


@pytest.fixture(scope="module")
def tuner_decks():
    return o2_decks()


def test_preset_decks_match_the_tuner(tuner_decks):
    g_e, g_n = tuner_decks
    p = presets.PRESETS["O2:tuner"]
    tg = p.ti_grid()
    for mine, theirs in zip(tuner_decks, tg.grids, strict=True):
        assert mine.n == theirs.n
        np.testing.assert_allclose(theirs.points, mine.points, rtol=0, atol=2e-5)
    assert g_e.n * g_n.n < 200_000
    # the carried nuclear deck is the truncated tuner mesh refined once
    n_cut = truncate_real(o2_grids()[1], NUCLEAR_R_CUT).n
    assert g_n.n > 1.8 * n_cut


def test_truncation_keeps_the_mesh_and_the_tail():
    _, g_n_full = o2_grids()
    g_n = truncate_real(g_n_full, NUCLEAR_R_CUT)
    d_full, d_cut = deck(g_n_full), deck(g_n)
    # every real segment kept is a segment of the full deck, in order
    assert d_cut["real_segments"] == d_full["real_segments"][: len(d_cut["real_segments"])]
    assert d_cut["real_segments"][-1][1] >= NUCLEAR_R_CUT
    # the tail's element LENGTHS are the tuner's, re-based at the cut
    ends_full = [x for _, x in d_full["complex_segments"]]
    ends_cut = [x for _, x in d_cut["complex_segments"]]
    shift = d_full["real_segments"][-1][1] - d_cut["real_segments"][-1][1]
    np.testing.assert_allclose(np.array(ends_full) - shift, ends_cut, atol=1e-9)
    assert d_cut["angle_deg"] == d_full["angle_deg"] and g_n.n < g_n_full.n


def test_o2_is_registered_for_ve_only():
    assert presets.MODELS["O2"] is O2
    assert presets.VALIDITY["O2"] == frozenset({"ve"})
