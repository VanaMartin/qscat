"""Tests for the eMoScat per-molecule nuclear DA grids (NO/F2) and σ_DA
well-posedness on them.

σ_DA does NOT converge on the shared N2-style nuclear grid (the K_R~58
dissociation wave is under-resolved there); eMoScat's per-molecule nuclear
decks (transcribed from `reference/eMoScat/input/{NO,F2}/grids.txt`) resolve
it. These tests confirm `MoleculeConfig.da_grid()` builds the fine nuclear
grid and that σ_DA on it is finite, non-negative, and softly unitary.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.dissociation import da_cross_section
from qscat.core.vibrational import vibrational_states

from validation.diatomic.config import CONFIGS


def test_da_grid_uses_emoscat_nuclear_resolution() -> None:
    # F2 nuclear real region ends at 10.7 (eMoScat deck), finely resolved.
    tg = CONFIGS["F2"].da_grid()
    assert tg.grids[1].R0 == pytest.approx(10.7)
    # the [2.7,10.7] region is tiled at ~0.2 bohr -> many nuclear elements
    assert tg.grids[1].n > 700  # ~960 for the F2 deck at quad 14


@pytest.mark.parametrize("name", ["F2", "NO"])
def test_diatomic_decks_match_presets(name: str) -> None:
    """Guard: the eMoScat deck here and `qscat_run.presets`' copy of it must
    stay byte-identical. The two exist separately because layering forbids
    `qscat_run` importing `validation` (and the tuner should not reach into the
    app's grid internals) -- this test is what keeps them from drifting."""
    import qscat_run.presets as presets

    here = CONFIGS[name].da_grid().grids[1]
    there = {"F2": presets._f2_nuc_grid, "NO": presets._no_nuc_grid}[name]()
    assert here.n == there.n
    assert here.R0 == there.R0
    assert np.allclose(here.points, there.points)
    assert np.allclose(here.weights, there.weights)
    assert np.allclose(here.real_points, there.real_points)


@pytest.mark.slow
def test_f2_sigma_da_wellposed_on_emoscat_grid() -> None:
    cfg = CONFIGS["F2"]
    tg = cfg.da_grid()
    eps, chi = vibrational_states(tg.grids[1], cfg.model.mu, cfg.n_vib, cfg.model.v0)
    E = np.array([0.03])
    s = da_cross_section(tg, cfg.model, eps, chi, 0, E)[:, 0]
    assert np.all(np.isfinite(s)) and np.all(s >= 0.0)
    assert s[0] > 0.0                      # exothermic -> open
    assert s[0] < 50.0 * np.pi / (2.0 * E[0])   # soft unitarity
