from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import nuclear_grid
from qscat.model import F2, N2
from qscat.tuning import probe_channel_representation, propose_grid


def test_propose_grid_nuclear_n2_is_a_valid_grid_in_a_sane_range():
    g = propose_grid(N2, "nuclear", (0.04, 0.18))
    assert np.all(np.isfinite(g.points))
    assert np.all(np.isfinite(g.weights))
    assert 10.0 <= g.R0 <= 22.0

    deck_n = nuclear_grid().n  # the committed N2 nuclear deck's point count
    assert g.n < 2 * deck_n


def test_propose_grid_rejects_empty_energy_range():
    with pytest.raises(ValueError):
        propose_grid(N2, "nuclear", (0.18, 0.04))


def test_propose_grid_electronic_n2_is_a_valid_grid():
    g = propose_grid(N2, "electronic", (0.04, 0.18))
    assert np.all(np.isfinite(g.points))
    assert g.R0 > 0.0
    assert g.n > 0


def test_propose_grid_nuclear_f2_resolves_the_da_channel_wave():
    # The known F2 dissociative-attachment wavenumber (docs/physics/
    # diatomic-ve-cross-sections.md: "F2: K_R ~ 58, wavelength ~0.107 bohr")
    # -- the coarse shared N2-style grid famously failed to resolve this
    # (sigma_DA off by ~36 orders); the a-priori proposed grid must not
    # repeat that failure.
    g = propose_grid(F2, "nuclear", (0.01, 0.05))
    result = probe_channel_representation(g, k=58.0, l=0, rtol=1e-3)
    assert result.converged, result.detail
