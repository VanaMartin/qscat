from __future__ import annotations

import math

import numpy as np
import pytest
from qscat.core.dissociation import anion_electronic_states
from qscat.core.grids import electronic_grid, nuclear_grid
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
    # F2's actual dissociative-attachment (DA) wavenumber at the top of this
    # energy range: K = sqrt(2*mu*(E_max - eps_e)), eps_e the anion bound
    # electronic state at the dissociation limit (docs/physics/
    # diatomic-ve-cross-sections.md quotes the rough "K_R ~ 58" at E~0.03;
    # this is the precise value at E_max=0.05, ~78 -- Task 8's
    # validation.tuning.calibrate calibrated the tuner's phase constant
    # against exactly this wave). The coarse shared N2-style grid famously
    # failed to resolve it (sigma_DA off by ~36 orders); the a-priori
    # proposed grid must not repeat that failure.
    e_max = 0.05
    elec = electronic_grid(r_max=16.0, order=8, n_complex=6)
    eps_e, _ = anion_electronic_states(elec, F2, R_inf=10.7, n_states=1)
    K = math.sqrt(2.0 * F2.mu * (e_max - float(eps_e[0])))

    g = propose_grid(F2, "nuclear", (0.01, e_max))
    result = probe_channel_representation(g, k=K, l=0, mass=F2.mu, rtol=1e-3)
    assert result.converged, result.detail
