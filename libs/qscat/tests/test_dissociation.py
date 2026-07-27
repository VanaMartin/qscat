from __future__ import annotations

import numpy as np
import pytest
from qscat.core.dissociation import anion_electronic_states
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.vibrational import vibrational_states
from qscat.model import F2, N2, NO


def _eps0(model):
    g_R = nuclear_grid(r_max=22.0, n_complex=8, quadrature=12)
    eps, _ = vibrational_states(g_R, model.mu, 3, model.v0)
    return eps[0], g_R.R0


@pytest.mark.parametrize("model", [N2, NO, F2], ids=["N2", "NO", "F2"])
def test_one_bound_anion_state_real(model):
    g_r = electronic_grid(r_max=16.0, order=7, n_complex=6)
    _, R0 = _eps0(model)
    eps_e, phi = anion_electronic_states(g_r, model, R0, n_states=1)
    assert eps_e.shape == (1,) and phi.shape == (1, g_r.n)
    # c-product self-normalized over the real region ~ 1
    real = g_r.real_points <= g_r.R0
    p = phi[0].copy()
    p[~real] = 0.0
    assert abs(complex(p @ p) - 1.0) < 1e-6


def test_thresholds_have_correct_signs():
    # threshold(E_coll) = eps_e - eps[0]; F2 exothermic (<0), N2 closed (>0.3),
    # NO opens above its resonance (~0.17). No independent data -> sign/band gate.
    g_r = electronic_grid(r_max=16.0, order=7, n_complex=6)
    thr = {}
    for name, model in (("N2", N2), ("NO", NO), ("F2", F2)):
        eps0, R0 = _eps0(model)
        eps_e, _ = anion_electronic_states(g_r, model, R0, 1)
        thr[name] = float(eps_e[0]) - eps0
    assert thr["F2"] < 0.0            # exothermic: DA open at all E>0
    assert thr["N2"] > 0.3            # closed in the measurement window
    assert 0.10 < thr["NO"] < 0.25    # opens above the resonance


def test_raises_when_too_many_states_requested():
    g_r = electronic_grid(r_max=16.0, order=7, n_complex=6)
    _, R0 = _eps0(F2)
    with pytest.raises(ValueError):
        anion_electronic_states(g_r, F2, R0, n_states=50)
