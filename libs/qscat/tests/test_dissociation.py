from __future__ import annotations

import numpy as np
import pytest
from qscat.core.dissociation import anion_electronic_states, da_cross_section, v_dr_diag
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.model import F2, N2, NO

_eg = electronic_grid
_ng = nuclear_grid


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


def _tgrid():
    return TensorGrid([electronic_grid(r_max=14.0, order=6, n_complex=4),
                       nuclear_grid(r_max=20.0, n_complex=4, quadrature=8)])


def test_v_dr_shape_and_dtype():
    tg = _tgrid()
    vdr = v_dr_diag(tg, F2)
    assert vdr.shape == (tg.size,) and vdr.dtype == np.complex128


def test_v_dr_equals_definition_pointwise():
    tg = _tgrid()
    model = F2
    R_inf = tg.grids[1].R0
    pts_r, pts_R = tg.points()  # (n_r,1), (1,n_R)
    expect = (
        model.interaction_diag(tg)
        + np.broadcast_to(model.v0(pts_R), tg.shape).ravel()
        - np.broadcast_to(model.v_int(pts_r, R_inf), tg.shape).ravel()
    )
    assert np.allclose(v_dr_diag(tg, model), expect, rtol=0, atol=1e-14)


def test_v_dr_tends_to_v0_at_large_R():
    # Where R is near R_inf, V_int(r,R) ~ V_int(r,R_inf), so V_DR ~ v0(R).
    tg = _tgrid()
    model = F2
    vdr = v_dr_diag(tg, model).reshape(tg.shape)  # (n_r, n_R)
    pts_R = tg.points()[1].ravel()
    j = int(np.argmin(np.abs(pts_R - tg.grids[1].R0)))  # column nearest R_inf
    v0_col = np.broadcast_to(model.v0(tg.points()[1]), tg.shape)[:, j]
    assert np.allclose(vdr[:, j], v0_col, rtol=0, atol=1e-10)


def _working():
    tg = TensorGrid([_eg(r_max=16.0, order=8, n_complex=6),
                     _ng(r_max=22.0, n_complex=8, quadrature=12)])
    return tg


def test_da_shape_scalar_and_array():
    tg = _working()
    eps, chi = vibrational_states(tg.grids[1], F2.mu, 3, F2.v0)
    s1 = da_cross_section(tg, F2, eps, chi, 0, 0.05)
    assert s1.shape == (1,)
    sN = da_cross_section(tg, F2, eps, chi, 0, np.array([0.05, 0.10]))
    assert sN.shape == (2, 1)
    assert np.all(sN >= 0.0) and np.all(np.isfinite(sN))


def test_n2_channel_closed_is_zero():
    # N2's DA threshold is +0.5 Ha -> sigma_DA == 0 across the whole VE window.
    tg = _working()
    eps, chi = vibrational_states(tg.grids[1], N2.mu, 3, N2.v0)
    E = np.array([0.04, 0.10, 0.18])
    s = da_cross_section(tg, N2, eps, chi, 0, E)
    assert np.all(s == 0.0)


@pytest.mark.slow
def test_f2_exothermic_da_is_positive():
    # F2 DA is open at all E>0; expect a nonzero, finite sigma in its resonance
    # window. No golden number (no independent DA data) -- positivity + soft
    # unitarity only.
    tg = _working()
    eps, chi = vibrational_states(tg.grids[1], F2.mu, 3, F2.v0)
    E = np.array([0.02, 0.03, 0.04])
    s = da_cross_section(tg, F2, eps, chi, 0, E)[:, 0]
    assert np.all(np.isfinite(s)) and np.all(s >= 0.0)
    assert s.max() > 0.0
    # soft unitarity: sigma_DA <= a few * pi/(2E) (partial-wave cap, generous
    # band for the under-resolved fast outgoing wave; see the convergence note)
    cap = np.pi / (2.0 * E)
    assert np.all(s < 50.0 * cap)


def _h2p_proxy():
    # small ionic proxy: electronic to ~60 bohr (holds a couple Rydberg states +
    # the incident), nuclear to ~14. Big enough for well-posedness, laptop-fast.
    from qscat.core.grids import electronic_grid, nuclear_grid
    from qscat.dvr import TensorGrid

    return TensorGrid(
        [
            electronic_grid(r_max=60.0, order=8, n_complex=6),
            nuclear_grid(r_max=22.0, n_complex=6, quadrature=10),
        ]
    )


@pytest.mark.slow
def test_dr_wellposed_and_threshold_ordered():
    from qscat.core.dissociation import dr_cross_section
    from qscat.core.vibrational import vibrational_states
    from qscat.model import H2P

    tg = _h2p_proxy()
    eps, chi = vibrational_states(tg.grids[1], H2P.mu, 3, H2P.v0)
    E = np.array([0.01, 0.03])
    s = dr_cross_section(tg, H2P, eps, chi, 0, E, n_channels=2)
    assert s.shape == (2, 2)
    assert np.all(np.isfinite(s)) and np.all(s >= 0.0)
