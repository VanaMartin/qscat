from __future__ import annotations

import numpy as np
import pytest
from qscat.core.dissociation import anion_electronic_states
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.lcp import lcp_da_cross_section, local_complex_potential
from qscat.core.vibrational import vibrational_states
from qscat.model import F2, N2


def _elec_grids():
    return (electronic_grid(r_max=16.0, order=7, n_complex=6, angle_deg=35.0),
            electronic_grid(r_max=16.0, order=7, n_complex=6, angle_deg=44.0))


def test_vd_gamma_shapes_and_gamma_nonneg():
    g_R = nuclear_grid(r_max=22.0, n_complex=6, quadrature=10)
    ga, gb = _elec_grids()
    Vd, Gamma = local_complex_potential(F2, g_R, ga, gb)
    assert Vd.shape == (g_R.n,) and Gamma.shape == (g_R.n,)
    assert Vd.dtype == np.complex128 and Gamma.dtype == np.float64
    assert np.all(Gamma >= 0.0)


def test_gamma_closes_and_vd_matches_anion_at_large_R():
    # As R -> R_inf the pole closes to the bound anion: Gamma -> ~0 and
    # V_d(R_inf) == eps_e (the exact-DA threshold from anion_electronic_states).
    g_R = nuclear_grid(r_max=22.0, n_complex=6, quadrature=10)
    ga, gb = _elec_grids()
    Vd, Gamma = local_complex_potential(F2, g_R, ga, gb)
    R = g_R.points
    real = R.imag == 0.0
    i_outer = np.flatnonzero(real)[np.argmax(R[real].real)]  # largest real R
    eps_e, _ = anion_electronic_states(ga, F2, g_R.R0, 1)
    assert Gamma[i_outer] < 1e-3                               # closed at the edge
    assert abs(Vd[i_outer].real - eps_e[0]) < 5e-3            # == anion asymptote


def test_gamma_positive_in_resonance_region():
    # At smaller R (inside the crossing) the anion is a real resonance: Gamma>0.
    g_R = nuclear_grid(r_max=22.0, n_complex=6, quadrature=10)
    ga, gb = _elec_grids()
    Vd, Gamma = local_complex_potential(F2, g_R, ga, gb)
    R = g_R.points.real
    band = (R > 1.5) & (R < 2.5)
    assert Gamma[band].max() > 1e-4                           # genuine width somewhere


def test_matches_n2_vres_oracle():
    from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
    from projects.n2_ti_cross_section.vres import vres_on_grid
    g_R = n2_nuclear_grid()
    ga, gb = _elec_grids()
    Vd, Gamma = local_complex_potential(N2, g_R, ga, gb)
    Vd_ref, Gamma_ref = vres_on_grid(g_R)
    real = g_R.points.imag == 0.0
    # compare on the resonance region where both are well-defined (R in [1.5,3.5])
    R = g_R.points.real
    band = real & (R > 1.5) & (R < 3.5)
    assert np.allclose(Vd[band].real, Vd_ref[band].real, atol=5e-3)
    assert np.allclose(Gamma[band], Gamma_ref[band], atol=5e-3)


def _f2_lcp_inputs():
    g_R = nuclear_grid(r_max=22.0, n_complex=6, quadrature=10)
    ga, gb = _elec_grids()
    Vd, Gamma = local_complex_potential(F2, g_R, ga, gb)
    eps, chi = vibrational_states(g_R, F2.mu, 3, F2.v0)
    real = g_R.points.imag == 0.0
    eps_e = float(Vd[np.flatnonzero(real)][np.argmax(g_R.points.real[real])].real)
    return g_R, Vd, Gamma, eps, chi, eps_e


def test_lcp_da_shape_and_nonneg():
    g_R, Vd, Gamma, eps, chi, eps_e = _f2_lcp_inputs()
    E = np.array([0.02, 0.03, 0.04])
    s = lcp_da_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, eps_e, E, dt=1.0, n_steps=400)
    assert s.shape == (3,)
    assert np.all(np.isfinite(s)) and np.all(s >= 0.0)


def test_lcp_da_closed_channel_is_zero():
    # F2 is exothermic (no closed region for E>0), so test the closed guard
    # directly: a threshold eps_e above E_tot forces E_DR<0 -> sigma must be 0.
    g_R, Vd, Gamma, eps, chi, _ = _f2_lcp_inputs()
    eps_e_high = float(eps[0]) + 1.0
    s = lcp_da_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, eps_e_high,
                             np.array([0.03]), dt=1.0, n_steps=400)
    assert s[0] == 0.0


@pytest.mark.slow
def test_lcp_da_f2_positive():
    g_R, Vd, Gamma, eps, chi, eps_e = _f2_lcp_inputs()
    E = np.array([0.02, 0.03, 0.04])
    s = lcp_da_cross_section(g_R, F2.mu, Vd, Gamma, eps, chi, 0, eps_e, E, dt=1.0, n_steps=2000)
    assert s.max() > 0.0                     # exothermic F2: DA is open and nonzero
