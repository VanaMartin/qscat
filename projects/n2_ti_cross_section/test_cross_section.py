"""Failing-first tests for the VE cross section via the resolvent/driven
equation (sub-project #3, Task 3 -- THE CRUX).

the eMoScat TI extraction sections 3-4: doorway
`d_v(R) = sqrt(Gamma(R)/(2*pi)) * chi_v(R)`; driven equation
`(E_tot - H_res) xi = d_{v_init}` with
`H_res = T_nuc(mu) + diag(V_d(R) - i*Gamma(R)/2)`; S-matrix
`S_{v'<-v_init} = <d_{v'} | xi>` via the DVR c-product (plain dot, no
conjugate -- the basis is already 1/sqrt(w)-normalized); cross section
`sigma = 4*pi**3*|S|**2/(2*E)`, zero if the final channel is energetically
closed (`E_tot - eps_{v'} <= 0`).

Two families of checks:

- INTERNAL (model-independent, non-negotiable): sigma real & >=0; a closed
  channel gives exactly 0; the v=0->1 cross section is resonance-enhanced
  in the ~2-3 eV region relative to near threshold.

the Houfek anchor comparison lives in `validation/n2/test_anchor_gate.py`
(validation may import projects; not the reverse).
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.model import N2

from projects.n2_ti_cross_section.cross_section import ve_cross_section
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states
from projects.n2_ti_cross_section.vres import vres_on_grid

MU = N2.mu  # N2 nuclear reduced mass (a.u.), 12766.36
N_VIB = 6  # v=0..5, enough to cover vprimes up to 3 used by the anchors


@pytest.fixture(scope="module")
def system():
    """Build the shared grid / vibrational states / V_d(R),Gamma(R) once.

    `vres_on_grid` walks a two-angle-matched electronic pole-finder
    continuation across ~300 nuclear grid points (~7s) -- module-scoped so
    the whole test file pays that cost exactly once.
    """
    grid = n2_nuclear_grid()
    eps, chi = vibrational_states(grid, MU, N_VIB)
    Vd, Gamma = vres_on_grid(grid)
    return grid, eps, chi, Vd, Gamma


def test_sigma_real_and_nonnegative(system):
    grid, eps, chi, Vd, Gamma = system
    # A handful of (E, v') spanning near-threshold, mid-range, resonance.
    for E in (0.02, 0.05, 0.1, 0.15, 0.2):
        sigma = ve_cross_section(grid, MU, Vd, Gamma, eps, chi, 0, [0, 1, 2, 3], E)
        assert sigma.shape == (4,)
        assert np.all(np.isreal(sigma)) or np.all(np.abs(sigma.imag) < 1e-10)
        assert np.all(sigma.real >= -1e-12)


def test_closed_channel_is_exactly_zero(system):
    grid, eps, chi, Vd, Gamma = system
    # eps[3]-eps[0] is several vibrational quanta (~0.037 Ha); a tiny
    # collision energy E=0.001 Ha leaves E_tot far below eps[3] -> closed.
    E = 0.001
    E_tot = E + eps[0]
    assert E_tot - eps[3] < 0  # sanity: channel is indeed closed
    sigma = ve_cross_section(grid, MU, Vd, Gamma, eps, chi, 0, [3], E)
    assert sigma[0] == 0.0


def test_v0_to_v1_resonance_enhancement(system):
    """sigma_{0->1}(E) should be larger in the ~2-3 eV (Pi_g resonance)
    region than near the v'=1 threshold (~0.0124 Ha above v=0)."""
    grid, eps, chi, Vd, Gamma = system

    e_near_threshold_ha = 0.02  # ~0.54 eV collision energy, just above the
    # v'=1 threshold (eps1-eps0 ~ 0.0124 Ha)
    e_resonance_ha = 0.1  # ~2.72 eV collision energy, inside the 2-3 eV
    # Pi_g shape-resonance window (literature ~2.3-2.5 eV)

    sigma_near = ve_cross_section(grid, MU, Vd, Gamma, eps, chi, 0, [1], e_near_threshold_ha)[0]
    sigma_res = ve_cross_section(grid, MU, Vd, Gamma, eps, chi, 0, [1], e_resonance_ha)[0]

    print(f"\nsigma_0->1(E={e_near_threshold_ha} Ha) = {sigma_near:.6e} bohr^2")
    print(f"sigma_0->1(E={e_resonance_ha} Ha) = {sigma_res:.6e} bohr^2")

    assert sigma_res > sigma_near
