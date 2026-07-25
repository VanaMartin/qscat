"""Incident Gaussian electron wavepacket and the 2-D initial state g(r) x chi_v(R)."""

from __future__ import annotations

import numpy as np
from qscat.dvr import TensorGrid

from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_2d_cross_section.hamiltonian2d import MU
from projects.n2_2d_td_cross_section.wavepacket import gaussian_coeffs, initial_state
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states

# A box big enough to hold a wavepacket launched near r0=20 (test-scale, not production).
TG = TensorGrid(
    [
        n2_electronic_grid(r_max=40.0, order=8, n_complex=6),
        n2_nuclear_grid(quadrature=10, r_max=22.0, n_complex=5),
    ]
)
EPS, CHI = vibrational_states(TG.grids[1], MU, 4)


def test_gaussian_is_localized_near_r0() -> None:
    g = n2_electronic_grid(r_max=40.0, order=8, n_complex=6)
    r0 = 20.0
    coeffs = gaussian_coeffs(g, r0=r0, p0=-0.35, sigma=4.0)
    # density peak (|coeff|^2, which is the DVR density weight) sits near r0
    r_real = g.real_points
    peak_r = r_real[np.argmax(np.abs(coeffs[: r_real.size]) ** 2)]
    assert abs(peak_r - r0) < 3.0


def test_gaussian_carries_inward_momentum() -> None:
    """<p> < 0 for p0 < 0 -- the wavepacket moves toward the molecule."""
    g = n2_electronic_grid(r_max=40.0, order=8, n_complex=6)
    coeffs = gaussian_coeffs(g, r0=20.0, p0=-0.35, sigma=4.0)
    # phase increments negatively in r: consecutive real-region coeffs rotate clockwise
    real = g.real_points <= g.R0
    c = coeffs[: g.real_points.size][real[: g.real_points.size]]
    phases = np.unwrap(np.angle(c[np.abs(c) > 1e-6]))
    assert phases[-1] - phases[0] < 0.0


def test_gaussian_masked_to_unscaled_region() -> None:
    g = n2_electronic_grid(r_max=40.0, order=8, n_complex=6)
    coeffs = gaussian_coeffs(g, r0=20.0, p0=-0.35, sigma=4.0)
    tail = g.real_points > g.R0
    assert np.all(coeffs[tail] == 0.0)


def test_initial_state_is_unit_hermitian_norm_and_separable() -> None:
    psi = initial_state(TG, CHI[0], r0=20.0, p0=-0.35, sigma=4.0)
    assert psi.shape == (TG.size,)
    # Hermitian L2 norm == 1 (the physical probability norm), NOT the c-product
    assert abs(float(np.linalg.norm(psi)) - 1.0) < 1e-10
    # separable: reshape factorizes to outer(g_coeff, chi) up to scale
    block = psi.reshape(TG.shape)
    u, s, vh = np.linalg.svd(block)
    assert s[0] / s.sum() > 0.999  # essentially rank-1


def test_initial_state_masked() -> None:
    psi = initial_state(TG, CHI[0], r0=20.0, p0=-0.35, sigma=4.0)
    assert np.all(psi[~TG.real_mask()] == 0.0)
