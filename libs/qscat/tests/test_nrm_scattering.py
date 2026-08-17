"""Tests for the NRM's real-energy electronic scattering primitive."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import electronic_grid
from qscat.core.nrm.scattering import (
    free_hamiltonian,
    incident_coefficients,
    scattering_state,
)
from qscat.dvr import kinetic
from qscat.linalg import c_product
from qscat.special import riccati_bessel_en


@pytest.fixture(scope="module")
def grid():
    return electronic_grid(r_max=16.0, order=8, n_complex=6)


def test_free_hamiltonian_is_kinetic_plus_centrifugal(grid):
    """H_free = T_r + l(l+1)/2r^2, nothing else."""
    h_free = free_hamiltonian(grid, ell=1)
    cent = 1.0 * 2.0 / (2.0 * grid.points**2)
    expected = kinetic(grid, 1.0) + np.diag(cent)
    assert np.allclose(h_free, expected, rtol=1e-12, atol=0.0)


def test_incident_coefficients_are_function_values_times_sqrt_w(grid):
    """The incident wave is a FUNCTION and picks up sqrt(w); the ECS tail is masked."""
    k = 0.5
    c = incident_coefficients(grid, k=k, ell=1)
    real = grid.points.imag == 0.0
    vals = riccati_bessel_en(grid.points[real].real, k, 1)
    expected = vals * np.sqrt(grid.weights[real])
    assert np.allclose(c[real], expected, rtol=1e-12, atol=0.0)
    assert np.all(c[~real] == 0.0)


def test_free_potential_gives_back_the_incident_wave(grid):
    """With h == H_free the scattered part vanishes: phi+ == the incident wave.

    NOTE this does NOT check the sign convention of the driven equation: with
    h == h_free the right-hand side `(h - h_free) @ inc` is identically zero
    regardless of the sign in front of either `(E*I - h)` or the source term,
    so `phi_sc = 0` (and this test passes) for every sign combination. What it
    does check is that `free_hamiltonian` and `incident_coefficients` are
    mutually consistent -- `H_free`'s regular solution really is
    `riccati_bessel_en` on this grid. See
    `test_scattering_is_unitary_for_a_real_potential` for the test that
    actually discriminates the sign.
    """
    h_free = free_hamiltonian(grid, ell=1)
    k = 0.5
    phi = scattering_state(h_free, grid, energy=0.5 * k**2, ell=1)
    inc = incident_coefficients(grid, k=k, ell=1)
    real = grid.points.imag == 0.0
    scale = np.max(np.abs(inc[real]))
    assert np.max(np.abs((phi - inc)[real])) < 1e-8 * scale


def test_attractive_well_scatters(grid):
    """A real attractive Gaussian must produce a non-trivial scattered wave.

    NOTE this does NOT discriminate the sign convention either: it thresholds
    `|phi - inc|`, and `|phi_sc|` is identical whichever sign phi_sc carries
    (`|x| == |-x|`). See `test_scattering_is_unitary_for_a_real_potential`.
    """
    h_free = free_hamiltonian(grid, ell=1)
    well = -5.0 * np.exp(-3.0 * grid.points**2)
    h = h_free + np.diag(well)
    k = 0.5
    phi = scattering_state(h, grid, energy=0.5 * k**2, ell=1)
    inc = incident_coefficients(grid, k=k, ell=1)
    real = grid.points.imag == 0.0
    scale = np.max(np.abs(inc[real]))
    assert np.max(np.abs((phi - inc)[real])) > 1e-3 * scale


def test_scattering_is_unitary_for_a_real_potential(grid):
    """|S| == 1 for a REAL potential is a sign-sensitive check the two tests
    above cannot make.

    A bug that flips the sign of either side of the driven equation
    (`(E*I - h)` or the RHS alone) leaves `|phi_sc|` unchanged everywhere --
    it just negates `phi_sc`, and `|x| == |-x|`, so no magnitude-threshold
    test on `phi` can catch it (confirmed empirically: both single-sided
    flips reproduce `test_attractive_well_scatters`'s scattered-wave
    magnitude exactly). Unitarity does not have this blind spot: with the
    elastic T-matrix `T = <J_k | (h - H_free) | phi+>` (c-product, the 1-D
    analogue of `qscat.core.driven`'s formula) and `S = 1 - 2 pi i T` (same
    convention as `qscat.core.driven`'s module docstring), a real potential
    must give `|S| = 1` exactly -- flipping the sign of `T` breaks that
    identity unless `T == 0`. Measured on this grid/well/k: `||S| - 1|` is
    ~6e-11 for the correct sign and ~1.4e-4 for either single-sided flip, a
    4-order-of-magnitude gap that `atol=1e-6` sits cleanly inside.
    """
    h_free = free_hamiltonian(grid, ell=1)
    well = -5.0 * np.exp(-3.0 * grid.points**2)
    h = h_free + np.diag(well)
    k = 0.5
    phi = scattering_state(h, grid, energy=0.5 * k**2, ell=1)
    inc = incident_coefficients(grid, k=k, ell=1)
    t = c_product(inc, (h - h_free) @ phi)
    s = 1.0 - 2.0j * np.pi * t
    assert abs(abs(s) - 1.0) < 1e-6


def test_rejects_non_positive_energy(grid):
    h_free = free_hamiltonian(grid, ell=1)
    with pytest.raises(ValueError, match="positive"):
        scattering_state(h_free, grid, energy=0.0, ell=1)
