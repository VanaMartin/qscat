"""Independent angular, hard-sphere, flux and decoupled pole controls."""

import mpmath as mp
import numpy as np

from projects.ukrmol_co.coupled_tail import (
    exact_tail,
    exact_threshold,
    exact_threshold_pole,
    tail_scattering,
    threshold_scattering,
)
from projects.ukrmol_co.physical_pole import exact_pole


def test_dipole_angular_element():
    x, w = np.polynomial.legendre.leggauss(16)
    theta1 = -np.sqrt(3 / 4) * np.sqrt(1 - x * x)
    theta2 = -np.sqrt(15 / 4) * x * np.sqrt(1 - x * x)
    np.testing.assert_allclose(np.dot(w, theta1 * x * theta2), 1 / np.sqrt(5), atol=1e-14, rtol=0)


def test_uncoupled_hard_sphere_and_matrix_tail():
    q = 1.3
    # l=1 outgoing/incoming Riccati functions give the hard-sphere phase directly.
    expected_l1 = -np.exp(-2j * q) * (1 - 1j / q) / (1 + 1j / q)
    result, _ = tail_scattering(q, 0, 2, 1e-11)
    np.testing.assert_allclose(result[0, 0], expected_l1, atol=1e-9, rtol=0)
    result, _ = tail_scattering(q, 0.4, 2, 1e-11)
    with mp.workdps(60):
        expected = np.array(exact_tail(q, 0.4).tolist(), complex)
    np.testing.assert_allclose(result, expected, atol=1e-9, rtol=0)
    np.testing.assert_allclose(result.conj().T @ result, np.eye(2), atol=1e-9, rtol=0)


def test_threshold_flux_and_decoupled_lambert_pole():
    for energy in (8, 18):
        result, _ = threshold_scattering(energy, 1, 1.5, 1e-11)
        with mp.workdps(60):
            expected = np.array(exact_threshold(energy, 1).tolist(), complex)
        np.testing.assert_allclose(result, expected, atol=1e-8, rtol=0)
        np.testing.assert_allclose(result.conj().T @ result, np.eye(len(result)), atol=1e-8, rtol=0)
    with mp.workdps(60):
        result = exact_threshold_pole(0, [[2.8, -0.2], [2.9, -0.1]])
        assert abs(result - exact_pole(10, 1)) < mp.mpf("1e-45")
