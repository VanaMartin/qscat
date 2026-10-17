"""Independent polynomial, analytic-sheet and operator controls."""

import mpmath as mp
import numpy as np
import pytest

from projects.ukrmol_co.combined_exterior import (
    chebyshev_grid,
    exact_map,
    riccati_map,
    spectral_map,
)


def test_polynomial_derivative_and_endpoint_orientation():
    t, derivative = chebyshev_grid(12, 4)
    np.testing.assert_allclose(derivative @ t**3, 3 * t**2, atol=1e-10, rtol=0)
    np.testing.assert_allclose(t[[0, -1]], [4, 1], atol=1e-14, rtol=0)


@pytest.mark.parametrize("energy,dipole,separation", [(8, 0, 16), (4 - 0.2j, 0.4, 0)])
def test_independent_hankel_limits(energy, dipole, separation):
    with mp.workdps(60):
        expected = np.array(exact_map(energy, dipole, separation).tolist(), complex)
    actual, _ = riccati_map(energy, dipole, 40, 35, 1e-11, separation)
    spectral, _, _ = spectral_map(energy, dipole, 40, 35, 192, separation)
    np.testing.assert_allclose(actual, expected, atol=1e-7, rtol=0)
    np.testing.assert_allclose(spectral, expected, atol=1e-7, rtol=0)


def test_invalid_outgoing_contour_is_rejected():
    with pytest.raises(ValueError, match="Contour"):
        riccati_map(8 - 0.2j, 1.2, 40, 0, 1e-8)
