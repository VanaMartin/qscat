"""Analytic local Chebyshev shape, electronic contour and no-extrapolation controls."""

import numpy as np
import pytest

from projects.ukrmol_co.local_gaussian_model import (
    local_gaussian_interaction,
    local_gaussian_parameters,
)


def test_known_quadratic_logs_and_two_coordinate_broadcast():
    coefficients = np.array([np.log(2.5), 0.1, -0.02, np.log(0.5), -0.03, 0.01])
    radius = np.array([1.9, 2.2, 2.5])
    x = (radius - 2.2) / 0.3
    lam, alpha = local_gaussian_parameters(coefficients, radius)
    expected_lam = 2.5 * np.exp(0.1 * x - 0.02 * (2 * x * x - 1))
    expected_alpha = 0.5 * np.exp(-0.03 * x + 0.01 * (2 * x * x - 1))
    np.testing.assert_allclose(lam, expected_lam, rtol=1e-14, atol=0)
    np.testing.assert_allclose(alpha, expected_alpha, rtol=1e-14, atol=0)
    r = np.array([0.0, 1.0, 2.0 + 0.3j])[:, None]
    actual = local_gaussian_interaction(coefficients, r, radius[None, :])
    np.testing.assert_allclose(
        actual, -expected_lam * np.exp(-expected_alpha * r * r), rtol=1e-14, atol=1e-14
    )


@pytest.mark.parametrize("radius", [1.89, 2.51, 2.2 + 0.1j, np.nan, True])
def test_rejects_extrapolation_complex_nuclear_R_and_malformed_inputs(radius):
    with pytest.raises(ValueError):
        local_gaussian_parameters(np.zeros(6), radius)


def test_no_implicit_neutral_zero_centrifugal_or_asymptote():
    c = np.array([np.log(2.0), 0.0, 0.0, np.log(0.4), 0.0, 0.0])
    np.testing.assert_allclose(local_gaussian_interaction(c, 0.0, 2.2), -2.0, rtol=0, atol=1e-14)
    for bad in [np.ones(5), np.full(6, np.nan), np.array([100.0, 0.0, 0.0, 100.0, 0.0, 0.0])]:
        with pytest.raises(ValueError):
            local_gaussian_parameters(bad, 2.2)
