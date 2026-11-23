"""Independent analytic controls for the synthetic radial boundary comparator."""

import numpy as np
import pytest
from scipy.special import spherical_jn, spherical_yn

from projects.potential_factory import radial_oracle as oracle


@pytest.mark.parametrize("energy", [0.07 - 0.04j, 0.12 - 0.02j])
def test_free_regular_solution_matches_spherical_bessel(energy):
    radius = 10.0
    k = np.sqrt(2 * energy)
    z = k * radius
    expected = (
        3 * spherical_jn(1, z) / z,
        3 * k * (spherical_jn(1, z, derivative=True) / z - spherical_jn(1, z) / z**2),
    )
    actual = oracle.regular_p_wave(energy, 0.0, 0.5, radius=radius, rtol=1e-11)
    np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=1e-10)


@pytest.mark.parametrize("energy", [0.07 - 0.04j, 0.12 - 0.02j])
def test_outgoing_boundary_matches_bessel_value_and_derivative(energy):
    radius = 10.0
    k = np.sqrt(2 * energy)
    z = k * radius
    h = spherical_jn(1, z) + 1j * spherical_yn(1, z)
    dh = spherical_jn(1, z, derivative=True) + 1j * spherical_yn(1, z, derivative=True)
    expected = (-k * radius * h, -k * (h + z * dh))
    actual = oracle.p_wave_boundary(energy, radius)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)


def test_finite_endpoint_boundary_vanishes_and_has_nonzero_derivative():
    energy = 0.07 - 0.04j
    endpoint = 55 + 20j
    value, derivative = oracle.p_wave_boundary(energy, endpoint, endpoint=endpoint)
    np.testing.assert_allclose(value, 0.0, rtol=0, atol=1e-12)
    np.testing.assert_allclose(derivative, 2j * np.sqrt(2 * energy), rtol=1e-12, atol=1e-12)


def test_free_regular_outgoing_wronskian_is_not_a_resonance():
    mismatch = oracle.gaussian_p_wave_mismatch(0.07 - 0.04j, 0.0, 0.5, radius=10, rtol=1e-11)
    assert abs(mismatch) > 1e-3


@pytest.mark.parametrize("kwargs", [dict(lam=-1), dict(alpha=0), dict(radius=0), dict(rtol=0)])
def test_regular_solution_rejects_invalid_inputs(kwargs):
    args = dict(lam=2.5, alpha=0.5, radius=10.0, rtol=1e-10) | kwargs
    with pytest.raises(ValueError):
        oracle.regular_p_wave(0.07 - 0.04j, **args)
