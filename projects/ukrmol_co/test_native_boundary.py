"""Free-box spectral normalization and complex symmetric resolvent controls."""

import numpy as np
import pytest

from projects.ukrmol_co.native_boundary import surface_matrix


def test_free_neumann_box_normalization_and_omitted_tail():
    radius, q = 5.0, 1.3
    exact = radius * np.tan(q) / q
    errors = []
    for count in [128, 512, 2048]:
        n = np.arange(count)
        poles = ((n + 0.5) * np.pi / radius) ** 2 / 2
        amplitude = ((-1.0) ** n / np.sqrt(radius))[:, None]
        value = surface_matrix(poles, amplitude, q * q / (2 * radius**2))[0, 0]
        error = exact - value
        assert error > 0
        np.testing.assert_allclose(error * count, 2 * radius / np.pi**2, atol=1e-5, rtol=0)
        errors.append(error)
        for wrong in [value / 2, 2 * value, value / radius]:
            assert abs(wrong - exact) > 1
    assert errors[2] < errors[1] / 3.9 < errors[0] / 3.9**2


def test_complex_total_energy_transpose_and_pole_rejection():
    poles = np.array([-2.0, 0.5, 3.0])
    amplitudes = np.array([[1.0, 0.2], [-0.4, 0.9], [0.7, -0.6]])
    energy = 0.1 - 0.2j
    actual = surface_matrix(poles, amplitudes, energy)
    expected = sum(
        np.outer(w, w) / (pole - energy) for pole, w in zip(poles, amplitudes, strict=True)
    )
    np.testing.assert_allclose(actual, expected, atol=1e-14, rtol=0)
    np.testing.assert_allclose(actual, actual.T, atol=1e-14, rtol=0)
    np.testing.assert_allclose(
        surface_matrix(poles, amplitudes, energy.conjugate()),
        actual.conjugate(),
        atol=1e-14,
        rtol=0,
    )
    with pytest.raises(ValueError, match="pole"):
        surface_matrix(poles, amplitudes, poles[1])
