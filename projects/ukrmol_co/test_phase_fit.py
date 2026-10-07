"""Analytic unitary S factors provide the independent phase-fit oracle."""

import numpy as np
import pytest

from projects.ukrmol_co.phase_fit import fit_phase, predict_phase


@pytest.mark.parametrize("points", [31, 61, 121])
@pytest.mark.parametrize("nback", [1, 2, 3, 4])
def test_analytic_unitary_s_factor(points, nback):
    energy = np.linspace(0.055, 0.135, points)
    position, width = 0.094, 0.012
    coefficients = np.array([0.1, -0.03, 0.006, -0.002])[:nback]
    x = (energy - 0.095) / 0.04
    background = sum(c * x**i for i, c in enumerate(coefficients))
    s = np.exp(2j * background) * (
        (position - energy + 0.5j * width) / (position - energy - 0.5j * width)
    )
    np.testing.assert_allclose(np.abs(s), 1, atol=3e-15, rtol=0)
    phase = np.angle(s) / 2
    result = fit_phase(energy, phase, nback, 0.09, 0.018)
    assert result["optimizer_success"] and not result["at_parameter_bound"]
    np.testing.assert_allclose(
        [result["position_hartree"], result["full_width_hartree"]],
        [position, width],
        atol=1e-9,
        rtol=0,
    )
    difference = np.remainder(predict_phase(energy, result) - phase + np.pi / 2, np.pi) - np.pi / 2
    np.testing.assert_allclose(difference, 0, atol=1e-9, rtol=0)
    branch = fit_phase(energy, phase + 3 * np.pi, nback, 0.09, 0.018)
    np.testing.assert_allclose(
        [branch["position_hartree"], branch["full_width_hartree"]],
        [position, width],
        atol=1e-9,
        rtol=0,
    )


def test_insufficient_background_leaves_resolved_prediction_error():
    energy = np.linspace(0.055, 0.135, 61)
    phase = (
        np.unwrap(np.angle((0.094 - energy + 0.006j) / (0.094 - energy - 0.006j))) / 2
        - 0.25 * (energy - 0.095) / 0.04
    )
    held = np.arange(len(energy)) % 5 == 0
    incomplete = fit_phase(energy[~held], phase[~held], 1, 0.09, 0.018)
    complete = fit_phase(energy[~held], phase[~held], 2, 0.09, 0.018)
    assert np.sqrt(np.mean((predict_phase(energy[held], incomplete) - phase[held]) ** 2)) > 0.02
    np.testing.assert_allclose(
        predict_phase(energy[held], complete), phase[held], atol=1e-9, rtol=0
    )


def test_rejects_duplicate_energy_points():
    with pytest.raises(ValueError, match="increasing grid"):
        fit_phase(np.array([0.05, 0.06, 0.06, 0.08]), np.zeros(4), 1, 0.065, 0.01)
