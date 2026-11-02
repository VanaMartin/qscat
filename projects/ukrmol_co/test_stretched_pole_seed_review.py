"""Analytic real phase controls only; no CO read, ODE, root search or width claim."""

import numpy as np
import pytest

from projects.ukrmol_co.stretched_pole_seed_review import phase_seed


def test_single_level_phase_seed_has_correct_rydberg_hartree_full_width_factors():
    energy = np.linspace(0.1, 0.3, 2001)
    center = 0.20003
    gamma = 0.02
    k = (gamma / (center - energy))[:, None, None]
    result = phase_seed(energy, k)
    assert result["energy_rydberg"] == pytest.approx(center, abs=1e-4, rel=0)
    assert result["diagnostic_full_width_rydberg"] == pytest.approx(2 * gamma, abs=1e-7, rel=1e-4)
    assert result["diagnostic_seed_incident_hartree"][1] == pytest.approx(
        -gamma / 2, abs=1e-7, rel=1e-4
    )


def test_constant_phase_does_not_manufacture_an_interior_resonance_seed():
    with pytest.raises(ValueError, match="No interior rising phase feature"):
        phase_seed(np.linspace(0.1, 0.3, 9), np.zeros((9, 1, 1)))
