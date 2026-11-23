"""Analytic attenuation and exposure controls, independent of pole searches."""

import cmath

import numpy as np
import pytest
from qscat.tuning import tune_ecs_tail

from projects.potential_factory.resonance_tail import tune_resonance_tail


@pytest.mark.parametrize("angle", [30.0, 35.0])
def test_complex_tail_accounts_for_real_region_growth(angle):
    k = cmath.sqrt(2 * (0.065 - 0.041j))
    lengths = tune_resonance_tail(k, 24.0, angle=angle, order=11, decay_target=1e-6, base=3.5)
    endpoint = 24 + sum(lengths) * cmath.exp(1j * np.deg2rad(angle))
    assert abs(cmath.exp(1j * k * endpoint)) <= 1e-6 * (1 + 1e-12)
    # The last element is needed for this precisely stated amplitude budget.
    previous_endpoint = 24 + sum(lengths[:-1]) * cmath.exp(1j * np.deg2rad(angle))
    assert abs(cmath.exp(1j * k * previous_endpoint)) > 1e-6


def test_real_momentum_retains_existing_tail_contract():
    actual = tune_resonance_tail(0.4 + 0j, 24, angle=30, order=11, decay_target=1e-6, base=3.5)
    expected = tune_ecs_tail(0.4, 24, angle=30, order=11, decay_target=1e-6, base=3.5)
    np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=1e-14)


def test_unexposed_resonance_is_rejected():
    with pytest.raises(ValueError, match="expos"):
        tune_resonance_tail(0.4 - 0.2j, 24, angle=10, order=11, decay_target=1e-6)


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(momentum=0j),
        dict(momentum=0.4 + 0.1j),
        dict(R0=-1),
        dict(angle=0),
        dict(decay_target=1),
        dict(decay_target=0),
        dict(R0=1e6),
        dict(tail_alpha=-1),
        dict(base=float("nan")),
    ],
)
def test_invalid_or_unrepresentable_budgets_reject(kwargs):
    args = dict(momentum=0.4 - 0.1j, R0=24, angle=30, order=11, decay_target=1e-6) | kwargs
    with pytest.raises(ValueError):
        tune_resonance_tail(**args)
