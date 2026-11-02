"""Analytic rate signs and saved-block localization; no exterior ODE or CO map."""

import numpy as np
import pytest

from projects.ukrmol_co.diagnose_stretched_real_guard import rate_separation, saved_blocks


def test_exact_two_momentum_critical_angle_separates_decay_from_amplification():
    # E=1 and T=2 give k=kappa=1, hence the exact crossover is45degrees.
    critical = rate_separation(1, 2, 45, 10)
    np.testing.assert_allclose(critical["relative_decay_rate_per_bohr"], 0, atol=1e-14, rtol=0)
    np.testing.assert_allclose(critical["finite_interval_damping_estimate"], 1, atol=1e-13, rtol=0)
    assert rate_separation(1, 2, 30, 10)["finite_interval_damping_estimate"] < 1
    assert rate_separation(1, 2, 60, 10)["finite_interval_damping_estimate"] > 1


def test_forbidden_closed_or_threshold_energies_are_not_assigned_a_decay_rate():
    for energy, threshold in [(0, 2), (2, 2), (3, 2)]:
        with pytest.raises(ValueError):
            rate_separation(energy, threshold, 40, 102)


def test_known_open_closed_flux_error_is_localized_without_a_bad_matching_solve():
    q = np.zeros((7, 7))
    a = np.column_stack([np.eye(6), np.zeros(6)])
    y = np.diag([1j] * 6 + [-1.0]).astype(complex)
    y[0, 6] = y[6, 0] = 2e-5j
    d = saved_blocks(q, y, a)
    np.testing.assert_allclose(d["open_closed_flux_defect"], 2e-5, atol=1e-15, rtol=0)
    np.testing.assert_allclose(d["open_open_flux_defect"], 0, atol=1e-15, rtol=0)
    np.testing.assert_allclose(d["matching_condition"], 1, atol=1e-15, rtol=0)
