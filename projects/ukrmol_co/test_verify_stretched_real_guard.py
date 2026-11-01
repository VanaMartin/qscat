"""Known-truth free Jost S and independent flux/normalization failure arithmetic."""

import numpy as np
import pytest

from projects.ukrmol_co.verify_stretched_real_guard import checks, measurements

POLICY = dict(
    real_s_absolute_budget=1e-6,
    real_phase_sum_budget_rad=1e-6,
    real_flux_absolute_budget=1e-7,
    real_invariant_absolute_budget=1e-6,
)


def test_free_jost_map_gives_identity_s_and_zero_invariant_errors():
    actual, metrics = measurements(np.zeros((6, 6)), 1j * np.eye(6), np.eye(6), np.eye(6))
    np.testing.assert_allclose(actual, np.eye(6), atol=1e-14, rtol=0)
    assert all(checks(metrics, POLICY).values())
    for value in metrics.values():
        np.testing.assert_allclose(value, 0, atol=1e-14, rtol=0)


def test_exact_s_agreement_cannot_hide_wrong_flux_normalization():
    actual, metrics = measurements(np.zeros((6, 6)), 2j * np.eye(6), np.eye(6), np.eye(6))
    np.testing.assert_allclose(actual, np.eye(6), atol=1e-14, rtol=0)
    np.testing.assert_allclose(metrics["flux_defect"], 1.0, atol=1e-14, rtol=0)
    gates = checks(metrics, POLICY)
    assert (
        gates["s_absolute_error"]
        and gates["phase_sum_error_mod_pi_rad"]
        and not gates["flux_defect"]
    )


def test_rank_deficient_jost_rows_cannot_produce_an_accepted_s():
    with pytest.raises(ValueError):
        measurements(np.zeros((6, 6)), 1j * np.eye(6), np.zeros((6, 6)), np.eye(6))


def test_frozen_flux_and_wronskian_budgets_classify_known_guard_stop():
    metrics = dict(
        s_absolute_error=1.7801e-7,
        phase_sum_error_mod_pi_rad=4.017e-7,
        flux_defect=6.9699e-6,
        wronskian_residual=1.4909e-5,
        unitarity_defect=9.1581e-8,
        reciprocity_defect=3.7829e-9,
    )
    gates = checks(metrics, POLICY)
    assert [n for n, ok in gates.items() if not ok] == ["flux_defect", "wronskian_residual"]
