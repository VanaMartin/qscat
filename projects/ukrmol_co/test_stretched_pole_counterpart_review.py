"""Physical comparison never borrows same-root or post-hoc identity thresholds."""

import numpy as np
import pytest
from qscat.units import ev_to_hartree

from projects.ukrmol_co.stretched_pole_counterpart_review import compare


def reading(position, width, profile):
    return dict(
        original_exits=[0] * 6,
        finite_reference_pole_verdict=True,
        physical_pi_support_verdict=True,
        roots=[dict(z_hartree=[ev_to_hartree(position), -ev_to_hartree(width) / 2])],
        character_witnesses=[{}, dict(profile=profile)],
    )


def test_full_width_factor_and_physical_budgets_not_same_root_tolerance():
    left = reading(1.0, 0.2, [1.0, 0.0])
    right = reading(1.02, 0.208, [0.9, 0.1])
    r = compare(left, right, [[1], [1]], [np.eye(2), np.eye(2)])
    assert r["candidate_insensitivity_verdict"] is True
    assert r["difference"]["full_width_ev"] == pytest.approx(0.008, abs=1e-14, rel=0)
    assert r["physical_budgets_dz_reference"]["full_width_ev"] == pytest.approx(
        0.01, abs=1e-14, rel=0
    )
    assert r["cross_basis_wavefunction_identity_verdict"] is None
    assert r["profile_or_residue_identity_threshold_applied"] is False


def test_matching_residue_does_not_waive_large_physical_shift_or_prove_identity():
    left = reading(0.9, 0.26, [1.0, 0.0])
    right = reading(1.01, 0.31, [1.0, 0.0])
    r = compare(left, right, [[1], [1]], [np.eye(2), np.eye(2)])
    assert r["candidate_budget_components"] == dict(position_ev="reject", full_width_ev="reject")
    assert r["candidate_insensitivity_verdict"] is False
    assert r["cross_basis_wavefunction_identity_verdict"] is None
    assert r["full_sweep_released"] is False


def test_unqualified_pole_or_mismatched_labels_block_comparison():
    left = reading(1.0, 0.2, [1.0, 0.0])
    right = reading(1.0, 0.2, [1.0, 0.0])
    with pytest.raises(ValueError):
        compare(left, right, [[1], [2]], [np.eye(2)] * 2)
    right["finite_reference_pole_verdict"] = None
    with pytest.raises(ValueError):
        compare(left, right, [[1], [1]], [np.eye(2)] * 2)


def test_residue_magnitude_diagnostic_is_invariant_to_channel_sign_conventions():
    data = reading(1.0, 0.2, [1.0, 0.0])
    residue = np.outer([1.0, 2.0], [1.0, 2.0]) * 1j
    signs = np.diag([-1.0, 1.0])
    result = compare(data, data, [[1], [1]], [residue, signs @ residue @ signs])
    assert result["normalized_residue_magnitude_distance"] == pytest.approx(0, abs=1e-15, rel=0)
    assert result["cross_basis_wavefunction_identity_verdict"] is None
