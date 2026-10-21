"""Physical subspace dominance must not depend on the retained root gauge."""

import numpy as np
import pytest

from projects.ukrmol_co.ci_overlap import frozen_core_overlaps
from projects.ukrmol_co.ensemble_frame_comparison import validate_cached
from projects.ukrmol_co.ensemble_frame_overlaps import diagnose_overlap, full_minor_overlaps
from projects.ukrmol_co.test_ci_overlap import full_determinant_oracle


@pytest.mark.parametrize("kind", ["permutation", "unitary"])
def test_retained_gauge_changes_pass(kind):
    overlap = np.eye(8, dtype=complex)
    if kind == "permutation":
        overlap[:5, :5] = np.eye(5)[[4, 0, 2, 1, 3]]
    else:
        overlap[:5, :5] = np.exp(2j * np.pi * np.outer(np.arange(5), np.arange(5)) / 5) / np.sqrt(5)
    result = diagnose_overlap(overlap)
    assert result["retained_subspace_dominance_passed"]
    np.testing.assert_allclose(result["retained_squared_principal_overlaps"], 1, atol=1e-14, rtol=0)
    if kind == "unitary":
        assert all(row["weak_assignment_trigger"] for row in result["old_to_new_assignments"])


def test_retained_buffer_exchange_rejects_both_directions():
    overlap = np.eye(8)[:, [0, 1, 2, 3, 5, 4, 6, 7]]
    result = diagnose_overlap(overlap)
    assert not result["retained_subspace_dominance_passed"]
    for direction in ("old_to_new_assignments", "new_to_old_assignments"):
        assert result[direction][4]["retained_to_extra_trigger"]


def test_lost_physical_overlap_is_not_normalized_away():
    result = diagnose_overlap(0.7 * np.eye(12, 8))
    assert not result["retained_subspace_dominance_passed"]
    np.testing.assert_allclose(
        result["retained_squared_principal_overlaps"], 0.49, atol=1e-14, rtol=0
    )
    np.testing.assert_allclose(
        result["old_retained_projection_into_new_buffer"], 0.49, atol=1e-14, rtol=0
    )


def test_buffer_rows_do_not_compete_as_mandatory_assignment_rows():
    overlap = np.zeros((12, 8))
    overlap[:5, :5] = 0.6 * np.eye(5)
    overlap[5:10, :5] = 0.8 * np.eye(5)
    result = diagnose_overlap(overlap)
    assert [r["assigned_opposite_root"] for r in result["old_to_new_assignments"]] == list(
        range(1, 6)
    )
    assert all(r["retained_to_extra_trigger"] for r in result["new_to_old_assignments"])


@pytest.mark.parametrize("matrix", [np.zeros((4, 8)), np.ones((8, 8)), np.full((8, 8), np.nan)])
def test_invalid_physical_overlap_rejects(matrix):
    with pytest.raises(ValueError):
        diagnose_overlap(matrix)


def test_explicit_full_minors_against_scalar_slater_oracle_and_schur():
    rng = np.random.default_rng(401)
    metric = np.eye(5) + 0.03 * (rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5)))
    left = rng.normal(size=(2, 3, 3)) + 1j * rng.normal(size=(2, 3, 3))
    right = rng.normal(size=(3, 3, 3)) + 1j * rng.normal(size=(3, 3, 3))
    expected = full_determinant_oracle(left, right, metric, 2, (2, 1))
    explicit = full_minor_overlaps(left, right, metric, 2, (2, 1))
    np.testing.assert_allclose(explicit, expected, atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(
        frozen_core_overlaps(left, right, metric, 2, (2, 1)), explicit, atol=1e-12, rtol=1e-12
    )
    reverse = full_minor_overlaps(right, left, metric.conj().T, 2, (2, 1))
    np.testing.assert_allclose(reverse, explicit.conj().T, atol=1e-12, rtol=1e-12)


@pytest.mark.parametrize("defect", [None, "checkpoint", "residual", "spin", "native", "norm"])
def test_cached_wave_gates_are_reconstructed_not_trusted_as_a_flag(defect):
    vectors = np.zeros((12, 462, 462))
    vectors[np.arange(12), 0, np.arange(12)] = 1
    energies = -100 + np.arange(12) * 0.1
    native = {f"singlet.A1.{i + 1}": float(e) for i, e in enumerate(energies[:5])}
    row = dict(
        sector="singlet.A1",
        checkpoint_sha256="pinned",
        gate_passed=True,
        ci_file="singlet.A1.ci.npy",
        native_comparison_root_count=5,
        ci_converged=[True] * 12,
        energy_hartree=energies.tolist(),
        physical_residual_norms_hartree=[1e-10] * 12,
        penalized_residual_norms_hartree=[1e-10] * 12,
        spin_square=[0] * 12,
        original_eight_energy_errors_hartree=[0] * 8,
        original_eight_energy_matched_projection_weights=[1] * 8,
        native_root_errors_hartree=[0] * 5,
    )
    if defect == "checkpoint":
        row["checkpoint_sha256"] = "another"
    elif defect == "residual":
        row["physical_residual_norms_hartree"][-1] = 2e-9
    elif defect == "spin":
        row["spin_square"][-1] = 2
    elif defect == "native":
        native["singlet.A1.5"] += 2e-7
    elif defect == "norm":
        vectors[-1] *= 0.9
    if defect is None:
        assert validate_cached(row, vectors, "pinned", native, 0, "A1")[
            "cached_ci_gates_reconstructed"
        ]
    else:
        with pytest.raises(AssertionError):
            validate_cached(row, vectors, "pinned", native, 0, "A1")
