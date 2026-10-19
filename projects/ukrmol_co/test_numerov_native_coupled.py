"""Mixed-solve integration, independent residual decisions and failure restoration."""

from functools import partial

import numpy as np
import pytest
from scipy.linalg import LinAlgWarning

import projects.ukrmol_co.numerov_native_coupled as trial
from projects.ukrmol_co.numerov_native_residual import numpy_residual


def stats():
    return dict(
        maximum_solve_backward_error=0.0, maximum_refinement_corrections=0, large_system_solves=0
    )


def test_algorithm_hook_preserves_refinement_with_numpy_residual(monkeypatch):
    monkeypatch.setattr(trial, "native_residual", numpy_residual)
    v = np.arange(1, 7, dtype=float) / 16
    a = np.eye(6) + np.outer(v, v)
    b = np.column_stack([v, np.ones(6)])
    old, new = stats(), stats()
    want = trial.REFERENCE_SOLVE(a, b, old)
    actual = trial.native_refined_solve(a, b, new)
    np.testing.assert_array_equal(actual, want)
    assert old == new and new["large_system_solves"] == 1


def test_failure_restores_oracle_binding(monkeypatch):
    monkeypatch.setattr(trial, "native_residual", numpy_residual)
    with pytest.raises(LinAlgWarning), trial.installed_solver(trial.native_refined_solve):
        trial.reference.refined_solve(np.zeros((4, 4)), np.eye(4), stats())
    assert trial.reference.refined_solve is trial.REFERENCE_SOLVE


def test_observer_detects_false_acceptance():
    dtype = np.longdouble
    a = x = np.eye(4, dtype=dtype)
    b = 2 * a
    audit = trial.ResidualAudit()
    audit(a, b, x, np.zeros_like(b), dtype(0), 0)
    assert audit.summary()["gate_disagreements"] == 1
    assert audit.summary()["maximum_scaled_difference"] > 32 * np.finfo(dtype).eps


def test_refinement_exhaustion_preserves_three_correction_bound(monkeypatch):
    monkeypatch.setattr(trial, "native_residual", lambda a, b, x: np.ones_like(b))
    iterations = []
    s = stats()
    with pytest.raises(np.linalg.LinAlgError, match="three corrections"):
        trial.native_refined_solve(
            np.eye(4),
            np.eye(4),
            s,
            observer=lambda a, b, x, residual, error, iteration: iterations.append(iteration),
        )
    assert iterations == [0, 1, 2, 3]
    assert s["large_system_solves"] == 0


def test_native_rank_one_solve_and_independent_acceptance():
    native = pytest.importorskip("qscat_extended_kernels")
    if not native.capability()[3] or np.finfo(np.longdouble).nmant + 1 != 64:
        pytest.skip("Requires qualified native extended arithmetic")
    v = np.arange(1, 13, dtype=np.longdouble) / 512
    a = np.eye(12, dtype=np.longdouble) + np.outer(v, v) / 512
    b = np.column_stack([v, np.ones(12, dtype=np.longdouble)])
    audit = trial.ResidualAudit()
    s = stats()
    actual = partial(trial.native_refined_solve, observer=audit)(a, b, s)
    expected = b - v[:, None] * (v @ b)[None, :] / (512 + v @ v)
    np.testing.assert_allclose(actual, expected, atol=2e-18, rtol=0)
    assert audit.summary()["evaluations"] > 0 and audit.summary()["gate_disagreements"] == 0
    assert audit.summary()["maximum_accepted_numpy_error"] <= 32 * np.finfo(np.longdouble).eps
