"""Recurrence-prefix hook termination, arithmetic and exceptional restoration."""

import numpy as np
import pytest

from projects.ukrmol_co import numerov_co_canary as canary


def arguments():
    q = np.array([[0.3, 0.12], [0.12, -0.2]])
    c2 = np.array([[2.0, -0.4], [-0.4, 6.0]])
    return q, np.array([0.0, 1.8]), np.zeros(2, int), c2, np.zeros((2, 2))


def test_prefix_is_exact_last_original_recurrence_and_skips_tail(monkeypatch):
    original_solve = canary.kernel.refined_solve
    original_taylor = canary.kernel.taylor_change
    captured = []

    def observe(a, b, stats):
        x = original_solve(a, b, stats)
        captured.append(x.copy())
        return x

    monkeypatch.setattr(canary, "REFERENCE_SOLVE", observe)
    monkeypatch.setattr(canary.kernel, "typed_tail", lambda *a: pytest.fail("Tail was evaluated"))
    row = canary.prefix(arguments(), 1.0, 8, 0.025, "reference", dtype=np.float64)
    assert len(captured) == 15
    assert row["final_delta_sha256"] == canary.canonical_hash(captured[-1])
    assert row["matching_evaluations"] == 0
    assert canary.kernel.refined_solve is original_solve
    assert canary.kernel.taylor_change is original_taylor


def test_exception_restores_both_hooks(monkeypatch):
    original = canary.kernel.refined_solve, canary.kernel.taylor_change

    def fail(*args):
        raise np.linalg.LinAlgError("test failure")

    monkeypatch.setattr(canary, "REFERENCE_SOLVE", fail)
    with pytest.raises(np.linalg.LinAlgError, match="test failure"):
        canary.prefix(arguments(), 1.0, 8, 0.025, "reference", dtype=np.float64)
    assert (canary.kernel.refined_solve, canary.kernel.taylor_change) == original


def test_condition_estimate_and_numerical_prefix_against_direct_recurrence(monkeypatch):
    np.testing.assert_allclose(
        canary.condition_estimate(np.diag([1.0, 4.0])), 0.25, rtol=1e-14, atol=0
    )
    saved = []
    solve = canary.REFERENCE_SOLVE

    def observer(a, b, stats):
        x = solve(a, b, stats)
        saved.append((a.copy(), b.copy(), x.copy()))
        return x

    monkeypatch.setattr(canary, "REFERENCE_SOLVE", observer)
    canary.prefix(arguments(), 1.0, 8, 0.025, "reference", dtype=np.float64)
    for a, b, x in saved:
        np.testing.assert_allclose(a @ x, b, rtol=1e-12, atol=1e-15)
