"""Embedding invariance and independent integration timing arithmetic."""

import json
from pathlib import Path

import numpy as np
import pytest

from projects.ukrmol_co import numerov_native_coupled as trial
from projects.ukrmol_co.numerov_dense_scaling import timed_propagation
from projects.ukrmol_co.numerov_native_integration import embed, timing_summary
from projects.ukrmol_co.numerov_native_residual import numpy_residual
from projects.ukrmol_co.numerov_typed_coupled import typed_coupled_k


def test_closed_embedding_preserves_open_observable_and_exercises_mixed_solves():
    q = np.array([[0.3, 0.12], [0.12, -0.2]])
    ts = [0.0, 1.8]
    c2 = np.array([[2.0, -0.4], [-0.4, 6.0]])
    c3 = np.array([[0.1, 0.3], [0.3, -0.2]])
    coordinate = [1.0, 2.0, 4.0, 0.025, 8]
    reference, small_stats = typed_coupled_k(q, ts, [0, 0], c2, c3, *coordinate)
    embedded = embed(q, ts, [0, 0], c2, c3, 1.0)
    result, stats = typed_coupled_k(*embedded, *coordinate)
    np.testing.assert_allclose(result, reference, rtol=0, atol=1e-13)
    assert small_stats["large_system_solves"] == 0
    assert stats["large_system_solves"] == 160
    assert np.all(embedded[1][2:] > 1.0)


def test_dense_timer_hook_nesting_and_exception_restoration(monkeypatch):
    monkeypatch.setattr(trial, "native_residual", numpy_residual)
    c = json.loads(Path(__file__).with_name("numerov-dense-scaling-contract.json").read_text())
    c.update(channels=8, open_channels=2, active_indices=[0, 1, 2], other_open_q=[])
    expected = timed_propagation(c, "varying", 8, "binary64")
    with trial.installed_solver(trial.native_refined_solve):
        actual = timed_propagation(c, "varying", 8, "binary64")
        assert trial.reference.refined_solve is trial.native_refined_solve
    assert trial.reference.refined_solve is trial.REFERENCE_SOLVE
    assert actual["k_decimal"] == expected["k_decimal"]
    assert actual["diagnostics"] == expected["diagnostics"]
    with (
        pytest.raises(RuntimeError, match="already instrumented"),
        trial.installed_solver(trial.native_refined_solve),
    ):
        with trial.installed_solver(trial.native_refined_solve):
            pass
    assert trial.reference.refined_solve is trial.REFERENCE_SOLVE


def test_interval_speedup_separates_setup_and_does_not_gate_accuracy():
    rows = [
        dict(model=model, backend=backend, intervals=n, kernel_wall_seconds=10 + n * rate)
        for model in ["free", "varying"]
        for backend, rate in [("reference", 2.0), ("native", 1.0)]
        for n in [32, 128]
        for _ in range(2)
    ]
    result = timing_summary(rows, 1.5)
    for row in result:
        np.testing.assert_allclose(row["marginal_speedup"], 2.0, rtol=0, atol=1e-14)
        assert row["hypothesis_native_marginal_speedup_at_least_minimum"]
    with pytest.raises(ValueError, match="Incomplete"):
        timing_summary(rows[:-1], 1.5)
