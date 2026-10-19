"""Tracer isolation, exact numerical preservation, and reusable operand capture."""

import json
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import LinAlgWarning

import projects.ukrmol_co.numerov_typed_coupled as kernel
from projects.ukrmol_co.numerov_arithmetic_control import decimal_matrix
from projects.ukrmol_co.numerov_dense_scaling import inputs
from projects.ukrmol_co.numerov_solve_profile import profiled_row, trace_solve


@pytest.mark.parametrize("model", ["free", "varying"])
def test_trace_preserves_propagation_and_captured_solves(model, tmp_path):
    c = json.loads(Path(__file__).with_name("numerov-dense-scaling-contract.json").read_text())
    c.update(
        channels=8, open_channels=2, active_indices=[0, 1, 2], other_open_q=[], profile_intervals=8
    )
    original = kernel.refined_solve
    measured = profiled_row(c, model, "binary64", True, tmp_path, 4)
    assert kernel.refined_solve is original and sys.gettrace() is None
    direct, diagnostics = kernel.typed_coupled_k(
        *inputs(c, model), 0.7, 18.0, 18.125, 0.015625, c["tail_orders"][model]
    )
    assert measured["k_decimal"] == decimal_matrix(direct)
    assert measured["diagnostics"] == diagnostics
    assert len(measured["captures"]) == (1 if model == "free" else 2)
    for capture in measured["captures"]:
        with np.load(tmp_path / capture["file"]) as data:
            np.testing.assert_allclose(data["a"] @ data["x"], data["b"], atol=2e-15, rtol=2e-14)
            assert data["a"].shape == data["b"].shape == data["x"].shape == (8, 8)
    for call in measured["calls"]:
        assert 0 <= sum(call["category_seconds"].values()) <= call["wall_seconds"]
        assert any(x["source"] == "residual = b - a @ x" for x in call["lines"])


def test_exception_restores_tracer():
    stats = dict(
        maximum_solve_backward_error=0.0, maximum_refinement_corrections=0, large_system_solves=0
    )
    with pytest.raises(LinAlgWarning):
        trace_solve(kernel.refined_solve, np.zeros((4, 4)), np.eye(4), stats)
    assert sys.gettrace() is None


def test_existing_tracer_is_preserved():
    def prior_trace(_frame, _event, _arg):
        return None

    try:
        sys.settrace(prior_trace)
        with pytest.raises(ValueError, match="Existing Python tracer"):
            trace_solve(kernel.refined_solve, np.eye(4), np.eye(4), {})
        assert sys.gettrace() is prior_trace
    finally:
        sys.settrace(None)
