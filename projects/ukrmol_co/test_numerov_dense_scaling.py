"""Timing wrappers preserve numerical outputs; affine fits recover known costs."""

import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

import projects.ukrmol_co.numerov_typed_coupled as kernel
from projects.ukrmol_co.numerov_arithmetic_control import decimal_matrix
from projects.ukrmol_co.numerov_dense_scaling import (
    inputs,
    oracle,
    solve_roles,
    timed_propagation,
    timing_analysis,
)


def contract():
    c = json.loads(Path(__file__).with_name("numerov-dense-scaling-contract.json").read_text())
    return c | dict(channels=8, open_channels=2, active_indices=[0, 1, 2], other_open_q=[])


@pytest.mark.parametrize("model", ["free", "varying"])
def test_timer_wrappers_preserve_arrays_and_restore_bindings(model):
    c = contract()
    original = {
        name: getattr(kernel, name) for name in ["refined_solve", "taylor_change", "typed_tail"]
    }
    measured = timed_propagation(c, model, 8, "binary64")
    for name, want in original.items():
        assert getattr(kernel, name) is want
    direct, _ = kernel.typed_coupled_k(
        *inputs(c, model), c["energy_rydberg"], 18.0, 18.125, 0.015625, c["tail_orders"][model]
    )
    assert measured["k_decimal"] == decimal_matrix(direct)
    reference = oracle(c, model, 8, 80)
    np.testing.assert_allclose(direct, np.array(reference.tolist(), float), atol=1e-12, rtol=0)
    assert [
        call["role"] for call in measured["calls"] if call["function"] == "refined_solve"
    ] == solve_roles(model, 8)
    assert measured["timed_function_seconds"] <= measured["kernel_wall_seconds"] + 1e-6


def test_exception_restores_all_numerical_bindings(monkeypatch):
    original = {
        name: getattr(kernel, name) for name in ["refined_solve", "taylor_change", "typed_tail"]
    }

    def failure(*_args, **_kwargs):
        raise ValueError("Deliberate worker failure")

    monkeypatch.setattr(kernel, "typed_coupled_k", failure)
    with pytest.raises(ValueError, match="Deliberate"):
        timed_propagation(contract(), "free", 8, "binary64")
    for name, want in original.items():
        assert getattr(kernel, name) is want


def test_affine_analysis_recovers_known_setup_and_recurring_rates():
    c = contract()
    rows = []
    for model in c["models"]:
        for backend in c["backends"]:
            slope = (0.1 if backend == "longdouble" else 0.001) * (2 if model == "varying" else 1)
            for n in c["interval_counts"]:
                for repetition in range(2):
                    rows.append(
                        dict(
                            model=model,
                            backend=backend,
                            intervals=n,
                            kernel_wall_seconds=2 + slope * n + (-0.01 if repetition else 0.01),
                            recurring_solve_seconds=slope * (n - 1),
                            timing_by_role={"taylor_change": 1.0, "typed_tail": 0.5},
                        )
                    )
    result = timing_analysis(rows, c)
    assert all(result["hypotheses"].values())
    for fit in result["fits"]:
        np.testing.assert_allclose(fit["affine_intercept_seconds"], 2.0, atol=1e-13, rtol=0)
    np.testing.assert_allclose(
        result["varying_to_free_extended_marginal_ratio"], 2.0, atol=1e-13, rtol=0
    )


def test_reference_precision_stabilizes_for_varying_active_block():
    c = contract()
    low, high = [oracle(c, "varying", 8, digits) for digits in [60, 90]]
    with mp.workdps(80):
        assert max(abs(a - b) for a, b in zip(low, high, strict=True)) < mp.mpf("1e-40")
