"""Analytic unitary phase factors independently test the comparator budgets."""

import json
from pathlib import Path

import numpy as np
import pytest

from projects.ukrmol_co.sparse_comparison import compare_phases


def contract():
    return json.loads(
        Path(__file__)
        .with_name("sparse-scattering-cas12-4096-comparator-contract.json")
        .read_text()
    )


def phases(position_ev=2.457, width_ev=0.5, branch=0):
    settings = contract()
    requested = 0.1 + 0.05 * np.arange(99)
    energy = requested * settings["native_rydberg_per_requested_ev"] / 2
    center, width = np.array([position_ev, width_ev]) / settings["hartree_to_ev"]
    s = (center - energy + 0.5j * width) / (center - energy - 0.5j * width)
    np.testing.assert_allclose(abs(s), 1, atol=3e-15, rtol=0)
    return np.column_stack((requested, np.angle(s) / 2 + branch * np.pi))


def test_branch_invariant_pair_and_all_background_orders():
    result = compare_phases(phases(), phases(branch=3), contract())
    assert result["root_pair_stability_verdict"] is True
    assert len(result["fixed_window_fits"]) == 4
    for row in result["fixed_window_fits"]:
        assert row["passed"] and row["reference_fit"]["points"] == 31
        np.testing.assert_allclose(
            row["reference_fit"]["position_hartree"] * contract()["hartree_to_ev"],
            2.457,
            atol=1e-7,
            rtol=0,
        )


def test_resolved_position_and_width_changes_fail_original_budgets():
    for shifted in (phases(position_ev=2.557), phases(width_ev=0.6)):
        result = compare_phases(phases(), shifted, contract())
        assert result["root_pair_stability_verdict"] is False
        assert all(row["passed"] is False for row in result["fixed_window_fits"])


def test_full_grid_gate_detects_change_outside_fit_window_and_rejects_missing_points():
    candidate = phases()
    candidate[0, 1] += 0.06
    result = compare_phases(phases(), candidate, contract())
    assert all(row["passed"] for row in result["fixed_window_fits"])
    assert result["root_pair_stability_verdict"] is False
    with pytest.raises(ValueError, match="Incomplete"):
        compare_phases(phases()[:-1], candidate, contract())


def test_flat_nonresonant_phase_does_not_become_a_passing_fit_gate():
    flat = phases()
    flat[:, 1] = 0
    result = compare_phases(flat, flat, contract())
    assert result["root_pair_stability_verdict"] is None
    assert any(row["passed"] is None for row in result["fixed_window_fits"])
