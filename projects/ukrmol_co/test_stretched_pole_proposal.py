"""Analytic sheet and finite allocation checks; no CO integration or pole search."""

import json
from pathlib import Path

import numpy as np
import pytest


def proposal():
    return json.loads(
        Path(__file__).with_name("stretched-dz-pole-proposal-20261102.json").read_text()
    )


def test_every_open_box_energy_is_exposed_by_both_prescribed_angles():
    c = proposal()
    box = c["search_box_incident_hartree"]
    theta_min = min(s[1] for s in c["settings"])
    assert np.arctan(max(abs(box[2]), abs(box[3])) / box[0]) / 2 < np.deg2rad(theta_min)
    assert 2 * box[1] < 0.370  # strictly below either authenticated first closed threshold
    assert c["execution_authorized"] is False and c["pole_search_released"] is False


def test_exact_map_cap_includes_near_feature_guards_and_all_witnesses():
    c = proposal()
    total = (
        len(c["near_feature_real_settings"])
        + sum(c["root_evaluation_caps"])
        + c["contour_nodes"]
        + 2 * len(c["derivative_steps_hartree"])
        + 1
    )
    assert total == c["resources"]["total_map_cap"] == 79
    assert c["same_root_position_full_width_budget_ev"] == pytest.approx(0.0001, abs=0, rel=0)
    assert c["tz_counterpart_released"] is False and c["native_calls"] == 0


def test_seed_and_initial_witness_stencil_stay_inside_the_declared_sheet_box():
    c = proposal()
    left, right, bottom, top = c["search_box_incident_hartree"]
    radius = max(c["contour_radius_hartree"], *c["derivative_steps_hartree"])
    for real, imag in c["seeds_incident_hartree"]:
        assert left + radius < real < right - radius
        assert bottom + radius < imag < top - radius
    # Final-root stencils must be checked again by the future adapter; no box expansion.
    np.testing.assert_allclose(
        c["reference_energy_rydberg"],
        2 * c["seeds_incident_hartree"][0][0],
        atol=1e-15,
        rtol=0,
    )


def test_near_feature_real_guards_keep_the_original_physical_budgets():
    c = proposal()
    original = json.loads(
        Path(__file__).with_name("stretched-real-guard-contract-20261101.json").read_text()
    )
    for name in [
        "real_s_absolute_budget",
        "real_phase_sum_budget_rad",
        "real_invariant_absolute_budget",
        "real_flux_absolute_budget",
    ]:
        assert c[name] == original[name]
    assert c["near_feature_real_index"] == 16
    assert [row[1] for row in c["near_feature_real_settings"]] == [15, 20]
    assert c["automatic_successor"] is c["full_sweep_released"] is False
