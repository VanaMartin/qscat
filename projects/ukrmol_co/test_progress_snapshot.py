"""Keep the current evidence overview separate from a CO scientific release."""

import hashlib
import json

import pytest

from projects.ukrmol_co.progress_figure import DEFAULT_SNAPSHOT, ROOT, _snapshot_data


def test_current_snapshot_pins_sources_and_does_not_requalify_the_co_stop():
    snapshot = _snapshot_data(DEFAULT_SNAPSHOT)
    for name, expected in snapshot["source_file_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
    live = snapshot["live_status"]
    assert live["toy_integrations_verified"] == 72 and live["toy_mechanism_verdict"] is True
    assert live["real_guard_original_terminal_exits"] == [2] * 6
    assert live["co_guard_repaired"] is False and live["real_guard_tz_started"] is False
    assert live["new_co_maps"] == 104 and live["route_holdout_s_matrices_analyzed"] == 0
    assert live["automatic_successor"] is False and live["full_sweep_released"] is False
    compact = json.loads((DEFAULT_SNAPSHOT.parent / live["verification"]).read_text())
    assert compact["integrations_verified"] == live["toy_integrations_verified"]
    assert compact["science_wall_seconds"] == pytest.approx(
        live["toy_science_wall_seconds"], abs=1e-12, rel=0
    )


def test_historical_proposal_is_not_retroactively_released_by_the_coupled_pass():
    live = _snapshot_data(DEFAULT_SNAPSHOT)["live_status"]
    assert live["coupled_integrations_verified"] == 54 and live["coupled_control_verdict"] is True
    result = json.loads((DEFAULT_SNAPSHOT.parent / live["coupled_verification"]).read_text())
    assert result["integrations_verified"] == live["coupled_integrations_verified"]
    assert result["science_wall_seconds"] == pytest.approx(
        live["coupled_science_wall_seconds"], abs=1e-12, rel=0
    )
    proposal = json.loads(
        DEFAULT_SNAPSHOT.with_name("stretched-contour-study-proposal-20261102.json").read_text()
    )
    assert proposal["execution_authorized"] is False
    assert proposal["maps_max"] == 3 and proposal["development_index"] == 80


def test_separate_authorized_study_pass_is_not_a_remaining_guard_or_pole_release():
    live = _snapshot_data(DEFAULT_SNAPSHOT)["live_status"]
    assert (
        live["contour_study_verdict"] is True and live["new_co_study_execution_authorized"] is True
    )
    result = json.loads(DEFAULT_SNAPSHOT.with_name(live["study_verification"]).read_text())
    assert result["contour_study_verdict"] is True and result["original_guard_scientific_exit"] == 2
    assert result["tz_started"] is False and result["original_guard_requalified"] is False
    proposal = json.loads(
        DEFAULT_SNAPSHOT.with_name("stretched-real-completion-proposal-20261102.json").read_text()
    )
    assert proposal["execution_authorized"] is False and len(proposal["maps"]) == 5


def test_failed_completion_is_preserved_separately_from_the_corrected_release():
    live = _snapshot_data(DEFAULT_SNAPSHOT)["live_status"]
    assert live["contour_study_verdict"] is True
    assert live["stopped_completion_original_exits"] == [2] * 6
    assert live["stopped_completion_maps_reserved"] == 0
    assert live["stopped_completion_verdict"] is None
    assert len(live["stopped_completion_missing_source_prerequisites"]) == 1
    assert live["capsule_repair_execution_authorized"] is True
    assert live["corrected_completion_original_exits"] == [0] * 6
    assert live["corrected_completion_maps_reserved"] == 5


def test_real_prerequisite_pass_does_not_authorize_complex_use_or_a_pole_search():
    live = _snapshot_data(DEFAULT_SNAPSHOT)["live_status"]
    result = json.loads(DEFAULT_SNAPSHOT.with_name(live["completion_verification"]).read_text())
    assert result["real_completion_verdict"] is live["corrected_completion_verdict"] is True
    assert result["combined_20degree_real_prerequisite_verdict"] is True
    assert live["combined_20degree_real_prerequisite_verdict"] is True
    assert result["previous_dz80_repeated"] is live["previous_dz80_repeated"] is False
    assert result["original_exits"] == live["corrected_completion_original_exits"]
    assert result["zero_map_failure_original_exits"] == live["stopped_completion_original_exits"]
    assert result["tz80_last_two_stability"] == live["tz80_last_two_stability"]
    assert result["profile_wall_seconds"] == pytest.approx(
        live["corrected_completion_profile_wall_seconds"], abs=1e-12, rel=0
    )
    proposal = json.loads(
        DEFAULT_SNAPSHOT.with_name("stretched-dz-pole-proposal-20261102.json").read_text()
    )
    assert result["complex_use_released"] is live["complex_use_released"] is False
    assert live["pole_search_execution_authorized"] is proposal["execution_authorized"] is False
    assert proposal["pole_search_released"] is proposal["tz_counterpart_released"] is False
    assert live["route_holdout_s_matrices_still_unanalyzed"] == 114


def test_zero_map_pole_adapter_stop_is_preserved_separately_from_corrected_success():
    live = _snapshot_data(DEFAULT_SNAPSHOT)["live_status"]
    result = json.loads(DEFAULT_SNAPSHOT.with_name(live["pole_stop_verification"]).read_text())
    assert result["original_exits"] == live["pole_attempt_original_exits"] == [2] * 6
    assert result["real_maps"] == result["complex_maps"] == 0
    assert live["pole_attempt_was_authorized"] is True
    assert live["pole_attempt_real_maps"] == live["pole_attempt_complex_maps"] == 0
    assert result["finite_reference_pole_verdict"] is result["physical_pi_support_verdict"] is None
    assert live["combined_20degree_real_prerequisite_verdict"] is True
    proposal = json.loads(
        DEFAULT_SNAPSHOT.with_name(
            "stretched-dz-pole-adapter-repair-proposal-20261102.json"
        ).read_text()
    )
    assert proposal["execution_authorized"] is False  # Historical proposal, not the later release.
    assert live["pole_adapter_repair_execution_authorized"] is True
    assert proposal["numerical_algorithms_changed"] is False


def test_finite_dz_success_remains_distinct_from_later_tz_or_physical_verdict():
    snapshot = _snapshot_data(DEFAULT_SNAPSHOT)
    live = snapshot["live_status"]
    result = json.loads(DEFAULT_SNAPSHOT.with_name(live["dz_pole_verification"]).read_text())
    assert result["original_exits"] == live["corrected_pole_original_exits"] == [0] * 6
    assert result["finite_reference_pole_verdict"] is live["dz_finite_pole_qualified"] is True
    assert result["physical_pi_support_verdict"] is live["dz_pi_support_qualified"] is True
    assert result["total_maps"] == live["corrected_pole_total_maps"] == 48
    assert result["prior_adapter_failure_original_exits"] == live["pole_attempt_original_exits"]
    point = snapshot["finite_poles"][0]
    for key in ["position_ev", "full_width_ev"]:
        assert point[key] == pytest.approx(result[key], abs=1e-12, rel=0)
    tz = json.loads(
        DEFAULT_SNAPSHOT.with_name("stretched-tz-pole-proposal-20261102.json").read_text()
    )
    assert tz["execution_authorized"] is False  # Historical proposal stays unchanged.
    assert live["tz_pole_execution_authorized"] is True
    assert live["basis_model_verdict"] is None and live["full_sweep_released"] is False


def test_two_numerical_pole_successes_do_not_waive_physical_candidate_budget_rejection():
    snapshot = _snapshot_data(DEFAULT_SNAPSHOT)
    live = snapshot["live_status"]
    tz = json.loads(DEFAULT_SNAPSHOT.with_name(live["tz_pole_verification"]).read_text())
    comparison = json.loads(
        DEFAULT_SNAPSHOT.with_name(live["counterpart_verification"]).read_text()
    )
    assert tz["original_exits"] == live["tz_pole_original_exits"] == [0] * 6
    assert tz["finite_reference_pole_verdict"] is live["tz_finite_pole_qualified"] is True
    assert tz["physical_pi_support_verdict"] is live["tz_pi_support_qualified"] is True
    assert tz["total_maps"] == live["tz_pole_total_maps"] == 48
    assert (
        comparison["candidate_insensitivity_verdict"]
        is live["candidate_insensitivity_verdict"]
        is False
    )
    assert comparison["cross_basis_wavefunction_identity_verdict"] is None
    assert comparison["profile_or_residue_identity_threshold_applied"] is False
    assert comparison["position_difference_ev"] > comparison["position_budget_ev"]
    assert comparison["full_width_difference_ev"] > comparison["full_width_budget_dz_reference_ev"]
    assert [p["basis"] for p in snapshot["finite_poles"]] == ["DZ", "TZ"]
    for key in ["position_ev", "full_width_ev"]:
        assert snapshot["finite_poles"][1][key] == pytest.approx(tz[key], abs=1e-12, rel=0)
    terminal = _snapshot_data(
        DEFAULT_SNAPSHOT.with_name("progress-snapshot-20261102T1635-basis-sensitivity.json")
    )["live_status"]
    assert terminal["active_numerical_job"] is None and terminal["active_monitor"] is False
    assert live["production_model_selected"] is live["full_sweep_released"] is False


def test_terminal_snapshot_does_not_rewrite_the_historical_pending_peer():
    historical = _snapshot_data(
        DEFAULT_SNAPSHOT.with_name("progress-snapshot-20261102T1555-dz-pole.json")
    )
    assert historical["live_status"]["tz_pole_execution_authorized"] is False
    assert len(historical["finite_poles"]) == 1


def test_target_correspondence_does_not_certify_resonance_identity_or_remove_sensitivity():
    snapshot = _snapshot_data(DEFAULT_SNAPSHOT)
    live = snapshot["live_status"]
    result = json.loads(DEFAULT_SNAPSHOT.with_name(live["target_overlap_verification"]).read_text())
    assert result["original_exits"] == [0] * 6
    assert (
        result["retained_eight_dominance_verdict"] is live["target_correspondence_verdict"] is True
    )
    assert result["all64_assignments_preserve_rank"] is True
    rows = snapshot["ensemble_frame_subspaces"]
    assert len(rows) == 8 and all(row["passed"] for row in rows)
    assert [row["minimum_squared_principal_overlap"] for row in rows] == result["sector_minima"]
    assert min(result["sector_minima"]) == pytest.approx(0.9780694190482758)
    assert result["global_resonance_identity_verdict"] is None
    assert live["cross_basis_wavefunction_identity_verdict"] is None
    assert live["candidate_insensitivity_verdict"] is False
    assert live["production_model_selected"] is live["full_sweep_released"] is False
    for name in [
        "native_calls",
        "new_ci_solves",
        "new_hamiltonian_actions",
        "new_maps",
        "route_holdout_s_matrices_analyzed",
    ]:
        assert result[name] == 0
    historical = _snapshot_data(
        DEFAULT_SNAPSHOT.with_name("progress-snapshot-20261102T1635-basis-sensitivity.json")
    )
    assert any(not row["passed"] for row in historical["ensemble_frame_subspaces"])


def test_active_aug_property_job_does_not_invent_terminal_results_or_release_scattering():
    snapshot = _snapshot_data(
        DEFAULT_SNAPSHOT.with_name("progress-snapshot-20261102T1755-augtz-property-launch.json")
    )
    live = snapshot["live_status"]
    launch = json.loads(
        DEFAULT_SNAPSHOT.with_name(live["augtz_property_launch_verification"]).read_text()
    )
    proof = json.loads(
        DEFAULT_SNAPSHOT.with_name(live["augtz_feasibility_verification"]).read_text()
    )
    assert proof["reuse_feasibility_verdict"] is True and proof["missing_states_for_all64"] == 24
    assert proof["native_root_count"] == 64 and proof["new_native_calls"] == 0
    assert launch["property_inventory_verdict"] is launch["original_terminal_exits"] is None
    assert launch["native_denprop_calls"] == 1 and launch["runtime_controls_passed"] == 10
    assert live["active_numerical_job"] == "stretched/augtz/DENPROP64"
    assert live["active_monitor"] is True and live["candidate_insensitivity_verdict"] is False
    assert live["production_model_selected"] is live["full_sweep_released"] is False
    assert live["route_holdout_s_matrices_still_unanalyzed"] == 114


def test_completed_aug_properties_keep_monitor_exception_and_physical_gates_separate():
    snapshot = _snapshot_data(DEFAULT_SNAPSHOT)
    live = snapshot["live_status"]
    result = json.loads(DEFAULT_SNAPSHOT.with_name(live["augtz_property_verification"]).read_text())
    assert (
        result["independent_property_inventory_verdict"]
        is live["augtz_property_inventory_verdict"]
        is True
    )
    assert result["target_count"] == 64 and result["roots_per_sector"] == [8] * 8
    assert result["multipole_records"] == 1928
    assert result["original_exits"] == live["augtz_property_terminal_exits"] == [0] * 8
    assert result["strict_300second_observation_cadence_passed"] is False
    assert result["maximum_observation_gap_seconds"] > 300
    assert result["original_independent_reader_exit"] == 1
    assert result["corrected_independent_reader_exit"] == 0
    assert live["active_numerical_job"] is None and live["active_monitor"] is False
    assert live["candidate_insensitivity_verdict"] is False and live["full_sweep_released"] is False
    assert result["new_ci_calls"] == result["new_qc_calls"] == result["new_maps"] == 0
