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
    assert live["new_co_maps"] == 8 and live["route_holdout_s_matrices_analyzed"] == 0
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


def test_zero_map_pole_adapter_stop_preserves_real_pass_and_requires_new_decision():
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
    assert (
        proposal["execution_authorized"]
        is live["pole_adapter_repair_execution_authorized"]
        is False
    )
    assert proposal["numerical_algorithms_changed"] is False
