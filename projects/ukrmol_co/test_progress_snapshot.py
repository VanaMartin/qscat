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
    assert live["new_co_maps"] == live["route_holdout_s_matrices_analyzed"] == 0
    assert live["automatic_successor"] is False and live["full_sweep_released"] is False
    compact = json.loads((DEFAULT_SNAPSHOT.parent / live["verification"]).read_text())
    assert compact["integrations_verified"] == live["toy_integrations_verified"]
    assert compact["science_wall_seconds"] == pytest.approx(
        live["toy_science_wall_seconds"], abs=1e-12, rel=0
    )


def test_coupled_pass_does_not_release_the_proposed_co_study():
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
