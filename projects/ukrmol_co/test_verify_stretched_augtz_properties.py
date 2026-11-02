"""Preserve observation overruns separately from native/scientific qualification."""

import pytest

from projects.ukrmol_co.verify_stretched_augtz_properties import observation_summary


def test_subsecond_monitor_overrun_is_not_rounded_into_a_strict_cadence_pass():
    rows = [
        dict(observed_utc="2026-11-02T17:52:39+00:00", owner=dict(state="running")),
        dict(observed_utc="2026-11-02T17:57:39.6+00:00", owner=dict(state="finished")),
    ]
    r = observation_summary(rows, "2026-11-02T17:52:04+00:00", "2026-11-02T17:57:30+00:00")
    assert r["maximum_observation_gap_seconds"] == pytest.approx(300.6, rel=0, abs=1e-12)
    assert r["strict_300second_observation_cadence_passed"] is False
    rows[1]["observed_utc"] = "2026-11-02T17:57:39+00:00"
    assert (
        observation_summary(rows, "2026-11-02T17:52:04+00:00", "2026-11-02T17:57:30+00:00")[
            "strict_300second_observation_cadence_passed"
        ]
        is True
    )


def test_missing_terminal_or_reordered_observation_is_not_completion_evidence():
    rows = [dict(observed_utc="2026-11-02T17:52:39+00:00", owner=dict(state="running"))] * 2
    with pytest.raises(ValueError):
        observation_summary(rows, "2026-11-02T17:52:04+00:00", "2026-11-02T17:57:30+00:00")
