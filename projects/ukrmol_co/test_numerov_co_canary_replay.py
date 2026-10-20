"""Independent declared coverage and invariant numerical-record extraction."""

import json
from pathlib import Path

import pytest

from projects.ukrmol_co.numerov_co_canary_replay import audit, coordinates, numerical_row


def test_contract_covers_fixed_six_cases_and_reversed_timing_order():
    c = json.loads(Path(__file__).with_name("numerov-co-canary-contract.json").read_text())
    rows = coordinates(c)
    assert len(rows) == len(set(rows)) == 36
    assert sum(2 * row[2] - 1 for row in rows) == 3996
    assert sum(row[2] - 1 for row in rows) == 1980
    assert rows[0] == ("accuracy", "compressed/nine/B1", 8, "reference", 0)
    assert rows[-1] == ("timing", "compressed/nine/B1", 32, "reference", 1)


def test_replay_preserves_samples_but_discards_fresh_timing_observations():
    old = dict(
        kernel_wall_seconds=1.0,
        recurring_seconds=0.5,
        checks={},
        calls=[dict(role="initialization", wall_seconds=0.3)],
        samples=[dict(x_sha256="abc")],
    )
    fresh = old | dict(
        kernel_wall_seconds=2.0,
        recurring_seconds=1.0,
        calls=[dict(role="initialization", wall_seconds=0.6)],
    )
    assert numerical_row(old) == numerical_row(fresh)
    fresh = fresh | dict(samples=[dict(x_sha256="def")])
    assert numerical_row(old) != numerical_row(fresh)
    assert old["calls"][0]["wall_seconds"] == 0.3


def test_saved_auditor_rejects_wrong_energy_order_before_reading_operands(tmp_path):
    c = json.loads(Path(__file__).with_name("numerov-co-canary-contract.json").read_text())
    row = dict(
        stage="accuracy",
        case="compressed/nine/B1",
        intervals=8,
        backend="reference",
        repetition=0,
        energy_index=88,
        energy_rydberg_hex=c["identities"][0]["energy_rydberg_hex"],
    )
    with pytest.raises(ValueError, match="Energy identity differs"):
        audit(dict(contract=c, rows=[row]), tmp_path)
