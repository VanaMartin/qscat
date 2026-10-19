"""Replay coordinate coverage, decimal precision and numeric/timing separation."""

import copy
import json
from pathlib import Path

from projects.ukrmol_co.numerov_native_integration_replay import (
    expected_coordinates,
    maximum_decimal_difference,
    numerical_records,
)


def test_replay_declares_all_fifty_coordinates_with_reversed_timing_pairs():
    root = Path(__file__).parent
    coupled = dict(
        contract=json.loads((root / "numerov-coupled-precision-contract.json").read_text())
    )
    dense = dict(contract=json.loads((root / "numerov-dense-scaling-contract.json").read_text()))
    coordinates = expected_coordinates(coupled, dense)
    assert len(coordinates) == 50
    assert coordinates[0] == ("short", "unequal_open", 0.1, "reference", 0)
    assert coordinates[-1] == ("timing", "free", 32, "reference", 1)
    assert len(set(coordinates)) == 50


def test_decimal_reference_error_does_not_round_through_binary64():
    error = maximum_decimal_difference([["1.0000000000000000001"]], [["1.0"]])
    assert 0.999e-19 < error < 1.001e-19


def test_numerical_replay_excludes_timing_but_retains_audit_evidence():
    old = dict(
        rows=[
            dict(
                k_decimal=[["1.0"]],
                residual_audit=dict(accepted_solves=2),
                kernel_wall_seconds=1.0,
                calls=[dict(role="matching", wall_seconds=0.1)],
                checks=dict(timer_partition=True, paired_k_exact=True),
            )
        ]
    )
    fresh = copy.deepcopy(old)
    fresh["rows"][0]["kernel_wall_seconds"] = 2.0
    fresh["rows"][0]["calls"][0]["wall_seconds"] = 0.2
    assert numerical_records(old) == numerical_records(fresh)
    fresh["rows"][0]["residual_audit"]["accepted_solves"] = 1
    assert numerical_records(old) != numerical_records(fresh)
    assert old["rows"][0]["calls"][0]["wall_seconds"] == 0.1
