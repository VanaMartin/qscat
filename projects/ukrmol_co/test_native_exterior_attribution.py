"""Deck fidelity and separate observable gates for the fixed attribution study."""

import json
from pathlib import Path

import numpy as np
import pytest

from projects.ukrmol_co.native_exterior_attribution import (
    attribution_checks,
    differences,
    subrange_deck,
)


def test_native_deck_mutation_changes_exactly_one_declared_line():
    original = b"&rslvin\n  nescat=99,\n/\n&bprop ! comment retained\n  raf=100.0,\n/\n"
    assert subrange_deck(original, 0) == original
    for count in [16, 32]:
        actual = subrange_deck(original, count)
        inserted = f"  nrange = {count},\n".encode()
        assert actual.count(inserted) == 1
        assert actual.replace(inserted, b"") == original
    with pytest.raises(ValueError):
        subrange_deck(original, 64)
    with pytest.raises(ValueError):
        subrange_deck(subrange_deck(original, 16), 32)


def test_matrix_sequence_comparison_checks_all_records():
    left = np.zeros((99, 3, 3))
    right = left.copy()
    right[-1] = 4e-7 * np.eye(3)
    errors = differences(left, right)
    assert errors["k"] == pytest.approx(4e-7, abs=1e-20, rel=0)
    assert errors["s"] < 1e-6
    assert errors["phase_rad_mod_pi"] > 1e-6


@pytest.mark.parametrize("offset, expected", [(0.0, True), (4e-7, False)])
def test_attribution_cannot_waive_phase_gate_with_passing_k_s(offset, expected):
    contract = json.loads(
        Path(__file__).with_name("native-exterior-attribution-contract.json").read_text()
    )
    contract["energy_index_zero_based"] = 0
    zero = np.zeros((3, 3))
    independent = [{"k": zero.tolist()}, {"k": zero.tolist()}]
    prior = {"k": (1e-8 * np.eye(3)).tolist()}
    native = [
        np.array([5e-6 * np.eye(3)]),
        np.array([(offset + 1e-8) * np.eye(3)]),
        np.array([offset * np.eye(3)]),
    ]
    checks, metrics = attribution_checks(native, independent, prior, contract)
    assert all(checks.values()) is expected
    assert metrics["native_subrange_attribution"] is expected
    assert checks["native_k_reproduction"] and checks["native_s_reproduction"]
    assert checks["native_phase_reproduction"] is expected
