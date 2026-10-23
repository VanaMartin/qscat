"""Known-support residues: angular character is not degeneracy or a single l."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co.pole_character import character, compare_pair, screen_witnesses

CONTRACT = json.loads(
    Path(__file__).with_name("paired-pole-character-contract-20261023.json").read_text()
)
POLICY = CONTRACT["character"]
B1 = CONTRACT["open_labels"]["B1"]
B2 = CONTRACT["open_labels"]["B2"]
Z = complex(*CONTRACT["b1_pole_incident_hartree"])


def residue(*entries):
    vector = np.zeros(6, complex)
    for index, value in entries:
        vector[index] = value
    return -0.1j * np.outer(vector, vector)


def screen(value, labels=B1):
    return screen_witnesses([value] * 4, labels, -2 * Z.imag, POLICY)


def test_pi_in_p_or_f_wave_passes_and_higher_m_degenerate_pair_rejects():
    for index in [0, 2]:
        result = screen(residue((index, 1)))
        assert result["disposition"] == "pass"
        np.testing.assert_allclose(result["reference"]["leakage"], 0, atol=1e-15, rtol=0)
    left = screen(residue((3, 1)))
    right = screen(residue((2, 1)), B2)
    paired = compare_pair(left, right, Z, Z, POLICY, 1e-4)
    assert left["disposition"] == right["disposition"] == paired["disposition"] == "reject"


def test_mixed_support_uses_both_sides_and_has_an_inconclusive_band():
    result = screen(residue((0, np.sqrt(0.999)), (3, np.sqrt(0.001))))
    np.testing.assert_allclose(result["reference"]["leakage"], 0.001, atol=1e-14, rtol=0)
    assert result["disposition"] == "inconclusive"
    result = screen(residue((0, 0.1), (3, np.sqrt(0.99))))
    assert result["disposition"] == "reject"  # A nonzero p wave does not suffice.


def test_permutation_complex_phases_and_scaling_do_not_change_character():
    base = residue((0, 1), (2, 0.3j), (3, 0.003))
    reference = character(base, B1, -2 * Z.imag, POLICY)
    order = np.array([3, 0, 5, 2, 1, 4])
    phases = np.exp(1j * np.arange(6) * 0.43)
    changed = (2 - 3j) * phases[:, None] * base * phases[None, :]
    changed = changed[np.ix_(order, order)]
    actual = character(changed, np.asarray(B1)[order], -2 * Z.imag, POLICY)
    for name in ["leakage", "incoming_leakage", "outgoing_leakage"]:
        np.testing.assert_allclose(actual[name], reference[name], atol=1e-14, rtol=0)
    np.testing.assert_allclose(actual["profile"], reference["profile"], atol=1e-14, rtol=0)


def test_subspace_unitary_mixing_preserves_support_not_l_profile():
    rng = np.random.default_rng(81)
    transform = np.zeros((6, 6), complex)
    for indices in [[0, 1, 2, 4], [3, 5]]:
        n = len(indices)
        block, _ = np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
        transform[np.ix_(indices, indices)] = block
    base = residue((0, 1), (3, 0.02j))
    left = character(base, B1, -2 * Z.imag, POLICY)
    right = character(transform @ base @ transform.T, B1, -2 * Z.imag, POLICY)
    np.testing.assert_allclose(right["leakage"], left["leakage"], atol=1e-14, rtol=0)
    assert not np.allclose(left["profile"], right["profile"], atol=1e-6, rtol=0)


def test_b2_label_order_and_signs_are_not_component_identity():
    left = screen(residue((0, 1), (2, 0.4j), (4, 0.2)))
    right = screen(residue((0, -1), (3, -0.4j), (5, -0.2)), B2)
    result = compare_pair(left, right, Z, Z, POLICY, 1e-4)
    assert result["disposition"] == "pass"
    np.testing.assert_allclose(result["profile_tv"], 0, atol=1e-14, rtol=0)
    wrong = screen(residue((0, 1), (1, 0.4j), (5, 0.2)), B2)
    assert compare_pair(left, wrong, Z, Z, POLICY, 1e-4)["disposition"] == "reject"


def test_signal_rank_reciprocity_and_witness_instability_are_inconclusive():
    for value in [np.zeros((6, 6)), residue((0, 1)) + residue((2, 1)), np.eye(6, k=1)]:
        assert screen(value)["disposition"] == "inconclusive"
    first = residue((0, np.sqrt(1 - 0.00009)), (3, np.sqrt(0.00009)))
    last = residue((0, np.sqrt(1 - 0.00011)), (3, np.sqrt(0.00011)))
    assert (
        screen_witnesses([first, first, first, last], B1, -2 * Z.imag, POLICY)["disposition"]
        == "inconclusive"
    )


def test_pair_position_full_width_and_strength_gates():
    value = screen(residue((0, 1)))
    for other_z in [Z + 1e-5, Z - 1e-5j]:
        assert compare_pair(value, value, Z, other_z, POLICY, 1e-4)["disposition"] == "reject"
    right = screen(1.02 * residue((0, 1)), B2)
    assert compare_pair(value, right, Z, Z, POLICY, 1e-4)["disposition"] == "reject"


def test_invalid_labels_width_and_nonfinite_input_cannot_qualify():
    case = unittest.TestCase()
    for labels in [B1[:-1], [B1[0]] * 6, [[i + 1, 2, 1, 1] for i in range(6)]]:
        with case.assertRaises(ValueError):
            character(residue((0, 1)), labels, -2 * Z.imag, POLICY)
    for width in [0, -1, float("nan")]:
        with case.assertRaises(ValueError):
            character(residue((0, 1)), B1, width, POLICY)
    with case.assertRaises(ValueError):
        character(np.full((6, 6), np.nan), B1, -2 * Z.imag, POLICY)


def test_saved_stage_does_not_release_b2_or_hide_after_integrity_failure():
    from projects.ukrmol_co import saved_pole_character

    calls = []

    def verify(*_):
        calls.append(1)
        if len(calls) == 2:
            raise ValueError("Changed input after scoring")
        return 123

    with (
        TemporaryDirectory() as temporary,
        patch.object(saved_pole_character, "verify_inputs", verify),
        patch.object(saved_pole_character, "verified_labels", lambda *_: {"B1": B1, "B2": B2}),
        patch.object(saved_pole_character, "load_b1", lambda *_: (Z, [residue((0, 1))] * 4, 86)),
    ):
        root = Path(temporary)
        output = root / "screen"
        assert saved_pole_character.run(root, output, CONTRACT) == 2
        report = json.loads((output / "result.json").read_text())
    assert report["b1_character_disposition"] == "pass"
    assert report["paired_verdict"] is None
    assert report["b2_calculation_eligible"] is False
    assert report["b2_execution_released"] is False
    assert report["new_propagations"] == report["route_holdout_s_matrices_analyzed"] == 0


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
