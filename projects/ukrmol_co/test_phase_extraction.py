"""Independent S-factor recovery, held-data isolation and constructive width ambiguity."""

import copy
import json
from pathlib import Path

import numpy as np

from projects.ukrmol_co.phase_extraction import estimate_phase, width_alternatives
from projects.ukrmol_co.phase_extraction_study import analytic_phase
from projects.ukrmol_co.phase_fit import fit_phase, predict_phase

CONTRACT = json.loads(Path(__file__).with_name("phase-extraction-contract.json").read_text())
SETTINGS = CONTRACT["estimator"] | {"hartree_to_ev": CONTRACT["hartree_to_ev"]}


def test_known_linear_background_recovers_pole():
    case = CONTRACT["calibration"]["cases"][0]
    energy = np.linspace(*case["range_hartree"], 61)
    phase, unitarity = analytic_phase(energy, case)
    assert unitarity < 3e-15
    result = estimate_phase(energy, phase, case["windows_hartree"], SETTINGS)
    assert result["estimator_verdict"] is True
    fit = result["reference_fit"]
    np.testing.assert_allclose(
        [fit["position_hartree"], fit["full_width_hartree"]],
        [case["position_hartree"], case["full_width_hartree"]],
        atol=1e-9,
        rtol=0,
    )
    np.testing.assert_allclose(
        np.exp(2j * predict_phase(energy, fit)),
        np.exp(2j * phase),
        atol=1e-8,
        rtol=0,
    )


def test_audit_values_cannot_select_background_or_profile_alternatives():
    case = CONTRACT["calibration"]["cases"][0]
    energy = np.linspace(*case["range_hartree"], 61)
    phase, _ = analytic_phase(energy, case)
    original = estimate_phase(energy, phase, case["windows_hartree"], SETTINGS)
    changed = phase.copy()
    changed[np.arange(len(energy)) % 5 == 0] += 0.1
    audit_failure = estimate_phase(energy, changed, case["windows_hartree"], SETTINGS)
    assert audit_failure["frozen_choices_before_audit"] == original["frozen_choices_before_audit"]
    assert audit_failure["estimator_verdict"] is False
    for a, b in zip(original["rows"], audit_failure["rows"], strict=True):
        for key in ("starts", "folds", "width_alternatives"):
            assert a.get(key) == b.get(key)
        audit = set(a["audit_indices"])
        for fold in a["folds"]:
            assert audit.isdisjoint(fold["training_indices"] + fold["validation_indices"])


def test_narrow_window_has_explicit_beyond_budget_width_counterexample():
    case = next(c for c in CONTRACT["calibration"]["cases"] if c["name"] == "underidentified")
    energy = np.linspace(*case["range_hartree"], 61)
    phase, _ = analytic_phase(energy, case)
    fit = fit_phase(energy, phase, 4, case["position_hartree"], case["full_width_hartree"])
    alternatives = width_alternatives(energy, phase, fit, SETTINGS)
    compatible = next(a for a in alternatives if a["compatible"])
    assert (
        abs(compatible["full_width_hartree"] - fit["full_width_hartree"])
        > 0.05 * fit["full_width_hartree"]
    )
    np.testing.assert_allclose(
        predict_phase(energy, compatible),
        compatible["development_prediction_rad"],
        atol=1e-13,
        rtol=0,
    )
    assert compatible["maximum_phase_error_rad_mod_pi"] <= 0.001


def test_three_start_failure_is_not_replaced_by_a_passing_subset():
    case = CONTRACT["calibration"]["cases"][0]
    energy = np.linspace(*case["range_hartree"], 61)
    phase, _ = analytic_phase(energy, case)
    settings = copy.deepcopy(SETTINGS)
    settings["grid_only_starts"][0][1] = 20
    result = estimate_phase(energy, phase, case["windows_hartree"], settings)
    assert result["estimator_verdict"] is None
    assert all(not row["valid"] for row in result["rows"])
