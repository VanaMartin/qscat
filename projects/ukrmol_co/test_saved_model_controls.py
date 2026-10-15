"""Analytic unitary S factors distinguish prediction quality from fit-parameter stability."""

import json
from pathlib import Path

import numpy as np
import pytest

from projects.ukrmol_co.saved_model_controls import extraction_control, read_model


def contract():
    return json.loads(Path(__file__).with_name("saved-model-controls-contract.json").read_text())


def analytic_table(slope=0.0, branch=0.0):
    requested = 0.1 + 0.05 * np.arange(99)
    energy = requested * 0.0735 / 2
    position, width = 0.036, 0.0105
    background = 0.1 + slope * (energy - position) / 0.02
    s = np.exp(2j * background) * (
        (position - energy + 0.5j * width) / (position - energy - 0.5j * width)
    )
    np.testing.assert_allclose(abs(s), 1, atol=3e-15, rtol=0)
    return np.column_stack((requested, np.angle(s) / 2 + branch)), {
        "energy_hartree": position,
        "width_hartree": width,
    }


@pytest.mark.parametrize("branch", [0.0, 3 * np.pi])
def test_exact_s_factor_passes_every_window_start_and_unused_fold(branch):
    table, seed = analytic_table(branch=branch)
    c = contract()
    result = extraction_control(table, seed, c["extraction"]["stretched_windows_requested_ev"], c)
    assert result["extraction_stability_verdict"] is True
    assert len(result["rows"]) == 12
    assert all(len(r["starts"]) == 3 and len(r["folds"]) == 5 for r in result["rows"])
    assert result["maximum_held_phase_error_rad_mod_pi"] < 1e-8
    assert result["full_width_span_ev"] < 1e-7
    assert result["position_span_ev"] < 1e-7
    for row in result["rows"]:
        np.testing.assert_allclose(
            [row["best_fit"]["position_hartree"], row["best_fit"]["full_width_hartree"]],
            [seed["energy_hartree"], seed["width_hartree"]],
            atol=1e-9,
            rtol=0,
        )


def test_inadequate_background_can_predict_well_but_fail_width_stability():
    table, seed = analytic_table(slope=-0.04)
    c = contract()
    result = extraction_control(table, seed, c["extraction"]["stretched_windows_requested_ev"], c)
    assert all(row["valid"] for row in result["rows"])
    assert result["maximum_held_phase_error_rad_mod_pi"] < 0.05
    assert result["full_width_span_fraction_of_median"] > 0.05
    assert result["extraction_stability_verdict"] is False


def test_one_invalid_declared_case_leaves_verdict_unset():
    table, seed = analytic_table()
    c = contract()
    c["extraction"]["starts"].append([10.0, 1.0])
    result = extraction_control(table, seed, c["extraction"]["stretched_windows_requested_ev"], c)
    assert result["extraction_stability_verdict"] is None
    assert len(result["rows"]) == 12
    assert all("error" in row for row in result["rows"])


def test_original_failed_cas11_record_cannot_be_read_as_success(tmp_path):
    (tmp_path / "comparison.json").write_text('{"all_stages_complete":true}')
    (tmp_path / "resources.json").write_text('{"exit_code":1}')
    with pytest.raises(ValueError, match="Original CAS11"):
        read_model(tmp_path, {"directory": ".", "record": "comparison.json"}, contract())
    assert not (tmp_path / "pilot.json").exists()
