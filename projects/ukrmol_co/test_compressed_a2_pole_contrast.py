"""Adapter allocation/guard and known-truth pole-budget controls; no CO propagation."""

import hashlib
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
from qscat.units import ev_to_hartree

from projects.ukrmol_co.compressed_a2_pole_contrast import (
    contrast_contract,
    load_pole_model,
    pole_insensitivity,
    require_real_guard,
)
from projects.ukrmol_co.pole_character import screen_witnesses

DIRECTORY = Path(__file__).parent
PROPOSAL = json.loads(
    (DIRECTORY / "compressed-a2-pole-contrast-proposal-20261025.json").read_text()
)
BASE = json.loads((DIRECTORY / "compressed-pole-reference-contract-20261023.json").read_text())
CHARACTER = json.loads((DIRECTORY / "paired-pole-character-contract-20261023.json").read_text())


def contract(sector="B1"):
    return contrast_contract(
        PROPOSAL,
        BASE,
        sector,
        "guard.json",
        "pinned",
        seed=0.14 - 0.04j if sector == "B2" else None,
    )


def guard():
    invariants = dict(flux_defect=0, wronskian_residual=0, unitarity_defect=0, reciprocity_defect=0)
    return dict(
        case="compressed/five/B1",
        scientific_exit=0,
        real_rows=[
            dict(
                index=i,
                completed=True,
                s_absolute_error=0,
                phase_sum_error_mod_pi_rad=0,
                invariants=invariants,
            )
            for i in [60, 80]
        ],
    )


def screen(index=0):
    residue = np.zeros((6, 6), complex)
    residue[index, index] = 0.1j
    return screen_witnesses(
        [residue] * 4, CHARACTER["open_labels"]["B1"], 0.08, CHARACTER["character"]
    )


def test_adapter_caps_and_numerics_preserve_the_proposal():
    before = json.dumps(BASE, sort_keys=True)
    for sector, cap in [("B1", 77), ("B2", 65)]:
        actual = contract(sector)
        assert actual["case"] == "compressed/five/" + sector
        assert actual["model_view"] == "compressed-five-property-ids-v1"
        assert actual["settings"] == BASE["settings"]
        assert actual["root_evaluation_caps"] == PROPOSAL["root_caps"][sector]
        assert actual["resources"]["total_map_cap"] == cap
        assert actual["input_sha256"]["guard.json"] == "pinned"
        assert (
            "prepared/saved-matrices-source-20261017/writer-evidence/compressed/five/B1/scatci.inp"
            in actual["input_sha256"]
        )
    assert json.dumps(BASE, sort_keys=True) == before
    with unittest.TestCase().assertRaises(ValueError):
        contrast_contract(PROPOSAL, BASE, "B2", "guard.json", "pinned")


def test_case_dispatch_needs_explicit_view_and_never_selects_equilibrium():
    from projects.ukrmol_co import compressed_a2_pole_contrast

    case = contract()
    with patch.object(
        compressed_a2_pole_contrast, "load_compressed_five", return_value="five"
    ) as load:
        assert load_pole_model(Path("raw/compressed/five/B1"), case) == "five"
        load.assert_called_once()
    for changed in [dict(case, model_view=None), dict(case, case="equilibrium/five/B1")]:
        with unittest.TestCase().assertRaises(ValueError):
            load_pole_model(Path("raw/compressed/five/B1"), changed)


def test_real_guard_requires_pinned_complete_measured_development_gates():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        case = contract()
        for change in ["valid", "holdout", "case", "phase", "incomplete", "flux"]:
            data = guard()
            if change == "holdout":
                data["real_rows"][0]["index"] = 64
            if change == "case":
                data["case"] = "equilibrium/five/B1"
            if change == "phase":
                data["real_rows"][0]["phase_sum_error_mod_pi_rad"] = 2e-6
            if change == "incomplete":
                data["real_rows"].pop()
            if change == "flux":
                data["real_rows"][0]["invariants"]["flux_defect"] = 2e-7
            payload = json.dumps(data).encode()
            (root / "guard.json").write_bytes(payload)
            case["input_sha256"]["guard.json"] = hashlib.sha256(payload).hexdigest()
            if change == "valid":
                require_real_guard(root, case)
            else:
                with unittest.TestCase().assertRaises(ValueError):
                    require_real_guard(root, case)
        case["input_sha256"]["guard.json"] = "wrong"
        with unittest.TestCase().assertRaises(ValueError):
            require_real_guard(root, case)


def test_position_and_full_width_boundary_bands_are_not_best_fit_selection():
    reference = 0.14 - 0.04j
    for quantity in ["position", "width"]:
        tolerance = 0.05 if quantity == "position" else max(0.001, 0.05 * (0.08 / ev_to_hartree(1)))
        band = PROPOSAL["model_budgets"][
            "position_boundary_band_ev" if quantity == "position" else "full_width_boundary_band_ev"
        ]
        for difference, expected in [
            (tolerance - 2 * band, "pass"),
            (tolerance, "inconclusive"),
            (tolerance + 2 * band, "reject"),
        ]:
            delta = ev_to_hartree(difference) * (1 if quantity == "position" else -0.5j)
            result = pole_insensitivity(reference + delta, reference, screen(), screen(), PROPOSAL)
            assert result["disposition"] == expected


def test_half_width_and_wrong_pi_profile_cannot_manufacture_passes():
    reference = 0.14 - 0.04j
    actual = reference - 0.5j * ev_to_hartree(0.15)
    assert (
        pole_insensitivity(actual, reference, screen(), screen(), PROPOSAL)["disposition"]
        == "reject"
    )
    narrow = 0.14 - 0.5j * ev_to_hartree(0.0005)
    result = pole_insensitivity(
        narrow - 0.5j * ev_to_hartree(0.0015), narrow, screen(), screen(), PROPOSAL
    )
    np.testing.assert_allclose(result["difference"]["full_width_ev"], 0.0015, atol=1e-15, rtol=0)
    np.testing.assert_allclose(result["tolerances"]["full_width_ev"], 0.001, atol=1e-15, rtol=0)
    assert result["disposition"] == "reject"  # A half-width error would admit this floor case.
    result = pole_insensitivity(reference, reference, screen(2), screen(0), PROPOSAL)
    assert result["disposition"] == "inconclusive"
    assert result["identity_profile_tv"] > 0.25


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
