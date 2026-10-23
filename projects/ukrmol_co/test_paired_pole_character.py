"""Adapter-only controls; no physical pole, propagation or actual residue analysis."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co.paired_pole_character import b2_contract, real_checks, run

DIRECTORY = Path(__file__).parent
POLICY = json.loads((DIRECTORY / "paired-pole-character-contract-20261023.json").read_text())
BASE = json.loads((DIRECTORY / "compressed-pole-reference-contract-20261023.json").read_text())
EXECUTION = json.loads((DIRECTORY / "paired-pole-execution-contract-20261023.json").read_text())


def test_b2_translation_keeps_numerics_and_exact_total_map_cap():
    before = json.dumps(BASE, sort_keys=True)
    result = b2_contract(POLICY, BASE, EXECUTION)
    assert json.dumps(BASE, sort_keys=True) == before
    assert result["case"] == "compressed/nine/B2"
    assert result["settings"] == BASE["settings"]
    assert result["root_evaluation_caps"] == [12, 8, 8]
    assert result["resources"]["total_map_cap"] == sum(result["root_evaluation_caps"]) + 37 == 65
    assert result["resources"]["total_map_cap"] + len(POLICY["b2_development_indices"]) == 67
    np.testing.assert_allclose(
        result["seeds_incident_hartree"][0], POLICY["b1_pole_incident_hartree"], atol=1e-15, rtol=0
    )
    for key in [
        "contour_radius_hartree",
        "contour_nodes",
        "derivative_steps_hartree",
        "same_root_position_full_width_budget_ev",
    ]:
        assert result[key] == BASE[key]


def test_real_phase_gate_is_not_hidden_by_small_individual_entries():
    invariants = dict(flux_defect=0, wronskian_residual=0, unitarity_defect=0, reciprocity_defect=0)
    reference = np.eye(6, dtype=complex)
    assert all(real_checks(reference, reference, invariants, POLICY)["checks"].values())
    result = real_checks(np.exp(0.4e-6j) * reference, reference, invariants, POLICY)
    assert result["checks"]["s"] is True
    assert result["checks"]["phase"] is False


def test_nonpassing_b1_never_loads_b2_or_calls_pole_runner():
    from projects.ukrmol_co import paired_pole_character

    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        left = root / EXECUTION["b1_screen"]
        owner = root / EXECUTION["b1_screen_execution"]
        left.parent.mkdir(parents=True)
        left.write_text(json.dumps(dict(scientific_exit=1, b1_character_disposition="reject")))
        owner.write_text("{}")
        with (
            patch.object(paired_pole_character, "digest", return_value="pinned"),
            patch.object(paired_pole_character, "verify_inputs", return_value=123),
            patch.object(paired_pole_character, "load_model") as load,
            patch.object(paired_pole_character, "run_pole") as pole,
        ):
            execution = dict(
                EXECUTION, policy_sha256="pinned", pole_contract_sha256="pinned", input_sha256={}
            )
            assert run(root, DIRECTORY.parent.parent, root / "out", execution) == 2
            load.assert_not_called()
            pole.assert_not_called()
        report = json.loads((root / "out/result.json").read_text())
        assert report["paired_verdict"] is None
        assert report["new_propagations"] == 0


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
