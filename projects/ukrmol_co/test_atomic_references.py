"""No-quantum atomic inventory, energy-zero, release and first-failure controls."""

import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from projects.ukrmol_co import atomic_references as atom
from projects.ukrmol_co.atomic_references_owner import container_command

ROOT = Path(__file__).resolve().parents[2]
INPUTS = "projects/ukrmol_co/atomic-asymptote-inputs-20261103.json"


def test_four_ground_term_occupations_have_correct_charge_spin_and_one_core():
    inputs = json.loads((ROOT / INPUTS).read_text())
    for s in inputs["states"]:
        occ = atom.occupations(s)
        assert occ["Ag"] == (2, 2)
        assert sum(sum(v) for v in occ.values()) == s["electrons"]
        assert sum(a - b for a, b in occ.values()) == s["spin_2s"]
    assert inputs["frozen_spatial_orbitals_per_atom"] == 1


def test_inconsistent_valence_is_rejected_before_quantum_import():
    s = json.loads((ROOT / INPUTS).read_text())["states"][0]
    for change in [dict(electrons=7), dict(spin_2s=0), dict(p_occupations=[3, -1, 0])]:
        with unittest.TestCase().assertRaises(ValueError):
            atom.occupations(s | change)


def test_known_atomic_zero_sums_eas_and_ordering():
    r = atom.endpoint_sums({"C": -37.0, "O": -74.0, "O-": -74.05, "C-": -37.04})
    for key, want in [("C+O", -111.0), ("C+O-", -111.05), ("C-+O", -111.04)]:
        assert math.isclose(r["fragment_sum_hartree"][key], want, abs_tol=1e-12, rel_tol=0)
    for key, want in [("O", 0.05), ("C", 0.04)]:
        assert math.isclose(r["electron_affinity_hartree"][key], want, abs_tol=1e-12, rel_tol=0)
    with unittest.TestCase().assertRaises(ValueError):
        atom.endpoint_sums({"C": -37.0})
    with unittest.TestCase().assertRaises(ValueError):
        atom.endpoint_sums({"C": -37.0, "O": -74.0, "O-": float("nan"), "C-": -37.04})


def test_contract_authenticates_and_no_expanded_basis_or_unauthorized_release():
    c = json.loads((ROOT / atom.CONTRACT).read_text())
    assert atom.authenticate(ROOT, c)["basis"] == "aug-cc-pVQZ"
    for change in [
        dict(execution_authorized=False),
        dict(maximum_atoms=8),
        dict(automatic_second_basis=True),
    ]:
        with unittest.TestCase().assertRaises(ValueError):
            atom.authenticate(ROOT, c | change)
    with unittest.TestCase().assertRaises(ValueError):
        atom.authenticate(ROOT, c | dict(source_sha256={INPUTS: "0" * 64}))


def test_first_failed_atom_preserves_exit_and_never_starts_peer():
    c = json.loads((ROOT / atom.CONTRACT).read_text())
    with (
        tempfile.TemporaryDirectory() as tmp,
        patch("projects.ukrmol_co.run.run_profiled", return_value=124) as run,
    ):
        assert atom.run_queue(ROOT, Path(tmp), c) == 2
        assert run.call_count == 1
        r = json.loads((Path(tmp) / "result.json").read_text())
        assert r["new_scf_calls_reserved"] == 1 and r["atoms"][0]["original_profile_exit"] == 124
        assert "total_energy_hartree" not in r


def test_owner_four_physical_core_readonly_networkless_runtime_and_controls():
    c = json.loads((ROOT / atom.CONTRACT).read_text())
    command = container_command(Path("/src"), Path("/out"), c)
    for token in [
        "--pull=never",
        "--cpuset-cpus=0-3",
        "--memory=16g",
        "--network=none",
        "/src:/opt/qmodeling:ro",
    ]:
        assert token in command
    assert "test_atomic_references" in command[-1]
    assert c["resources"]["observation_interval_seconds"] == 600
    assert not c["automatic_second_basis"] and not c["full_sweep_released"]
    compile(command[-1], "<atomic-owner>", "exec")


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
