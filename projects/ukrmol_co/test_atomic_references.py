"""No-quantum atomic inventory, energy-zero, release and first-failure controls."""

import hashlib
import inspect
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
RELEASE = "projects/ukrmol_co/atomic-asymptote-5z-source-compat-release-20261103.json"


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
    c = json.loads((ROOT / RELEASE).read_text())
    assert atom.authenticate(ROOT, c)["basis"] == "aug-cc-pV5Z"
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
    c = json.loads((ROOT / RELEASE).read_text())
    with (
        tempfile.TemporaryDirectory() as tmp,
        patch("projects.ukrmol_co.run.run_profiled", return_value=124) as run,
    ):
        assert atom.run_queue(ROOT, Path(tmp), c, RELEASE) == 2
        assert run.call_count == 1
        r = json.loads((Path(tmp) / "result.json").read_text())
        assert r["new_scf_calls_reserved"] == 1 and r["atoms"][0]["original_profile_exit"] == 124
        assert "total_energy_hartree" not in r
        assert str(ROOT / RELEASE) in run.call_args.args[2]


def test_owner_four_physical_core_readonly_networkless_runtime_and_controls():
    c = json.loads((ROOT / RELEASE).read_text())
    command = container_command(Path("/src"), Path("/out"), c, RELEASE)
    for token in [
        "--pull=never",
        "--cpuset-cpus=0-3",
        "--memory=16g",
        "--network=none",
        "/src:/opt/qmodeling:ro",
    ]:
        assert token in command
    assert "test_atomic_references" in command[-1]
    assert RELEASE in command[-1]
    assert c["resources"]["observation_interval_seconds"] == 600
    assert not c["automatic_second_basis"] and not c["full_sweep_released"]
    compile(command[-1], "<atomic-owner>", "exec")


def test_calculation_algorithm_is_identical_to_the_executed_qz_implementation():
    c = json.loads((ROOT / RELEASE).read_text())
    assert (
        hashlib.sha256(inspect.getsource(atom.calculate).encode()).hexdigest()
        == c["qz_calculation_source_sha256"]
    )


def test_five_z_rejects_failed_predecessor_and_unsupported_basis_before_quantum_import():
    c = json.loads((ROOT / RELEASE).read_text())
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        paths = {
            key: root / c[key]
            for key in [
                "inputs",
                "predecessor_result",
                "predecessor_inputs",
                "predecessor_public_proof",
            ]
        }
        for key, path in paths.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes((ROOT / c[key]).read_bytes())
        previous = json.loads(paths["predecessor_result"].read_text())
        paths["predecessor_result"].write_text(json.dumps(previous | dict(original_exits=[2] * 8)))
        with unittest.TestCase().assertRaises(ValueError):
            atom.authenticate(root, c | dict(source_sha256={}))
    with unittest.TestCase().assertRaises(ValueError):
        atom.authenticate(ROOT, c | dict(basis="aug-cc-pV6Z"))


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
