"""Routing/stop controls for the missing TZ input; no CO solver or package install."""

import json
import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch

from projects.ukrmol_co.compressed_cas12_tz_qc import (
    CONTRACT,
    container_command,
    execute_commands,
    inside,
    qc_command,
)


def test_exact_recipe_changes_basis_not_ensemble_or_solver():
    contract = json.loads(
        Path(__file__).with_name("compressed-cas12-tz-qc-contract-20261026.json").read_text()
    )
    command = qc_command(Path("/old"), Path("/new"), contract)
    assert "--qc-only" in command and "--target-only" not in command
    assert command[command.index("--basis") + 1] == "cc-pVTZ"
    assert command[
        command.index("--active-orbitals") + 1 : command.index("--active-orbitals") + 5
    ] == ["6", "3", "3", "0"]
    assert "--sa-singlet-roots" in command and "--sa-triplet-roots" in command
    assert command[command.index("--target-singlet-a1-driver") + 1] == "spin0"
    assert command[command.index("--target-ci-max-space") + 1] == "160"
    assert command[-1].endswith("co.casscf.chk")


def test_runner_failure_stops_before_analyzer_without_repair():
    with TemporaryDirectory() as temporary:
        path = Path(temporary)
        report = {}
        with patch(
            "projects.ukrmol_co.compressed_cas12_tz_qc.subprocess.run",
            return_value=Mock(returncode=1),
        ) as run:
            with unittest.TestCase().assertRaisesRegex(ValueError, "runner failed"):
                execute_commands(["runner"], path, report)
        assert (
            run.call_count == 1
            and report["original_runner_exit"] == 1
            and "original_analyzer_exit" not in report
        )


def test_analyzer_failure_stays_separate_from_successful_qc():
    with TemporaryDirectory() as temporary:
        path = Path(temporary)
        report = {}
        with patch(
            "projects.ukrmol_co.compressed_cas12_tz_qc.subprocess.run",
            side_effect=[Mock(returncode=0), Mock(returncode=1)],
        ) as run:
            with unittest.TestCase().assertRaisesRegex(ValueError, "analyzer failed"):
                execute_commands(["runner"], path, report)
        assert run.call_count == 2 and report == dict(
            original_runner_exit=0, original_analyzer_exit=1
        )


def test_success_is_exactly_runner_then_existing_analyzer():
    with TemporaryDirectory() as temporary:
        path = Path(temporary)
        report = {}
        with patch(
            "projects.ukrmol_co.compressed_cas12_tz_qc.subprocess.run",
            return_value=Mock(returncode=0),
        ) as run:
            execute_commands(["runner"], path, report)
        assert run.call_count == 2 and report == dict(
            original_runner_exit=0, original_analyzer_exit=0
        )
        assert run.call_args_list[1].args[0][2] == "projects.ukrmol_co.analyze"


def test_owner_has_no_network_or_writable_old_input_and_four_cores():
    command = container_command(
        Path("/old"), Path("/source"), Path("/new"), dict(image="sha256:frozen")
    )
    assert (
        "--cpuset-cpus=0-3" in command and "--memory=24g" in command and "--network=none" in command
    )
    assert [command[i + 1] for i, value in enumerate(command) if value == "-v"] == [
        "/old:/work:ro",
        "/source:/opt/qmodeling:ro",
        "/new:/out",
    ]
    compile(command[-1], "<tz-qc-owner>", "exec")


def test_native_mode_or_wrong_basis_cannot_be_dispatched():
    for args in [
        ["--qc-only", "--target-only", "--basis", "cc-pVTZ"],
        ["--basis", "cc-pVTZ"],
        ["--qc-only", "--basis", "cc-pVDZ"],
    ]:
        with unittest.TestCase().assertRaises(ValueError):
            qc_command(Path("/old"), Path("/new"), dict(runner_args=args))


def test_changed_prerequisite_starts_no_qc_or_dependency_import():
    with TemporaryDirectory() as temporary:
        path = Path(temporary)
        contract = path / CONTRACT
        contract.parent.mkdir(parents=True)
        contract.write_text(json.dumps(dict(input_sha256={"seed": "wanted"}, source_sha256={})))
        with (
            patch("projects.ukrmol_co.compressed_cas12_tz_qc.digest", return_value="wrong"),
            patch("projects.ukrmol_co.compressed_cas12_tz_qc.subprocess.run") as run,
        ):
            code = inside(path, path, path)
        report = json.loads((path / "result.json").read_text())
        run.assert_not_called()
        assert code == 2 and report["qc_calls"] == 0 and report["qc_candidate_verdict"] is None


def test_pinned_runtime_singlet_driver_matches_analytic_hubbard_energy():
    try:
        from pyscf import fci, gto, lib
    except ImportError:
        raise unittest.SkipTest(
            "PySCF is required in the pinned compute image, optional locally"
        ) from None
    import numpy as np

    from projects.ukrmol_co.target import _sector_ci_driver

    lib.num_threads(1)
    mol = gto.M(atom="H 0 0 0; H 0 0 1.4", unit="Bohr", basis="sto-3g", symmetry="C2v", verbose=0)
    solver = _sector_ci_driver(fci, 0, "A1", {"target_singlet_a1_driver": "spin0"})(mol)
    solver.orbsym, solver.wfnsym = np.array([0, 0]), "A1"
    solver = fci.addons.fix_spin_(solver, shift=1.0, ss=0.0)
    solver.conv_tol, solver.conv_tol_residual = 1e-14, 1e-12
    h1 = np.array([[0.0, -1.0], [-1.0, 0.0]])
    h2 = np.zeros((2, 2, 2, 2))
    h2[0, 0, 0, 0] = h2[1, 1, 1, 1] = 4.0
    energy, ci = solver.kernel(h1, h2, 2, (1, 1), ecore=3.0)
    np.testing.assert_allclose(energy, 3 + (4 - np.sqrt(32)) / 2, atol=1e-12, rtol=0)
    effective = solver.absorb_h1e(h1, h2, 2, (1, 1), 0.5)
    action = fci.direct_spin1.contract_2e(effective, ci, 2, (1, 1))
    assert np.linalg.norm(action - (energy - 3) * ci) / np.linalg.norm(ci) <= 1e-12
    np.testing.assert_allclose(solver.spin_square(ci, 2, (1, 1))[0], 0, atol=1e-10, rtol=0)


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(
        0
        if result.wasSuccessful()
        and not (os.environ.get("QSCAT_REQUIRE_PYSCF_CONTROLS") == "1" and result.skipped)
        else 1
    )
