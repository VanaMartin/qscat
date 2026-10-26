"""Staged driver diagnostics do not release QC or promote a failed checkpoint."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from projects.ukrmol_co import singlet_a2_pilot as pilot
from projects.ukrmol_co import singlet_a2_reference as reference


def staged(failure):
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        (root / "old.json").write_text(json.dumps(dict(state_energies_hartree=[0] * 40)))
        (root / "failed.json").write_text(
            json.dumps(dict(scientific_exit=2, original_runner_exit=1, qc_candidate_verdict=None))
        )
        contract = dict(
            input_sha256={},
            source_sha256={},
            checkpoint="failed.chk",
            failed_diagnostics="old.json",
            failed_result="failed.json",
        )
        calls = []

        def exact():
            if failure == "exact":
                raise ValueError("failed exact control")
            return dict(passed=True)

        def mixed():
            if failure == "mixed":
                raise ValueError("failed mixed-spin control")
            return dict(passed=True)

        def fixed(path, roots, output):
            calls.append(roots)
            return dict(
                passed=failure != "fixed",
                states=[dict(energy_hartree=float(i)) for i in range(roots)],
            )

        with (
            patch.object(pilot, "authenticate_runtime"),
            patch.object(reference, "exact_controls", side_effect=exact),
            patch.object(reference, "mixed_optimizer_control", side_effect=mixed),
            patch.object(pilot, "fixed_a2", side_effect=fixed),
        ):
            status = pilot.compute(root, root, root, contract)
        return status, json.loads((root / "result.json").read_text()), calls


def test_failed_exact_control_blocks_all_molecular_use():
    code, report, calls = staged("exact")
    assert code == 2 and calls == [] and report["molecular_ci_calls_reserved"] == 0


def test_failed_non_a1_mixed_optimizer_control_blocks_molecular_use():
    code, report, calls = staged("mixed")
    assert code == 2 and calls == [] and report["pilot_verdict"] is None


def test_failed_first_molecular_request_stops_without_larger_request_or_qc():
    code, report, calls = staged("fixed")
    assert code == 2 and calls == [5] and report["molecular_ci_calls_reserved"] == 1
    assert report["new_qc_calls"] == 0 and not report["failed_optimization_qualified"]


def test_known_truth_success_still_releases_no_qc_or_failed_checkpoint():
    code, report, calls = staged("pass")
    assert code == 0 and report["pilot_verdict"] is True and calls == [5, 8]
    assert not report["new_qc_released"] and not report["failed_optimization_qualified"]
    assert report["new_qc_calls"] == report["native_calls"] == report["new_maps"] == 0


def test_singlet_guard_rejects_an_exact_triplet_and_bad_residual():
    row = {
        name: 0.0
        for name in [
            "spin_square",
            "full_residual_hartree",
            "penalized_residual_hartree",
            "spin_penalty_action_hartree",
            "restricted_full_action_difference_hartree",
            "spin_action_norm",
            "alpha_beta_exchange_error",
            "spatial_leakage",
        ]
    }
    assert pilot.state_passes(row)
    assert not pilot.state_passes(row | dict(spin_square=2.0))
    assert not pilot.state_passes(row | dict(full_residual_hartree=1e-3))


def test_small_owner_has_read_only_input_and_no_network():
    command = pilot.container_command(
        Path("/old"), Path("/source"), Path("/new"), dict(image="sha256:pinned")
    )
    assert (
        "--cpuset-cpus=0-3" in command and "--memory=8g" in command and "--network=none" in command
    )
    assert [command[i + 1] for i, value in enumerate(command) if value == "-v"] == [
        "/old:/work:ro",
        "/source:/opt/qmodeling:ro",
        "/new:/out",
    ]


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
