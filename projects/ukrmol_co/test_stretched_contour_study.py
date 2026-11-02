"""Synthetic authorization,receipt,routing,failure and owner controls; no CO ODE."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co import stretched_contour_study as study
from projects.ukrmol_co.compressed_cas12_pole_owner import container_command


def fixture(root):
    contract = json.loads(
        Path(__file__).with_name("stretched-contour-study-contract-20261102.json").read_text()
    )
    receipt = dict(
        case="dz",
        inner_export_verdict=True,
        public_reader_exit=0,
        original_exits=[0] * 6,
        contract_sha256=contract["model"]["inner_contract_sha256"],
        roots=2048,
        target_count=64,
        boundary_channels=400,
        spectral_truncation_verdict=None,
        full_sweep_released=False,
    )
    data = {
        "receipt.json": receipt,
        "coupled.json": dict(
            reader_exit=0,
            result=dict(
                coupled_control_verdict=True,
                integrations_verified=54,
                original_exits=[0] * 6,
                co_guard_repaired=False,
            ),
        ),
        "original.json": dict(scientific_exit=2, real_guard_verdict=None),
    }
    contract["model"]["receipt"] = "receipt.json"
    contract["coupled_reading"] = "coupled.json"
    contract["original_guard_result"] = "original.json"
    contract["source_sha256"] = {}
    for name, value in data.items():
        (root / name).write_text(json.dumps(value))
    contract["input_sha256"] = {name: study.digest(root / name) for name in data}
    return contract


def test_explicit_authorization_and_pinned_public_receipts_required():
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        c = fixture(root)
        assert study.authenticate(root, root, c)["case"] == "dz"
        for key, value in [
            ("execution_authorized", False),
            ("angle_degrees", 40),
            ("development_index", 60),
            ("endpoints_bohr", [60, 120]),
        ]:
            with unittest.TestCase().assertRaises(ValueError):
                study.authenticate(root, root, c | {key: value})
        (root / "coupled.json").write_text("{}")
        with unittest.TestCase().assertRaises(ValueError):
            study.authenticate(root, root, c)


def model():
    return dict(
        energies=np.arange(99, dtype=float),
        poles=np.asarray([1.0]),
        amplitudes=np.asarray([[1.0]]),
        ground=0.0,
        saved=np.zeros((99, 1, 1)),
        thresholds=np.asarray([0.0]),
        angular=np.asarray([0]),
        c2=np.zeros((1, 1)),
        c3=np.zeros((1, 1)),
    )


def mocks():
    return [
        patch.object(study, "own_model", return_value=model()),
        patch.object(study, "surface_matrix", return_value=np.zeros((1, 1))),
        patch.object(study, "scattering_matrix", return_value=np.eye(1)),
        patch.object(
            study,
            "physical_s",
            return_value=(
                np.eye(1),
                dict(
                    flux_defect=0.0,
                    wronskian_residual=0.0,
                    unitarity_defect=0.0,
                    reciprocity_defect=0.0,
                ),
            ),
        ),
    ]


def test_fixed_three_maps_preserve_short_diagnostic_failures_and_gate_last():
    from contextlib import ExitStack

    with TemporaryDirectory() as tmp, ExitStack() as stack:
        root = Path(tmp)
        c = fixture(root)
        for mock in mocks():
            stack.enter_context(mock)
        maps = stack.enter_context(
            patch.object(
                study, "outgoing_map", return_value=(np.zeros((1, 1)), np.zeros((1, 1)), {})
            )
        )
        stack.enter_context(
            patch.object(
                study,
                "real_checks",
                side_effect=[
                    dict(checks=dict(flux=False)),
                    dict(checks=dict(flux=False)),
                    dict(checks=dict(flux=True)),
                ],
            )
        )
        assert study.compute(root, root, root, c) == 0
        assert [call.args[6] for call in maps.call_args_list] == [60, 90, 120]
        assert all(call.args[7:9] == (20, 10) for call in maps.call_args_list)
        result = json.loads((root / "result.json").read_text())
        assert result["real_maps_reserved"] == 3 and result["original_guard_requalified"] is False


def test_integration_failure_is_preserved_without_remaining_maps():
    from contextlib import ExitStack

    with TemporaryDirectory() as tmp, ExitStack() as stack:
        root = Path(tmp)
        c = fixture(root)
        for mock in mocks():
            stack.enter_context(mock)
        maps = stack.enter_context(
            patch.object(study, "outgoing_map", side_effect=TimeoutError("synthetic cap"))
        )
        assert study.compute(root, root, root, c) == 2 and maps.call_count == 1
        result = json.loads((root / "result.json").read_text())
        assert result["contour_study_verdict"] is None and result["real_maps_reserved"] == 1


def test_failed_longest_gate_or_last_two_stability_cannot_pass():
    from contextlib import ExitStack

    with TemporaryDirectory() as tmp, ExitStack() as stack:
        root = Path(tmp)
        c = fixture(root)
        for mock in mocks():
            stack.enter_context(mock)
        stack.enter_context(
            patch.object(
                study,
                "outgoing_map",
                side_effect=[(np.full((1, 1), v), np.zeros((1, 1)), {}) for v in [0, 0, 1e-4]],
            )
        )
        assert study.compute(root, root, root, c) == 2
        assert (
            json.loads((root / "result.json").read_text())["last_two_stability"]["y_absolute"]
            > 1e-7
        )


def test_exact_owner_command_has_original_caps_readonly_inputs_and_no_successor():
    c = json.loads(
        Path(__file__).with_name("stretched-contour-study-contract-20261102.json").read_text()
    )
    command = container_command(
        Path("/root"),
        Path("/source"),
        Path("/output"),
        c,
        "projects/ukrmol_co/stretched-contour-study-contract-20261102.json",
    )
    assert (
        "--cpuset-cpus=0-3" in command and "--memory=8g" in command and "--network=none" in command
    )
    assert [command[i + 1] for i, v in enumerate(command) if v == "-v"] == [
        "/root:/work:ro",
        "/source:/opt/qmodeling:ro",
        "/output:/out",
    ]
    assert "projects.ukrmol_co.stretched_contour_study" in command[-1]
    compile(command[-1], "<contour-study-owner>", "exec")


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
