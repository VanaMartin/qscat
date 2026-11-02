"""Synthetic five-map routing,authorization,receipt and failure controls; no CO ODE."""

import json
import unittest
from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co import stretched_real_completion as completion
from projects.ukrmol_co.compressed_cas12_pole_owner import container_command


def resolved():
    source = Path(__file__).resolve().parents[2]
    execution = json.loads(
        Path(__file__).with_name("stretched-real-completion-contract-20261102.json").read_text()
    )
    return completion.resolve(source, execution)


def fixture(root):
    c = resolved()
    c["input_sha256"] = {}
    c["source_sha256"] = {}
    for case, spec in c["models"].items():
        receipt = dict(
            case=case,
            inner_export_verdict=True,
            public_reader_exit=0,
            original_exits=[0] * 6,
            contract_sha256=spec["inner_contract_sha256"],
            roots=2048,
            target_count=64,
            boundary_channels=400,
            spectral_truncation_verdict=None,
            full_sweep_released=False,
        )
        spec["receipt"] = case + ".json"
        (root / spec["receipt"]).write_text(json.dumps(receipt))
    c["previous_public_reading"] = "previous.json"
    c["original_guard_result"] = "original.json"
    (root / "previous.json").write_text(
        json.dumps(
            dict(
                reader_exit=0,
                result=dict(
                    contour_study_verdict=True,
                    original_exits=[0] * 6,
                    original_guard_scientific_exit=2,
                    original_guard_requalified=False,
                    tz_started=False,
                ),
            )
        )
    )
    (root / "original.json").write_text(
        json.dumps(dict(scientific_exit=2, real_guard_verdict=None))
    )
    c["input_sha256"] = {
        n: completion.digest(root / n)
        for n in ["dz.json", "tz.json", "previous.json", "original.json"]
    }
    return c


def install_mocks(stack):
    model = dict(
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
    stack.enter_context(patch.object(completion, "own_model", return_value=model))
    stack.enter_context(patch.object(completion, "surface_matrix", return_value=np.zeros((1, 1))))
    stack.enter_context(patch.object(completion, "scattering_matrix", return_value=np.eye(1)))
    stack.enter_context(
        patch.object(
            completion,
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
        )
    )


def test_original_input_and_method_contract_resolves_without_relaxed_gates():
    c = resolved()
    assert len(c["input_sha256"]) == 12 and len(c["source_sha256"]) == 9
    assert c["real_flux_absolute_budget"] == 1e-7 and c["real_s_absolute_budget"] == 1e-6
    assert c["execution_authorized"] is True and c["maps"] == [list(row) for row in completion.MAPS]


def test_unapproved_route_or_changed_public_receipt_starts_no_map():
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        c = fixture(root)
        with (
            patch.object(completion, "own_model") as model,
            patch.object(completion, "outgoing_map") as maps,
        ):
            assert completion.compute(root, root, root, c | dict(execution_authorized=False)) == 2
            (root / "previous.json").write_text("{}")
            assert completion.compute(root, root, root, c) == 2
        model.assert_not_called()
        maps.assert_not_called()


def test_exact_five_maps_never_repeat_passed_dz80():
    with TemporaryDirectory() as tmp, ExitStack() as stack:
        root = Path(tmp)
        c = fixture(root)
        install_mocks(stack)
        maps = stack.enter_context(
            patch.object(
                completion, "outgoing_map", return_value=(np.zeros((1, 1)), np.zeros((1, 1)), {})
            )
        )
        assert completion.compute(root, root, root, c) == 0
        r = json.loads((root / "result.json").read_text())
        assert [(v["case"], v["index"], v["endpoint"]) for v in r["real_rows"]] == completion.MAPS
        assert maps.call_count == 5 and r["previous_dz80_repeated"] is False


def test_failed_first_standalone_gate_leaves_tz_unstarted():
    with TemporaryDirectory() as tmp, ExitStack() as stack:
        root = Path(tmp)
        c = fixture(root)
        install_mocks(stack)
        maps = stack.enter_context(
            patch.object(
                completion, "outgoing_map", return_value=(np.zeros((1, 1)), np.zeros((1, 1)), {})
            )
        )
        stack.enter_context(
            patch.object(completion, "real_checks", return_value=dict(checks=dict(flux=False)))
        )
        assert completion.compute(root, root, root, c) == 2 and maps.call_count == 1


def test_integration_failure_has_no_retry_or_tz_successor():
    with TemporaryDirectory() as tmp, ExitStack() as stack:
        root = Path(tmp)
        c = fixture(root)
        install_mocks(stack)
        maps = stack.enter_context(
            patch.object(completion, "outgoing_map", side_effect=TimeoutError("synthetic cap"))
        )
        assert completion.compute(root, root, root, c) == 2 and maps.call_count == 1


def test_tz_endpoint_stability_and_exact_owner_caps_are_required():
    with TemporaryDirectory() as tmp, ExitStack() as stack:
        root = Path(tmp)
        c = fixture(root)
        install_mocks(stack)
        stack.enter_context(
            patch.object(
                completion,
                "outgoing_map",
                side_effect=[
                    (np.full((1, 1), v), np.zeros((1, 1)), {}) for v in [0, 0, 0, 0, 1e-4]
                ],
            )
        )
        assert completion.compute(root, root, root, c) == 2
        assert (
            json.loads((root / "result.json").read_text())["tz80_last_two_stability"]["y_absolute"]
            > 1e-7
        )
    c = resolved()
    command = container_command(
        Path("/root"),
        Path("/source"),
        Path("/out"),
        c,
        "projects/ukrmol_co/stretched-real-completion-contract-20261102.json",
    )
    assert (
        "--cpuset-cpus=0-3" in command and "--memory=8g" in command and "--network=none" in command
    )
    assert "projects.ukrmol_co.stretched_real_completion" in command[-1]
    compile(command[-1], "<real-completion-owner>", "exec")


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
