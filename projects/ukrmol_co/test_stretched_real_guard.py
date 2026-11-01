"""Known-truth receipt/layout/routing/owner controls; never propagates a CO model."""

import hashlib
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co import stretched_real_guard as guard
from projects.ukrmol_co.compressed_cas12_pole_owner import container_command


def receipt(case):
    return dict(
        case=case,
        inner_export_verdict=True,
        public_reader_exit=0,
        original_exits=[0] * 6,
        contract_sha256=case,
        roots=2048,
        target_count=64,
        boundary_channels=400,
        spectral_truncation_verdict=None,
        full_sweep_released=False,
    )


def fixture(root, changes=None):
    contract = dict(input_sha256={}, source_sha256={}, models={})
    for case in ["dz", "tz"]:
        path = root / (case + ".json")
        path.write_text(json.dumps(receipt(case) | (changes or {})))
        contract["input_sha256"][path.name] = guard.digest(path)
        contract["models"][case] = dict(
            receipt=path.name, inner_contract_sha256=case, native_directory="native/" + case
        )
    return contract


def test_only_hash_pinned_public_success_receipts_allow_model_access():
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        c = fixture(root)
        assert list(guard.authenticate(root, root, c)) == ["dz", "tz"]
        for changes in [
            dict(public_reader_exit=1),
            dict(case="tz"),
            dict(original_exits=[0] * 5 + [1]),
            dict(roots=4096),
            dict(target_count=40),
            dict(inner_export_verdict=False),
        ]:
            c = fixture(root, changes)
            with unittest.TestCase().assertRaises(ValueError):
                guard.authenticate(root, root, c)
        c = fixture(root)
        (root / "dz.json").write_text("{}")
        with unittest.TestCase().assertRaises(ValueError):
            guard.authenticate(root, root, c)
        c = fixture(root)
        (root / "tz.json").unlink()
        with unittest.TestCase().assertRaises(FileNotFoundError):
            guard.authenticate(root, root, c)


def test_invalid_receipt_starts_no_model_or_map():
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        c = fixture(root, dict(public_reader_exit=2))
        with (
            patch.object(guard, "own_model") as model,
            patch.object(guard, "real_guard") as propagation,
        ):
            assert guard.compute(root, root, root, c) == 2
        model.assert_not_called()
        propagation.assert_not_called()


def test_strict64_400_assembly_and_independent_array_hash_are_required():
    model = dict(poles=np.asarray([1.0, 2.0]))
    want = hashlib.sha256(model["poles"].astype("<f8").tobytes()).hexdigest()
    with (
        patch.object(guard, "read_interface", return_value={}),
        patch.object(guard, "_model_from_interface", return_value=model) as assembly,
    ):
        assert guard.own_model(Path("/synthetic"), dict(array_sha256=dict(poles=want))) is model
        assert assembly.call_args.args[-1] == (64, 400)
        with unittest.TestCase().assertRaises(ValueError):
            guard.own_model(Path("/synthetic"), dict(array_sha256=dict(poles="0" * 64)))


def test_exactly_four_real_maps_no_complex_or_model_verdict():
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        c = fixture(root)

        def real(model, directory, case, proposal, report):
            report["real_maps_reserved"] += 2
            path = directory / "guard.json"
            path.write_text("{}")
            return path

        with (
            patch.object(guard, "own_model", return_value={}),
            patch.object(guard, "real_guard", side_effect=real) as maps,
        ):
            assert guard.compute(root, root, root, c) == 0
        assert [call.args[2] for call in maps.call_args_list] == [
            "stretched/cas12/dz/B1",
            "stretched/cas12/tz/B1",
        ]
        report = json.loads((root / "result.json").read_text())
        assert report["real_maps_reserved"] == 4 and report["real_guard_verdict"] is True
        assert report["new_complex_maps"] == report["native_calls"] == 0
        assert report["pole_verdict"] is report["model_insensitivity_verdict"] is None


def test_failed_first_real_map_leaves_tz_unstarted():
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        c = fixture(root)

        def fail(model, directory, case, proposal, report):
            report["real_maps_reserved"] += 1
            raise ValueError("known failed real-axis check")

        with (
            patch.object(guard, "own_model", return_value={}) as model,
            patch.object(guard, "real_guard", side_effect=fail),
        ):
            assert guard.compute(root, root, root, c) == 2
        assert model.call_count == 1
        report = json.loads((root / "result.json").read_text())
        assert report["real_maps_reserved"] == 1 and report["real_guard_verdict"] is None


def test_four_core_readonly_owner_is_real_only_and_has_no_successor():
    c = dict(
        image="sha256:pinned",
        container_name="co-stretched-real-guard-20261101",
        runtime_controls_module="projects.ukrmol_co.test_stretched_real_guard",
        runner_module="projects.ukrmol_co.stretched_real_guard",
    )
    command = container_command(
        Path("/root"),
        Path("/source"),
        Path("/out"),
        c,
        "projects/ukrmol_co/stretched-real-guard-contract-20261101.json",
    )
    assert (
        "--cpuset-cpus=0-3" in command and "--memory=8g" in command and "--network=none" in command
    )
    assert c["container_name"] in command and "stretched_real_guard" in command[-1]
    assert [command[i + 1] for i, v in enumerate(command) if v == "-v"] == [
        "/root:/work:ro",
        "/source:/opt/qmodeling:ro",
        "/out:/out",
    ]
    compile(command[-1], "<real-guard-owner>", "exec")


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
