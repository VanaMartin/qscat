"""Strict own-data routing and stop controls plus existing analytic pole oracles."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co import stretched_dz_pole as pole
from projects.ukrmol_co import test_outgoing_map, test_pole_character, test_pole_reference
from projects.ukrmol_co.stretched_dz_pole_owner import CONTRACT, command


def resolved():
    source = Path(__file__).resolve().parents[2]
    return pole.resolve(source, json.loads((source / CONTRACT).read_text()))


def test_actual_source_prerequisites_and_authorized_proposal_match():
    c = resolved()
    source = Path(__file__).resolve().parents[2]
    for name, want in c["source_sha256"].items():
        assert pole.core.digest(source / name) == want, name
    assert c["execution_authorized"] and c["resources"]["total_map_cap"] == 79
    assert len(c["input_sha256"]) == 7 and not c["tz_counterpart_released"]


def test_loader_rejects_a_sliced_or_wrong_open_model():
    with patch.object(
        pole,
        "own_model",
        return_value=dict(
            poles=np.zeros(2048), amplitudes=np.zeros((2048, 296)), thresholds=np.zeros(296)
        ),
    ):
        with unittest.TestCase().assertRaisesRegex(ValueError, "Strict own64"):
            pole.load_own(Path("/unused"), resolved(), {})


def test_positive_actual_return_schema_uses_native_metadata_not_a_model_interface():
    c = resolved()
    model = dict(
        poles=np.zeros(2048),
        amplitudes=np.zeros((2048, 400)),
        thresholds=np.r_[np.zeros(6), np.ones(394)],
    )
    assert "interface" not in model
    with (
        patch.object(pole, "own_model", return_value=model),
        patch.object(
            pole, "read_interface", return_value=dict(channel_labels=c["open_labels"])
        ) as metadata,
    ):
        assert pole.load_own(Path("/native-root"), c, {}) is model
        assert metadata.call_args.args[0].name == "fort.10"
        assert metadata.call_args.args[1].name == "fort.21"
    with (
        patch.object(pole, "own_model", return_value=model),
        patch.object(pole, "read_interface", return_value=dict(channel_labels=[[1, 2, 1, 1]] * 6)),
    ):
        with unittest.TestCase().assertRaises(ValueError):
            pole.load_own(Path("/native-root"), c, {})


def test_corrected_release_keeps_all_numerical_parameters_and_binds_original_stop():
    source = Path(__file__).resolve().parents[2]
    execution = json.loads(
        (source / "projects/ukrmol_co/stretched-dz-pole-adapter-release-20261102.json").read_text()
    )
    original = resolved()
    corrected = pole.resolve(source, execution)
    changed = {k for k in set(original) | set(corrected) if original.get(k) != corrected.get(k)}
    assert changed == {
        "status",
        "container_name",
        "input_sha256",
        "failed_stop_reading",
        "base_execution_contract",
        "base_execution_sha256",
    }
    assert corrected["input_sha256"].items() >= original["input_sha256"].items()
    assert (
        corrected["input_sha256"][corrected["failed_stop_reading"]]
        == "d3ccd7d14b8576d407f03bdfae40153538e6f9309fa25a05d1154b7643827d7f"
    )
    cmd = command(
        Path("/root"),
        Path("/source"),
        Path("/out"),
        corrected,
        "projects/ukrmol_co/stretched-dz-pole-adapter-release-20261102.json",
    )
    assert "stretched-dz-pole-adapter-release-20261102.json" in cmd[-1]


def test_new_authorization_and_missing_input_failure_start_no_model_or_maps():
    with (
        TemporaryDirectory() as tmp,
        patch.object(pole, "load_own") as model,
        patch.object(pole.core, "calculate") as calc,
    ):
        root = Path(tmp)
        c = resolved()
        c["execution_authorized"] = False
        assert pole.run(root, root, root, c) == 2
        r = json.loads((root / "result.json").read_text())
        assert r["new_maps_reserved"] == 0 and r["finite_reference_pole_verdict"] is None
        model.assert_not_called()
        calc.assert_not_called()


def test_failed_near_feature_real_gate_prevents_all_complex_work():
    c = resolved()
    model = dict(
        energies=np.ones(99) * 0.06615,
        poles=np.array([1.0]),
        amplitudes=np.array([[1.0]]),
        ground=-1.0,
        thresholds=np.array([0.0]),
        angular=np.array([0]),
        c2=np.zeros((1, 1)),
        c3=np.zeros((1, 1)),
        saved=np.zeros((99, 1, 1)),
    )
    with (
        TemporaryDirectory() as tmp,
        patch.object(pole, "authenticate", return_value={}),
        patch.object(pole, "load_own", return_value=model),
        patch.object(
            pole, "outgoing_map", return_value=(np.zeros((1, 1)), np.ones((1, 1)), {})
        ) as maps,
        patch.object(pole, "real_checks", return_value=dict(checks=dict(flux=False))),
        patch.object(pole.core, "calculate") as calc,
    ):
        root = Path(tmp)
        assert pole.run(root, root, root, c) == 2
        assert maps.call_count == 1
        calc.assert_not_called()
        r = json.loads((root / "result.json").read_text())
        assert r["evaluations"] == [] and len(r["real_rows"]) == 1


def test_adapter_populates_unchanged_calculator_without_compressed_loader():
    c = resolved()
    model = dict(
        ground=-1.0,
        poles=np.array([0.1]),
        amplitudes=np.array([[0.2]]),
        thresholds=np.array([0.0]),
        angular=np.array([0]),
        c2=np.zeros((1, 1)),
        c3=np.zeros((1, 1)),
    )
    with (
        TemporaryDirectory() as tmp,
        patch.object(pole, "load_own", return_value=model),
        patch.object(pole.core, "load_pole_model") as compressed,
    ):
        root = Path(tmp)
        np.savez(root / "probes.npz", x=np.ones((1, 1)))
        pole.initialize(root, root, c, {})
        with (
            patch.object(
                pole.core, "outgoing_map", return_value=(np.array([[1j]]), np.array([[1.0]]), {})
            ) as maps,
            patch.object(pole.core, "surface_matrix", return_value=np.array([[0.1]])) as surface,
        ):
            z = 0.033075 - 0.00525j
            pole.core.calculate(
                dict(
                    z_hartree=[z.real, z.imag],
                    setting=c["settings"][0],
                    amplitudes=True,
                    kind="derivative",
                    file="capture.npz",
                )
            )
            np.testing.assert_allclose(maps.call_args.args[4], 2 * z, atol=1e-15, rtol=0)
            np.testing.assert_allclose(surface.call_args.args[2], -1 + z, atol=1e-15, rtol=0)
            assert maps.call_args.kwargs["continued_amplitudes"] is True
        compressed.assert_not_called()


def test_box_and_owner_keep_exact_caps_no_network_and_readonly_inputs():
    c = resolved()
    with unittest.TestCase().assertRaisesRegex(ValueError, "search box"):
        pole.require_box(0.1 - 0.005j, c)
    cmd = command(Path("/root"), Path("/source"), Path("/out"), c)
    assert "--cpuset-cpus=0-7" in cmd and "--memory=8g" in cmd and "--network=none" in cmd
    assert [cmd[i + 1] for i, v in enumerate(cmd) if v == "-v"] == [
        "/root:/work:ro",
        "/source:/opt/qmodeling:ro",
        "/out:/out",
    ]
    assert "projects.ukrmol_co.stretched_dz_pole" in cmd[-1]
    compile(cmd[-1], "<dz-pole-owner>", "exec")


if __name__ == "__main__":
    functions = [
        fn
        for module in [test_outgoing_map, test_pole_reference, test_pole_character]
        for name, fn in vars(module).items()
        if name.startswith("test_") and callable(fn)
    ]
    functions += [
        fn for name, fn in list(globals().items()) if name.startswith("test_") and callable(fn)
    ]
    suite = unittest.TestSuite(unittest.FunctionTestCase(fn) for fn in functions)
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
