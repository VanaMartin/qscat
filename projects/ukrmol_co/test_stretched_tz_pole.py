"""Own-TZ receipt/seed routing and unit controls, plus retained analytic oracles."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co import stretched_dz_pole as pole
from projects.ukrmol_co import test_stretched_dz_pole as dz_controls

CONTRACT = "projects/ukrmol_co/stretched-tz-pole-execution-20261102.json"


def resolved():
    source = Path(__file__).resolve().parents[2]
    return pole.resolve(source, json.loads((source / CONTRACT).read_text()))


def test_tz_uses_own_seed_and_exact_shared_numerical_allocation_with_complete_source():
    c = resolved()
    source = Path(__file__).resolve().parents[2]
    for name, want in c["source_sha256"].items():
        assert pole.core.digest(source / name) == want, name
    dz = dz_controls.resolved()
    for name in [
        "settings",
        "near_feature_real_settings",
        "root_evaluation_caps",
        "character",
        "resources",
        "same_root_position_full_width_budget_ev",
        "contour_radius_hartree",
        "derivative_steps_hartree",
    ]:
        assert c[name] == dz[name], name
    assert c["case"] == "stretched/cas12/tz/B1" and c["near_feature_real_index"] == 18
    assert c["seeds_incident_hartree"] == [[0.03675, -0.00606], [0.03677, -0.00606]]
    assert len(c["input_sha256"]) == 8 and len(c["source_sha256"]) == 10


def test_mismatched_own_receipt_and_failed_public_dz_stop_before_model_or_maps():
    source = Path(__file__).resolve().parents[2]
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        c = resolved()
        c["source_sha256"] = {}
        c["input_sha256"] = {}
        c["model"]["receipt"] = "own.json"
        c["completion_reading"] = "real.json"
        c["seed_reading"] = "seed.json"
        c["dz_pole_reading"] = "dz.json"
        receipt = dict(
            case="dz",
            inner_export_verdict=True,
            public_reader_exit=0,
            original_exits=[0] * 6,
            contract_sha256=c["model"]["inner_contract_sha256"],
            roots=2048,
            target_count=64,
            boundary_channels=400,
            spectral_truncation_verdict=None,
            full_sweep_released=False,
        )
        real = dict(
            reader_exit=0,
            result=dict(
                combined_20degree_real_prerequisite_verdict=True,
                original_exits=[0] * 6,
                zero_map_failure_original_exits=[2] * 6,
                original40degree_scientific_exit=2,
                complex_use_released=False,
            ),
        )
        seed = dict(
            result=dict(
                real_prerequisite_verdict=True,
                pole_measurement=False,
                cases=dict(tz=dict(index=18)),
            )
        )
        for name, value in [
            ("own.json", receipt),
            ("real.json", real),
            ("seed.json", seed),
            ("dz.json", dict(reader_exit=2, result={})),
        ]:
            (root / name).write_text(json.dumps(value))
        with patch.object(pole, "load_own") as model, patch.object(pole.core, "calculate") as calc:
            assert pole.run(root, source, root, c) == 2
            receipt["case"] = "tz"
            (root / "own.json").write_text(json.dumps(receipt))
            assert pole.run(root, source, root, c) == 2
            model.assert_not_called()
            calc.assert_not_called()


def test_tz_continues_own_incident_rydberg_and_neutral_hartree_reference():
    c = resolved()
    model = dict(
        ground=-2.0,
        poles=np.array([0.1]),
        amplitudes=np.array([[0.2]]),
        thresholds=np.array([0.0]),
        angular=np.array([0]),
        c2=np.zeros((1, 1)),
        c3=np.zeros((1, 1)),
    )
    with TemporaryDirectory() as tmp, patch.object(pole, "load_own", return_value=model):
        root = Path(tmp)
        np.savez(root / "probes.npz", x=np.ones((1, 1)))
        pole.initialize(root, root, c, {})
        with (
            patch.object(
                pole.core, "outgoing_map", return_value=(np.array([[1j]]), np.array([[1.0]]), {})
            ) as maps,
            patch.object(pole.core, "surface_matrix", return_value=np.array([[0.1]])) as surface,
        ):
            z = complex(*c["seeds_incident_hartree"][0])
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
            np.testing.assert_allclose(
                maps.call_args.kwargs["reference"], 0.0735, atol=1e-15, rtol=0
            )
            np.testing.assert_allclose(surface.call_args.args[2], -2 + z, atol=1e-15, rtol=0)


if __name__ == "__main__":
    functions = [
        fn
        for module in [
            dz_controls.test_outgoing_map,
            dz_controls.test_pole_reference,
            dz_controls.test_pole_character,
            dz_controls,
        ]
        for name, fn in vars(module).items()
        if name.startswith("test_") and callable(fn)
    ]
    functions += [
        fn for name, fn in list(globals().items()) if name.startswith("test_") and callable(fn)
    ]
    suite = unittest.TestSuite(unittest.FunctionTestCase(fn) for fn in functions)
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
