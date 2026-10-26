"""CAS12 model dispatch and staged known-truth contrast controls, no CO propagation."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch

import numpy as np

from projects.ukrmol_co import compressed_a2_pole_contrast as adapter
from projects.ukrmol_co import compressed_cas12_pole_contrast as owner
from projects.ukrmol_co import compressed_pole_reference as pole
from projects.ukrmol_co.compressed_cas12_pole_owner import container_command
from projects.ukrmol_co.pole_character import screen_witnesses
from projects.ukrmol_co.test_compressed_a2_pole_contrast import BASE, CHARACTER, PROPOSAL


def screen():
    r = np.zeros((6, 6), complex)
    r[0, 0] = 0.1j
    return screen_witnesses([r] * 4, CHARACTER["open_labels"]["B1"], 0.08, CHARACTER["character"])


def staged(kind):
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        (root / "outputs").mkdir()
        output = root / "outputs/check"
        execution = dict(
            resources=dict(
                cpuset_cpus="0-3", threads_per_worker=1, contour_workers=4, memory_gib=8
            ),
            raw_directory="raw",
            source_sha256={},
            input_sha256={},
            output_relative="outputs",
            per_case_hard_cap_seconds=1800,
            total_map_cap=79,
            inner_result="native.json",
            cas11_five_reference="reference.json",
        )
        for name, value in [
            ("base_pole_contract", BASE),
            ("numerical_proposal", PROPOSAL),
            ("character_policy", CHARACTER),
            ("inner_verification", dict(inner_export_verdict=True, array_sha256={})),
        ]:
            execution[name] = name + ".json"
            (root / execution[name]).write_text(json.dumps(value))
        (root / "native.json").write_text(
            json.dumps(dict(inner_export_verdict=True, scientific_exit=0))
        )
        z = 0.13859147623159174 - 0.03838640724468573j
        (root / "reference.json").write_text(
            json.dumps(
                dict(
                    pole_incident_hartree=[z.real, z.imag],
                    disposition="pass",
                    completed=True,
                    screen=screen(),
                )
            )
        )

        def guard(model, directory, case, proposal, report):
            report["real_maps_reserved"] = 1 if kind == "guard-failure" else 2
            if kind == "guard-failure":
                raise ValueError("failed real guard")
            path = directory / "guard.json"
            path.write_text("{}")
            return path

        def child(*args, **kwargs):
            directory = output / "pole"
            directory.mkdir()
            actual = z + (0.004 if kind == "model-reject" else 0)
            (directory / "result.json").write_text(
                json.dumps(
                    dict(
                        scientific_exit=1 if kind == "pole-failure" else 0,
                        finite_reference_pole_verdict=kind != "pole-failure",
                        pole_incident_hartree=[actual.real, actual.imag],
                        evaluations=[{}] * 40,
                    )
                )
            )
            for name in ["moments-0", "moments-1", "derivative-0", "derivative-1"]:
                r = np.zeros((6, 6), complex)
                r[0, 0] = 0.1j
                np.savez(directory / (name + ".npz"), residue=r)
            return Mock(returncode=1 if kind == "pole-failure" else 0)

        with (
            patch.object(owner, "verify_inputs", return_value=123),
            patch.object(
                owner,
                "load_pole_model",
                return_value=dict(interface=dict(channel_labels=CHARACTER["open_labels"]["B1"])),
            ),
            patch.object(owner, "real_guard", side_effect=guard),
            patch.object(owner, "require_real_guard"),
            patch.object(owner.subprocess, "Popen", side_effect=child) as process,
        ):
            code = owner.run(root, root, output, execution)
        return code, json.loads((output / "result.json").read_text()), process.call_count


def test_failed_real_guard_starts_no_complex_process():
    code, report, calls = staged("guard-failure")
    assert calls == 0 and code == 2 and report["model_insensitivity_verdict"] is None
    assert report["new_maps_reserved"] == 1


def test_native_pole_failure_is_not_model_sensitivity_rejection():
    code, report, calls = staged("pole-failure")
    assert calls == 1 and code == 2 and report["model_insensitivity_verdict"] is None
    assert report["original_pole_command_exit"] == report["original_pole_scientific_exit"] == 1


def test_known_truth_pass_and_resolved_rejection_preserve_unset_pi_pair():
    for kind, expected_code, verdict in [("pass", 0, True), ("model-reject", 1, False)]:
        code, report, calls = staged(kind)
        assert code == expected_code and report["model_insensitivity_verdict"] is verdict
        assert calls == 1 and report["pi_pair_verdict"] is None and report["native_calls"] == 0
        assert report["new_maps_reserved"] == 42


def test_resolved_seed_is_saved_cas11_pole_and_numerics_are_unchanged():
    ref = dict(pole_incident_hartree=[0.13859147623159174, -0.03838640724468573])
    c = owner.resolved_contract(
        dict(
            raw_directory="overlay",
            input_sha256={},
            resources=dict(
                cpuset_cpus="0-3", threads_per_worker=1, contour_workers=4, memory_gib=8
            ),
        ),
        BASE,
        PROPOSAL,
        ref,
        "guard",
        "sha",
    )
    assert (
        c["case"] == "compressed/cas12/B1"
        and c["seeds_incident_hartree"][0] == ref["pole_incident_hartree"]
    )
    assert (
        c["settings"] == PROPOSAL["settings"]
        and c["root_evaluation_caps"] == PROPOSAL["root_caps"]["B1"]
    )
    assert (
        c["resources"]["total_map_cap"] == 77 and c["model_view"] == "compressed-cas12-native40-v1"
    )


def test_cas12_dispatch_is_explicit_strict40_not_sparse_id_mapping():
    source = Path("/saved/compressed/cas12/B1")
    with (
        patch.object(adapter, "read_interface", return_value={}) as strict,
        patch.object(adapter, "_model_from_interface", return_value={}) as assembly,
    ):
        model = adapter.load_pole_model(
            source, dict(case="compressed/cas12/B1", model_view="compressed-cas12-native40-v1")
        )
    strict.assert_called_once_with(source / "fort.10", source / "fort.21")
    assert assembly.call_args.args[-1] == (40, 250)
    assert "interface" in model
    with unittest.TestCase().assertRaises(ValueError):
        adapter.load_pole_model(source, dict(case="compressed/cas12/B1"))


def test_existing_pole_cli_requires_cas12_guard_before_loading_any_model():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        (root / "raw").mkdir()
        (root / "raw/input-manifest.json").write_text(
            json.dumps({str(i): "sha" for i in range(123)})
        )
        (root / "canary.json").write_text(
            json.dumps(dict(canary_verdict=True, rows=[dict(checks=dict(pass_gate=True))] * 11))
        )
        contract = dict(
            raw_directory="raw",
            case="compressed/cas12/B1",
            claim_boundary="test",
            input_sha256={},
            canary_result="canary.json",
        )
        with (
            patch.object(pole, "digest", return_value="sha"),
            patch.object(pole.subprocess, "run", return_value=Mock(returncode=0)),
            patch.object(
                pole, "require_real_guard", side_effect=ValueError("missing measured guard")
            ) as guard,
            patch.object(pole, "load_pole_model") as load,
        ):
            code = pole.run(root, root / "out", contract)
        guard.assert_called_once()
        load.assert_not_called()
        assert code == 2


def test_small_job_owner_uses_four_cores_and_no_writable_input_mount():
    command = container_command(
        Path("/campaign"), Path("/frozen"), Path("/new"), dict(image="sha256:pinned")
    )
    assert "--cpuset-cpus=0-3" in command and "--memory=8g" in command
    assert [command[i + 1] for i, value in enumerate(command) if value == "-v"] == [
        "/campaign:/work:ro",
        "/frozen:/opt/qmodeling:ro",
        "/new:/out",
    ]
    compile(command[-1], "<cas12-pole-owner>", "exec")
    assert (
        "test_compressed_cas12_pole_contrast" in command[-1]
        and "compressed_cas12_inner" not in command[-1]
    )


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
