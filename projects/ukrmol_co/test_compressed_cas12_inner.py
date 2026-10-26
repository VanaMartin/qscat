"""Bounded B1 routing/strict-export controls without a quantum or native call."""

import argparse
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co import compressed_cas12_inner as inner
from projects.ukrmol_co import sparse_pilot
from projects.ukrmol_co.compressed_cas12_inner_owner import container_command


def test_invalid_single_sector_requests_stop_before_any_file_or_native_access():
    for sectors in [("B2",), ("B2", "B1"), ("B1", "B1")]:
        with unittest.TestCase().assertRaises(ValueError):
            sparse_pilot.compute(argparse.Namespace(), sectors=sectors)


def test_existing_pilot_b1_only_never_runs_b2_or_claims_pi_pairing():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        geometry = root / "preflight/output/CO/geom1"
        geometry.mkdir(parents=True)
        (geometry / "inputs").mkdir()
        (geometry / "doublet.B1").mkdir()
        work = root / "work"
        work.mkdir()
        (root / "preflight/preflight-result.json").write_text(
            json.dumps(
                dict(
                    original_exit=0,
                    sectors=dict(
                        B1=dict(contracted_dimension=86352, contracted_continuum_dimension=1410)
                    ),
                )
            )
        )
        (root / "preflight/target.json").write_text(json.dumps(dict(states=[])))
        (geometry / "prop.out").write_text("mocked property map")
        for name in ["fort.24", "fort.26"]:
            (geometry / name).write_bytes(b"native saved input")
        for name in ["fort.16", "fort.721"]:
            (geometry / "doublet.B1" / name).write_bytes(b"same-frame native input")
        template = (
            Path(__file__).parent
            / "compressed-cas12-preparation-inputs-20261026/scattering.scatci.doublet.B1.inp"
        )
        (geometry / "inputs/scattering.scatci.doublet.B1.inp").write_bytes(template.read_bytes())
        outer = root / "outer"
        outer.mkdir()
        for name in ["rsolve", "eigenp", "reson"]:
            (outer / f"scattering.{name}.doublet.B1.inp").write_text("&input\n/\n")
        args = argparse.Namespace(
            preflight=root / "preflight",
            workdir=work,
            outer_inputs=outer,
            roots=2048,
            ranks=12,
            stage_timeout=43200,
            telemetry=root / "telemetry.so",
        )
        with (
            patch.object(sparse_pilot, "target_map", return_value=list(range(1, 41))),
            patch.object(
                sparse_pilot, "audit_sector", return_value=dict(raw_export_passed=True)
            ) as audit,
            patch.object(sparse_pilot.subprocess, "Popen") as native,
        ):
            native.return_value.returncode = 0
            sparse_pilot.compute(args, sectors=("B1",))
        assert native.call_count == 4 and audit.call_count == 1
        first = native.call_args_list[0].args[0]
        assert first[:3] == ["mpiexec", "-n", "12"]
        report = json.loads((work / "pilot.json").read_text())
        assert list(report["sectors"]) == ["B1"] and report["all_stages_complete"]
        assert report["pi_check_passed"] is None and not (work / "B2").exists()


def test_failed_preparation_never_calls_runtime_or_native():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        with (
            patch.object(inner, "verify", side_effect=ValueError("preparation unqualified")),
            patch.object(inner, "runtime_gate") as runtime,
            patch.object(inner, "native_pilot") as native,
        ):
            exit_code = inner.compute(root, root, root, {})
        runtime.assert_not_called()
        native.assert_not_called()
        report = json.loads((root / "result.json").read_text())
        assert exit_code == 2 and report["inner_export_verdict"] is None
        assert report["pi_pair_verdict"] is None and report["pole_verdict"] is None


def test_saved_view_changes_paths_only_and_contains_no_b2_input():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        saved = root / "saved"
        saved.mkdir()
        for name in [
            "target.json",
            "fort.24",
            "fort.26",
            "fort.16",
            "fort.721",
            "scattering.scatci.doublet.B1.inp",
        ]:
            (saved / name).write_bytes(b"saved immutable input")
        before = {p.name: p.read_bytes() for p in saved.iterdir()}
        view = inner.preflight_view(
            root / "view", saved, dict(sectors=dict(B1=dict(contracted_dimension=86352)))
        )
        assert (view / "target.json").is_symlink()
        assert not (view / "output/CO/geom1/doublet.B2").exists()
        assert before == {p.name: p.read_bytes() for p in saved.iterdir()}


def test_export_gate_reuses_strict_reader_without_sparse_id_selection():
    with TemporaryDirectory() as temporary:
        root = Path(temporary)
        own = root / "fort.24"
        own.write_text(
            "\n".join(f"5 {i + 1} 0 0 0 0 0 0 {-100 + i * 0.1:.15g} Singlet A1" for i in range(40))
        )
        interface = dict(
            target_labels=[[i + 1, 0, 0, 0] for i in range(40)],
            target_energies_hartree=[-100 + i * 0.1 for i in range(40)],
            channel_labels=[[1, 1, 1, 1]] * 6,
            boundary_poles=2048,
        )
        model = {
            name: np.ones((2, 2))
            for name in ["poles", "amplitudes", "c2", "c3", "thresholds", "saved"]
        }
        model["ground"] = -100
        with (
            patch.object(inner, "read_interface", return_value=interface) as strict,
            patch.object(inner, "_model_from_interface", return_value=model) as assemble,
        ):
            result = inner.export_gate(root, own)
        strict.assert_called_once_with(root / "fort.10", root / "fort.21")
        assert assemble.call_args.args[-1] == (40, 250)
        assert result["strict40_target_250_channel_assembly_verdict"]


def test_owner_uses_twelve_physical_cpus_and_read_only_audited_engine():
    command = container_command(
        Path("/campaign"),
        Path("/frozen"),
        Path("/new"),
        dict(engine_source_directory="audited", image_id="sha256:pinned"),
    )
    assert (
        "--cpuset-cpus=0-11" in command and "--memory=64g" in command and "--shm-size=8g" in command
    )
    mounts = [command[i + 1] for i, word in enumerate(command) if word == "-v"]
    assert mounts == [
        "/campaign:/work:ro",
        "/frozen:/opt/qmodeling:ro",
        "/new:/out",
        "/campaign/audited/build/install/bin:/opt/ukrmolp/bin.double:ro",
        "/campaign/audited/build/install/lib:/opt/ukrmolp/lib.double:ro",
    ]
    compile(command[-1], "<inner-owner>", "exec")
    assert "test_compressed_cas12_inner" in command[-1]


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
