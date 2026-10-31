"""Known-truth64/400 adapter controls; never invokes a native or quantum program."""

import argparse
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co import compressed_cas12_inner as inner
from projects.ukrmol_co import native_full_reconstruction as reconstruction
from projects.ukrmol_co import sparse_pilot
from projects.ukrmol_co.compressed_cas12_inner_owner import container_command


def test_native_64_target_map_keeps_own_energy_order_and_eight_per_sector():
    labels = [(s, i) for i in ("A1", "B1", "B2", "A2") for s in ("singlet", "triplet")]
    states = [
        dict(spin=s, irrep=i, energy_hartree=-113 + 0.01 * (8 * r + k))
        for r in range(8)
        for k, (s, i) in enumerate(labels)
    ]
    prop = "\n".join(
        f"5 {n + 1} 0 0 0 0 0 0 {s['energy_hartree']:.12f} {s['spin']} {s['irrep']}"
        for n, s in enumerate(states)
    )
    assert sparse_pilot.target_map(prop, states, roots=8) == [
        8 * r + k + 1 for k in range(8) for r in range(8)
    ]
    with unittest.TestCase().assertRaises(ValueError):
        sparse_pilot.target_map(prop, states)


def test_64_400_model_assembly_has_known_shapes_thresholds_and_angular_diagonal():
    targets = -113 + 0.3 * np.arange(64)
    target_rows = np.r_[np.ones(6, dtype=int), 2 + np.arange(394) % 63]
    angular = np.arange(400) % 5
    interface = dict(
        target_components=64,
        boundary_channels=400,
        symmetry_charge=[0, 0, 0, 0],
        geometry_mass=[0, 0],
        boundary_poles=2048,
        multipole_orders=2,
        radius_bohr=18.0,
        channel_labels=np.column_stack([np.arange(400), target_rows, angular]).tolist(),
        channel_thresholds_rydberg=(2 * (targets[target_rows - 1] - targets[0])).tolist(),
        target_energies_hartree=targets.tolist(),
    )
    header = {
        k: interface[k]
        for k in ["target_components", "boundary_channels", "symmetry_charge", "geometry_mass"]
    }
    raw = [b""] * 5 + [
        np.zeros(400 * 401).tobytes(),
        np.linspace(-114, -110, 2048).tobytes(),
        np.zeros(2048 * 400).tobytes(),
    ]
    with (
        patch.object(reconstruction, "unit", return_value="fort.921"),
        patch.object(
            reconstruction,
            "read_k",
            return_value=(0.0735 * (0.1 + 0.05 * np.arange(99)), np.zeros((99, 6, 6)), header),
        ),
        patch.object(reconstruction, "records", return_value=raw),
    ):
        model = reconstruction._model_from_interface(Path("/synthetic"), interface, (64, 400))
        assert model["amplitudes"].shape == (2048, 400) and model["c2"].shape == (400, 400)
        np.testing.assert_allclose(np.diag(model["c2"]), angular * (angular + 1), atol=0, rtol=0)
        np.testing.assert_allclose(model["c3"], 0, atol=0, rtol=0)
        interface["channel_thresholds_rydberg"][6] = 0.01
        with unittest.TestCase().assertRaises((ValueError, AssertionError)):
            reconstruction._model_from_interface(Path("/synthetic"), interface, (64, 400))


def test_undeclared_layout_fails_before_any_native_or_file_access():
    for retained, boundary in [(8, 250), (5, 400), (9, 450)]:
        with unittest.TestCase().assertRaises(ValueError):
            sparse_pilot.compute(
                argparse.Namespace(),
                sectors=("B1",),
                target_roots_per_sector=retained,
                boundary_channels=boundary,
            )


def test_64_400_export_uses_strict_parser_and_own_thresholds():
    with TemporaryDirectory() as tmp:
        p = Path(tmp)
        own = p / "fort.24"
        own.write_text(
            "\n".join(f"5 {i + 1} 0 0 0 0 0 0 {-113 + 0.01 * i:.12f} Singlet A1" for i in range(64))
        )
        interface = dict(
            target_labels=[[i + 1, 0, 0, 0] for i in range(64)],
            target_energies_hartree=[-113 + 0.01 * i for i in range(64)],
            channel_labels=[[1, 1, 1, 1]] * 6,
            boundary_poles=2048,
        )
        model = {
            n: np.zeros((2, 2)) for n in ["poles", "amplitudes", "c2", "c3", "thresholds", "saved"]
        }
        model["ground"] = -113.0
        with (
            patch.object(inner, "read_interface", return_value=interface),
            patch.object(inner, "_model_from_interface", return_value=model) as assemble,
        ):
            result = inner.export_gate(p, own, (64, 400))
        assert result["strict64_target_400_channel_assembly_verdict"]
        assert assemble.call_args.args[-1] == (64, 400)


def test_stretched_saved_view_retains_actual_b1_dimensions_and_unchanged_inputs():
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        saved = root / "saved"
        saved.mkdir()
        for n in [
            "target.json",
            "fort.24",
            "fort.26",
            "fort.16",
            "fort.721",
            "scattering.scatci.doublet.B1.inp",
        ]:
            (saved / n).write_bytes(b"original")
        view = inner.preflight_view(
            root / "view",
            saved,
            dict(B1=dict(contracted_dimension=87198, contracted_continuum_dimension=2256)),
        )
        assert (
            json.loads((view / "preflight-result.json").read_text())["sectors"]["B1"][
                "contracted_dimension"
            ]
            == 87198
        )
        assert not (view / "output/CO/geom1/doublet.B2").exists()
        assert all(p.read_bytes() == b"original" for p in saved.iterdir())


def test_64_target_native_adapter_runs_b1_only_and_audits_400_channels():
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        view = root / "view"
        geometry = view / "output/CO/geom1"
        (geometry / "inputs").mkdir(parents=True)
        (geometry / "doublet.B1").mkdir()
        (view / "target.json").write_text(json.dumps(dict(states=[])))
        (view / "preflight-result.json").write_text(
            json.dumps(
                dict(
                    original_exit=0,
                    sectors=dict(
                        B1=dict(contracted_dimension=87198, contracted_continuum_dimension=2256)
                    ),
                )
            )
        )
        for n in ["prop.out", "fort.24", "fort.26"]:
            (geometry / n).write_text("saved input")
        for n in ["fort.16", "fort.721"]:
            (geometry / "doublet.B1" / n).write_text("saved input")
        template = (
            Path(__file__).with_name("compressed-cas12-preparation-inputs-20261026")
            / "scattering.scatci.doublet.B1.inp"
        )
        (geometry / "inputs/scattering.scatci.doublet.B1.inp").write_bytes(template.read_bytes())
        outer = root / "outer"
        outer.mkdir()
        for stage in ["rsolve", "eigenp", "reson"]:
            (outer / f"scattering.{stage}.doublet.B1.inp").write_text("&input\n/\n")
        work = root / "work"
        work.mkdir()
        args = argparse.Namespace(
            preflight=view,
            workdir=work,
            outer_inputs=outer,
            roots=2048,
            ranks=12,
            stage_timeout=43200,
            telemetry=root / "telemetry.so",
        )
        with (
            patch.object(sparse_pilot, "target_map", return_value=list(range(1, 65))) as mapping,
            patch.object(
                sparse_pilot, "audit_sector", return_value=dict(export_checks_passed=True)
            ) as audit,
            patch.object(sparse_pilot.subprocess, "Popen") as native,
        ):
            native.return_value.returncode = 0
            sparse_pilot.compute(
                args, sectors=("B1",), target_roots_per_sector=8, boundary_channels=400
            )
        assert native.call_count == 4 and mapping.call_args.kwargs["roots"] == 8
        assert audit.call_args.kwargs["boundary_channels"] == 400
        assert "ntarg = 64" in (work / "B1/scatci.inp").read_text()
        report = json.loads((work / "pilot.json").read_text())
        assert report["pi_check_passed"] is None and not (work / "B2").exists()


def test_separate_stretched_owner_changes_contract_identity_not_engine_or_resources():
    contract = dict(
        engine_source_directory="audited",
        image_id="sha256:pinned",
        container_name="co-stretched-dz-b1-2048-20261031",
        runtime_controls_module="projects.ukrmol_co.test_stretched_basis_inner",
    )
    command = container_command(
        Path("/root"),
        Path("/source"),
        Path("/out"),
        contract,
        "projects/ukrmol_co/stretched-dz-inner-contract-20261031.json",
    )
    assert "--cpuset-cpus=0-11" in command and "--memory=64g" in command
    assert contract["container_name"] in command and "stretched-dz-inner-contract" in command[-1]
    assert contract["runtime_controls_module"] in command[-1]
    compile(command[-1], "<stretched-inner-owner>", "exec")


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
