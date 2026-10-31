"""Saved-preparation adapter for a declared finite audited-engine CAS12 B1/2048 solve."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

from projects.ukrmol_co.compressed_cas12_properties import digest
from projects.ukrmol_co.engine_baseline import runtime_gate
from projects.ukrmol_co.native_full_reconstruction import _model_from_interface
from projects.ukrmol_co.saved_matrices import read_interface
from projects.ukrmol_co.sparse_pilot import compute as native_pilot
from projects.ukrmol_co.sparse_pilot import target_map
from projects.ukrmol_co.sparse_scattering import mpi_boundary_deck, selected_deck


def verify(root, source, contract):
    for base, key in [(root, "input_sha256"), (source, "source_prerequisite_sha256")]:
        for name, want in contract[key].items():
            if digest(base / name) != want:
                raise ValueError(f"Changed inner prerequisite: {name}")
    prepared = root / contract["preparation_directory"]
    saved = json.loads((prepared / "result.json").read_text())
    public_path = (
        root / contract["preparation_receipt"]
        if "preparation_receipt" in contract
        else source / contract["preparation_verification"]
    )
    public = json.loads(public_path.read_text())
    if saved["preparation_verdict"] is not True or public["preparation_verdict"] is not True:
        raise ValueError("Independent completed preparation required")
    if (
        saved["scientific_exit"] != 0
        or saved.get("sectors", {"B1": saved.get("B1")})["B1"]["contracted_dimension"]
        != contract["dimension"]
    ):
        raise ValueError("Preparation dimension/exit differs from the costed case")
    baseline = json.loads((root / contract["engine_baseline_result"]).read_text())
    baseline_owner = json.loads((root / contract["engine_baseline_owner"]).read_text())
    if (
        baseline["matched_engine_baseline_verdict"] is not True
        or baseline["independent_raw_reading_passed"] is not True
        or baseline_owner["original_controller_exit"] != 0
    ):
        raise ValueError("Audited engine scientific baseline has not qualified")
    target = json.loads((prepared / "target.json").read_text())
    retained = contract.get("target_roots_per_sector", 5)
    if retained == 8 and not (
        public["case"] == saved["case"] == contract["case"]
        and public["public_reader_exit"] == 0
        and public["original_scientific_command_container_owner_exits"] == [0] * 4
        and public["roots_per_sector"] == [8] * 8
        and public["contract_sha256"] == contract["preparation_contract_sha256"]
        and public["moints_sha256"] == digest(prepared / "moints")
        and public["moints_sha256"] == digest(prepared / "fort.16")
        and public["csf_sha256"] == digest(prepared / "fort.721")
    ):
        raise ValueError("Own public64-state preparation receipt differs")
    if retained == 8:
        peer = json.loads((root / contract["peer_preparation_receipt"]).read_text())
        if not (
            peer["case"] == {"dz": "tz", "tz": "dz"}[contract["case"]]
            and peer["preparation_verdict"] is True
            and peer["public_reader_exit"] == 0
            and peer["contract_sha256"] == contract["preparation_contract_sha256"]
            and peer["original_scientific_command_container_owner_exits"] == [0] * 4
            and peer["roots_per_sector"] == [8] * 8
        ):
            raise ValueError("Both finite basis preparations must independently qualify")
    mapping = target_map((prepared / "fort.24").read_text(), target["states"], roots=retained)
    if mapping != contract["own_target_map"]:
        raise ValueError("Own matched contraction-to-property IDs changed")
    text = selected_deck(
        (prepared / "scattering.scatci.doublet.B1.inp").read_text(), 2048, 1e-12, 1000
    )
    expected = mpi_boundary_deck(
        text, f"ntarg = {8 * retained}, idtarg = {','.join(map(str, mapping))}, rmatr = 18,"
    )
    if expected != (source / contract["frozen_native_deck"]).read_text():
        raise ValueError("Resolved native2048-root input differs before execution")
    return prepared, saved


def preflight_view(directory, prepared, saved):
    """Translate paths only; original physical files remain read-only and unchanged."""
    directory.mkdir()
    geometry = directory / "output/CO/geom1"
    geometry.mkdir(parents=True)
    (geometry / "inputs").mkdir()
    (geometry / "doublet.B1").mkdir()
    (directory / "target.json").symlink_to(prepared / "target.json")
    (directory / "preflight-result.json").write_text(
        json.dumps(dict(original_exit=0, sectors=saved.get("sectors", {"B1": saved.get("B1")})))
        + "\n"
    )
    for name in ["fort.24", "fort.26"]:
        (geometry / name).symlink_to(prepared / name)
    (geometry / "prop.out").symlink_to(prepared / "fort.24")
    (geometry / "inputs/scattering.scatci.doublet.B1.inp").symlink_to(
        prepared / "scattering.scatci.doublet.B1.inp"
    )
    for name in ["fort.16", "fort.721"]:
        (geometry / "doublet.B1" / name).symlink_to(prepared / name)
    return directory


def export_gate(directory, own, expected_layout=(40, 250)):
    """Strict declared native assembly, not a sparse-ID remapping framework."""
    if expected_layout not in [(40, 250), (64, 400)]:
        raise ValueError("Undeclared CAS12 target/channel layout")
    targets, channels = expected_layout
    interface = read_interface(directory / "fort.10", directory / "fort.21")
    if [label[0] for label in interface["target_labels"]] != list(range(1, targets + 1)):
        raise ValueError("Unexpected native selected-target row IDs")
    expected = [
        float(line.split()[8].replace("D", "E"))
        for line in own.read_text().splitlines()
        if line.split()[:1] == ["5"]
    ]
    np.testing.assert_allclose(interface["target_energies_hartree"], expected, atol=1e-9, rtol=0)
    model = _model_from_interface(directory, interface, expected_layout)
    arrays = {
        name: dict(shape=list(model[name].shape), finite=bool(np.isfinite(model[name]).all()))
        for name in ["poles", "amplitudes", "c2", "c3", "thresholds", "saved"]
    }
    if not all(row["finite"] for row in arrays.values()):
        raise ValueError("Nonfinite own model arrays")
    return dict(
        **{f"strict{targets}_target_{channels}_channel_assembly_verdict": True},
        arrays=arrays,
        ground_hartree=model["ground"],
        open_labels=interface["channel_labels"][:6],
        native_boundary_poles=interface["boundary_poles"],
        own_model_interface=interface,
        scope="Own W/D/CF/K assembly; no S construction, propagation, pole or Pi pairing",
    )


def compute(root, source, work, contract):
    report = dict(
        inner_export_verdict=None,
        scientific_exit=2,
        sectors=["B1"],
        roots=2048,
        mpi_ranks=12,
        pi_pair_verdict=None,
        pole_verdict=None,
        electronic_model_verdict=None,
        spectral_truncation_verdict=None,
        full_sweep_released=False,
        new_qc_calls=0,
        new_target_diagonalizations=0,
        new_property_exports=0,
        new_preparation_calls=0,
        new_propagations=0,
        propagation_counter_scope=(
            "Complex/Python maps only; native RSOLVE is a separate recorded stage"
        ),
        route_holdout_s_matrices_analyzed=0,
    )

    def save():
        (work / "result.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    try:
        prepared, saved = verify(root, source, contract)
        report["inputs_verified_before"] = len(contract["input_sha256"])
        report["runtime_gate"] = runtime_gate(contract)
        save()
        view = preflight_view(work / "preflight-view", prepared, saved)
        native = work / "native"
        native.mkdir()
        args = argparse.Namespace(
            preflight=view,
            outer_inputs=source / contract["outer_templates"],
            workdir=native,
            roots=2048,
            ranks=12,
            stage_timeout=contract["resources"]["native_timeout_seconds"],
            outer_stage_timeout=contract["resources"]["outer_timeout_seconds"],
            telemetry=root / contract["telemetry_directory"] / "scattering_telemetry.so",
        )
        native_pilot(
            args,
            sectors=("B1",),
            target_roots_per_sector=contract.get("target_roots_per_sector", 5),
            boundary_channels=contract["boundary_channels"],
        )
        report["native_report"] = json.loads((native / "pilot.json").read_text())
        for name, want in report["native_report"]["sectors"]["B1"]["input_sha256"].items():
            if digest(native / "B1" / name) != want:
                raise ValueError(f"Native copied input changed: {name}")
        if digest(native / "B1/scatci.inp") != digest(source / contract["frozen_native_deck"]):
            raise ValueError("Actual native solver deck differs from the frozen input")
        report["copied_native_inputs_unchanged"] = True
        report["export_gate"] = export_gate(
            native / "B1",
            prepared / "fort.24",
            (8 * contract.get("target_roots_per_sector", 5), contract["boundary_channels"]),
        )
        report.update(inner_export_verdict=True, scientific_exit=0)
    except Exception as error:
        report["error"] = repr(error)
        path = work / "native/pilot.json"
        if path.is_file():
            report["native_report"] = json.loads(path.read_text())
    finally:
        try:
            verify(root, source, contract)
            report["inputs_verified_after"] = len(contract["input_sha256"])
        except Exception as error:
            report.update(integrity_error=repr(error), scientific_exit=2, inner_export_verdict=None)
        save()
    return report["scientific_exit"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["root", "source", "work", "contract"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--profiled", action="store_true")
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text())
    if args.profiled:
        return compute(args.root, args.source, args.work, contract)
    from projects.ukrmol_co.run import run_profiled

    args.work.mkdir(exist_ok=False)
    (args.work / "scratch").mkdir()
    env = os.environ | dict(
        TMPDIR=str(args.work / "scratch"), PSI_SCRATCH=str(args.work / "scratch")
    )
    return run_profiled(
        args.work,
        env,
        [
            sys.executable,
            "-m",
            "projects.ukrmol_co.compressed_cas12_inner",
            *sys.argv[1:],
            "--profiled",
        ],
    )


if __name__ == "__main__":
    raise SystemExit(main())
