"""Export explicit five-per-sector properties on the qualified fifty-component frame."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path

from projects.ukrmol_co.a2_channel_preflight import SECTORS, property_inventory, stage
from projects.ukrmol_co.fifty_import_verification import covered_sectors, digest

OLD_SELECTION = "  ntgtl = 7,6,6,8,6,8,3,6,"
NEW_SELECTION = "  ntgtl = 5,5,5,5,5,5,5,5,"


def selected_deck(text):
    """Change exactly the recorded selection; leave all other native tokens intact."""
    if text.count(OLD_SELECTION) != 1 or text.count("ntgtl") != 1:
        raise ValueError("Unexpected original property selection")
    return text.replace(OLD_SELECTION, NEW_SELECTION)


def run(root, work, contract):
    def verify():
        for name, want in contract["input_sha256"].items():
            if digest(root / name) != want:
                raise ValueError(f"Changed frozen property input: {name}")

    verify()
    frame = json.loads((root / contract["frame_result"]).read_text())
    owner = json.loads((root / contract["frame_execution"]).read_text())
    assert owner["source_commit"] == contract["frame_source_commit"]
    assert owner["original_owner_exit"] == 0
    assert [owner["steps"][0][k] for k in ("original_command_exit", "original_container_exit")] == [
        0,
        0,
    ]
    assert frame["new_frame_import_verdict"] is True
    assert frame["original_seed_preservation_passed"] is False
    reference = covered_sectors(frame, frame["checkpoint_sha256"])
    expected = {
        f"{spin}.{irrep}": reference[(0 if spin == "singlet" else 2, irrep)][:5]
        for spin, irrep in SECTORS
    }
    geometry = root / contract["geometry"]
    work.mkdir(exist_ok=False)
    copied = {}
    for unit in (16, 26, 710, 730, 711, 731, 712, 732, 713, 733):
        source = geometry / f"fort.{unit}"
        shutil.copy2(source, work / source.name)
        copied[source.name] = digest(source)
    original = (geometry / "inputs/target.denprop.inp").read_text()
    deck = work / "target.denprop.inp"
    deck.write_text(selected_deck(original))
    report = dict(
        stages=[],
        new_property_export_verdict=None,
        original_seed_preservation_passed=False,
        original_seed_verifier_and_owner_exits=[1, 1, 1],
        frame_checkpoint_sha256=frame["checkpoint_sha256"],
        reused_input_sha256=copied,
        scattering_ensemble_verdict=None,
    )
    stage(
        work,
        report,
        "target.denprop",
        ["/opt/ukrmolp/bin.double/denprop"],
        deck,
        contract["resources"]["native_timeout_seconds"],
        os.environ.copy(),
    )
    if "has completed execution" not in (work / "target.denprop.out").read_text():
        raise ValueError("DENPROP completion marker missing")
    inventory = property_inventory((work / "fort.24").read_text(), expected)
    assert inventory["target_count"] == 40
    dipole_error = (
        inventory["ground_dipole_z_au"] - frame["native_properties"]["ground_state_dipole_au"][2]
    )
    if abs(dipole_error) > 1e-5:
        raise ValueError("Selected40 ground dipole differs from the qualified frame")
    for name, want in copied.items():
        if digest(work / name) != want:
            raise ValueError(f"Copied native input changed: {name}")
    verify()
    report.update(
        new_property_export_verdict=True,
        property_inventory=inventory,
        ground_dipole_difference_au=dipole_error,
        selected_property_sha256=digest(work / "fort.24"),
        original_deck_sha256=digest(geometry / "inputs/target.denprop.inp"),
        selected_deck_sha256=digest(deck),
        native_engine_invocations=1,
        new_orbital_optimizations=0,
        new_target_diagonalizations=0,
        new_scattering_inner_solves=0,
        physical_cross_frame_channel_correspondence_verdict=None,
    )
    (work / "result.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    args = parser.parse_args()
    run(args.root, args.output, json.loads(args.contract.read_text()))


if __name__ == "__main__":
    main()
