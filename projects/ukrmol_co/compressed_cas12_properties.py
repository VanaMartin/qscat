"""One matched first-five property export from saved compressed CAS12 eigenvectors.

Only DENPROP runs: no QC, target diagonalization, integral generation or scattering.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
from pathlib import Path


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def matched_deck(text):
    """Change only NTGTL; preserve the proved native CIDATA/congen set order."""
    from projects.ukrmol_co.a2_channel_preflight import replace_assignment

    active = "\n".join(line.split("!", 1)[0] for line in text.splitlines())
    expected = dict(
        ntgt=[8],
        nftsor=[710, 730, 711, 731, 712, 732, 713, 733],
        ntgtf=list(range(1, 9)),
        ntgtl=[6, 4, 5, 7, 5, 7, 3, 3],
    )
    for key, want in expected.items():
        values = re.findall(rf"\b{key}\s*=\s*([\d,\s]+)", active, re.I)
        if len(values) != 1 or list(map(int, re.findall(r"\d+", values[0]))) != want:
            raise ValueError(f"Changed original DENPROP assignment: {key}")
    return replace_assignment(text, "ntgtl", "5,5,5,5,5,5,5,5")


def inputs(root, source, contract):
    """Authenticate the completed64-root import instead of repeating it."""
    path = source / contract["input_inventory"]
    if digest(path) != contract["input_inventory_sha256"]:
        raise ValueError("Changed saved-input inventory")
    inventory = json.loads(path.read_text())
    for name, want in inventory["input_sha256"].items():
        if digest(root / name) != want["sha256"]:
            raise ValueError(f"Changed saved input: {name}")
    for name, want in contract["scientific_prerequisite_sha256"].items():
        if digest(source / name) != want:
            raise ValueError(f"Changed scientific prerequisite: {name}")
    contrast = json.loads((source / contract["contrast_result"]).read_text())
    verification = json.loads((source / contract["contrast_verification"]).read_text())
    target = root / inventory["target_directory"]
    original = json.loads((target / "result.json").read_text())
    imported = json.loads((root / contract["native_import_verification"]).read_text())
    if (
        contrast["model_insensitivity_verdict"] is not True
        or verification["model_insensitivity_verdict"] is not True
        or original["resources"]["exit_code"] != 0
        or imported["all_64_roots_pass"] is not True
    ):
        raise ValueError("Qualified contrast/native64-root prerequisites are required")
    config = json.loads((target / "config.json").read_text())
    for key, expected in contract["fixed_target_settings"].items():
        if config.get(key) != expected:
            raise ValueError(f"Changed target setting: {key}")
    return inventory, target


def compute(root, source, work, contract):
    from projects.ukrmol_co.a2_channel_preflight import SECTORS, property_inventory, stage

    report = dict(
        property_inventory_verdict=None,
        scientific_exit=2,
        stages=[],
        native_denprop_calls=0,
        qc_calls=0,
        target_diagonalizations=0,
        scattering_calls=0,
        route_holdout_s_matrices_analyzed=0,
        full_sweep_released=False,
        physical_co_accuracy_verdict=None,
    )
    try:
        inventory, target = inputs(root, source, contract)
        report["saved_inputs_verified_before"] = len(inventory["input_sha256"])
        geometry = target / "output/CO/geom1"
        names = ["fort.16", "fort.26", "co.molden", "co.casscf.chk"] + [
            f"fort.{unit}" for unit in [710, 730, 711, 731, 712, 732, 713, 733]
        ]
        copied = {}
        for name in names:
            shutil.copy2(geometry / name, work / name)
            copied[name] = digest(work / name)
        report["copied_inputs_sha256"] = copied
        deck = work / "target.denprop.inp"
        deck.write_text(matched_deck((geometry / "inputs/target.denprop.inp").read_text()))
        report["matched_deck_sha256"] = digest(deck)
        original_target = json.loads((target / "target.json").read_text())
        expected = {
            f"{spin}.{irrep}": [
                state["energy_hartree"]
                for state in original_target["states"]
                if state["spin"] == spin and state["irrep"] == irrep
            ]
            for spin, irrep in SECTORS
        }
        if not original_target["converged"] or any(
            len(values) != 5 for values in expected.values()
        ):
            raise ValueError("Expected exactly five already qualified roots per target sector")
        executable = Path("/opt/ukrmolp/bin.double/denprop")
        report["native_denprop_executable_sha256"] = digest(executable)
        report["native_denprop_calls"] = 1
        (work / "result.json").write_text(json.dumps(report, indent=2) + "\n")
        stage(
            work,
            report,
            "target.denprop",
            [str(executable)],
            deck,
            contract["resources"]["native_timeout_seconds"],
            os.environ.copy(),
        )
        if "has completed execution" not in (work / "target.denprop.out").read_text():
            raise ValueError("DENPROP completion marker missing")
        properties = property_inventory((work / "fort.24").read_text(), expected)
        difference = properties["ground_dipole_z_au"] - original_target["ground_dipole_au"][2]
        if abs(difference) > contract["dipole_budget_au"]:
            raise ValueError("Native ground dipole prerequisite failed")
        if any(digest(work / name) != want for name, want in copied.items()):
            raise ValueError("Copied original orbital/CI/integral bytes changed during DENPROP")
        report.update(
            property_inventory_verdict=True,
            scientific_exit=0,
            property_inventory=properties,
            ground_dipole_difference_au=difference,
            output_property_sha256=digest(work / "fort.24"),
        )
    except Exception as error:
        report.update(error=repr(error), scientific_exit=2)
    finally:
        try:
            inventory, _ = inputs(root, source, contract)
            report["saved_inputs_verified_after"] = len(inventory["input_sha256"])
        except Exception as error:
            report.update(
                integrity_error=repr(error), scientific_exit=2, property_inventory_verdict=None
            )
        (work / "result.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
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
    env = os.environ | {
        "TMPDIR": str(args.work / "scratch"),
        "PSI_SCRATCH": str(args.work / "scratch"),
    }
    return run_profiled(
        args.work,
        env,
        [
            sys.executable,
            "-m",
            "projects.ukrmol_co.compressed_cas12_properties",
            *sys.argv[1:],
            "--profiled",
        ],
    )


if __name__ == "__main__":
    raise SystemExit(main())
