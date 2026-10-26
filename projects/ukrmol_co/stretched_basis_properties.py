"""One missing64-state DENPROP table on a pinned stretched target; never run QC."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
from pathlib import Path

CONTRACT = "projects/ukrmol_co/stretched-basis-properties-contract-20261026.json"
COPY = ["fort.16", "fort.26", "co.molden", "co.casscf.chk"] + [
    f"fort.{unit}" for unit in [710, 730, 711, 731, 712, 732, 713, 733]
]


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path, value):
    partial = path.with_suffix(".partial")
    partial.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    partial.replace(path)


def matched_deck(text, original_counts):
    """Change only NTGTL in the exact saved eight-set property prescription."""
    from projects.ukrmol_co.a2_channel_preflight import replace_assignment

    active = "\n".join(line.split("!", 1)[0] for line in text.splitlines())
    required = dict(
        ntgt=[8],
        nftsor=[710, 730, 711, 731, 712, 732, 713, 733],
        ntgtf=list(range(1, 9)),
        ntgtl=original_counts,
        nftint=[16],
        nftg=[26],
        nftmt=[24],
    )
    for key, want in required.items():
        values = re.findall(rf"\b{key}\s*=\s*([\d,\s]+)", active, re.I)
        if len(values) != 1 or list(map(int, re.findall(r"\d+", values[0]))) != want:
            raise ValueError("Changed original DENPROP assignment: " + key)
    return replace_assignment(text, "ntgtl", "8,8,8,8,8,8,8,8")


def predecessor(root, case, contract_sha256, receipt=None, receipt_sha256=None):
    """TZ is a separate launch after an independently verified DZ export, not a queue."""
    if case == "dz":
        if receipt is not None or receipt_sha256 is not None:
            raise ValueError("DZ does not accept a predecessor override")
        return None
    if not receipt or not receipt_sha256:
        raise ValueError("TZ requires a separately pinned independent DZ receipt")
    path = root / receipt
    if not path.resolve().is_relative_to(root.resolve()) or digest(path) != receipt_sha256:
        raise ValueError("Changed or out-of-root DZ receipt")
    value = json.loads(path.read_text())
    if not (
        value.get("case") == "dz"
        and value.get("contract_sha256") == contract_sha256
        and value.get("independent_property_inventory_verdict") is True
        and value.get("target_count") == 64
        and value.get("roots_per_sector") == [8] * 8
        and value.get("original_scientific_command_container_owner_exits") == [0] * 4
        and value.get("public_reader_exit") == 0
        and re.fullmatch(r"[0-9a-f]{64}", value.get("property_sha256", ""))
    ):
        raise ValueError("DZ independent receipt does not qualify TZ launch")
    return dict(path=receipt, sha256=receipt_sha256, **value)


def authenticate(root, source, case, contract, receipt=None, receipt_sha256=None):
    model = contract["models"][case]
    target = root / model["target_directory"]
    for name, want in contract["source_sha256"].items():
        if digest(source / name) != want:
            raise ValueError("Changed frozen source prerequisite: " + name)
    proof = json.loads((source / contract["verification"]).read_text())
    if not (
        proof["saved_proof_reconstruction_pass"]
        and proof["native_root_counts_dz_tz"] == [64, 64]
        and proof["original_import_seed_preservation_verdict"] is False
    ):
        raise ValueError("Saved distinct-frame proof prerequisite failed")
    for name, want in model["input_sha256"].items():
        if digest(target / name) != want:
            raise ValueError("Changed frozen native input: " + name)
    cfg = json.loads((target / "config.json").read_text())
    if not (
        cfg["basis"] == model["basis"]
        and cfg["bond_length"] == 2.5
        and cfg["active_orbitals"] == [6, 3, 3, 0]
        and cfg["sa_singlet_roots"] == cfg["sa_triplet_roots"] == [5] * 4
    ):
        raise ValueError("Wrong full-ensemble target construction")
    previous = predecessor(root, case, digest(source / CONTRACT), receipt, receipt_sha256)
    return target, previous


def compute(root, source, work, case, contract, receipt=None, receipt_sha256=None):
    from projects.ukrmol_co.a2_channel_preflight import property_inventory, stage
    from projects.ukrmol_co.stretched_basis_feasibility import read_cidata

    report = dict(
        case=case,
        scientific_exit=2,
        property_inventory_verdict=None,
        stages=[],
        native_denprop_calls=0,
        qc_calls=0,
        target_diagonalizations=0,
        scattering_calls=0,
        new_maps=0,
        route_holdout_s_matrices_analyzed=0,
        full_sweep_released=False,
    )
    try:
        target, previous = authenticate(root, source, case, contract, receipt, receipt_sha256)
        report["predecessor"] = previous
        geometry = target / "output/CO/geom1"
        native = read_cidata((geometry / "fort.26").read_bytes())
        expected = {key: value["energies"].tolist() for key, value in native.items()}
        copied = {}
        for name in COPY:
            shutil.copy2(geometry / name, work / name)
            copied[name] = digest(work / name)
        report["copied_input_sha256"] = copied
        deck = work / "target.denprop.inp"
        deck.write_text(
            matched_deck(
                (geometry / "inputs/target.denprop.inp").read_text(),
                contract["models"][case]["original_property_counts"],
            )
        )
        report["deck_sha256"] = digest(deck)
        if report["deck_sha256"] != contract["matched_deck_sha256"]:
            raise ValueError("Exact prospective property deck differs")
        executable = Path("/opt/ukrmolp/bin.double/denprop")
        if digest(executable) != contract["denprop_sha256"]:
            raise ValueError("Changed pinned native DENPROP executable")
        report["native_denprop_calls"] = 1
        save(work / "result.json", report)
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
        inventory = property_inventory((work / "fort.24").read_text(), expected)
        original = json.loads((target / "target.json").read_text())
        difference = inventory["ground_dipole_z_au"] - original["ground_dipole_au"][2]
        if inventory["target_count"] != 64 or abs(difference) > 1e-5:
            raise ValueError("64-state property/dipole gate failed")
        if any(digest(work / name) != want for name, want in copied.items()):
            raise ValueError("Copied native input changed during DENPROP")
        report.update(
            scientific_exit=0,
            property_inventory_verdict=True,
            property_inventory=inventory,
            ground_dipole_difference_au=difference,
            property_sha256=digest(work / "fort.24"),
        )
    except Exception as error:
        report["error"] = repr(error)
    finally:
        try:
            authenticate(root, source, case, contract, receipt, receipt_sha256)
        except Exception as error:
            report.update(
                integrity_error=repr(error), scientific_exit=2, property_inventory_verdict=None
            )
        save(work / "result.json", report)
    return report["scientific_exit"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=["dz", "tz"], required=True)
    parser.add_argument("--profiled", action="store_true")
    parser.add_argument("--dz-receipt")
    parser.add_argument("--dz-receipt-sha256")
    args = parser.parse_args()
    root, source, work = Path("/work"), Path("/opt/qmodeling"), Path("/out/properties")
    contract = json.loads((source / CONTRACT).read_text())
    if args.profiled:
        return compute(
            root, source, work, args.case, contract, args.dz_receipt, args.dz_receipt_sha256
        )
    from projects.ukrmol_co.run import run_profiled

    work.mkdir()
    (work / "scratch").mkdir()
    env = os.environ | {"TMPDIR": str(work / "scratch"), "PSI_SCRATCH": str(work / "scratch")}
    return run_profiled(
        work,
        env,
        [
            sys.executable,
            "-m",
            "projects.ukrmol_co.stretched_basis_properties",
            *sys.argv[1:],
            "--profiled",
        ],
    )


if __name__ == "__main__":
    raise SystemExit(main())
