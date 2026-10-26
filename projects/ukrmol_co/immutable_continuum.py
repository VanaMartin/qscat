"""Construct continuum integrals and target CSFs without reoptimizing CO orbitals."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import struct
from pathlib import Path

from projects.ukrmol_co.a2_channel_preflight import stage
from projects.ukrmol_co.fifty_frame_coverage import checkpoint_export_overlap
from projects.ukrmol_co.fifty_import_verification import digest

COUNTS = ("target", "PCOs", "TGT+PCOs", "continuum", "total")
IRREPS = ("A1", "B1", "B2", "A2")


def native_inventory(text, *, target_counts=(7, 3, 3, 0, 0, 0, 0, 0)):
    """Parse complete native inventories, rejecting duplicates and inconsistent sums."""
    inventory = {}
    for kind in COUNTS:
        rows = re.findall(
            r"Final number of orbitals for " + re.escape(kind) + r"[^\S\n]*\n([^\n]+)", text
        )
        if len(rows) != 1 or not re.fullmatch(r"\s*\d+(?:\s+\d+){7}\s*", rows[0]):
            raise ValueError(f"Incomplete or duplicate native orbital inventory: {kind}")
        inventory[kind] = list(map(int, rows[0].split()))
    target = list(target_counts)
    if target not in ([7, 3, 3, 0, 0, 0, 0, 0], [8, 3, 3, 0, 0, 0, 0, 0]):
        raise ValueError("Only the fixed CAS11/CAS12 target inventories are supported")
    if inventory["target"] != target or inventory["TGT+PCOs"] != target or any(inventory["PCOs"]):
        raise ValueError("Native target/PCO space differs from the fixed CAS")
    continuum = inventory["continuum"]
    if min(continuum[:4]) < 1 or any(continuum[4:]) or continuum[1] != continuum[2]:
        raise ValueError("Missing or asymmetric C2v continuum")
    if inventory["total"] != [a + b for a, b in zip(target, continuum, strict=True)]:
        raise ValueError("Native total inventory differs from target plus continuum")
    errors = re.findall(r"NUMBER OF ERROR MESSAGES:\s*(\d+)", text)
    if len(errors) != 4 or any(map(int, errors)):
        raise ValueError("Missing four-rank completion or native integral errors")
    if "Two_electron_integrals took" not in text:
        raise ValueError("Two-electron integral completion is missing")
    return inventory


def congen_deck(original, total):
    """Only update full orbital addressing, never nob0 or the target CSF groups."""
    pattern = r"(?m)^(  nob\s*=\s*)58,37,37,22,(.*)$"
    if len(re.findall(pattern, original)) != 1:
        raise ValueError("Unknown baseline target-CSF total-orbital card")
    if len(total) != 8 or any(total[4:]) or min(total[:4]) < 1:
        raise ValueError("Invalid full C2v orbital count")
    return re.sub(pattern, lambda m: m[1] + ",".join(map(str, total[:4])) + "," + m[2], original)


def congen_completion(text, unit, dimension):
    """Use CONGEN's own final projection/write markers, not DENPROP's trailer."""
    if text.count("Wavefn projection (projec()) - completed") != 1:
        raise ValueError("Missing or duplicate final CONGEN projection marker")
    units = re.findall(r"Data on CSFs has been written to file \(MEGUL\)\s*=\s*(\d+)", text)
    counts = re.findall(r"Number of CSFs\s*\(nocsf\)\s*=\s*(\d+)", text)
    if units != [str(unit)] or counts != [str(dimension)]:
        raise ValueError("Wrong or incomplete CONGEN final unit/CSF count")


def csf_records(path, total, sector, dimension):
    """Inspect the pinned native header and record framing, not CSF phase meaning."""
    import hashlib

    records = []
    spin, irrep = sector.split(".")
    expected_spin = 0.0 if spin == "singlet" else 1.0
    with path.open("rb") as stream:
        while marker := stream.read(4):
            if len(marker) != 4:
                raise ValueError("Truncated CSF record marker")
            length = struct.unpack("<I", marker)[0]
            if length < 1 or length > path.stat().st_size - stream.tell() - 4:
                raise ValueError("Invalid CSF record length")
            data = stream.read(length)
            if stream.read(4) != marker:
                raise ValueError("Truncated or mismatched CSF record")
            if not records:
                if length != 304:
                    raise ValueError("Unexpected pinned target CSF header")
                header = struct.unpack("<23q", data[120:])
                if data[:120].decode().strip() != f"CO - {spin} {irrep}":
                    raise ValueError("Unexpected target-CSF sector title")
                if header[0] != IRREPS.index(irrep) or struct.unpack("<dd", data[128:144]) != (
                    expected_spin,
                    expected_spin,
                ):
                    raise ValueError("Native target-CSF spin/symmetry mismatch")
                if (header[5], header[6], header[7], header[8]) != (
                    sum(total),
                    2 * sum(total),
                    dimension,
                    14,
                ):
                    raise ValueError("Native target-CSF orbital/electron/spin/dimension mismatch")
            records.append(dict(bytes=length, sha256=hashlib.sha256(data).hexdigest()))
    if len(records) != 6:
        raise ValueError("Incomplete target-CSF record inventory")
    return records


def prepare(root, output, contract):
    def verify():
        for name, want in contract["input_sha256"].items():
            if digest(root / name) != want:
                raise ValueError(f"Pinned preparation input changed: {name}")

    verify()
    frame = json.loads((root / contract["frame_result"]).read_text())
    rejected = json.loads((root / contract["overlap_result"]).read_text())
    owner = json.loads((root / contract["overlap_execution"]).read_text())
    assert frame["new_frame_import_verdict"] is True
    assert rejected["retained_subspace_dominance_verdict"] is False
    assert rejected["all_eight_sectors_complete"] is True and owner["original_owner_exit"] == 1
    assert frame["checkpoint_sha256"] == digest(root / contract["checkpoint"])
    qc = json.loads((root / contract["target_report"]).read_text())
    assert qc["molden_sha256"] == digest(root / contract["molden"])
    output.mkdir(exist_ok=False)
    report = dict(
        stages=[],
        native_layout_generation_verdict=None,
        retained_subspace_dominance_verdict=False,
        original_overlap_exits=[1, 1, 1],
        native_target_vector_reuse_verdict=None,
        scattering_ensemble_verdict=None,
        new_orbital_optimizations=0,
        target_diagonalizations=0,
        property_exports=0,
        scattering_inner_solves=0,
    )
    try:
        report["checkpoint_export_overlap"] = checkpoint_export_overlap(
            root / contract["checkpoint"], root / contract["molden"]
        )
        shutil.copy2(root / contract["molden"], output / "co.molden")
        deck = output / "integrals.inp"
        shutil.copy2(root / contract["integral_deck"], deck)
        (output / "inp").symlink_to(deck.name)
        stage(
            output,
            report,
            "integrals",
            ["mpirun", "-np", "4", "/opt/ukrmolp/bin.double/scatci_integrals"],
            deck,
            120,
            os.environ.copy(),
        )
        log = output / "log_file.0"
        native_text = log.read_text() + (output / "integrals.out").read_text()
        (output / "integrals.native.out").write_text(native_text)
        inventory = native_inventory(native_text)
        report["native_inventory"] = inventory
        assert (output / "moints").is_file() and (output / "moints").stat().st_size > 0
        (output / "fort.16").symlink_to("moints")
        dimensions = contract["target_csf_dimensions"]
        report["target_csf_records"] = {}
        for irrep in IRREPS:
            for spin in ("singlet", "triplet"):
                label = f"{spin}.{irrep}"
                original = root / contract["baseline_inputs"] / f"target.congen.{label}.inp"
                card = output / f"target.congen.{label}.inp"
                card.write_text(congen_deck(original.read_text(), inventory["total"]))
                stage(
                    output,
                    report,
                    f"congen.{label}",
                    ["/opt/ukrmolp/bin.double/congen"],
                    card,
                    30,
                    os.environ.copy(),
                )
                unit = (710 if spin == "singlet" else 730) + IRREPS.index(irrep)
                congen_completion(
                    (output / f"congen.{label}.out").read_text(), unit, dimensions[label]
                )
                report["target_csf_records"][label] = csf_records(
                    output / f"fort.{unit}",
                    inventory["total"],
                    label,
                    dimensions[label],
                )
        assert digest(output / "co.molden") == digest(root / contract["molden"])
        assert digest(deck) == digest(root / contract["integral_deck"])
        report["native_layout_generation_verdict"] = True
        report["scope"] = (
            "Full-continuum integrals and target CSFs generated only. "
            "No cached-vector conversion, target Hamiltonian/import, properties "
            "or scattering qualification."
        )
        verify()
        (output / "result.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    except BaseException as error:
        report["error"] = repr(error)
        (output / "failure.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        raise
    finally:
        verify()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.root, args.output, json.loads(args.contract.read_text()))


if __name__ == "__main__":
    main()
