"""One continuum-integral/two-CONGEN preparation; no target or scattering solve."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
from itertools import pairwise
from pathlib import Path

from projects.ukrmol_co.a2_channel_preflight import property_inventory, replace_assignment, stage
from projects.ukrmol_co.compressed_cas12_properties import digest
from projects.ukrmol_co.immutable_continuum import congen_completion, native_inventory

TARGET = [8, 3, 3, 0, 0, 0, 0, 0]


def values(text, name):
    """Read one active native integer assignment, never a commented example."""
    clean = "\n".join(line.split("!", 1)[0] for line in text.splitlines())
    found = re.findall(rf"\b{re.escape(name)}\s*=\s*([\d,\s]+)", clean, re.I)
    if len(found) != 1:
        raise ValueError(f"No unique native assignment: {name}")
    return list(map(int, re.findall(r"\d+", found[0])))


def resolved_cards(congen, scatci, inventory, sector):
    """Reuse geometry-free CAS12 algebra; only update full continuum addressing."""
    irrep = {"B1": 1, "B2": 2}[sector]
    mcont = [irrep ^ (i // 2) for i in range(8)]
    checks = [
        (congen, "nob", [59, 37, 37, 22]),
        (congen, "nob0", TARGET[:4]),
        (congen, "nelect", [15]),
        (congen, "qntot", [2, irrep, 0]),
        (congen, "iscat", [2]),
        (congen, "megul", [720 + irrep]),
        (scatci, "mcont", mcont),
        (scatci, "numtgt", [5] * 8),
        (scatci, "ntgsym", [8]),
        (scatci, "icitg", [1]),
        (scatci, "iexpc", [1]),
        (scatci, "ntgtf", [i for i in range(1, 9) for _ in range(5)]),
        (scatci, "ntgts", list(range(1, 6)) * 8),
    ]
    if any(values(text, name) != want for text, name, want in checks):
        raise ValueError("Changed geometry-free CAS12 contraction/CONGEN cards")
    total = inventory["total"]
    radial = [inventory["continuum"][i] for i in mcont]
    if inventory["target"] != TARGET or len(total) != 8 or any(total[4:]):
        raise ValueError("Not the full CAS12 C2v target/continuum inventory")
    return (
        replace_assignment(congen, "nob", ",".join(map(str, total[:4]))),
        replace_assignment(scatci, "notgt", ",".join(map(str, radial))),
    )


def dimensions(text, scatci):
    """Count L2 CSFs plus five targets times the actual radial inventories."""
    groups = re.split(r"WFN GROUP\s+\d+", text)[1:]
    if len(groups) != 9 or "L^2:" not in groups[-1]:
        raise ValueError("No complete eight-target/one-L2 CONGEN grouping")
    for group, label in zip(
        groups[:8],
        [s + " " + i for i in ["A1", "B1", "B2", "A2"] for s in ["singlet", "triplet"]],
        strict=True,
    ):
        if label not in group or "x continuum^1" not in group:
            raise ValueError("Wrong native target-group order")
    ranges = [
        tuple(map(int, row))
        for row in re.findall(r"CSF NUMBERS\s+(\d+) TO\s+(\d+)\s*\(\s*(\d+) CSFS\)", groups[-1])
    ]
    totals = re.findall(r"TOTAL NUMBER OF CSF'S GENERATED IS\s+(\d+)", text)
    if not ranges or len(totals) != 1:
        raise ValueError("Missing or duplicate final CONGEN dimensions")
    first, last = ranges[0][0], ranges[-1][1]
    if (
        int(totals[0]) != last
        or any(b - a + 1 != count for a, b, count in ranges)
        or any(left[1] + 1 != right[0] for left, right in pairwise(ranges))
    ):
        raise ValueError("Inconsistent or noncontiguous L2 CSFs")
    radial, roots = values(scatci, "notgt"), values(scatci, "numtgt")
    if len(radial) != 8 or roots != [5] * 8 or min(radial) < 1:
        raise ValueError("Not the matched five-per-sector contraction")
    bound = last - first + 1
    continuum = sum(a * b for a, b in zip(radial, roots, strict=True))
    return dict(
        raw_congen_csfs=last,
        last_continuum_prototype_csf=first - 1,
        uncontracted_l2_csfs=bound,
        contracted_continuum_dimension=continuum,
        contracted_dimension=bound + continuum,
        radial_continuum_per_target_sector=radial,
        retained_roots_per_sector=roots,
    )


def authenticate(root, source, contract):
    for base, key in [(root, "input_sha256"), (source, "source_prerequisite_sha256")]:
        for name, want in contract[key].items():
            if digest(base / name) != want:
                raise ValueError(f"Changed preparation prerequisite: {name}")
    props = root / contract["property_directory"]
    producer = json.loads((props / "result.json").read_text())
    owner = json.loads((props.parent.parent / "execution.json").read_text())
    public = json.loads((source / contract["property_verification"]).read_text())
    if (
        producer["property_inventory_verdict"] is not True
        or public["property_inventory_verdict"] is not True
        or [
            producer["scientific_exit"],
            owner["original_command_exit"],
            owner["original_container_exit"],
            owner["original_owner_exit"],
        ]
        != [0] * 4
        or owner["timed_out"]
        or owner["container_oom_killed"]
        or not owner["source_integrity_after"]
    ):
        raise ValueError("Completed independently authenticated properties required")
    for name, want in contract["engine_sha256"].items():
        if digest(Path(name)) != want:
            raise ValueError(f"Changed pinned native engine: {name}")
    target = root / contract["target_directory"]
    config = json.loads((target / "config.json").read_text())
    if (
        config["bond_length"] != 1.9
        or config["basis"] != "cc-pVDZ"
        or config["active_orbitals"] != [6, 3, 3, 0]
    ):
        raise ValueError("Wrong physical compressed CAS12 frame")
    saved = json.loads((target / "target.json").read_text())
    expected = {
        spin + "." + irrep: [
            s["energy_hartree"] for s in saved["states"] if (s["spin"], s["irrep"]) == (spin, irrep)
        ]
        for irrep in ["A1", "B1", "B2", "A2"]
        for spin in ["singlet", "triplet"]
    }
    own = property_inventory((props / "fort.24").read_text(), expected)
    if own != producer["property_inventory"]:
        # JSON stringifies the label dictionary's integer IDs.
        if json.loads(json.dumps(own)) != producer["property_inventory"]:
            raise ValueError("Changed own property map")
    return target, props, own


def compute(root, source, work, contract):
    report = dict(
        preparation_verdict=None,
        scientific_exit=2,
        stages=[],
        sectors={},
        integral_calls=0,
        congen_calls=0,
        qc_calls=0,
        target_diagonalizations=0,
        property_exports=0,
        scattering_inner_solves=0,
        new_propagations=0,
        route_holdout_s_matrices_analyzed=0,
        full_sweep_released=False,
    )

    def save():
        (work / "result.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    try:
        target, props, own = authenticate(root, source, contract)
        report["own_properties"] = own
        report["inputs_verified_before"] = len(contract["input_sha256"])
        geometry = target / "output/CO/geom1"
        for name in ["co.molden", "co.casscf.chk", "fort.26"] + [
            f"fort.{n}" for n in [710, 730, 711, 731, 712, 732, 713, 733]
        ]:
            shutil.copy2(geometry / name, work / name)
        shutil.copy2(props / "fort.24", work / "fort.24")
        shutil.copy2(target / "target.json", work / "target.json")
        templates = source / contract["templates"]
        deck = work / "integrals.inp"
        shutil.copy2(templates / "target.scatci_integrals.inp", deck)
        if values(deck.read_text(), "nob") != TARGET[:4]:
            raise ValueError("Full-continuum deck target orbitals changed")
        (work / "inp").symlink_to(deck.name)
        report["integral_calls"] = 1
        save()
        stage(
            work,
            report,
            "integrals",
            ["mpirun", "-np", "4", "/opt/ukrmolp/bin.double/scatci_integrals"],
            deck,
            contract["resources"]["integral_timeout_seconds"],
            os.environ.copy(),
        )
        text = (work / "log_file.0").read_text() + (work / "integrals.out").read_text()
        (work / "integrals.native.out").write_text(text)
        inventory = native_inventory(text, target_counts=TARGET)
        report["native_inventory"] = inventory
        if not (work / "moints").is_file() or (work / "moints").stat().st_size == 0:
            raise ValueError("No completed full-continuum integrals")
        (work / "fort.16").symlink_to("moints")
        for sector, unit in [("B1", 721), ("B2", 722)]:
            congen, scatci = resolved_cards(
                (templates / f"scattering.congen.doublet.{sector}.inp").read_text(),
                (templates / f"scattering.scatci.doublet.{sector}.inp").read_text(),
                inventory,
                sector,
            )
            card = work / f"scattering.congen.doublet.{sector}.inp"
            card.write_text(congen)
            (work / f"scattering.scatci.doublet.{sector}.inp").write_text(scatci)
            report["congen_calls"] += 1
            save()
            stage(
                work,
                report,
                "congen." + sector,
                ["/opt/ukrmolp/bin.double/congen"],
                card,
                contract["resources"]["congen_timeout_seconds"],
                os.environ.copy(),
            )
            log = (work / f"congen.{sector}.out").read_text()
            counts = dimensions(log, scatci)
            congen_completion(log, unit, counts["raw_congen_csfs"])
            csf = work / f"fort.{unit}"
            if not csf.is_file() or csf.stat().st_size == 0:
                raise ValueError("No complete native scattering CSFs")
            report["sectors"][sector] = counts | dict(
                csf_sha256=digest(csf), csf_bytes=csf.stat().st_size
            )
            if counts["contracted_dimension"] > contract["maximum_contracted_dimension"]:
                raise ValueError("Contracted dimension exceeds preparation planning cap")
        if (
            report["sectors"]["B1"]["contracted_dimension"]
            != report["sectors"]["B2"]["contracted_dimension"]
        ):
            raise ValueError("Broken B1/B2 contracted dimensions")
        report["output_sha256"] = {
            str(p.relative_to(work)): digest(p)
            for p in sorted(work.iterdir())
            if p.is_file() and p.name not in ["result.json", "resources.csv", "run.log"]
        }
        report.update(preparation_verdict=True, scientific_exit=0)
    except Exception as error:
        report["error"] = repr(error)
    finally:
        try:
            authenticate(root, source, contract)
            report["inputs_verified_after"] = len(contract["input_sha256"])
            for name in ["co.molden", "co.casscf.chk", "fort.26"] + [
                f"fort.{n}" for n in [710, 730, 711, 731, 712, 732, 713, 733]
            ]:
                if digest(work / name) != digest(
                    root / contract["target_directory"] / "output/CO/geom1" / name
                ):
                    raise ValueError(f"Copied native target/frame changed: {name}")
            if digest(work / "fort.24") != digest(
                root / contract["property_directory"] / "fort.24"
            ):
                raise ValueError("Copied matched properties changed")
            report["copied_target_and_properties_unchanged"] = True
        except Exception as error:
            report.update(integrity_error=repr(error), preparation_verdict=None, scientific_exit=2)
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
            "projects.ukrmol_co.compressed_cas12_preparation",
            *sys.argv[1:],
            "--profiled",
        ],
    )


if __name__ == "__main__":
    raise SystemExit(main())
