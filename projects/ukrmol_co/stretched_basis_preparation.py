"""One missing stretched64-state full-continuum/B1 preparation; no inner solve."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
from pathlib import Path

from projects.ukrmol_co.stretched_basis_properties import COPY, digest, save

CONTRACT = "projects/ukrmol_co/stretched-basis-preparation-contract-20261031.json"
PROPERTY_CONTRACT = "projects/ukrmol_co/stretched-basis-properties-contract-20261026.json"
TEMPLATES = "projects/ukrmol_co/compressed-cas12-preparation-inputs-20261026"
REUSE = "projects/ukrmol_co/stretched-tz-preparation-reuse-release-20261031.json"


def prior_preparation(root, source, case, receipt=None, receipt_sha256=None):
    if case == "dz":
        if receipt is not None or receipt_sha256 is not None:
            raise ValueError("DZ has no preparation predecessor override")
        return None
    if not receipt or not receipt_sha256:
        raise ValueError("TZ preparation requires independent public DZ preparation success")
    path = root / receipt
    if not path.resolve().is_relative_to(root.resolve()) or digest(path) != receipt_sha256:
        raise ValueError("Changed/out-of-root DZ preparation receipt")
    value = json.loads(path.read_text())
    if not (
        value.get("case") == "dz"
        and value.get("preparation_verdict") is True
        and value.get("contract_sha256") == digest(source / CONTRACT)
        and value.get("original_scientific_command_container_owner_exits") == [0] * 4
        and value.get("public_reader_exit") == 0
        and value.get("roots_per_sector") == [8] * 8
        and re.fullmatch(r"[0-9a-f]{64}", value.get("moints_sha256", ""))
        and re.fullmatch(r"[0-9a-f]{64}", value.get("csf_sha256", ""))
    ):
        raise ValueError("DZ preparation receipt does not qualify TZ")
    return dict(path=receipt, sha256=receipt_sha256)


def authenticate(root, source, case, contract, receipt=None, receipt_sha256=None):
    for name, want in contract["source_sha256"].items():
        if digest(source / name) != want:
            raise ValueError("Changed frozen source prerequisite: " + name)
    properties = json.loads((source / PROPERTY_CONTRACT).read_text())
    for selected, model in contract["models"].items():
        public = json.loads((source / model["property_receipt"]).read_text())
        if not (
            public["case"] == selected
            and public["contract_sha256"] == digest(source / PROPERTY_CONTRACT)
            and public["independent_property_inventory_verdict"] is True
            and public["public_reader_exit"] == 0
            and public["target_count"] == 64
            and public["roots_per_sector"] == [8] * 8
            and public["original_scientific_command_container_owner_exits"] == [0] * 4
        ):
            raise ValueError("Both independent64-state property inputs must qualify")
    model = contract["models"][case]
    target = root / properties["models"][case]["target_directory"]
    for name, want in properties["models"][case]["input_sha256"].items():
        if digest(target / name) != want:
            raise ValueError("Changed original target: " + name)
    for name, want in model["input_sha256"].items():
        if digest(root / name) != want:
            raise ValueError("Changed completed property prerequisite: " + name)
    props = root / model["property_directory"]
    actual = json.loads((props / "result.json").read_text())
    public = json.loads((source / model["property_receipt"]).read_text())
    if not (
        actual["scientific_exit"] == 0
        and actual["property_inventory_verdict"] is True
        and digest(props / "fort.24") == public["property_sha256"]
    ):
        raise ValueError("Original property output no longer qualifies")
    previous = prior_preparation(root, source, case, receipt, receipt_sha256)
    return target, props, previous


def cards(inventory, source):
    from projects.ukrmol_co.a2_channel_preflight import replace_assignment
    from projects.ukrmol_co.compressed_cas12_preparation import resolved_cards

    folder = source / TEMPLATES
    congen, scatci = resolved_cards(
        (folder / "scattering.congen.doublet.B1.inp").read_text(),
        (folder / "scattering.scatci.doublet.B1.inp").read_text(),
        inventory,
        "B1",
    )
    for name, values in dict(
        numtgt=[8] * 8,
        ntgtf=[i for i in range(1, 9) for _ in range(8)],
        ntgts=list(range(1, 9)) * 8,
    ).items():
        scatci = replace_assignment(scatci, name, ",".join(map(str, values)))
    return congen, scatci


def reuse_prerequisites(root, source, case, previous):
    """Authenticate the separate opt-in release; the executed DZ contract stays unchanged."""
    release = json.loads((source / REUSE).read_text())
    if (
        case != "tz"
        or previous != release["predecessor"]
        or digest(source / CONTRACT) != release["contract_sha256"]
    ):
        raise ValueError("CSF reuse requires its exact TZ-only release and public DZ predecessor")
    prepared = root / release["preparation_directory"]
    if not prepared.resolve().is_relative_to(root.resolve()):
        raise ValueError("Out-of-root preparation reuse")
    for name, want in release["input_sha256"].items():
        if digest(prepared / name) != want:
            raise ValueError("Changed qualified DZ reuse input: " + name)
    receipt = json.loads((root / previous["path"]).read_text())
    if receipt["csf_sha256"] != release["input_sha256"]["fort.721"]:
        raise ValueError("Reuse CSF is not the public qualified DZ output")
    return prepared, release


def reuse_csf(work, prepared, release, inventory, congen):
    """Copy algebraic CSFs, never basis-dependent integrals, CIDATA or properties."""
    if (
        inventory != release["native_inventory"]
        or (prepared / "scattering.congen.doublet.B1.inp").read_text() != congen
    ):
        raise ValueError("TZ resolved inventory/deck differs; CSF reuse not released")
    for name, want in release["target_csf_sha256"].items():
        if digest(work / name) != want:
            raise ValueError("TZ target-CSF bytes differ: " + name)
    for name in ["fort.721", "congen.B1.out"]:
        shutil.copy2(prepared / name, work / name)
        if digest(work / name) != release["input_sha256"][name]:
            raise ValueError("Copied qualified CSF/log changed")
    return dict(
        source_directory=release["preparation_directory"],
        csf_sha256=release["input_sha256"]["fort.721"],
        original_congen_exit=0,
        reason="Identical resolved deck and eight target-CSF hashes; geometry-free algebra only",
    )


def compute(
    root, source, work, case, contract, receipt=None, receipt_sha256=None, reuse_dz_csf=False
):
    from projects.ukrmol_co.a2_channel_preflight import property_inventory, stage
    from projects.ukrmol_co.compressed_cas12_preparation import TARGET, dimensions, values
    from projects.ukrmol_co.immutable_continuum import congen_completion, native_inventory
    from projects.ukrmol_co.stretched_basis_feasibility import read_cidata

    report = dict(
        case=case,
        scientific_exit=2,
        preparation_verdict=None,
        stages=[],
        integral_calls=0,
        congen_calls=0,
        qc_calls=0,
        target_diagonalizations=0,
        property_exports=0,
        inner_solves=0,
        maps=0,
        automatic_successor=False,
        route_holdout_s_matrices_analyzed=0,
        full_sweep_released=False,
    )
    copied = {}
    try:
        target, props, previous = authenticate(
            root, source, case, contract, receipt, receipt_sha256
        )
        report["predecessor"] = previous
        reuse = reuse_prerequisites(root, source, case, previous) if reuse_dz_csf else None
        for name, want in contract["engine_sha256"].items():
            if digest(Path(name)) != want:
                raise ValueError("Changed preparation engine: " + name)
        geometry = target / "output/CO/geom1"
        native = read_cidata((geometry / "fort.26").read_bytes())
        expected = {name: value["energies"].tolist() for name, value in native.items()}
        own = property_inventory((props / "fort.24").read_text(), expected)
        report["own_property_inventory"] = own
        target_meta = json.loads((target / "target.json").read_text())
        target_meta.update(
            orbital_ensemble_components=40,
            retained_target_components=64,
            states=[
                dict(spin=s.split(".")[0], irrep=s.split(".")[1], energy_hartree=e, root=i)
                for s, energies in expected.items()
                for i, e in enumerate(energies, 1)
            ],
        )
        save(work / "target.json", target_meta)
        for name in [n for n in COPY if n != "fort.16"]:
            shutil.copy2(geometry / name, work / name)
            copied[name] = digest(work / name)
        shutil.copy2(props / "fort.24", work / "fort.24")
        copied["fort.24"] = digest(work / "fort.24")
        report["copied_input_sha256"] = copied
        deck = work / "integrals.inp"
        shutil.copy2(source / TEMPLATES / "target.scatci_integrals.inp", deck)
        if values(deck.read_text(), "nob") != TARGET[:4]:
            raise ValueError("Changed target orbital inventory in full-continuum deck")
        (work / "inp").symlink_to(deck.name)
        report["integral_calls"] = 1
        save(work / "result.json", report)
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
        if not (work / "moints").is_file() or not (work / "moints").stat().st_size:
            raise ValueError("No completed full-continuum integrals")
        (work / "fort.16").symlink_to("moints")
        congen, scatci = cards(inventory, source)
        card = work / "scattering.congen.doublet.B1.inp"
        card.write_text(congen)
        (work / "scattering.scatci.doublet.B1.inp").write_text(scatci)
        if reuse is not None:
            report["csf_reuse"] = reuse_csf(work, *reuse, inventory, congen)
            report["csf_reuse"]["release_sha256"] = digest(source / REUSE)
        else:
            report["congen_calls"] = 1
            save(work / "result.json", report)
            stage(
                work,
                report,
                "congen.B1",
                ["/opt/ukrmolp/bin.double/congen"],
                card,
                contract["resources"]["congen_timeout_seconds"],
                os.environ.copy(),
            )
        log = (work / "congen.B1.out").read_text()
        counts = dimensions(log, scatci, roots_per_sector=8)
        congen_completion(log, 721, counts["raw_congen_csfs"])
        if counts["contracted_dimension"] > contract["maximum_contracted_dimension"]:
            raise ValueError("Actual B1 dimension exceeds preparation cap")
        report["B1"] = counts
        report["output_sha256"] = {
            name: digest(work / name)
            for name in [
                "moints",
                "fort.721",
                "target.json",
                "integrals.inp",
                "scattering.congen.doublet.B1.inp",
                "scattering.scatci.doublet.B1.inp",
                "integrals.native.out",
                "congen.B1.out",
            ]
        }
        report.update(scientific_exit=0, preparation_verdict=True)
    except Exception as error:
        report["error"] = repr(error)
    finally:
        try:
            _, _, previous = authenticate(root, source, case, contract, receipt, receipt_sha256)
            if reuse_dz_csf:
                reuse_prerequisites(root, source, case, previous)
            if any(digest(work / n) != want for n, want in copied.items()):
                raise ValueError("Copied target/property changed during preparation")
            report["input_integrity_after"] = True
        except Exception as error:
            report.update(scientific_exit=2, preparation_verdict=None, integrity_error=repr(error))
        save(work / "result.json", report)
    return report["scientific_exit"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=["dz", "tz"], required=True)
    parser.add_argument("--profiled", action="store_true")
    parser.add_argument("--previous-receipt")
    parser.add_argument("--previous-receipt-sha256")
    parser.add_argument("--reuse-dz-csf", action="store_true")
    args = parser.parse_args()
    root, source, work = Path("/work"), Path("/opt/qmodeling"), Path("/out/preparation")
    contract = json.loads((source / CONTRACT).read_text())
    if args.profiled:
        return compute(
            root,
            source,
            work,
            args.case,
            contract,
            args.previous_receipt,
            args.previous_receipt_sha256,
            args.reuse_dz_csf,
        )
    from projects.ukrmol_co.run import run_profiled

    work.mkdir()
    (work / "scratch").mkdir()
    return run_profiled(
        work,
        os.environ | {"TMPDIR": str(work / "scratch"), "PSI_SCRATCH": str(work / "scratch")},
        [
            sys.executable,
            "-m",
            "projects.ukrmol_co.stretched_basis_preparation",
            *sys.argv[1:],
            "--profiled",
        ],
    )


if __name__ == "__main__":
    raise SystemExit(main())
