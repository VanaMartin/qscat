"""Stream saved B1 preparation bytes and reconstruct counts independently of the producer."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import itertools
import json
import re
import struct
import tarfile
from pathlib import Path

from projects.ukrmol_co.verify_stretched_properties import (
    files_in,
    native_energies,
    require,
    sha,
)
from projects.ukrmol_co.verify_stretched_properties import inventory as property_inventory

PREFIX = "outputs/preparation/"
CONTRACT = "projects/ukrmol_co/stretched-basis-preparation-contract-20261031.json"
PROPERTY_CONTRACT = "projects/ukrmol_co/stretched-basis-properties-contract-20261026.json"
CONTRACT_SHA = "d00a24a308e3b8da289e631ceba937f7086a61d37e8984ac79cb6da151cf6f4f"
REUSE = "projects/ukrmol_co/stretched-tz-preparation-reuse-release-20261031.json"


def values(text, key):
    active = "\n".join(line.split("!", 1)[0] for line in text.splitlines())
    matches = re.findall(rf"\b{key}\s*=\s*([\d,\s]+)", active, re.I)
    require(len(matches) == 1, "No unique assignment: " + key)
    return list(map(int, re.findall(r"\d+", matches[0])))


def orbital_inventory(text):
    result = {}
    for kind in ["target", "PCOs", "TGT+PCOs", "continuum", "total"]:
        rows = re.findall(
            r"Final number of orbitals for " + re.escape(kind) + r"[^\S\n]*\n([^\n]+)", text
        )
        require(len(rows) == 1, "Incomplete/duplicate native inventory: " + kind)
        result[kind] = list(map(int, rows[0].split()))
        require(len(result[kind]) == 8, "Wrong native inventory width")
    target = [8, 3, 3, 0, 0, 0, 0, 0]
    cont = result["continuum"]
    require(
        result["target"] == result["TGT+PCOs"] == target and result["PCOs"] == [0] * 8,
        "Changed physical target/PCO space",
    )
    require(min(cont[:4]) > 0 and cont[1] == cont[2] and cont[4:] == [0] * 4, "Bad continuum")
    require(
        result["total"] == [a + b for a, b in zip(target, cont, strict=True)], "Bad orbital sum"
    )
    require(
        re.findall(r"NUMBER OF ERROR MESSAGES:\s*(\d+)", text) == ["0"] * 4
        and "Two_electron_integrals took" in text,
        "Incomplete four-rank integral result",
    )
    return result


def contraction(log, scatci, orbitals):
    groups = re.split(r"WFN GROUP\s+\d+", log)[1:]
    require(len(groups) == 9 and "L^2:" in groups[-1], "Wrong CONGEN grouping")
    labels = [s + " " + i for i in ["A1", "B1", "B2", "A2"] for s in ["singlet", "triplet"]]
    require(
        all(
            label in group and "x continuum^1" in group
            for group, label in zip(groups[:8], labels, strict=True)
        ),
        "Changed target groups",
    )
    ranges = [
        tuple(map(int, row))
        for row in re.findall(r"CSF NUMBERS\s+(\d+) TO\s+(\d+)\s*\(\s*(\d+) CSFS\)", groups[-1])
    ]
    require(bool(ranges), "No L2 ranges")
    first, last = ranges[0][0], ranges[-1][1]
    require(
        all(b - a + 1 == n for a, b, n in ranges)
        and all(left[1] + 1 == right[0] for left, right in itertools.pairwise(ranges))
        and sum(n for _, _, n in ranges) == last - first + 1,
        "Broken L2 ranges",
    )
    require(
        re.findall(r"TOTAL NUMBER OF CSF'S GENERATED IS\s+(\d+)", log) == [str(last)]
        and log.count("Wavefn projection (projec()) - completed") == 1
        and re.findall(r"Data on CSFs has been written to file \(MEGUL\)\s*=\s*(\d+)", log)
        == ["721"]
        and re.findall(r"Number of CSFs\s*\(nocsf\)\s*=\s*(\d+)", log) == [str(last)],
        "Incomplete B1 CONGEN result",
    )
    mcont = [1 ^ (i // 2) for i in range(8)]
    radial = [orbitals["continuum"][i] for i in mcont]
    require(
        values(scatci, "mcont") == mcont
        and values(scatci, "notgt") == radial
        and values(scatci, "numtgt") == [8] * 8
        and values(scatci, "ntgtf") == [i for i in range(1, 9) for _ in range(8)]
        and values(scatci, "ntgts") == list(range(1, 9)) * 8,
        "Wrong64-state contraction",
    )
    continuum = 8 * sum(radial)
    return dict(
        raw_congen_csfs=last,
        last_continuum_prototype_csf=first - 1,
        uncontracted_l2_csfs=last - first + 1,
        contracted_continuum_dimension=continuum,
        contracted_dimension=last - first + 1 + continuum,
        radial_continuum_per_target_sector=radial,
        retained_roots_per_sector=[8] * 8,
    )


def csf_header(data, total, raw):
    require(
        len(data) >= 312 and struct.unpack_from("<I", data)[0] == 304 and data[:4] == data[308:312],
        "Invalid native CSF header record",
    )
    header = data[4:308]
    require(header[:120].decode().strip() == "e + CO - doublet B1", "Wrong CSF title")
    integers = struct.unpack("<23q", header[120:])
    require(
        integers[0] == 1
        and struct.unpack_from("<2d", header, 128) == (0.5, 0.5)
        and [integers[i] for i in [5, 6, 7, 8]] == [sum(total), 2 * sum(total), raw, 15],
        "Wrong native CSF spin/orbital/electron/dimension",
    )


def stream_files(path, collection):
    with path.open("rb") as stream:
        require(
            path.stat().st_size == collection["archive"]["bytes"]
            and hashlib.file_digest(stream, "sha256").hexdigest()
            == collection["archive"]["sha256"],
            "Changed full evidence archive",
        )
    kept, digests, headers = {}, {}, {}
    with tarfile.open(path) as archive:
        for entry in archive:
            require(entry.isfile() and entry.name not in digests, "Unexpected archive member")
            pin = collection["files"][entry.name]
            require(entry.size == pin["bytes"], "Wrong archived file size")
            capture = Path(entry.name).suffix in (".json", ".inp", ".out", ".log", ".csv")
            capture = (
                capture
                or entry.name == "source.tar.gz"
                or entry.name in (PREFIX + "fort.24", PREFIX + "fort.26")
            )
            chunks, header, digest = [], b"", hashlib.sha256()
            with archive.extractfile(entry) as stream:
                while chunk := stream.read(1024 * 1024):
                    if not header:
                        header = chunk[:312]
                    digest.update(chunk)
                    if capture:
                        chunks.append(chunk)
            require(digest.hexdigest() == pin["sha256"], "Changed archived data: " + entry.name)
            digests[entry.name], headers[entry.name] = digest.hexdigest(), header
            if capture:
                kept[entry.name] = b"".join(chunks)
    require(set(digests) == set(collection["files"]), "Incomplete evidence archive")
    return kept, digests, headers


def verify(path, collection):
    kept, digests, headers = stream_files(path, collection)

    def load(name):
        return json.loads(kept[name])

    result, owner = load(PREFIX + "result.json"), load("execution.json")
    service, resources = load("service-return.json"), load(PREFIX + "resources.json")
    source = files_in(kept["source.tar.gz"])
    pins = load("source-manifest.json")
    reused = result.get("csf_reuse")
    require(
        set(source) == set(pins) and len(pins) == (31 if reused else 30), "Wrong frozen source set"
    )
    for name, want in pins.items():
        require(sha(source[name]) == want, "Changed frozen source: " + name)
    require(
        digests["source.tar.gz"] == owner["source_archive_sha256"]
        and sha(source["projects/ukrmol_co/stretched_basis_preparation_owner.py"])
        == owner["owner_sha256"],
        "Wrong executing source/owner",
    )
    require(sha(source[CONTRACT]) == CONTRACT_SHA, "Changed executed preparation contract")
    contract = json.loads(source[CONTRACT])
    for name, want in contract["source_sha256"].items():
        require(sha(source[name]) == want, "Changed preparation source prerequisite")
    container = load("container-inspect.json")
    require(
        container["Image"] == contract["image"]
        and not container["State"]["Running"]
        and container["State"]["ExitCode"] == 0
        and not container["State"]["OOMKilled"],
        "Wrong image/container state",
    )
    exits = [
        result["scientific_exit"],
        owner["original_command_exit"],
        owner["original_container_exit"],
        owner["original_owner_exit"],
        service["original_service_exit"],
        resources["exit_code"],
    ]
    require(
        exits == [0] * 6
        and result["preparation_verdict"] is True
        and result["input_integrity_after"]
        and owner["source_integrity_after"]
        and not owner["timed_out"]
        and not owner["container_oom_killed"],
        "Failed original preparation",
    )
    if reused:
        release = json.loads(source[REUSE])
        predecessor = load("independent-dz-preparation-receipt.json")
        require(
            result["case"] == "tz"
            and sha(source[REUSE]) == reused["release_sha256"]
            and release["contract_sha256"] == CONTRACT_SHA
            and result["predecessor"] == release["predecessor"]
            and digests["independent-dz-preparation-receipt.json"]
            == release["predecessor"]["sha256"]
            and predecessor["case"] == "dz"
            and predecessor["preparation_verdict"] is True
            and predecessor["public_reader_exit"] == 0
            and predecessor["original_scientific_command_container_owner_exits"] == [0] * 4
            and predecessor["contract_sha256"] == CONTRACT_SHA
            and predecessor["csf_sha256"] == release["input_sha256"]["fort.721"]
            and predecessor["csf_sha256"] == digests[PREFIX + "fort.721"]
            and reused["original_congen_exit"] == 0,
            "Unqualified conditional CSF reuse",
        )
        for name in ("fort.721", "congen.B1.out", "scattering.congen.doublet.B1.inp"):
            require(
                digests[PREFIX + name] == release["input_sha256"][name], "Changed reused algebra"
            )
        require(
            digests["reuse-source/result.json"] == release["input_sha256"]["result.json"],
            "Changed DZ source result",
        )
        for name, want in release["target_csf_sha256"].items():
            require(digests[PREFIX + name] == want, "Changed own target CSF")
    require(
        result["integral_calls"] == 1
        and result["congen_calls"] == (0 if reused else 1)
        and [s["name"] for s in result["stages"]]
        == (["integrals"] if reused else ["integrals", "congen.B1"])
        and all(s["original_exit"] == 0 and not s["timed_out"] for s in result["stages"]),
        "Wrong native execution scope",
    )
    require(
        all(
            result[k] == 0
            for k in [
                "qc_calls",
                "target_diagonalizations",
                "property_exports",
                "inner_solves",
                "maps",
                "route_holdout_s_matrices_analyzed",
            ]
        ),
        "Unexpected scientific successor",
    )
    props = json.loads(source[PROPERTY_CONTRACT])["models"][result["case"]]
    require(collection["original_input_sha256"] == props["input_sha256"], "Wrong original pins")
    model = contract["models"][result["case"]]
    require(
        collection["property_prerequisite_sha256"] == model["input_sha256"], "Wrong property pins"
    )
    for name, want in model["input_sha256"].items():
        require(digests["prerequisites/" + name] == want, "Changed saved property prerequisite")
    for name in ("config.json", "target.json"):
        require(
            digests["original/" + name] == props["input_sha256"][name], "Wrong original metadata"
        )
    expected_copies = {"fort.26", "fort.24", "co.molden", "co.casscf.chk"} | {
        f"fort.{n}" for n in (710, 730, 711, 731, 712, 732, 713, 733)
    }
    require(set(result["copied_input_sha256"]) == expected_copies, "Wrong copied input set")
    for name, want in result["copied_input_sha256"].items():
        require(digests[PREFIX + name] == want, "Changed copied input")
        if name != "fort.24":
            require(
                want == props["input_sha256"]["output/CO/geom1/" + name], "Wrong native target copy"
            )
    receipt = json.loads(source[contract["models"][result["case"]]["property_receipt"]])
    require(
        digests[PREFIX + "fort.24"] == receipt["property_sha256"]
        and receipt["independent_property_inventory_verdict"] is True
        and receipt["public_reader_exit"] == 0,
        "Unqualified property copy",
    )
    for name, want in result["output_sha256"].items():
        require(digests[PREFIX + name] == want, "Changed preparation output")
    orbitals = orbital_inventory(kept[PREFIX + "integrals.native.out"].decode())
    require(orbitals == result["native_inventory"], "Wrong reported orbital inventory")
    if reused:
        require(
            orbitals == release["native_inventory"] == predecessor["native_inventory"],
            "Changed reused continuum addresses",
        )
    card = kept[PREFIX + "scattering.congen.doublet.B1.inp"].decode()
    require(
        values(card, "nob") == orbitals["total"][:4] and values(card, "nob0") == [8, 3, 3, 0],
        "Changed CONGEN target/full addresses",
    )
    original = source[
        "projects/ukrmol_co/compressed-cas12-preparation-inputs-20261026/scattering.congen.doublet.B1.inp"
    ].decode()
    require(card.split("&wfngrp") == original.split("&wfngrp"), "Changed CONGEN algebraic groups")
    counts = contraction(
        kept[PREFIX + "congen.B1.out"].decode(),
        kept[PREFIX + "scattering.scatci.doublet.B1.inp"].decode(),
        orbitals,
    )
    require(
        counts == result["B1"]
        and counts["contracted_dimension"] <= contract["maximum_contracted_dimension"],
        "Wrong contracted dimension",
    )
    csf_header(headers[PREFIX + "fort.721"], orbitals["total"], counts["raw_congen_csfs"])
    energies, norm = native_energies(kept[PREFIX + "fort.26"])
    target = load(PREFIX + "target.json")
    require(
        target["orbital_ensemble_components"] == 40 and target["retained_target_components"] == 64,
        "Changed orbital objective/retention",
    )
    own = property_inventory(
        kept[PREFIX + "fort.24"].decode(),
        energies,
        load("original/target.json")["ground_dipole_au"][2],
    )
    require(
        load("outputs/controls.json")["original_exit"] == 0
        and (b"Ran 10 tests" if reused else b"Ran 7 tests") in kept["outputs/controls.log"],
        "Runtime controls failed",
    )
    samples = list(csv.DictReader(io.StringIO(kept[PREFIX + "resources.csv"].decode())))
    for column, key in [
        ("container_memory_bytes", "sampled_peak_container_memory_bytes"),
        ("run_disk_bytes", "sampled_peak_run_allocated_bytes"),
        ("scratch_bytes", "sampled_peak_scratch_allocated_bytes"),
    ]:
        require(
            samples and max(int(s[column]) for s in samples) == resources[key],
            "Resource peak mismatch: " + column,
        )
    return dict(
        case=result["case"],
        contract_sha256=sha(source[CONTRACT]),
        preparation_verdict=True,
        original_scientific_command_container_owner_exits=exits[:4],
        original_service_exit=exits[4],
        original_profile_exit=exits[5],
        roots_per_sector=[8] * 8,
        moints_sha256=digests[PREFIX + "moints"],
        csf_sha256=digests[PREFIX + "fort.721"],
        source_entries_verified=len(pins),
        archived_files_verified=len(digests),
        resource_samples_verified=len(samples),
        native_inventory=orbitals,
        B1=counts,
        own_properties=own,
        maximum_native_norm_error=norm,
        resources=resources,
        finished_utc=owner["finished_utc"],
        native_stages=result["stages"],
        **({"csf_reuse": reused} if reused else {}),
        new_scientific_calls=0,
        full_sweep_released=False,
        independence=(
            "Full payload streaming hashes; separate native inventory, group/range/contraction "
            "arithmetic and CSF header checks. No native executable, CIRMAT phase audit, "
            "scattering interface or electronic-accuracy claim."
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("collection", type=Path)
    args = parser.parse_args()
    print(
        json.dumps(
            verify(args.archive, json.loads(args.collection.read_text())), indent=2, allow_nan=False
        )
    )
