"""Read saved property evidence independently, without importing the producer or running it."""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import io
import itertools
import json
import math
import struct
import tarfile
from pathlib import Path

CONTRACT = "projects/ukrmol_co/stretched-basis-properties-contract-20261026.json"
CONTRACT_SHA = "7bcf27226b67d7a54b74de4f84423b96e832122be289a2ae1921cc09512ea712"
SECTORS = tuple(
    f"{spin}.{irrep}" for irrep in ("A1", "B1", "B2", "A2") for spin in ("singlet", "triplet")
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def files_in(data):
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as archive:
        files = {}
        for entry in archive:
            if entry.isdir():
                continue
            require(
                entry.isfile() and entry.name not in files, "Nonregular/duplicate evidence member"
            )
            files[entry.name] = archive.extractfile(entry).read()
        return files


def native_energies(data):
    """Independent little-endian record/energy/vector-norm reader; no Hamiltonian actions."""
    cursor = 0

    def record():
        nonlocal cursor
        require(cursor + 8 <= len(data), "Truncated CIDATA record")
        size = struct.unpack_from("<I", data, cursor)[0]
        end = cursor + 4 + size
        require(
            end + 4 <= len(data) and struct.unpack_from("<I", data, end)[0] == size,
            "Broken CIDATA framing",
        )
        answer = memoryview(data)[cursor + 4 : end]
        cursor = end + 4
        return answer

    answer, maximum_norm_error = {}, 0.0
    for index, sector in enumerate(SECTORS, 1):
        require(bytes(record()) == b"********                CIDATA  ", "Wrong CIDATA marker")
        header = record()
        require(len(header) == 208, "Wrong CIDATA header")
        number, count = struct.unpack_from("<qq", header)
        centers, dimension, roots, symmetry = struct.unpack_from("<4q", header, 136)
        label = bytes(header[16:136]).decode().strip()
        require(
            number == index
            and count == 11
            and centers == 2
            and roots == 8
            and dimension >= 8
            and symmetry == (index - 1) // 2
            and label == "CO CI calculation - " + sector.replace(".", " "),
            "Wrong native sector",
        )
        spin = float(index % 2 == 0)
        require(struct.unpack_from("<ddq", header, 168) == (spin, spin, 14), "Wrong spin/electrons")
        nuclear = struct.unpack_from("<d", header, 192)[0]
        for charge in (6, 8):
            center = record()
            require(
                len(center) == 40 and struct.unpack_from("<d", center, 32)[0] == charge,
                "Wrong nuclear center",
            )
        mixed = record()
        require(len(mixed) == 8 * (2 * dimension + roots), "Wrong native mixed record")
        require(
            all(v[0] in (-1, 1) for v in struct.iter_unpack("<q", mixed[: 8 * dimension])),
            "Wrong native CSF phases",
        )
        energies = [x + nuclear for x in struct.unpack_from("<8d", mixed, 8 * dimension)]
        require(
            all(map(math.isfinite, energies)) and energies == sorted(energies), "Invalid energies"
        )
        answer[sector] = energies
        for root in range(1, 9):
            vector = record()
            require(
                len(vector) == 8 * (dimension + 1) and struct.unpack_from("<q", vector)[0] == root,
                "Wrong native vector",
            )
            norm = math.fsum(v[0] * v[0] for v in struct.iter_unpack("<d", vector[8:]))
            require(math.isfinite(norm) and abs(norm - 1) <= 1e-9, "Unnormalized native vector")
            maximum_norm_error = max(maximum_norm_error, abs(norm - 1))
    require(cursor == len(data), "Trailing native data")
    return answer, maximum_norm_error


def inventory(text, energies, ground_dipole):
    states, moments, headers = [], {}, []
    for line in text.splitlines():
        tokens = line.split()
        if not tokens:
            continue
        if tokens[0] == "6":
            headers.append(int(tokens[4]))
        elif tokens[0] == "5":
            states.append(
                (
                    int(tokens[1]),
                    tokens[-2].lower() + "." + tokens[-1],
                    float(tokens[8].replace("D", "E")),
                )
            )
        elif tokens[0] == "1":
            key = tuple(map(int, tokens[1:8]))
            value = float(tokens[8].replace("D", "E"))
            require(key not in moments and math.isfinite(value), "Invalid/duplicate moment")
            moments[key] = value
    require(
        headers == [64] and [s[0] for s in states] == list(range(1, 65)), "Wrong state IDs/count"
    )
    require(
        all(math.isfinite(s[2]) for s in states)
        and all(b[2] >= a[2] - 1e-10 for a, b in itertools.pairwise(states)),
        "Unordered energies",
    )
    counts = collections.Counter(s[1] for s in states)
    require(set(counts) == set(SECTORS) and set(counts.values()) == {8}, "Wrong retained sectors")
    error = max(
        abs(a - b)
        for sector in SECTORS
        for a, b in zip([s[2] for s in states if s[1] == sector], energies[sector], strict=True)
    )
    require(error <= 1e-7, "Own-native energy mismatch")
    identifiers, covered, ranks = {s[0] for s in states}, set(), set()
    for left, _, right, _, kind, rank, projection in moments:
        require(
            left in identifiers
            and right in identifiers
            and kind == 3
            and rank in (1, 2)
            and abs(projection) <= rank,
            "Invalid moment reference/rank",
        )
        covered.update((left, right))
        ranks.add(rank)
    require(covered == identifiers and ranks == {1, 2}, "Missing moment coverage")
    ground = next(s[0] for s in states if s[1] == "singlet.A1")
    dipole_error = abs(-moments[(ground, 0, ground, 0, 3, 1, 0)] - ground_dipole)
    require(dipole_error <= 1e-5, "Ground dipole mismatch")
    return dict(
        target_count=64,
        roots_per_sector=[counts[s] for s in SECTORS],
        multipole_records=len(moments),
        maximum_energy_error_hartree=error,
        ground_dipole_error_au=dipole_error,
    )


def verify(archive):
    files = files_in(archive)

    def load(name):
        return json.loads(files[name])

    collection = load("collection.json")
    require(
        set(files) == set(collection["files_sha256"]) | {"collection.json"}, "Wrong evidence set"
    )
    for name, want in collection["files_sha256"].items():
        require(sha(files[name]) == want, "Changed archived evidence: " + name)
    source, pins = files_in(files["source.tar.gz"]), load("source-manifest.json")
    require(set(source) == set(pins) and len(pins) == 12, "Wrong capsule source set")
    for name, want in pins.items():
        require(sha(source[name]) == want, "Changed capsule source: " + name)
    require(sha(source[CONTRACT]) == CONTRACT_SHA, "Changed exact contract")
    contract = json.loads(source[CONTRACT])
    result, owner = load("outputs/properties/result.json"), load("execution.json")
    require(sha(files["source.tar.gz"]) == owner["source_archive_sha256"], "Changed source archive")
    require(
        sha(source["projects/ukrmol_co/stretched_basis_properties_owner.py"])
        == owner["owner_sha256"],
        "Changed executing owner",
    )
    case = result["case"]
    model = contract["models"][case]
    require(
        collection["original_input_sha256"] == model["input_sha256"], "Original input pins differ"
    )
    for leaf in ("config.json", "target.json", "inputs/target.denprop.inp"):
        key = "output/CO/geom1/" + leaf if leaf.startswith("inputs/") else leaf
        require(
            sha(files["original/" + leaf]) == model["input_sha256"][key],
            "Changed original metadata",
        )
    for name, want in result["copied_input_sha256"].items():
        require(want == model["input_sha256"]["output/CO/geom1/" + name], "Wrong copy pin")
    require(
        collection["copied_input_sha256"] == result["copied_input_sha256"],
        "Changed copy attestation",
    )
    require(
        sha(files["outputs/properties/fort.26"]) == result["copied_input_sha256"]["fort.26"],
        "Changed CIDATA",
    )
    require(
        sha(files["outputs/properties/target.denprop.inp"]) == contract["matched_deck_sha256"],
        "Changed property deck",
    )
    require(
        sha(files["outputs/properties/fort.24"]) == result["property_sha256"], "Changed properties"
    )
    resources, service = load("outputs/properties/resources.json"), load("service-return.json")
    container = load("container-inspect.json")
    require(
        container["Image"] == contract["image"]
        and not container["State"]["Running"]
        and container["State"]["ExitCode"] == 0
        and not container["State"]["OOMKilled"],
        "Wrong image/container state",
    )
    samples = list(csv.DictReader(io.StringIO(files["outputs/properties/resources.csv"].decode())))
    for column, key in [
        ("container_memory_bytes", "sampled_peak_container_memory_bytes"),
        ("run_disk_bytes", "sampled_peak_run_allocated_bytes"),
        ("scratch_bytes", "sampled_peak_scratch_allocated_bytes"),
    ]:
        require(
            samples and max(int(row[column]) for row in samples) == resources[key],
            "Resource sample peak differs: " + column,
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
        exits == [0] * 6 and result["property_inventory_verdict"] is True, "Failed original exits"
    )
    require(
        len(result["stages"]) == 1
        and result["stages"][0]["original_exit"] == 0
        and not result["stages"][0]["timed_out"]
        and not owner["timed_out"]
        and not owner["container_oom_killed"]
        and owner["source_integrity_after"],
        "Native/integrity stop",
    )
    require(
        load("outputs/controls.json")["original_exit"] == 0
        and b"Ran 8 tests" in files["outputs/controls.log"],
        "Runtime controls failed",
    )
    require(
        b"has completed execution" in files["outputs/properties/target.denprop.out"],
        "No completion",
    )
    energies, norm_error = native_energies(files["outputs/properties/fort.26"])
    checked = inventory(
        files["outputs/properties/fort.24"].decode(),
        energies,
        load("original/target.json")["ground_dipole_au"][2],
    )
    return dict(
        case=case,
        contract_sha256=CONTRACT_SHA,
        independent_property_inventory_verdict=True,
        original_scientific_command_container_owner_exits=exits[:4],
        original_service_exit=exits[4],
        original_profile_exit=exits[5],
        property_sha256=result["property_sha256"],
        maximum_native_norm_error=norm_error,
        finished_utc=owner["finished_utc"],
        resources=resources,
        archived_files_verified=len(files),
        resource_samples_verified=len(samples),
        source_entries_verified=12,
        original_input_pins_authenticated=len(model["input_sha256"]),
        independence=(
            "Separate stdlib parser on saved native energies/vectors and properties; "
            "no producer imports, native calls, new CI or Hamiltonian actions. "
            "Original hash attestations authenticated, not regenerated integrals."
        ),
        **checked,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.archive.read_bytes()), indent=2, allow_nan=False))
