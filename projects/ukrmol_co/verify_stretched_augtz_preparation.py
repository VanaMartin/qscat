"""Stream and independently reconstruct saved aug-TZ own-integral/CSF-reuse evidence."""

from __future__ import annotations

import argparse
import datetime
import json
from itertools import pairwise
from pathlib import Path

from projects.ukrmol_co.verify_stretched_preparation import stream_files
from projects.ukrmol_co.verify_stretched_preparation import verify as preparation_verify
from projects.ukrmol_co.verify_stretched_properties import files_in, require

SOURCE = "bc3b33a2faeb2a6320196907c1239a2f17ad73fb"
ARCHIVE = "f9b0123339798c2f9e949f7e97a107ed72a888f3e33c5229ec7ce3d4242ef9ba"
CONTRACT = "projects/ukrmol_co/stretched-augtz-preparation-contract-20261103.json"
CONTRACT_SHA = "32f88b36febe06f7bb6a6e1e42c3846a0a180b047b42013d1029c77f80062f54"
PROPERTY = "projects/ukrmol_co/stretched-augtz-properties-contract-20261102.json"
REUSE = "projects/ukrmol_co/stretched-augtz-preparation-reuse-20261103.json"


def original_exits(owner, result, resources, service, process):
    exits = [
        result["stages"][0]["original_exit"],
        result["scientific_exit"],
        resources["exit_code"],
        owner["original_command_exit"],
        owner["original_container_exit"],
        owner["original_owner_exit"],
        service["original_service_exit"],
        process["original_exit"],
    ]
    require(
        exits == [0] * 8,
        "Original native/science/profile/command/container/owner/service/process exits",
    )
    return exits


def verify(path, collection):
    checked = preparation_verify(
        path,
        collection,
        contract_path=CONTRACT,
        contract_sha256=CONTRACT_SHA,
        property_contract_path=PROPERTY,
        reuse_path=REUSE,
        source_entries=37,
        runtime_controls=14,
        source_commit=SOURCE,
    )
    kept, digests, _ = stream_files(path, collection)

    def load(n):
        return json.loads(kept[n])

    source = files_in(kept["source.tar.gz"])
    c = json.loads(source[CONTRACT])
    owner = load("execution.json")
    result = load("outputs/preparation/result.json")
    container = load("container-inspect.json")
    require(
        digests["source.tar.gz"] == ARCHIVE
        and owner["source_commit"] == SOURCE
        and result["case"] == "augtz"
        and result["integral_calls"] == 1
        and result["congen_calls"] == 0
        and owner["max_congen_calls"] == 0
        and not owner["automatic_successor"],
        "Exact aug-TZ-only scope",
    )
    exits = original_exits(
        owner,
        result,
        checked["resources"],
        load("service-return.json"),
        load("owner-process-exit.json"),
    )
    host = container["HostConfig"]
    mounts = {m["Destination"]: m for m in container["Mounts"]}
    require(
        host["CpusetCpus"] == "0-3"
        and host["Memory"] == host["MemorySwap"] == 32 * 2**30
        and host["NetworkMode"] == "none"
        and set(mounts) == {"/work", "/opt/qmodeling", "/out"}
        and not mounts["/work"]["RW"]
        and not mounts["/opt/qmodeling"]["RW"],
        "Original immutable four-core allocation",
    )
    require(
        checked["native_stages"][0]["wall_seconds"] < c["resources"]["integral_timeout_seconds"]
        and checked["resources"]["wall_seconds"] < c["resources"]["container_timeout_seconds"]
        and owner["wall_seconds"] < c["resources"]["owner_timeout_seconds"],
        "Native/container/owner caps",
    )
    observations = [json.loads(s) for s in kept["monitor.jsonl"].decode().splitlines()]
    started = datetime.datetime.fromisoformat(owner["started_utc"])
    times = [datetime.datetime.fromisoformat(row["observed_utc"]) for row in observations]
    gaps = [(b - a).total_seconds() for a, b in pairwise([started, *times])]
    require(
        times
        and all(g > 0 for g in gaps)
        and observations[-1]["owner_process_exit"] == 0
        and observations[-1]["owner"]["state"] == "finished",
        "Terminal observation and ordered cadence",
    )
    previous = load("independent-dz-preparation-receipt.json")
    return checked | dict(
        source_commit=SOURCE,
        executed_git_sources_verified=37,
        source_prerequisites_verified=len(c["source_sha256"]),
        original_exits=exits,
        runtime_controls_passed=14,
        original_input_hash_gate_attested=18,
        property_prerequisite_pins_verified=5,
        signed_csf_reused=True,
        new_integral_calls=1,
        new_congen_calls=0,
        new_ci_calls=0,
        new_qc_calls=0,
        new_property_exports=0,
        new_inner_solves=0,
        new_maps=0,
        route_holdout_s_matrices_analyzed=0,
        automatic_successor=False,
        own_integrals_differ_from_dz=checked["moints_sha256"] != previous["moints_sha256"],
        observation_count=len(observations),
        maximum_observation_gap_seconds=max(gaps),
        strict_30second_observation_cadence_passed=max(gaps) <= 30,
        inherited_60second_preparation_observation_bound_passed=max(gaps) <= 60,
        owner_wall_seconds=owner["wall_seconds"],
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("archive", type=Path)
    p.add_argument("collection", type=Path)
    a = p.parse_args()
    print(
        json.dumps(
            verify(a.archive, json.loads(a.collection.read_text())), indent=2, allow_nan=False
        )
    )
