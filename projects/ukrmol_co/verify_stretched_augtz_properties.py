"""Saved-only independent qualification of the completed aug-TZ64 property table."""

from __future__ import annotations

import argparse
import datetime
import json
from itertools import pairwise
from pathlib import Path

from projects.ukrmol_co.verify_stretched_augtz_properties_launch import ARCHIVE, CONTRACT, SOURCE
from projects.ukrmol_co.verify_stretched_properties import files_in, require, sha
from projects.ukrmol_co.verify_stretched_properties import verify as property_verify

CONTRACT_SHA = "3f82307b03773ef276897a64e62eae2516a9eb0075debdf2cf7235cfd7e7266b"


def observation_summary(observations, started, finished):
    times = [datetime.datetime.fromisoformat(r["observed_utc"]) for r in observations]
    require(len(times) >= 2 and all(b > a for a, b in pairwise(times)), "Ordered observations")
    require(
        observations[0]["owner"]["state"] == "running"
        and observations[-1]["owner"]["state"] == "finished",
        "Original launch/terminal witnesses",
    )
    gaps = [(b - a).total_seconds() for a, b in pairwise(times)]
    first_delay = (times[0] - datetime.datetime.fromisoformat(started)).total_seconds()
    last_delay = (times[-1] - datetime.datetime.fromisoformat(finished)).total_seconds()
    require(first_delay >= 0 and last_delay >= 0, "Observation timing")
    return dict(
        observation_count=len(times),
        maximum_observation_gap_seconds=max(gaps),
        first_observation_delay_seconds=first_delay,
        terminal_observation_delay_seconds=last_delay,
        strict_300second_observation_cadence_passed=max([first_delay, *gaps, last_delay]) <= 300,
    )


def verify(archive):
    checked = property_verify(
        archive,
        contract_path=CONTRACT,
        contract_sha256=CONTRACT_SHA,
        source_entries=20,
        runtime_controls=10,
        source_commit=SOURCE,
    )
    files = files_in(archive)

    def load(name):
        return json.loads(files[name])

    source = files_in(files["source.tar.gz"])
    contract = json.loads(source[CONTRACT])
    owner, result, container = [
        load(n)
        for n in ["execution.json", "outputs/properties/result.json", "container-inspect.json"]
    ]
    require(
        sha(files["source.tar.gz"]) == ARCHIVE and owner["source_commit"] == SOURCE,
        "Exact original aug-TZ source freeze",
    )
    require(
        result["case"] == "augtz"
        and result["predecessor"] is None
        and result["native_denprop_calls"] == 1
        and all(
            result[n] == 0
            for n in [
                "qc_calls",
                "target_diagonalizations",
                "scattering_calls",
                "new_maps",
                "route_holdout_s_matrices_analyzed",
            ]
        )
        and not result["full_sweep_released"]
        and not contract["automatic_successor"],
        "One property call with no scientific successor",
    )
    require(len(contract["source_sha256"]) == 6, "Six exact source prerequisites")
    for name, pin in contract["source_sha256"].items():
        require(sha(source[name]) == pin, "Source prerequisite " + name)
    after = load("integrity-after.json")
    model = contract["models"]["augtz"]
    require(
        after["original_input_sha256"] == model["input_sha256"]
        and after["copied_input_sha256"] == result["copied_input_sha256"]
        and after["source_sha256"] == load("source-manifest.json"),
        "After-run immutable input attestations",
    )
    for name in ["fort.16", "fort.26", "co.casscf.chk", "co.molden"]:
        require(
            sha(files["outputs/properties/" + name]) == result["copied_input_sha256"][name],
            "Retained copied raw input " + name,
        )
    host = container["HostConfig"]
    mounts = {m["Destination"]: m for m in container["Mounts"]}
    require(
        host["CpusetCpus"] == "4-11"
        and host["Memory"] == host["MemorySwap"] == 8 * 2**30
        and host["NetworkMode"] == "none"
        and set(mounts) == {"/work", "/opt/qmodeling", "/out"}
        and not mounts["/work"]["RW"]
        and not mounts["/opt/qmodeling"]["RW"],
        "Original immutable allocation",
    )
    limits = contract["resources"]
    require(
        result["stages"][0]["wall_seconds"] < limits["native_timeout_seconds"]
        and checked["resources"]["wall_seconds"] < limits["container_timeout_seconds"]
        and owner["wall_seconds"] < limits["owner_timeout_seconds"],
        "Original hard time caps",
    )
    require(load("owner-process-exit.json")["original_exit"] == 0, "Original launcher process exit")
    observations = [json.loads(line) for line in files["monitor.jsonl"].decode().splitlines()]
    observation = observation_summary(observations, owner["started_utc"], owner["finished_utc"])
    return checked | dict(
        source_commit=SOURCE,
        executed_git_sources_verified=20,
        source_prerequisites_verified=6,
        original_exits=[0] * 8,
        exit_order=[
            "native",
            "scientific",
            "profile",
            "command",
            "container",
            "owner",
            "service",
            "launcher_process",
        ],
        original_input_pins_attested=18,
        copied_input_pins_attested=12,
        native_denprop_calls=1,
        new_qc_calls=0,
        new_ci_calls=0,
        new_hamiltonian_actions=0,
        new_maps=0,
        route_holdout_s_matrices_analyzed=0,
        automatic_successor=False,
        full_sweep_released=False,
        owner_wall_seconds=owner["wall_seconds"],
        **observation,
        independence=(
            "Separate stdlib native-energy/vector-norm/property parser,not producer "
            "inventory;raw CIDATA/property/checkpoint/Molden/one-electron-file hashes checked;"
            "large CSF and original input after-hashes are authenticated owner/collection "
            "attestations,not independent raw CSF replay"
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("archive", type=Path)
    a = p.parse_args()
    print(json.dumps(verify(a.archive.read_bytes()), indent=2, allow_nan=False))
