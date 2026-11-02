"""Reconstruct the zero-map own-DZ adapter stop from immutable saved evidence."""

from __future__ import annotations

import argparse
import ast
import csv
import io
import json
import subprocess
from pathlib import Path

import numpy as np

from projects.ukrmol_co.verify_stretched_inner import numbers, records
from projects.ukrmol_co.verify_stretched_preparation import stream_files
from projects.ukrmol_co.verify_stretched_properties import files_in, require, sha

SOURCE_COMMIT = "d04d3b0392ac3e78e4b70a7332cf55bb5cb5b3a7"
SOURCE_ARCHIVE_SHA = "45676b0080b1c457609df70b7ead1437e9d09a45e7aaf0e6f254eeb3ebd99ee6"
ARCHIVE_SHA = "7a33fc3888e44ee376df313531b4bbfbf1611775c4755c60dd7f7853a937e1a6"
CONTRACT = "projects/ukrmol_co/stretched-dz-pole-execution-20261102.json"
PROPOSAL = "projects/ukrmol_co/stretched-dz-pole-proposal-20261102.json"


def returned_dict_keys(raw, function):
    tree = ast.parse(raw)
    body = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == function)
    result = next(n.value for n in body.body if isinstance(n, ast.Return))
    require(isinstance(result, ast.Dict), "Explicit model dictionary required")
    return [ast.literal_eval(n) for n in result.keys]


def verify(archive, collection):
    require(collection["archive"]["sha256"] == ARCHIVE_SHA, "Original allocation archive")
    capture = [n for n in collection["files"] if Path(n).name.endswith(("fort.10", "fort.21"))]
    files, hashes, _ = stream_files(archive, collection, capture_names=capture)

    def load(n):
        return json.loads(files[n])

    source = files_in(files["source.tar.gz"])
    pins = load("source-manifest.json")
    require(
        set(source) == set(pins)
        and len(pins) == 47
        and hashes["source.tar.gz"] == SOURCE_ARCHIVE_SHA,
        "Frozen source inventory",
    )
    for n, pin in pins.items():
        require(
            sha(source[n]) == pin
            and source[n] == subprocess.check_output(["git", "show", SOURCE_COMMIT + ":" + n]),
            "Executed Git source " + n,
        )
    execution = json.loads(source[CONTRACT])
    proposal = source[PROPOSAL]
    require(sha(proposal) == execution["proposal_sha256"], "Frozen proposal")
    c = json.loads(proposal) | execution
    refs = collection["reference_files"]
    require(
        c["execution_authorized"] and set(refs) == set(c["input_sha256"]) and len(refs) == 7,
        "Authorized raw input inventory",
    )
    for n, pin in c["input_sha256"].items():
        require(hashes[refs[n]] == pin, "Original input " + n)
    for n, pin in c["source_sha256"].items():
        require(sha(source[n]) == pin, "Method prerequisite " + n)
    owner, result, resources, container = [
        load(n)
        for n in [
            "execution.json",
            "outputs/check/result.json",
            "outputs/check/resources.json",
            "container-inspect.json",
        ]
    ]
    require(
        owner["source_commit"] == SOURCE_COMMIT
        and owner["source_archive_sha256"] == SOURCE_ARCHIVE_SHA
        and owner["owner_sha256"] == sha(source["projects/ukrmol_co/stretched_dz_pole_owner.py"])
        and owner["source_integrity_after"],
        "Immutable owner/source integrity",
    )
    exits = [
        result["scientific_exit"],
        resources["exit_code"],
        owner["original_command_exit"],
        owner["original_container_exit"],
        owner["original_owner_exit"],
        load("owner-process-exit.json")["original_exit"],
    ]
    require(
        exits == [2] * 6 and not owner["timed_out"] and not owner["container_oom_killed"],
        "Original scientific/wrapper stop exits",
    )
    require(
        load("outputs/controls.json")["original_exit"] == 0
        and b"Ran 26 tests" in files["outputs/controls.log"],
        "26 original runtime controls",
    )
    require(
        container["Image"] == c["image"] and not container["State"]["Running"],
        "Pinned terminal container",
    )
    host = container["HostConfig"]
    mounts = {m["Destination"]: m for m in container["Mounts"]}
    require(
        host["CpusetCpus"] == "0-7"
        and host["Memory"] == host["MemorySwap"] == 8 * 2**30
        and host["NetworkMode"] == "none"
        and set(mounts) == {"/work", "/opt/qmodeling", "/out"}
        and not mounts["/work"]["RW"]
        and not mounts["/opt/qmodeling"]["RW"],
        "Eight-core caps/mounts",
    )
    for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
        require(name + "=1" in container["Config"]["Env"], "Library threads")
    require(
        result["real_rows"] == result["roots"] == result["evaluations"] == []
        and result["new_maps_reserved"] == 0
        and not any(n.endswith(".npz") for n in files)
        and result["error"] == "KeyError('interface')",
        "Zero-map adapter failure",
    )
    require(
        result["finite_reference_pole_verdict"] is result["physical_pi_support_verdict"] is None
        and not result["tz_started"]
        and result["native_calls"] == result["route_holdout_s_matrices_analyzed"] == 0
        and result["inputs_unchanged"]
        and not result["automatic_successor"]
        and not result["full_sweep_released"],
        "Protected scientific scope",
    )
    keys = returned_dict_keys(
        source["projects/ukrmol_co/native_full_reconstruction.py"], "_model_from_interface"
    )
    require(
        "interface" not in keys
        and b'model["interface"]["channel_labels"]'
        in source["projects/ukrmol_co/stretched_dz_pole.py"],
        "Exact source contract mismatch",
    )
    # Authenticate the model ingredients independently, without invoking the failed adapter.
    prefix = c["model"]["native_directory"] + "/"
    boundary = records(files[refs[prefix + "fort.21"]])
    channel = records(files[refs[prefix + "fort.10"]])
    receipt = load(refs[c["model"]["receipt"]])
    require(
        receipt["inner_export_verdict"] is True and receipt["original_exits"] == [0] * 6,
        "Own native prerequisite",
    )
    w = numbers(boundary[7], "<f8", 2048 * 400).reshape(2048, 400)
    poles = numbers(boundary[6], "<f8", 2048)
    for name, value in [("amplitudes", w), ("poles", poles)]:
        require(
            sha(np.asarray(value, dtype="<f8").tobytes()) == receipt["array_sha256"][name],
            "Independent own array " + name,
        )
    import struct

    labels = [list(struct.unpack("<4q", row[:32])) for row in channel[68:74]]
    require(labels == c["open_labels"], "Pinned native labels exist despite omitted model metadata")
    previous = load(refs[c["completion_reading"]])
    require(
        previous["reader_exit"] == 0
        and previous["result"]["combined_20degree_real_prerequisite_verdict"]
        and previous["result"]["original_exits"] == [0] * 6,
        "Prior real prerequisite unchanged",
    )
    samples = list(csv.DictReader(io.StringIO(files["outputs/check/resources.csv"].decode())))
    require(
        max(int(r["container_memory_bytes"]) for r in samples)
        == resources["sampled_peak_container_memory_bytes"],
        "Original resource sample peak",
    )
    monitor = [json.loads(line) for line in files["monitor.jsonl"].decode().splitlines()]
    require(
        monitor[-1]["owner_process_exit"] == 2 and monitor[-1]["science"]["complex_reserved"] == 0,
        "Original terminal observation",
    )
    return dict(
        adapter_stop_verified=True,
        original_exits=exits,
        scientific_exit=2,
        real_maps=0,
        complex_maps=0,
        finite_reference_pole_verdict=None,
        physical_pi_support_verdict=None,
        error=result["error"],
        model_returned_keys=keys,
        own_open_labels=labels,
        prior_real_prerequisite_verdict=True,
        source_commit=SOURCE_COMMIT,
        raw_files_verified=len(hashes),
        executed_git_sources_verified=len(pins),
        raw_input_pins_verified=7,
        source_prerequisites_verified=9,
        resources=resources,
        resource_samples_verified=len(samples),
        owner_wall_seconds=owner["wall_seconds"],
        finished_utc=owner["finished_utc"],
        native_calls=0,
        new_ode_maps=0,
        automatic_successor=False,
        full_sweep_released=False,
        interpretation=(
            "Adapter requests absent nested interface metadata after own arrays validate;"
            "not an ODE,pole,character or electronic-model rejection;no retry released"
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
