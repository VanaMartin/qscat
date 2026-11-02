"""Authenticate the zero-map completion stop and its missing source pointer; no ODE."""

from __future__ import annotations

import argparse
import csv
import io
import json
import subprocess
from pathlib import Path

from projects.ukrmol_co.verify_stretched_preparation import stream_files
from projects.ukrmol_co.verify_stretched_properties import files_in, require, sha

SOURCE_COMMIT = "048bb941a785ff4276a8153270c6dd7eac2ec903"
ARCHIVE_SHA = "0f51c413b9ff7841563a7b3b92bc000c9417ce25f6acb015e9170c9eaf9ac407"
CONTRACT = "projects/ukrmol_co/stretched-real-completion-contract-20261102.json"
MISSING = "projects/ukrmol_co/qualification-evidence/stretched-tz-inner-receipt/artifacts.json"


def source_prerequisites(source, pins):
    missing = []
    for name, pin in pins.items():
        if name not in source:
            missing.append(name)
        else:
            require(sha(source[name]) == pin, "Changed present prerequisite: " + name)
    return missing


def verify(archive, collection):
    require(collection["archive"]["sha256"] == ARCHIVE_SHA, "Original failed allocation archive")
    files, hashes, _ = stream_files(archive, collection)

    def load(name):
        return json.loads(files[name])

    source = files_in(files["source.tar.gz"])
    pins = load("source-manifest.json")
    require(set(source) == set(pins) and len(pins) == 34, "Original34-file capsule")
    for name, pin in pins.items():
        require(sha(source[name]) == pin, "Frozen source pin " + name)
        require(
            source[name] == subprocess.check_output(["git", "show", SOURCE_COMMIT + ":" + name]),
            "Original Git bytes " + name,
        )
    execution = json.loads(source[CONTRACT])
    raw = source[execution["base_guard_contract"]]
    require(sha(raw) == execution["base_guard_sha256"], "Original base contract")
    base = json.loads(raw)
    inputs = base["input_sha256"] | execution["input_sha256"]
    methods = base["source_sha256"] | execution["source_sha256"]
    refs = collection["reference_files"]
    require(set(refs) == set(inputs) and len(inputs) == 12, "Original input inventory")
    for name, pin in inputs.items():
        require(hashes[refs[name]] == pin, "Present raw input " + name)
    missing = source_prerequisites(source, methods)
    require(missing == [MISSING] and len(methods) == 9, "Unique omitted TZ source pointer")
    result, owner, resources, container = [
        load(name)
        for name in [
            "outputs/check/result.json",
            "execution.json",
            "outputs/check/resources.json",
            "container-inspect.json",
        ]
    ]
    exits = [
        result["scientific_exit"],
        resources["exit_code"],
        owner["original_command_exit"],
        owner["original_container_exit"],
        owner["original_owner_exit"],
        load("owner-process-exit.json")["original_exit"],
    ]
    require(
        exits == [2] * 6
        and load("outputs/controls.json")["original_exit"] == 0
        and b"Ran 6 tests" in files["outputs/controls.log"],
        "Original stop/controls exits",
    )
    require(
        owner["source_commit"] == SOURCE_COMMIT
        and owner["source_archive_sha256"] == hashes["source.tar.gz"]
        and owner["owner_sha256"]
        == sha(source["projects/ukrmol_co/compressed_cas12_pole_owner.py"])
        and owner["source_integrity_after"]
        and not owner["timed_out"]
        and not owner["container_oom_killed"],
        "Original owner and source integrity",
    )
    require(
        container["Image"] == execution["image"]
        and not container["State"]["Running"]
        and container["State"]["ExitCode"] == 2,
        "Original container terminal state",
    )
    require(
        result["real_maps_reserved"] == 0
        and result["real_rows"] == []
        and result["real_completion_verdict"] is None
        and "FileNotFoundError" in result["error"]
        and not any(name.endswith(".npz") for name in files),
        "Zero-map prerequisite failure",
    )
    require(
        result["native_calls"]
        == result["new_complex_maps"]
        == result["route_holdout_s_matrices_analyzed"]
        == 0
        and not result["previous_dz80_repeated"]
        and not result["original_guard_requalified"]
        and not result["automatic_successor"]
        and not result["full_sweep_released"],
        "No successor/actions",
    )
    previous = load(refs[execution["previous_public_reading"]])
    old = load(refs[execution["original_guard_result"]])
    require(
        previous["reader_exit"] == 0
        and previous["result"]["contour_study_verdict"] is True
        and old["scientific_exit"] == 2
        and old["real_guard_verdict"] is None,
        "Prior qualifications unchanged",
    )
    samples = list(csv.DictReader(io.StringIO(files["outputs/check/resources.csv"].decode())))
    require(
        max(int(row["container_memory_bytes"]) for row in samples)
        == resources["sampled_peak_container_memory_bytes"],
        "Original resource sampling",
    )
    return dict(
        allocation_stop_verified=True,
        real_completion_verdict=None,
        original_exits=exits,
        source_commit=SOURCE_COMMIT,
        source_entries_verified=34,
        raw_files_verified=len(hashes),
        raw_input_pins_verified=12,
        present_source_prerequisites_verified=8,
        missing_source_prerequisites=missing,
        maps_reserved=0,
        resource_samples_verified=len(samples),
        resources=resources,
        owner_wall_seconds=owner["wall_seconds"],
        finished_utc=owner["finished_utc"],
        previous_dz80_contour_verdict=True,
        original40degree_scientific_exit=2,
        new_maps=0,
        native_calls=0,
        automatic_successor=False,
        full_sweep_released=False,
        interpretation=(
            "Source capsule omitted a pinned TZ receipt pointer;input gate stops before "
            "model/map access. Evidence reconstruction is not scientific success "
            "or permission to retry."
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
