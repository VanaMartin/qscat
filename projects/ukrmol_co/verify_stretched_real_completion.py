"""Read the saved five-map completion with own native surface/K algebra; no ODE."""

from __future__ import annotations

import argparse
import ast
import csv
import datetime
import io
import itertools
import json
import struct
import subprocess
from pathlib import Path

import numpy as np

from projects.ukrmol_co.verify_stretched_inner import numbers, records, symmetric
from projects.ukrmol_co.verify_stretched_preparation import stream_files
from projects.ukrmol_co.verify_stretched_properties import files_in, require, sha
from projects.ukrmol_co.verify_stretched_real_completion_stop import source_prerequisites
from projects.ukrmol_co.verify_stretched_real_guard import checks, measurements

SOURCE_COMMIT = "a4f86e2886792831e73903efca69d448dbfd5c3a"
ARCHIVE_SHA = "126165bf5a6bcaf2e5e31950b43bd6bac5f7a0e073b8569f1cbda01de53af359"
CONTRACT = "projects/ukrmol_co/stretched-real-completion-capsule-release-20261102.json"
MAPS = [("dz", 60, 120), ("tz", 60, 120), ("tz", 80, 60), ("tz", 80, 90), ("tz", 80, 120)]


def native_at(files, refs, spec, index):
    prefix = spec["native_directory"] + "/"
    boundary = records(files[refs[prefix + "fort.21"]])
    channel = records(files[refs[prefix + "fort.10"]])
    km = records(files[refs[prefix + "fort.921"]])
    poles = numbers(boundary[6], "<f8", 2048)
    w = numbers(boundary[7], "<f8", 2048 * 400).reshape(2048, 400)
    ground = struct.unpack_from("<d", channel[4], 32)[0]
    values = numbers(km[5 + index][24:], "<f8", 22)
    energy = values[0]
    k = symmetric(values[1:], 6)
    q = w.T @ (w / (poles - ground - energy / 2)[:, None])
    s = np.linalg.solve(np.eye(6) - 1j * k, np.eye(6) + 1j * k)
    return energy, q, s


def verify(archive, collection):
    require(collection["archive"]["sha256"] == ARCHIVE_SHA, "Original completion archive")
    capture = [
        n
        for n in collection["files"]
        if n.endswith(".npz") or Path(n).name.endswith(("fort.10", "fort.21", "fort.921"))
    ]
    files, hashes, _ = stream_files(archive, collection, capture_names=capture)

    def load(name):
        return json.loads(files[name])

    source = files_in(files["source.tar.gz"])
    pins = load("source-manifest.json")
    require(set(source) == set(pins) and len(pins) == 36, "Complete frozen source capsule")
    for name, pin in pins.items():
        require(sha(source[name]) == pin, "Frozen source pin " + name)
        require(
            source[name] == subprocess.check_output(["git", "show", SOURCE_COMMIT + ":" + name]),
            "Executed Git bytes " + name,
        )
    execution = json.loads(source[CONTRACT])
    base_raw = source[execution["base_guard_contract"]]
    require(sha(base_raw) == execution["base_guard_sha256"], "Original base contract")
    base = json.loads(base_raw)
    inputs = base["input_sha256"] | execution["input_sha256"]
    method = base["source_sha256"] | execution["source_sha256"]
    refs = collection["reference_files"]
    require(set(refs) == set(inputs) and len(inputs) == 13, "Exact input inventory")
    for name, pin in inputs.items():
        require(hashes[refs[name]] == pin, "Original input " + name)
    require(
        source_prerequisites(source, method) == [] and len(method) == 9, "All source prerequisites"
    )
    require(
        sha(source["projects/ukrmol_co/stretched_real_completion.py"])
        == "e79d7885650da54cfc8f00a6c77352a7e09aa37f9aa18a86d99eca2f1bd71a24",
        "Unchanged numerical runner",
    )
    original = json.loads(
        source["projects/ukrmol_co/stretched-real-completion-contract-20261102.json"]
    )
    changed = {n for n in set(original) | set(execution) if original.get(n) != execution.get(n)}
    require(
        changed == {"status", "container_name", "input_sha256", "failed_stop_reading"},
        "Metadata-only release",
    )
    policy = {
        n: base[n]
        for n in [
            "real_s_absolute_budget",
            "real_invariant_absolute_budget",
            "real_phase_sum_budget_rad",
            "real_flux_absolute_budget",
        ]
    }
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
        and owner["source_archive_sha256"] == hashes["source.tar.gz"]
        and owner["owner_sha256"]
        == sha(source["projects/ukrmol_co/compressed_cas12_pole_owner.py"])
        and owner["source_integrity_after"]
        and not owner["timed_out"]
        and not owner["container_oom_killed"],
        "Original owner/source integrity",
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
        exits == [0] * 6
        and container["Image"] == execution["image"]
        and not container["State"]["Running"],
        "Original successful termination",
    )
    host = container["HostConfig"]
    mounts = {m["Destination"]: m for m in container["Mounts"]}
    require(
        host["CpusetCpus"] == "0-3"
        and host["Memory"] == host["MemorySwap"] == 8 * 2**30
        and host["NetworkMode"] == "none"
        and set(mounts) == {"/work", "/opt/qmodeling", "/out"}
        and not mounts["/work"]["RW"]
        and not mounts["/opt/qmodeling"]["RW"],
        "Runtime caps/mounts",
    )
    for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
        require(name + "=1" in container["Config"]["Env"], "Library threads")
    controls = ast.parse(source["projects/ukrmol_co/test_stretched_real_completion.py"])
    count = sum(
        isinstance(n, ast.FunctionDef) and n.name.startswith("test_") for n in controls.body
    )
    require(
        count == 7
        and f"Ran {count} tests".encode() in files["outputs/controls.log"]
        and load("outputs/controls.json")["original_exit"] == 0,
        "Runtime source-closure and synthetic controls",
    )
    previous = load(refs[execution["previous_public_reading"]])
    stopped = load(refs[execution["failed_stop_reading"]])
    old = load(refs[execution["original_guard_result"]])
    require(
        previous["reader_exit"] == 0
        and previous["result"]["contour_study_verdict"] is True
        and previous["result"]["original_exits"] == [0] * 6
        and stopped["reader_exit"] == 0
        and stopped["result"]["original_exits"] == [2] * 6
        and stopped["result"]["maps_reserved"] == 0
        and old["scientific_exit"] == 2
        and old["real_guard_verdict"] is None,
        "Preserved prior passes/stops",
    )
    for case, spec in base["models"].items():
        receipt = load(refs[spec["receipt"]])
        require(
            receipt["case"] == case
            and receipt["inner_export_verdict"] is True
            and receipt["public_reader_exit"] == 0
            and receipt["original_exits"] == [0] * 6
            and (receipt["roots"], receipt["target_count"], receipt["boundary_channels"])
            == (2048, 64, 400),
            "Own-model public native export",
        )
    require(len(result["real_rows"]) == result["real_maps_reserved"] == 5, "Five new maps exactly")
    rows = []
    tz80 = []
    for row, (case, index, endpoint) in zip(result["real_rows"], MAPS, strict=True):
        require(
            (row["case"], row["index"], row["endpoint"], row["completed"])
            == (case, index, endpoint, True),
            "Fixed real identity/endpoint",
        )
        energy, q_ref, native = native_at(files, refs, base["models"][case], index)
        with np.load(
            io.BytesIO(files[f"outputs/check/{case}-{index}-{endpoint}.npz"]), allow_pickle=False
        ) as z:
            require(all(np.isfinite(z[n]).all() for n in z.files), "Finite saved matrices")
            q, y, a, s, reference = [z[n] for n in ["q", "y", "a", "s", "reference"]]
            np.testing.assert_allclose(float(z["energy_rydberg"]), energy, atol=1e-14, rtol=0)
        require(
            q.shape == y.shape == (400, 400)
            and a.shape == (6, 400)
            and s.shape == reference.shape == (6, 6),
            "Saved shapes",
        )
        np.testing.assert_allclose(q, q_ref, atol=1e-9, rtol=1e-10)
        np.testing.assert_allclose(reference, native, atol=1e-12, rtol=1e-12)
        actual, metrics = measurements(q, y, a, native)
        np.testing.assert_allclose(actual, s, atol=1e-10, rtol=1e-10)
        for name, value in metrics.items():
            np.testing.assert_allclose(
                value, row.get(name, row["invariants"].get(name)), atol=1e-13, rtol=1e-5
            )
        require(
            all(checks(metrics, policy).values()) == all(row["checks"].values()), "Saved real gates"
        )
        if index == 60 or endpoint == 120:
            require(all(checks(metrics, policy).values()), "Standalone/final physical gates")
        d = row["diagnostics"]
        require(
            (d["endpoint"], d["angle"], d["order"], d["rtol"]) == (endpoint, 20, 10, 1e-12)
            and d["rhs_evaluations"] <= 50000
            and d["wall_seconds"] < 90,
            "Exact prescription/caps",
        )
        rows.append(
            dict(
                case=case,
                index=index,
                endpoint=endpoint,
                metrics=metrics,
                checks=checks(metrics, policy),
            )
        )
        if case == "tz" and index == 80:
            tz80.append((y, a))
    stability = dict(
        y_absolute=float(abs(tz80[-1][0] - tz80[-2][0]).max()),
        a_absolute=float(abs(tz80[-1][1] - tz80[-2][1]).max()),
    )
    for name, value in stability.items():
        np.testing.assert_allclose(
            value, result["tz80_last_two_stability"][name], atol=1e-15, rtol=1e-10
        )
    require(
        max(stability.values()) <= 1e-7 and result["real_completion_verdict"] is True,
        "TZ80 endpoint stability",
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
        "Protected scientific scope",
    )
    samples = list(csv.DictReader(io.StringIO(files["outputs/check/resources.csv"].decode())))
    require(
        max(int(r["container_memory_bytes"]) for r in samples)
        == resources["sampled_peak_container_memory_bytes"],
        "Resource sample peak",
    )
    require(
        result["science_wall_seconds"] < 600
        and resources["wall_seconds"] < 600
        and owner["wall_seconds"] < 900,
        "Science/owner caps",
    )
    monitor = [json.loads(line) for line in files["monitor.jsonl"].decode().splitlines()]
    times = [datetime.datetime.fromisoformat(owner["started_utc"])] + [
        datetime.datetime.fromisoformat(r["observed_utc"]) for r in monitor
    ]
    gap = max((b - a).total_seconds() for a, b in itertools.pairwise(times))
    require(
        gap <= 60 and monitor[-1]["owner_process_exit"] == 0,
        "Observation cadence and terminal witness",
    )
    return dict(
        real_completion_verdict=True,
        combined_20degree_real_prerequisite_verdict=True,
        original_exits=exits,
        rows=rows,
        tz80_last_two_stability=stability,
        source_commit=SOURCE_COMMIT,
        source_entries_verified=36,
        raw_files_verified=len(hashes),
        reference_inputs_verified=13,
        source_prerequisites_verified=9,
        runtime_controls_verified=count,
        resource_samples_verified=len(samples),
        observation_count=len(monitor),
        maximum_observation_gap_seconds=gap,
        resources=resources,
        finished_utc=owner["finished_utc"],
        previous_dz80_repeated=False,
        zero_map_failure_original_exits=[2] * 6,
        original40degree_scientific_exit=2,
        new_maps=0,
        native_calls=0,
        complex_use_released=False,
        automatic_successor=False,
        full_sweep_released=False,
        scope=(
            "Same-model real-axis development60/80 for each finite own-DZ/TZ export "
            "plus endpoint controls only;not complex,pole,basis sensitivity or electronic accuracy"
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
