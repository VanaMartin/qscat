"""Authenticate saved authorized contour-study maps and native K; no new ODE."""

from __future__ import annotations

import argparse
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
from projects.ukrmol_co.verify_stretched_real_guard import checks, measurements

CONTRACT = "projects/ukrmol_co/stretched-contour-study-contract-20261102.json"
CONTRACT_SHA = "88dadad7934cf4b062d278ded22d8dc96f0bad11e7be055b3234f5be4bece5c6"
SOURCE_COMMIT = "102b5235ac24747887409d05d9745678723f7a08"
ARCHIVE_SHA = "a561209ab4d784c45dc878c95c1754e52d062739230574d6e565100f9a479f51"


def verify(archive, collection):
    require(collection["archive"]["sha256"] == ARCHIVE_SHA, "Unexpected allocation archive")
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
    require(set(source) == set(pins) and len(pins) == 33, "Source inventory")
    for name, pin in pins.items():
        require(sha(source[name]) == pin, "Source pin " + name)
        require(
            source[name] == subprocess.check_output(["git", "show", SOURCE_COMMIT + ":" + name]),
            "Executed Git bytes " + name,
        )
    require(sha(source[CONTRACT]) == CONTRACT_SHA, "Executed contract")
    c = json.loads(source[CONTRACT])
    refs = collection["reference_files"]
    require(c["execution_authorized"] and set(refs) == set(c["input_sha256"]), "Authorized inputs")
    for name, pin in c["input_sha256"].items():
        require(hashes[refs[name]] == pin, "Original input " + name)
    for name, pin in c["source_sha256"].items():
        require(sha(source[name]) == pin, "Method pin " + name)
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
        and owner["source_integrity_after"],
        "Executed owner/source",
    )
    require(
        not owner["timed_out"]
        and not owner["container_oom_killed"]
        and not container["State"]["Running"]
        and container["Image"] == c["image"],
        "Container termination/image",
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
        require(name + "=1" in container["Config"]["Env"], "Library thread cap")
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
        and load("outputs/controls.json")["original_exit"] == 0
        and b"Ran 5 tests" in files["outputs/controls.log"],
        "Original exits/controls",
    )
    receipt = load(refs[c["model"]["receipt"]])
    coupled = load(refs[c["coupled_reading"]])
    original = load(refs[c["original_guard_result"]])
    require(
        receipt["case"] == "dz"
        and receipt["inner_export_verdict"] is True
        and receipt["public_reader_exit"] == 0
        and receipt["original_exits"] == [0] * 6
        and (receipt["roots"], receipt["target_count"], receipt["boundary_channels"])
        == (2048, 64, 400),
        "Own public native receipt",
    )
    require(
        coupled["reader_exit"] == 0
        and coupled["result"]["coupled_control_verdict"] is True
        and original["scientific_exit"] == 2
        and original["real_guard_verdict"] is None,
        "Coupled prerequisite/original failure",
    )
    prefix = c["model"]["native_directory"] + "/"
    boundary = records(files[refs[prefix + "fort.21"]])
    channel = records(files[refs[prefix + "fort.10"]])
    km = records(files[refs[prefix + "fort.921"]])
    poles = numbers(boundary[6], "<f8", 2048)
    w = numbers(boundary[7], "<f8", 2048 * 400).reshape(2048, 400)
    ground = struct.unpack_from("<d", channel[4], 32)[0]
    native_values = numbers(km[85][24:], "<f8", 22)
    energy = native_values[0]
    k = symmetric(native_values[1:], 6)
    native = np.linalg.solve(np.eye(6) - 1j * k, np.eye(6) + 1j * k)
    q_reference = w.T @ (w / (poles - ground - energy / 2)[:, None])
    rows = []
    captures = []
    require(len(result["real_rows"]) == 3 and result["real_maps_reserved"] == 3, "Fixed map count")
    for row, endpoint in zip(result["real_rows"], [60, 90, 120], strict=True):
        require(
            row["index"] == 80 and row["endpoint"] == endpoint and row["completed"], "Frozen row"
        )
        with np.load(
            io.BytesIO(files[f"outputs/check/endpoint-{endpoint}.npz"]), allow_pickle=False
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
        np.testing.assert_allclose(q, q_reference, atol=1e-9, rtol=1e-10)
        np.testing.assert_allclose(reference, native, atol=1e-12, rtol=1e-12)
        actual, metrics = measurements(q, y, a, native)
        np.testing.assert_allclose(actual, s, atol=1e-10, rtol=1e-10)
        for name, value in metrics.items():
            expected = row.get(name, row["invariants"].get(name))
            np.testing.assert_allclose(value, expected, atol=1e-13, rtol=1e-5)
        require(
            all(checks(metrics, c).values()) == all(row["checks"].values()), "Saved gate arithmetic"
        )
        d = row["diagnostics"]
        require(
            (d["endpoint"], d["angle"], d["order"], d["rtol"]) == (endpoint, 20, 10, 1e-12)
            and d["rhs_evaluations"] <= 50000
            and d["wall_seconds"] < 90,
            "Exact numerical settings/caps",
        )
        rows.append(dict(endpoint=endpoint, metrics=metrics, checks=checks(metrics, c)))
        captures.append((y, a))
    stability = dict(
        y_absolute=float(abs(captures[-1][0] - captures[-2][0]).max()),
        a_absolute=float(abs(captures[-1][1] - captures[-2][1]).max()),
    )
    for name, value in stability.items():
        np.testing.assert_allclose(
            value, result["last_two_stability"][name], atol=1e-15, rtol=1e-10
        )
    require(
        all(checks(rows[-1]["metrics"], c).values())
        and max(stability.values()) <= 1e-7
        and result["contour_study_verdict"] is True,
        "Final physical/stability gates",
    )
    require(
        result["native_calls"]
        == result["new_complex_maps"]
        == result["route_holdout_s_matrices_analyzed"]
        == 0
        and not result["tz_started"]
        and not result["original_guard_requalified"]
        and not result["full_sweep_released"],
        "Protected scope",
    )
    samples = list(csv.DictReader(io.StringIO(files["outputs/check/resources.csv"].decode())))
    require(
        max(int(r["container_memory_bytes"]) for r in samples)
        == resources["sampled_peak_container_memory_bytes"],
        "Measured memory samples",
    )
    require(
        result["science_wall_seconds"] < 300
        and resources["wall_seconds"] < 300
        and owner["wall_seconds"] < 600,
        "Science/owner caps",
    )
    monitor = [json.loads(line) for line in files["monitor.jsonl"].decode().splitlines()]
    times = [datetime.datetime.fromisoformat(owner["started_utc"])] + [
        datetime.datetime.fromisoformat(r["observed_utc"]) for r in monitor
    ]
    maximum_gap = max((b - a).total_seconds() for a, b in itertools.pairwise(times))
    require(
        maximum_gap <= 60 and monitor[-1]["owner_process_exit"] == 0,
        "Observation cadence/terminal witness",
    )
    return dict(
        contour_study_verdict=True,
        original_exits=exits,
        rows=rows,
        last_two_stability=stability,
        source_commit=SOURCE_COMMIT,
        source_entries_verified=len(pins),
        raw_files_verified=len(hashes),
        reference_inputs_verified=len(refs),
        resource_samples_verified=len(samples),
        observation_count=len(monitor),
        maximum_observation_gap_seconds=maximum_gap,
        resources=resources,
        finished_utc=owner["finished_utc"],
        original_guard_scientific_exit=2,
        original_guard_requalified=False,
        native_calls=0,
        new_maps=0,
        tz_started=False,
        automatic_successor=False,
        full_sweep_released=False,
        scope=(
            "Own-DZ development80 real-axis point and endpoint convergence only; "
            "no independent exterior solution or basis/pole/complex/model verdict"
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
