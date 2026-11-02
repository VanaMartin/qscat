"""Reconstruct the exact toy from saved rows with scalar formulas, never an ODE."""

from __future__ import annotations

import argparse
import hashlib
import io
import itertools
import json
import math
import subprocess
import tarfile
from pathlib import Path

SOURCE_COMMIT = "708d3eb10e32a71f809681bb0e71fb773562654f"
ARCHIVE_SHA = "cb015867b097191fadc7cae7f7b67a92310597b5eb9745897eab173daf5f532d"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def members(data):
    with tarfile.open(fileobj=io.BytesIO(data)) as archive:
        rows = archive.getmembers()
        require(all(row.isfile() for row in rows), "Non-file archive entry")
        result = {row.name: archive.extractfile(row).read() for row in rows}
        require(len(result) == len(rows), "Duplicate archive entries")
        return result


def scalar_oracle(angle, length, epsilon):
    theta = math.radians(angle)
    magnitude = epsilon * math.exp(length * (math.sin(theta) - math.cos(theta)))
    phase = -length * (math.sin(theta) + math.cos(theta))
    return [complex(1), magnitude * complex(math.cos(phase), math.sin(phase))]


def measurements(actual, expected, floor):
    ao, ac = actual
    s = (ao * ao + ac * ac) / (abs(ao) ** 2 + abs(ac) ** 2)
    flux = max(abs(1 - abs(ao) ** 2), abs(ao) * abs(ac), abs(ac) ** 2)
    return dict(
        scaled_error=max(
            abs(a - e) / (floor + abs(e)) for a, e in zip(actual, expected, strict=True)
        ),
        s_identity_error=abs(s - 1),
        flux_defect=flux,
    )


def verify(archive, collection):
    data = archive.read_bytes()
    require(digest(data) == ARCHIVE_SHA == collection["archive"]["sha256"], "Archive hash")
    require(len(data) == collection["archive"]["bytes"], "Archive size")
    files = members(data)
    require(set(files) == set(collection["files"]) and len(files) == 12, "Raw inventory")
    for name, raw in files.items():
        pin = collection["files"][name]
        require(digest(raw) == pin["sha256"] and len(raw) == pin["bytes"], "Raw pin " + name)
    source = members(files["source.tar.gz"])
    pins = json.loads(files["source-manifest.json"])
    require(set(source) == set(pins) and len(pins) == 7, "Minimal source inventory")
    for name, raw in source.items():
        require(digest(raw) == pins[name], "Source pin " + name)
        committed = subprocess.check_output(["git", "show", SOURCE_COMMIT + ":" + name])
        require(raw == committed, "Original Git bytes " + name)
    execution = json.loads(files["execution.json"])
    frozen = json.loads(files["freeze.json"])
    require(execution["source_commit"] == frozen["source_commit"] == SOURCE_COMMIT, "Commit")
    require(execution["source_integrity_after"], "Source changed")
    require(
        execution["source_archive_sha256"]
        == frozen["source_archive_sha256"]
        == digest(files["source.tar.gz"]),
        "Frozen source archive",
    )
    require(
        execution["owner_sha256"] == pins["projects/ukrmol_co/dual_transport_toy_owner.py"],
        "Executed owner",
    )
    require(not execution["campaign_input_mounted"] and not execution["timed_out"], "Scope/cap")
    container = json.loads(files["container-inspect.json"])
    contract = json.loads(source["projects/ukrmol_co/dual-transport-toy-contract-20261102.json"])
    require(
        contract["open_energy_rydberg"] == 1 and contract["closed_threshold_rydberg"] == 2,
        "Exact scalar momentum convention",
    )
    require(container["Image"] == contract["image"], "Pinned image")
    host = container["HostConfig"]
    require(
        host["CpusetCpus"] == "0-3" and host["Memory"] == host["MemorySwap"] == 2**31,
        "CPU/memory caps",
    )
    require(host["NetworkMode"] == "none", "Network")
    mounts = {m["Destination"]: m for m in container["Mounts"]}
    require(set(mounts) == {"/source", "/out"} and not mounts["/source"]["RW"], "No CO mount")
    for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
        require(name + "=1" in container["Config"]["Env"], "Library threads")
    controls = json.loads(files["outputs/controls.json"])
    result = json.loads(files["outputs/toy/result.json"])
    exits = [
        result["scientific_exit"],
        controls["original_exit"],
        execution["original_command_exit"],
        execution["original_container_exit"],
        execution["original_owner_exit"],
        json.loads(files["owner-process-exit.json"])["original_exit"],
    ]
    require(exits == [0] * 6 and not container["State"]["OOMKilled"], "Original exits/OOM")
    require(
        result["co_maps"] == result["native_calls"] == result["pole_searches"] == 0,
        "Scope counters",
    )
    grid = list(
        itertools.product(
            contract["angles_degrees"],
            contract["lengths_bohr"],
            contract["endpoint_contamination"],
            contract["rtols"],
        )
    )
    require(len(result["rows"]) == len(grid) == 72, "Exact integration count")
    maxima = {str(tol): 0.0 for tol in contract["rtols"]}
    finest = []
    for row, (angle, length, epsilon, rtol) in zip(result["rows"], grid, strict=True):
        require(
            [row[n] for n in ["angle_degrees", "length_bohr", "epsilon", "rtol"]]
            == [angle, length, epsilon, rtol],
            "Ordered case grid",
        )
        actual = [complex(*v) for v in row["actual"]]
        expected = scalar_oracle(angle, length, epsilon)
        require(all(math.isfinite(v.real) and math.isfinite(v.imag) for v in actual), "Nonfinite")
        measured = measurements(actual, expected, contract["scaled_error_floor"])
        for value, saved in zip(expected, row["exact"], strict=True):
            require(abs(value - complex(*saved)) <= 1e-14 * (1e-12 + abs(value)), "Oracle")
        for name in ["scaled_error", "s_identity_error", "flux_defect"]:
            # Independent scalar arithmetic may differ by ordinary absolute roundoff.
            require(abs(measured[name] - row[name]) <= 2e-14, "Saved metric " + name)
        require(0 < row["rhs_evaluations"] <= contract["rhs_cap"], "RHS cap")
        maxima[str(rtol)] = max(maxima[str(rtol)], measured["scaled_error"])
        if rtol == contract["rtols"][-1]:
            require(measured["scaled_error"] <= contract["finest_scaled_error_budget"], "Fine gate")
            if epsilon == 0:
                require(
                    abs(actual[1]) <= contract["zero_injection_closed_budget"], "Zero injection"
                )
            finest.append(
                dict(angle_degrees=angle, length_bohr=length, epsilon=epsilon, **measured)
            )
    negative = [
        r
        for r in finest
        if r["epsilon"]
        and r["s_identity_error"] <= contract["s_negative_control_budget"]
        and r["flux_defect"] > contract["flux_negative_control_budget"]
    ]
    require(len(negative) == result["negative_control_cases"] and negative, "Negative control")
    require(
        result["toy_verdict"] is True
        and result["wall_seconds"] < 300
        and execution["wall_seconds"] < 360,
        "Verdict/time caps",
    )
    return dict(
        source_commit=SOURCE_COMMIT,
        archive_sha256=ARCHIVE_SHA,
        raw_files_verified=12,
        git_sources_verified=7,
        original_exits=exits,
        integrations_verified=72,
        maximum_scaled_error_by_rtol=maxima,
        negative_control_cases=len(negative),
        science_wall_seconds=result["wall_seconds"],
        owner_wall_seconds=execution["wall_seconds"],
        recorded_cgroup_memory_peak_bytes=int(result["cgroup_memory.peak"]),
        finest_cases=finest,
        toy_mechanism_verdict=True,
        co_guard_repaired=False,
        co_maps=0,
        native_calls=0,
        automatic_successor=False,
        full_sweep_released=False,
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
