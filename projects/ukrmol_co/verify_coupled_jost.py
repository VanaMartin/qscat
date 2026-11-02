"""Independent fundamental-seed reading of saved coupled maps; never integrates."""

from __future__ import annotations

import argparse
import itertools
import json
import subprocess
from pathlib import Path

import numpy as np

from projects.ukrmol_co.verify_dual_transport_toy import digest, members, require

SOURCE_COMMIT = "16f76fb78d95c9fe1055328f2bd361dcf5a4dd0d"
ARCHIVE_SHA = "bd830ce069eac7bb9030e95fca6dbb9eb586e330607683731bc550a7ac154066"


def fundamental_maps(z, g):
    k = np.asarray([1, np.sqrt(3)])
    c = np.diag(k) + g * np.ones((2, 2))
    u = np.diag(np.cosh(k * z)) + np.diag(np.sinh(k * z) / k) @ c
    du = np.diag(k * np.sinh(k * z)) + np.diag(np.cosh(k * z)) @ c
    w = du @ np.linalg.inv(u)
    dw = np.diag([1, 3]) - w @ w
    s = np.diag([1j, -1])
    normalization = np.diag(1 / (np.diag(s) - k))
    exponential = np.diag(np.exp(np.diag(s) * z))
    f = (s - w) @ exponential @ normalization
    df = (-dw + (s - w) @ s) @ exponential @ normalization
    inverse = np.linalg.inv(f)
    return df @ inverse, inverse[:1]


def decoded(value):
    array = np.asarray(value)
    return array[..., 0] + 1j * array[..., 1]


def measurements(y, a, reference):
    # One open row admits scalar least squares, independent of scipy/numpy lstsq.
    wall = np.sum(a * a) / np.sum(abs(a) ** 2)
    return dict(
        map_error=float(max(abs(y - reference[0]).max(), abs(a - reference[1]).max())),
        flux_defect=float(abs((y - y.conj().T) / (2j) - a.conj().T @ a).max()),
        wronskian_residual=float(abs(a.conj().T * wall - a.T).max()),
        s_identity_error=float(abs(wall - 1)),
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
    source, pins = members(files["source.tar.gz"]), json.loads(files["source-manifest.json"])
    require(set(source) == set(pins) and len(pins) == 7, "Minimal executed source")
    for name, raw in source.items():
        require(digest(raw) == pins[name], "Source pin " + name)
        require(
            raw == subprocess.check_output(["git", "show", SOURCE_COMMIT + ":" + name]),
            "Executed Git bytes " + name,
        )
    frozen, execution = json.loads(files["freeze.json"]), json.loads(files["execution.json"])
    require(frozen["source_commit"] == execution["source_commit"] == SOURCE_COMMIT, "Commit")
    require(
        frozen["source_archive_sha256"]
        == execution["source_archive_sha256"]
        == digest(files["source.tar.gz"])
        and execution["source_integrity_after"],
        "Source freeze",
    )
    contract = json.loads(source["projects/ukrmol_co/coupled-jost-contract-20261102.json"])
    container = json.loads(files["container-inspect.json"])
    require(container["Image"] == contract["image"], "Pinned image")
    host = container["HostConfig"]
    require(
        host["CpusetCpus"] == "0-3" and host["Memory"] == host["MemorySwap"] == 2**31,
        "Physical CPU/memory caps",
    )
    require(host["NetworkMode"] == "none", "No network")
    mounts = {m["Destination"]: m for m in container["Mounts"]}
    require(set(mounts) == {"/source", "/out"} and not mounts["/source"]["RW"], "No CO mount")
    for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
        require(name + "=1" in container["Config"]["Env"], "Library thread cap")
    result = json.loads(files["outputs/control/result.json"])
    exits = [
        result["scientific_exit"],
        json.loads(files["outputs/controls.json"])["original_exit"],
        execution["original_command_exit"],
        execution["original_container_exit"],
        execution["original_owner_exit"],
        json.loads(files["owner-process-exit.json"])["original_exit"],
    ]
    require(
        exits == [0] * 6 and not container["State"]["OOMKilled"] and not execution["timed_out"],
        "Original exits/OOM/timeout",
    )
    reference = fundamental_maps(contract["radius_bohr"], contract["g"])
    grid = list(
        itertools.product(
            contract["angles_degrees"],
            contract["lengths_bohr"],
            contract["rtols"],
            contract["initializations"],
        )
    )
    require(len(result["rows"]) == len(grid) == 54, "Exact integration count")
    exact_errors = {str(tol): 0.0 for tol in contract["rtols"]}
    free_errors, finest = {}, []
    for row, case in zip(result["rows"], grid, strict=True):
        require(
            tuple(row[n] for n in ["angle_degrees", "length_bohr", "rtol", "initialization"])
            == case,
            "Ordered frozen case grid",
        )
        angle, length, rtol, initialization = case
        y, a = decoded(row["y"]), decoded(row["a"])
        require(
            y.shape == (2, 2)
            and a.shape == (1, 2)
            and np.isfinite(y).all()
            and np.isfinite(a).all(),
            "Finite map shapes",
        )
        measured = measurements(y, a, reference)
        for name, value in measured.items():
            require(abs(value - row[name]) <= 3e-14, "Saved metric " + name)
        require(0 < row["rhs_evaluations"] <= contract["rhs_cap"], "RHS cap")
        if initialization == "exact":
            exact_errors[str(rtol)] = max(exact_errors[str(rtol)], measured["map_error"])
        if rtol == contract["rtols"][-1]:
            finest.append(
                dict(
                    angle_degrees=angle,
                    length_bohr=length,
                    initialization=initialization,
                    **measured,
                )
            )
            if initialization == "free":
                free_errors.setdefault(str(angle), []).append(measured["map_error"])
            if initialization == "exact" or length == contract["lengths_bohr"][-1]:
                require(
                    all(
                        measured[n] <= contract[b]
                        for n, b in [
                            ("map_error", "map_error_budget"),
                            ("flux_defect", "flux_budget"),
                            ("wronskian_residual", "wronskian_budget"),
                        ]
                    ),
                    "Frozen finest gate",
                )
    for values in free_errors.values():
        require(
            values[0] > contract["minimum_resolved_tail_error"]
            and values[0] >= contract["minimum_tail_error_reduction"] * values[-1],
            "Tail separation",
        )
    require(
        result["coupled_control_verdict"] is True
        and result["science_wall_seconds"] < 300
        and execution["owner_wall_seconds"] < 360,
        "Verdict/time caps",
    )
    require(
        result["co_maps"] == result["native_calls"] == 0 and not result["automatic_successor"],
        "Scope",
    )
    return dict(
        source_commit=SOURCE_COMMIT,
        archive_sha256=ARCHIVE_SHA,
        raw_files_verified=12,
        executed_git_sources_verified=7,
        original_exits=exits,
        integrations_verified=54,
        exact_start_maximum_map_error_by_rtol=exact_errors,
        finest_free_start_map_error_by_angle_and_length=free_errors,
        finest_rows=finest,
        science_wall_seconds=result["science_wall_seconds"],
        owner_wall_seconds=execution["owner_wall_seconds"],
        recorded_cgroup_memory_peak_bytes=int(result["cgroup_memory.peak"]),
        coupled_control_verdict=True,
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
