"""Saved-only full-Slater-minor reading of the finite DZ/TZ target screen."""

from __future__ import annotations

import argparse
import csv
import datetime
import io
import itertools
import json
import subprocess
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from projects.ukrmol_co.verify_stretched_preparation import stream_files
from projects.ukrmol_co.verify_stretched_properties import files_in, require, sha

SOURCE = "d4301b3ff012a37b7e29c0d7467941a49d6ff177"
SOURCE_ARCHIVE = "a8eb8e0391db041e5877eb77c285bd2c74a2d14aaefe9fda0ca8fd407e6c354b"
CONTRACT = "projects/ukrmol_co/stretched-target-overlap-contract-20261102.json"


def full_overlap(left, right, metric, electrons):
    """Full occupied determinants; no frozen-core Schur reduction or normalization."""
    minors = []
    for count in electrons:
        masks = sorted(i for i in range(1 << 12) if i.bit_count() == count)
        occupations = np.array(
            [[0, 1] + [j + 2 for j in range(12) if mask & (1 << j)] for mask in masks]
        )
        minors.append(
            np.linalg.det(metric[occupations[:, None, :, None], occupations[None, :, None, :]])
        )
    transformed = np.stack([minors[0] @ ci @ minors[1].T for ci in right])
    return left.reshape(8, -1).conj() @ transformed.reshape(8, -1).T


def verify(archive, collection):
    captures = [n for n in collection["files"] if n.endswith(".npz") or n.endswith("ci.npy")]
    files, hashes, _ = stream_files(archive, collection, capture_names=captures)

    def load(name):
        return json.loads(files[name])

    source = files_in(files["source.tar.gz"])
    pins = load("source-manifest.json")
    require(
        len(pins) == 62 and set(source) == set(pins) and hashes["source.tar.gz"] == SOURCE_ARCHIVE,
        "Frozen source inventory",
    )
    for name, pin in pins.items():
        require(
            sha(source[name]) == pin
            and source[name] == subprocess.check_output(["git", "show", SOURCE + ":" + name]),
            "Executed Git source " + name,
        )
    c = json.loads(source[CONTRACT])
    refs = collection["reference_files"]
    manifest = load("outputs/check/input-manifest.json")
    require(
        c["execution_authorized"]
        and len(refs) == len(manifest) == 41
        and set(refs) == set(manifest),
        "41 exact saved inputs",
    )
    require(
        hashes[refs[c["feasibility_collection"]]] == c["feasibility_collection_sha256"],
        "Original public input inventory",
    )
    original = load(refs[c["feasibility_collection"]])["files"]
    for name, pin in c["source_sha256"].items():
        require(sha(source[name]) == pin, "Source prerequisite " + name)
    models = json.loads(source[c["property_contract"]])["models"]
    expected = {c["feasibility_collection"]: c["feasibility_collection_sha256"]}
    for model in models.values():
        for name, pin in model["input_sha256"].items():
            if name == "target.json" or (
                name.startswith("output/CO/geom1/")
                and Path(name).name
                in [
                    "co.casscf.chk",
                    "co.molden",
                    "fort.26",
                    *[f"fort.{u}" for u in [710, 730, 711, 731, 712, 732, 713, 733]],
                ]
            ):
                expected[model["target_directory"] + "/" + name] = pin
    for spin in [0, 2]:
        for irrep in ["A1", "B1", "B2", "A2"]:
            for filename in ["ci.npy", "result.json"]:
                name = c["tz_reference_directory"] + f"/spin{spin}-{irrep}-space240/" + filename
                expected[name] = original[name]["sha256"]
    require(expected == manifest, "Original pin derivation")
    for name, pin in expected.items():
        require(hashes[refs[name]] == pin, "Saved raw input " + name)
    owner, result, resources, container = [
        load(n)
        for n in [
            "execution.json",
            "outputs/check/result.json",
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
        exits == [exits[0]] * 6
        and exits[0] in [0, 1, 2]
        and owner["source_commit"] == SOURCE
        and owner["source_archive_sha256"] == SOURCE_ARCHIVE
        and owner["source_integrity_after"]
        and owner["owner_sha256"] == sha(source["projects/ukrmol_co/stretched_dz_pole_owner.py"])
        and not owner["timed_out"]
        and not owner["container_oom_killed"],
        "Original allocation provenance/exits",
    )
    host = container["HostConfig"]
    mounts = {m["Destination"]: m for m in container["Mounts"]}
    require(
        container["Image"] == c["image"]
        and not container["State"]["Running"]
        and host["CpusetCpus"] == "0-3"
        and host["Memory"] == host["MemorySwap"] == 8 * 2**30
        and host["NetworkMode"] == "none"
        and set(mounts) == {"/work", "/opt/qmodeling", "/out"}
        and not mounts["/work"]["RW"]
        and not mounts["/opt/qmodeling"]["RW"],
        "Four physical-core immutable allocation",
    )
    require(
        load("outputs/controls.json")["original_exit"] == 0
        and b"Ran 10 tests" in files["outputs/controls.log"],
        "Original10 runtime controls",
    )
    for name in [
        "native_calls",
        "new_ci_solves",
        "new_hamiltonian_actions",
        "new_maps",
        "route_holdout_s_matrices_analyzed",
    ]:
        require(result[name] == 0, "Protected scope " + name)
    require(
        not result["automatic_successor"] and not result["full_sweep_released"],
        "No successor or production release",
    )

    def arrays(name):
        with np.load(io.BytesIO(files["outputs/check/" + name]), allow_pickle=False) as z:
            return {n: z[n] for n in z.files}

    metrics = [arrays(case + "-metric.npz") for case in ["dz", "tz"]]
    cross = arrays("cross-metric.npz")
    for data in metrics:
        np.testing.assert_allclose(
            data["mo"].conj().T @ data["ao"] @ data["mo"], np.eye(14), atol=1e-9, rtol=0
        )
    np.testing.assert_allclose(cross["ao"], cross["reverse_ao"].conj().T, atol=1e-12, rtol=0)
    metric = metrics[0]["mo"].conj().T @ cross["ao"] @ metrics[1]["mo"]
    np.testing.assert_allclose(metric, cross["mo"], atol=1e-12, rtol=0)
    rows = []
    for row in result["sectors"]:
        label = row["sector"]
        spin = 0 if label.startswith("singlet") else 1
        data = arrays(label + ".npz")
        left, right = data["dz"], data["tz"]
        for vector in [left, right]:
            np.testing.assert_allclose(
                vector.reshape(8, -1).conj() @ vector.reshape(8, -1).T, np.eye(8), atol=1e-9, rtol=0
            )
        overlap = full_overlap(left, right, metric, (5 + spin, 5 - spin))
        for name in ["overlap", "full_minors"]:
            np.testing.assert_allclose(overlap, data[name], atol=1e-10, rtol=0)
        np.testing.assert_allclose(overlap, data["reverse"].conj().T, atol=1e-10, rtol=0)
        weights = np.linalg.eigvalsh(overlap.conj().T @ overlap)[::-1]
        np.testing.assert_allclose(weights, row["squared_principal_overlaps"], atol=1e-10, rtol=0)
        require(
            bool(weights[-1] >= 0.5) == row["retained_eight_dominance_passed"],
            "Predeclared finite dominance screen",
        )
        rr, cc = linear_sum_assignment(-(abs(overlap) ** 2))
        for i, j, saved in zip(rr, cc, row["assignment"], strict=True):
            require(
                (int(i + 1), int(j + 1)) == (saved["dz_root"], saved["tz_root"]),
                "Root assignment reconstruction",
            )
            np.testing.assert_allclose(
                abs(overlap[i, j]) ** 2, saved["squared_overlap"], atol=1e-10, rtol=0
            )
        # Same-frame independent TZ CI is a representation oracle, not a new solve.
        order = metrics[1]["original_order"]
        permutation = np.eye(14)[order]
        name = (
            c["tz_reference_directory"] + f"/spin{2 * spin}-{label.split('.')[1]}-space240/ci.npy"
        )
        reference = np.load(io.BytesIO(files[refs[name]]), allow_pickle=False)
        conversion = full_overlap(right, reference, permutation, (5 + spin, 5 - spin))
        np.testing.assert_allclose(conversion, data["conversion_overlap"], atol=1e-10, rtol=0)
        require(
            min(abs(np.diag(conversion)) ** 2) >= 1 - 1e-8,
            "Saved independent TZ coefficient correspondence",
        )
        rows.append(
            dict(
                sector=label,
                minimum_squared_principal_overlap=float(weights[-1]),
                retained_eight_dominance_passed=bool(weights[-1] >= 0.5),
                maximum_overlap_reconstruction_error=float(abs(overlap - data["overlap"]).max()),
                dz_capture_in_tz=(abs(overlap) ** 2).sum(axis=1).tolist(),
                tz_capture_in_dz=(abs(overlap) ** 2).sum(axis=0).tolist(),
                minimum_assigned_squared_overlap=float(min(abs(overlap[rr, cc]) ** 2)),
                assignments_preserve_rank=bool(np.array_equal(rr, cc)),
            )
        )
    verdict = None
    if exits[0] in [0, 1]:
        require(len(rows) == 8, "Complete64-state/eight-sector screen")
        verdict = all(row["retained_eight_dominance_passed"] for row in rows)
        require(
            verdict == result["retained_eight_dominance_verdict"]
            and exits[0] == (0 if verdict else 1),
            "Original physical disposition",
        )
    samples = list(csv.DictReader(io.StringIO(files["outputs/check/resources.csv"].decode())))
    require(
        max(int(s["container_memory_bytes"]) for s in samples)
        == resources["sampled_peak_container_memory_bytes"],
        "Original resource samples",
    )
    monitor = [json.loads(line) for line in files["monitor.jsonl"].decode().splitlines()]
    times = [datetime.datetime.fromisoformat(owner["started_utc"])] + [
        datetime.datetime.fromisoformat(m["observed_utc"]) for m in monitor
    ]
    gap = max((b - a).total_seconds() for a, b in itertools.pairwise(times))
    require(
        gap <= 60
        and monitor[-1]["owner_process_exit"] == exits[0]
        and resources["wall_seconds"] < 1800
        and owner["wall_seconds"] < 2100,
        "Caps and observation cadence",
    )
    return dict(
        original_exits=exits,
        source_commit=SOURCE,
        raw_files_verified=len(hashes),
        executed_git_sources_verified=62,
        saved_input_pins_verified=41,
        source_prerequisites_verified=3,
        sectors=rows,
        retained_eight_dominance_verdict=verdict,
        error=result.get("error"),
        resources=resources,
        resource_samples_verified=len(samples),
        observation_count=len(monitor),
        maximum_observation_gap_seconds=gap,
        finished_utc=owner["finished_utc"],
        native_calls=0,
        new_ci_solves=0,
        new_hamiltonian_actions=0,
        new_maps=0,
        full_sweep_released=False,
        scope=(
            "Full occupied Slater minors/eigvalsh replay from authenticated saved "
            "converted CI and AO matrices;native-CSF expansion and AO integration are not "
            "independently regenerated;no N+1 identity or accuracy verdict"
        ),
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
