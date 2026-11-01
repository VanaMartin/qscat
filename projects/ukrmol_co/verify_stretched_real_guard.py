"""Reconstruct saved real S/flux checks with NumPy only; never integrates an exterior map."""

from __future__ import annotations

import argparse
import csv
import io
import json
import struct
from pathlib import Path

import numpy as np

from projects.ukrmol_co.verify_stretched_inner import disk_peak, numbers, records, symmetric
from projects.ukrmol_co.verify_stretched_preparation import stream_files
from projects.ukrmol_co.verify_stretched_properties import files_in, require, sha

CONTRACT = "projects/ukrmol_co/stretched-real-guard-contract-20261101.json"
CONTRACT_SHA = "bb9f9e38920c66125d1645799706b65a885fa916f649805f97d62bb400fd5abe"


def measurements(q, y, a, reference):
    wall, _, rank, _ = np.linalg.lstsq(a.conj().T, a.T, rcond=None)
    require(rank == 6, "Rank-deficient saved amplitude map")
    actual = wall + 2j * a @ np.linalg.solve(np.eye(len(q)) - q @ y, q) @ a.T
    metrics = dict(
        s_absolute_error=float(np.max(abs(actual - reference))),
        phase_sum_error_mod_pi_rad=float(
            abs(np.angle(np.linalg.slogdet(actual)[0] / np.linalg.slogdet(reference)[0])) / 2
        ),
        flux_defect=float(np.max(abs((y - y.conj().T) / (2j) - a.conj().T @ a))),
        wronskian_residual=float(np.max(abs(a.conj().T @ wall - a.T))),
        unitarity_defect=float(np.max(abs(actual.conj().T @ actual - np.eye(6)))),
        reciprocity_defect=float(np.max(abs(actual - actual.T))),
    )
    return actual, metrics


def checks(metrics, contract):
    keys = [
        ("s_absolute_error", "real_s_absolute_budget"),
        ("phase_sum_error_mod_pi_rad", "real_phase_sum_budget_rad"),
        ("flux_defect", "real_flux_absolute_budget"),
    ]
    keys += [
        (n, "real_invariant_absolute_budget")
        for n in ["wronskian_residual", "unitarity_defect", "reciprocity_defect"]
    ]
    return {name: metrics[name] <= contract[budget] for name, budget in keys}


def verify(path, collection, allocation):
    capture = [
        n
        for n in collection["files"]
        if n.endswith(".npz") or Path(n).name in ["fort.10", "fort.21", "fort.921"]
    ]
    files, hashes, _ = stream_files(path, collection, capture_names=capture)

    def load(n):
        return json.loads(files[n])

    source = files_in(files["source.tar.gz"])
    pins = load("source-manifest.json")
    require(set(source) == set(pins) and len(pins) == 31, "Wrong source capsule")
    for n, w in pins.items():
        require(sha(source[n]) == w, "Changed source: " + n)
    require(sha(source[CONTRACT]) == CONTRACT_SHA, "Changed executed real-guard contract")
    c = json.loads(source[CONTRACT])
    refs = collection["reference_files"]
    require(set(refs) == set(c["input_sha256"]), "Wrong input set")
    for n, w in c["input_sha256"].items():
        require(hashes[refs[n]] == w, "Changed original input: " + n)
    for n, w in c["source_sha256"].items():
        require(sha(source[n]) == w, "Changed method/policy prerequisite")
    owner = load("execution.json")
    result = load("outputs/check/result.json")
    resources = load("outputs/check/resources.json")
    container = load("container-inspect.json")
    exits = [
        result["scientific_exit"],
        resources["exit_code"],
        owner["original_command_exit"],
        owner["original_container_exit"],
        owner["original_owner_exit"],
        load("service-return.json")["original_service_exit"],
    ]
    require(
        owner["source_integrity_after"]
        and not owner["timed_out"]
        and not owner["container_oom_killed"]
        and container["Image"] == c["image"]
        and not container["State"]["Running"],
        "Invalid source/container termination",
    )
    require(
        owner["source_archive_sha256"] == hashes["source.tar.gz"]
        and owner["owner_sha256"]
        == sha(source["projects/ukrmol_co/compressed_cas12_pole_owner.py"]),
        "Wrong executed source",
    )
    require(
        load("outputs/controls.json")["original_exit"] == 0
        and b"Ran 6 tests" in files["outputs/controls.log"],
        "Failed controls",
    )
    rows = []
    for case, spec in c["models"].items():
        guard_path = "outputs/check/" + case + "/guard.json"
        if guard_path not in files:
            continue
        guard = load(guard_path)
        require(guard["case"] == "stretched/cas12/" + case + "/B1", "Wrong guard case")
        receipt = load(refs[spec["receipt"]])
        prefix = spec["native_directory"] + "/"
        require(
            receipt["case"] == case
            and receipt["inner_export_verdict"] is True
            and receipt["public_reader_exit"] == 0
            and receipt["original_exits"] == [0] * 6,
            "Unqualified own public input",
        )
        boundary = records(files[refs[prefix + "fort.21"]])
        channel = records(files[refs[prefix + "fort.10"]])
        km = records(files[refs[prefix + "fort.921"]])
        poles = numbers(boundary[6], "<f8", 2048)
        w = numbers(boundary[7], "<f8", 2048 * 400).reshape(2048, 400)
        ground = struct.unpack_from("<d", channel[4], 32)[0]
        for row in guard["real_rows"]:
            index = row["index"]
            require(index in [60, 80], "Not a frozen development identity")
            if not row["completed"]:
                rows.append(dict(case=case, index=index, completed=False))
                continue
            with np.load(
                io.BytesIO(files[f"outputs/check/{case}/real-{index}.npz"]), allow_pickle=False
            ) as z:
                require(all(np.isfinite(z[n]).all() for n in z.files), "Nonfinite saved matrix")
                q, y, a, s, reference = [z[n] for n in ["q", "y", "a", "s", "reference"]]
                energy = float(z["energy_rydberg"])
            require(
                q.shape == y.shape == (400, 400)
                and a.shape == (6, 400)
                and s.shape == reference.shape == (6, 6),
                "Wrong map shapes",
            )
            values = numbers(km[5 + index][24:], "<f8", 22)
            require(energy == values[0], "Wrong own native energy")
            np.testing.assert_allclose(
                q, w.T @ (w / (poles - ground - energy / 2)[:, None]), atol=1e-9, rtol=1e-10
            )
            k = symmetric(values[1:], 6)
            native = np.linalg.solve(np.eye(6) - 1j * k, np.eye(6) + 1j * k)
            np.testing.assert_allclose(reference, native, atol=1e-12, rtol=1e-12)
            actual, metrics = measurements(q, y, a, native)
            np.testing.assert_allclose(actual, s, atol=1e-10, rtol=1e-10)
            for n, v in metrics.items():
                expected = row.get(n, row["invariants"].get(n))
                np.testing.assert_allclose(v, expected, atol=1e-10, rtol=1e-5)
            passed = all(checks(metrics, c).values())
            require(passed == all(row["checks"].values()), "Stored gate arithmetic differs")
            rows.append(
                dict(
                    case=case,
                    index=index,
                    completed=True,
                    metrics=metrics,
                    checks=checks(metrics, c),
                )
            )
    passed = len(rows) == 4 and all(r["completed"] and all(r["checks"].values()) for r in rows)
    require(
        result["real_guard_verdict"] is (True if passed else None)
        and exits == ([0] * 6 if passed else [2] * 6),
        "Terminal verdict/exit mismatch",
    )
    require(
        result["pole_verdict"] is result["model_insensitivity_verdict"] is None
        and result["native_calls"]
        == result["new_complex_maps"]
        == result["route_holdout_s_matrices_analyzed"]
        == 0
        and not result["full_sweep_released"],
        "Unsupported scientific successor/verdict",
    )
    samples = list(csv.DictReader(io.StringIO(files["outputs/check/resources.csv"].decode())))
    for col, key in [
        ("container_memory_bytes", "sampled_peak_container_memory_bytes"),
        ("scratch_bytes", "sampled_peak_scratch_allocated_bytes"),
    ]:
        require(max(int(r[col]) for r in samples) == resources[key], "Wrong resource peak")
    require(
        {"outputs/check/" + n for n in allocation["files"]}
        == {n for n in hashes if n.startswith("outputs/check/")},
        "Wrong allocation witness file set",
    )
    for n, v in allocation["files"].items():
        require(
            v["sha256"] == hashes["outputs/check/" + n]
            and v["bytes"] == collection["files"]["outputs/check/" + n]["bytes"]
            and v["allocated_bytes"] >= 0
            and v["allocated_bytes"] % 512 == 0,
            "Changed allocation witness payload",
        )
    peak, final = disk_peak(max(int(r["run_disk_bytes"]) for r in samples), allocation)
    require(peak == resources["sampled_peak_run_allocated_bytes"], "Wrong combined final disk peak")
    return dict(
        real_guard_verdict=True if passed else None,
        original_exits=exits,
        rows=rows,
        source_commit=owner["source_commit"],
        source_entries_verified=len(pins),
        archived_files_verified=len(hashes),
        reference_inputs_verified=len(refs),
        resource_samples_verified=len(samples),
        resources=resources,
        final_pre_summary_allocated_bytes=final,
        finished_utc=owner["finished_utc"],
        native_calls=0,
        new_maps=0,
        pole_verdict=None,
        model_insensitivity_verdict=None,
        route_holdout_s_matrices_analyzed=0,
        full_sweep_released=False,
        independence=(
            "Full saved-byte/source/input authentication, NumPy native K and surface sum, "
            "Jost S algebra and flux/invariant arithmetic. No new ODE/native call, independent "
            "exterior solution, complex pole or basis-model verdict."
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ["archive", "collection", "allocation"]:
        p.add_argument(name, type=Path)
    args = p.parse_args()
    print(
        json.dumps(
            verify(
                args.archive,
                json.loads(args.collection.read_text()),
                json.loads(args.allocation.read_text()),
            ),
            indent=2,
            allow_nan=False,
        )
    )
