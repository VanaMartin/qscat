"""Saved-only native-sum, root, contour and character reconstruction for own DZ."""

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
from qscat.units import hartree_to_ev

from projects.ukrmol_co.verify_stretched_dz_pole_stop import returned_dict_keys
from projects.ukrmol_co.verify_stretched_inner import numbers, records
from projects.ukrmol_co.verify_stretched_preparation import stream_files
from projects.ukrmol_co.verify_stretched_properties import files_in, require, sha
from projects.ukrmol_co.verify_stretched_real_completion import native_at
from projects.ukrmol_co.verify_stretched_real_guard import checks, measurements

SOURCE_COMMIT = "fbd16ee059f0c5fa682444a0c3af2dba4fe090fc"
SOURCE_ARCHIVE_SHA = "701215d9927a67b5c3b39dc2b06f9cf89c46a31b58de9676beb6e7f3674c132a"
RELEASE = "projects/ukrmol_co/stretched-dz-pole-adapter-release-20261102.json"
PROPOSAL = "projects/ukrmol_co/stretched-dz-pole-proposal-20261102.json"


def difference(left, right):
    return dict(
        position_ev=float(hartree_to_ev(abs((left - right).real))),
        full_width_ev=float(hartree_to_ev(2 * abs((left - right).imag))),
    )


def relative(left, right):
    norm = np.linalg.norm(right)
    return float(np.linalg.norm(left - right) / norm) if norm else 1e300


def schur(poles, w, total, q, y, a, x):
    p = 1 / (poles - total)
    b = w.T @ (p[:, None] * x)
    solved = np.linalg.solve(np.eye(len(q)) - q @ y, np.column_stack([b, q @ a.T]))
    return -(x.T @ (p[:, None] * x) + b.T @ y @ solved[:, : x.shape[1]]), a @ solved[
        :, x.shape[1] :
    ]


def moments(center, points, interior, outgoing, phases):
    result = []
    for stride in [2, 1]:
        z = points[::stride]
        delta = z - center
        c0 = np.mean(delta[:, None, None] * interior[::stride], axis=0)
        c1 = np.mean((delta * z)[:, None, None] * interior[::stride], axis=0)
        residue = 2j * np.mean(delta[:, None, None] * outgoing[::stride], axis=0)
        pole = np.vdot(c0, c1) / np.vdot(c0, c0) if np.linalg.norm(c0) > 1e-20 else None
        phase = phases[::stride]
        jumps = np.angle(np.exp(1j * (np.roll(phase, -1) - phase)))
        result.append(
            dict(
                c0=c0,
                c1=c1,
                residue=residue,
                pole=pole,
                winding=float(jumps.sum() / (2 * np.pi)),
                max_phase_jump=float(abs(jumps).max()),
            )
        )
    return result


def character(residue, labels, width, policy):
    scale = float(abs(residue).max())
    r = residue / scale if scale else residue
    squared = abs(r) ** 2
    total = squared.sum()
    singular = np.linalg.svd(r, compute_uv=False)
    norm = float(np.sqrt(total) * scale)
    rank = float(singular[1] / singular[0]) if scale else None
    transpose = float(np.linalg.norm(r - r.T) / np.sqrt(total)) if scale else None
    complement = np.asarray([abs(row[3]) != 1 for row in labels])
    leakage = (
        float(max(squared[complement, :].sum(), squared[:, complement].sum()) / total)
        if scale
        else None
    )
    guards = dict(
        signal=norm / width > policy["signal_over_full_width_min"],
        rank=rank is not None and rank <= policy["rank_ratio_max"],
        transpose=transpose is not None and transpose <= policy["transpose_relative_max"],
    )
    disposition = (
        (
            "pass"
            if leakage <= policy["leakage_pass"]
            else "reject"
            if leakage >= policy["leakage_reject"]
            else "inconclusive"
        )
        if all(guards.values())
        else "inconclusive"
    )
    profile = (squared.sum(axis=0) + squared.sum(axis=1)) / (2 * total) if scale else np.zeros(6)
    order = sorted(range(6), key=lambda i: (labels[i][1], labels[i][2], abs(labels[i][3])))
    return dict(
        disposition=disposition,
        guards=guards,
        residue_norm=norm,
        rank_ratio=rank,
        transpose_relative_error=transpose,
        leakage=leakage,
        profile=profile[order].tolist(),
    )


def verify(archive, collection):
    capture = [
        n
        for n in collection["files"]
        if n.endswith(".npz") or Path(n).name.endswith(("fort.10", "fort.21", "fort.921"))
    ]
    files, hashes, _ = stream_files(archive, collection, capture_names=capture)

    def load(n):
        return json.loads(files[n])

    source = files_in(files["source.tar.gz"])
    pins = load("source-manifest.json")
    require(
        set(source) == set(pins)
        and len(pins) == 48
        and hashes["source.tar.gz"] == SOURCE_ARCHIVE_SHA,
        "Frozen source inventory",
    )
    for n, pin in pins.items():
        require(
            sha(source[n]) == pin
            and source[n] == subprocess.check_output(["git", "show", SOURCE_COMMIT + ":" + n]),
            "Executed Git source " + n,
        )
    release = json.loads(source[RELEASE])
    base_raw = source[release["base_execution_contract"]]
    require(sha(base_raw) == release["base_execution_sha256"], "Original execution contract")
    base = json.loads(base_raw)
    execution = base | release
    execution["input_sha256"] = base["input_sha256"] | release["input_sha256"]
    require(sha(source[PROPOSAL]) == execution["proposal_sha256"], "Unchanged numerical proposal")
    c = json.loads(source[PROPOSAL]) | execution
    refs = collection["reference_files"]
    require(
        c["execution_authorized"] and set(refs) == set(c["input_sha256"]) and len(refs) == 8,
        "Authorized raw inventory",
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
        exits == [exits[0]] * 6
        and exits[0] in [0, 1, 2]
        and not owner["timed_out"]
        and not owner["container_oom_killed"],
        "Original scientific/wrapper exits",
    )
    require(
        load("outputs/controls.json")["original_exit"] == 0
        and b"Ran 28 tests" in files["outputs/controls.log"],
        "28 original controls",
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
        "Eight-core runtime caps/mounts",
    )
    for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
        require(name + "=1" in container["Config"]["Env"], "Library threads")
    check = load("adapter-check.json")
    require(
        check["saved_native_adapter_verdict"]
        and check["new_ode_maps"] == check["native_calls"] == 0
        and check["adapter_source_sha256"]
        == sha(source["projects/ukrmol_co/stretched_dz_pole.py"]),
        "Positive saved adapter check",
    )
    require(
        sha(files["adapter-check.json"])
        == load("freeze.json")["positive_saved_adapter_check_sha256"]
        and check["owner_source_sha256"]
        == sha(source["projects/ukrmol_co/stretched_dz_pole_owner.py"])
        and check["model_keys"]
        == returned_dict_keys(
            source["projects/ukrmol_co/native_full_reconstruction.py"], "_model_from_interface"
        ),
        "Actual saved-load source and return-schema provenance",
    )
    old = load(refs[c["failed_stop_reading"]])
    previous = load(refs[c["completion_reading"]])
    require(
        old["reader_exit"] == 0
        and old["result"]["original_exits"] == [2] * 6
        and old["result"]["real_maps"] == old["result"]["complex_maps"] == 0
        and previous["reader_exit"] == 0
        and previous["result"]["combined_20degree_real_prerequisite_verdict"],
        "Preserved prior stop/real pass",
    )
    require(
        not result["tz_started"]
        and result["native_calls"] == result["route_holdout_s_matrices_analyzed"] == 0
        and not result["automatic_successor"]
        and not result["full_sweep_released"],
        "Protected scope",
    )
    real = []
    for row, setting in zip(result["real_rows"], c["near_feature_real_settings"], strict=False):
        require(row["index"] == 16 and row["setting"] == setting, "Near-feature real identity")
        if not row["completed"]:
            continue
        energy, q_ref, native = native_at(files, refs, c["model"], 16)
        with np.load(
            io.BytesIO(files[f"outputs/check/real-{setting[1]}.npz"]), allow_pickle=False
        ) as z:
            q, y, a, s, reference = [z[n] for n in ["q", "y", "a", "s", "reference"]]
            np.testing.assert_allclose(float(z["energy_rydberg"]), energy, atol=1e-14, rtol=0)
        np.testing.assert_allclose(q, q_ref, atol=1e-9, rtol=1e-10)
        np.testing.assert_allclose(reference, native, atol=1e-12, rtol=1e-12)
        actual, metrics = measurements(q, y, a, native)
        np.testing.assert_allclose(actual, s, atol=1e-10, rtol=1e-10)
        for n, v in metrics.items():
            np.testing.assert_allclose(
                v, row.get(n, row["invariants"].get(n)), atol=1e-13, rtol=1e-5
            )
        require(all(checks(metrics, c).values()) == all(row["checks"].values()), "Real gates")
        real.append(dict(setting=setting, metrics=metrics, checks=checks(metrics, c)))
    evaluations = result["evaluations"]
    require(
        len(evaluations) + len(result["real_rows"]) == result["new_maps_reserved"] <= 79,
        "Total map reservations",
    )
    require(
        not evaluations or (len(real) == 2 and all(all(r["checks"].values()) for r in real)),
        "Real gates before complex maps",
    )
    prefix = c["model"]["native_directory"] + "/"
    boundary = records(files[refs[prefix + "fort.21"]])
    channel = records(files[refs[prefix + "fort.10"]])
    poles = numbers(boundary[6], "<f8", 2048)
    w = numbers(boundary[7], "<f8", 2048 * 400).reshape(2048, 400)
    ground = struct.unpack_from("<d", channel[4], 32)[0]
    captures = {}
    by_point = {}
    probes = None
    if "outputs/check/probes.npz" in files:
        with np.load(io.BytesIO(files["outputs/check/probes.npz"]), allow_pickle=False) as z:
            probes = z["x"]
        expected, _ = np.linalg.qr(np.random.default_rng(314159).normal(size=(2048, 8)))
        np.testing.assert_allclose(probes, expected, atol=1e-14, rtol=0)
    for row in evaluations:
        if not row["completed"]:
            continue
        with np.load(io.BytesIO(files["outputs/check/" + row["file"]]), allow_pickle=False) as z:
            data = {n: z[n] for n in z.files}
        z = complex(*row["z_hartree"])
        q, y, a = [data[n] for n in ["q", "y", "a"]]
        require(
            q.shape == y.shape == (400, 400)
            and a.shape == ((6, 400) if row["amplitudes"] else (0, 400)),
            "Complex capture shape",
        )
        left, right, bottom, top = c["search_box_incident_hartree"]
        require(left <= z.real <= right and bottom <= z.imag <= top, "Box")
        np.testing.assert_allclose(data["z_hartree"], z, atol=1e-15, rtol=0)
        np.testing.assert_allclose(data["total_hartree"], ground + z, atol=1e-13, rtol=0)
        np.testing.assert_allclose(
            q, w.T @ (w / (poles - (ground + z))[:, None]), atol=1e-9, rtol=1e-10
        )
        d = row["diagnostics"]
        require(d["rhs_evaluations"] <= 50000 and d["wall_seconds"] < 90, "Per-map caps")
        require(
            [d["endpoint"], d["angle"], d["order"], d["rtol"]] == row["setting"], "Map settings"
        )
        if row["kind"] == "contour":
            interior, outgoing = schur(poles, w, ground + z, q, y, a, probes)
            np.testing.assert_allclose(interior, data["interior"], atol=1e-8, rtol=1e-8)
            np.testing.assert_allclose(outgoing, data["outgoing"], atol=1e-8, rtol=1e-8)
            sign, _ = np.linalg.slogdet(np.eye(400) - q @ y)
            require(
                abs(np.angle(np.exp(1j * (np.angle(sign) - row["det_phase"])))) < 1e-8,
                "Contour determinant phase",
            )
        captures[row["index"]] = data
        by_point[(z, tuple(row["setting"]), row["amplitudes"])] = data
    roots = []
    for stage, row in enumerate(result["roots"]):
        require(
            row["stage"] == stage and row["setting"] == c["settings"][stage], "Root stage routing"
        )
        reserved = [
            r for r in evaluations if r["kind"] == "root" and r["setting"] == row["setting"]
        ]
        require(len(reserved) <= c["root_evaluation_caps"][stage], "Stage root evaluation cap")
        if not row["completed"]:
            continue
        require(row["search"]["evaluations"] == len(reserved), "Saved root evaluation count")
        require(row["search"]["trace"][-1]["alpha_abs"] <= 1e-10, "Original bordered-root gate")
        z = complex(*row["z_hartree"])
        data = by_point[(z, tuple(row["setting"]), False)]
        singular = np.linalg.svd(np.eye(400) - data["q"] @ data["y"], compute_uv=False)
        require(singular[-1] <= 1e-8 and singular[-2] >= 1e-4, "Saved root null/isolation gates")
        roots.append(
            dict(
                z_hartree=row["z_hartree"],
                position_ev=float(hartree_to_ev(z.real)),
                full_width_ev=float(hartree_to_ev(-2 * z.imag)),
                sigma_min=float(singular[-1]),
                sigma_second=float(singular[-2]),
            )
        )
    root_errors = []
    if len(roots) == 3:
        center = complex(*roots[-1]["z_hartree"])
        root_errors = [difference(complex(*r["z_hartree"]), center) for r in roots]
        for actual, saved in zip(root_errors, result["root_setting_differences"], strict=True):
            for n, v in actual.items():
                np.testing.assert_allclose(v, saved[n], atol=1e-12, rtol=0)
    numerical = None
    contour_report = None
    derivative_errors = []
    stability = None
    screen_rows = None
    if result["finite_reference_pole_verdict"] is True:
        require(
            len(roots) == 3 and max(max(e.values()) for e in root_errors) <= 1e-4,
            "Root setting agreement",
        )
        center = complex(*result["pole_incident_hartree"])
        rows = [r for r in evaluations if r["kind"] == "contour"]
        require(len(rows) == 32 and all(r["completed"] for r in rows), "Complete contour")
        points = np.asarray([complex(*r["z_hartree"]) for r in rows])
        np.testing.assert_allclose(
            points, center + 0.001 * np.exp(2j * np.pi * np.arange(32) / 32), atol=1e-15, rtol=0
        )
        mm = moments(
            center,
            points,
            np.asarray([captures[r["index"]]["interior"] for r in rows]),
            np.asarray([captures[r["index"]]["outgoing"] for r in rows]),
            np.asarray([r["det_phase"] for r in rows]),
        )
        for ordinal, m in enumerate(mm):
            with np.load(
                io.BytesIO(files[f"outputs/check/moments-{ordinal}.npz"]), allow_pickle=False
            ) as z:
                for n in ["c0", "c1", "residue"]:
                    np.testing.assert_allclose(m[n], z[n], atol=1e-9, rtol=1e-8)
        high, low = mm[1], mm[0]
        sr = np.linalg.svd(high["c0"], compute_uv=False)
        ss = np.linalg.svd(high["residue"], compute_uv=False)
        require(
            np.linalg.norm(high["c0"]) > 1e-8
            and sr[1] / sr[0] <= 1e-3
            and ss[1] / ss[0] <= 1e-3
            and np.linalg.norm(high["residue"]) / (-2 * center.imag) > 1e-6,
            "Contour signal/rank",
        )
        require(
            relative(low["c0"], high["c0"]) <= 1e-5
            and relative(low["residue"], high["residue"]) <= 1e-5
            and max(difference(high["pole"], center).values()) <= 1e-4
            and max(difference(low["pole"], high["pole"]).values()) <= 1e-4
            and all(abs(m["winding"] - 1) <= 1e-6 and m["max_phase_jump"] <= np.pi / 2 for m in mm),
            "Contour moments/winding gates",
        )
        contour_report = dict(
            winding=[m["winding"] for m in mm],
            max_phase_jump=[m["max_phase_jump"] for m in mm],
            interior_rank_ratio=float(sr[1] / sr[0]),
            open_rank_ratio=float(ss[1] / ss[0]),
            interior_rule_change=relative(low["c0"], high["c0"]),
            open_rule_change=relative(low["residue"], high["residue"]),
            root_difference=difference(high["pole"], center),
            rule_pole_difference=difference(low["pole"], high["pole"]),
        )
        cc = by_point[(center, tuple(c["settings"][-1]), True)]
        matrix = np.eye(400) - cc["q"] @ cc["y"]
        u, s, vh = np.linalg.svd(matrix)
        v = vh[-1].conj()
        ell = u[:, -1].conj()
        residues = []
        for ordinal, h in enumerate(c["derivative_steps_hartree"]):
            plus = by_point[(center + h, tuple(c["settings"][-1]), False)]
            minus = by_point[(center - h, tuple(c["settings"][-1]), False)]
            derivative = (
                (np.eye(400) - plus["q"] @ plus["y"]) - (np.eye(400) - minus["q"] @ minus["y"])
            ) / (2 * h)
            boundary_residue = np.outer(v, ell @ cc["q"]) / (ell @ derivative @ v)
            residue = 2j * cc["a"] @ boundary_residue @ cc["a"].T
            residues.append(residue)
            with np.load(
                io.BytesIO(files[f"outputs/check/derivative-{ordinal}.npz"]), allow_pickle=False
            ) as z:
                np.testing.assert_allclose(derivative, z["derivative"], atol=1e-8, rtol=1e-8)
                np.testing.assert_allclose(residue, z["residue"], atol=1e-8, rtol=1e-8)
            require(
                s[-1] <= 1e-8
                and relative(residue, high["residue"]) <= 1e-3
                and relative(boundary_residue.T, boundary_residue) <= 1e-3,
                "Derivative residue gates",
            )
            derivative_errors.append(relative(residue, high["residue"]))
        screen_rows = [
            character(r, c["open_labels"], -2 * center.imag, c["character"])
            for r in [low["residue"], high["residue"], *residues]
        ]
        for actual, saved in zip(screen_rows, result["screen"]["witnesses"], strict=True):
            require(
                actual["disposition"] == saved["disposition"]
                and actual["guards"] == saved["guards"],
                "Character gates",
            )
            for n in ["leakage", "rank_ratio", "residue_norm", "transpose_relative_error"]:
                np.testing.assert_allclose(actual[n], saved[n], atol=1e-9, rtol=1e-6)
            np.testing.assert_allclose(actual["profile"], saved["profile"], atol=1e-9, rtol=0)
        stability = None
        disposition = "inconclusive"
        if all(all(r["guards"].values()) for r in screen_rows):
            reference = screen_rows[1]
            stability = dict(
                leakage_change=max(abs(r["leakage"] - reference["leakage"]) for r in screen_rows),
                profile_tv_change=max(
                    float(np.sum(abs(np.asarray(r["profile"]) - reference["profile"])) / 2)
                    for r in screen_rows
                ),
                norm_relative_change=max(
                    abs(r["residue_norm"] - reference["residue_norm"])
                    / max(r["residue_norm"], reference["residue_norm"])
                    for r in screen_rows
                ),
            )
            p = c["character"]
            if (
                len({r["disposition"] for r in screen_rows}) == 1
                and stability["leakage_change"] <= p["witness_leakage_change_max"]
                and stability["profile_tv_change"] <= p["witness_profile_tv_max"]
                and stability["norm_relative_change"] <= p["witness_norm_relative_max"]
            ):
                disposition = reference["disposition"]
            for name, value in stability.items():
                np.testing.assert_allclose(
                    value, result["screen"]["stability"][name], atol=1e-9, rtol=1e-5
                )
        require(
            disposition == result["disposition"], "Four-witness character stability/disposition"
        )
        numerical = True
    require(
        (exits[0] != 0)
        or (
            numerical is True
            and result["physical_pi_support_verdict"] is True
            and all(r["disposition"] == "pass" for r in screen_rows)
        ),
        "No successful verdict without all witnesses",
    )
    samples = list(csv.DictReader(io.StringIO(files["outputs/check/resources.csv"].decode())))
    require(
        max(int(r["container_memory_bytes"]) for r in samples)
        == resources["sampled_peak_container_memory_bytes"],
        "Resource sample peak",
    )
    require(
        result["wall_seconds"] < 7200
        and resources["wall_seconds"] < 7200
        and owner["wall_seconds"] < 7800,
        "Wall caps",
    )
    monitor = [json.loads(line) for line in files["monitor.jsonl"].decode().splitlines()]
    times = [datetime.datetime.fromisoformat(owner["started_utc"])] + [
        datetime.datetime.fromisoformat(r["observed_utc"]) for r in monitor
    ]
    gap = max((b - a).total_seconds() for a, b in itertools.pairwise(times))
    require(
        gap <= 60 and monitor[-1]["owner_process_exit"] == exits[0],
        "Observation cadence/terminal witness",
    )
    return dict(
        original_exits=exits,
        scientific_exit=result["scientific_exit"],
        disposition=result["disposition"],
        finite_reference_pole_verdict=numerical,
        physical_pi_support_verdict=result["physical_pi_support_verdict"],
        real_rows=real,
        complex_maps_reserved=len(evaluations),
        complex_maps_completed=len(captures),
        roots=roots,
        root_setting_differences=root_errors,
        contour_witness=contour_report,
        derivative_relative_errors=derivative_errors,
        character_witnesses=screen_rows,
        character_stability=stability,
        error=result.get("error"),
        source_commit=SOURCE_COMMIT,
        raw_files_verified=len(hashes),
        executed_git_sources_verified=len(pins),
        raw_input_pins_verified=8,
        source_prerequisites_verified=9,
        resources=resources,
        resource_samples_verified=len(samples),
        observation_count=len(monitor),
        maximum_observation_gap_seconds=gap,
        finished_utc=owner["finished_utc"],
        native_calls=0,
        new_ode_maps=0,
        automatic_successor=False,
        full_sweep_released=False,
        scope=(
            "Finite own-DZ pole-or-stop evidence using saved native sums/Jost matrices only;"
            "no independent exterior engine or electronic-model/basis/geometry verdict"
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
