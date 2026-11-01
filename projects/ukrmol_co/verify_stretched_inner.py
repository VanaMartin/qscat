"""Independent saved64/400 native-export reader; no producer imports or new scientific calls."""

from __future__ import annotations

import argparse
import csv
import io
import json
import struct
from pathlib import Path

import numpy as np

from projects.ukrmol_co.verify_stretched_preparation import stream_files
from projects.ukrmol_co.verify_stretched_properties import (
    SECTORS,
    files_in,
    inventory,
    native_energies,
    require,
    sha,
)

PREFIX = "outputs/inner/native/B1/"
CONTRACT = "projects/ukrmol_co/stretched-dz-inner-contract-20261031.json"
CONTRACT_SHA = "47f150582e19e3072d2c9c1d7ddea0fe9aba9071958f6991ac715f089bad95cc"
CASES = {
    "dz": (CONTRACT, CONTRACT_SHA, 34, 7),
    "tz": (
        "projects/ukrmol_co/stretched-tz-inner-contract-20261101.json",
        "e107ee7aaf256d82c9603727be6336e402a87bfd6d52f4702452c448b4c62c96",
        36,
        8,
    ),
}


def select_contract(source):
    matched = [(case, spec) for case, spec in CASES.items() if spec[0] in source]
    require(len(matched) == 1, "Missing/ambiguous frozen inner case")
    case, spec = matched[0]
    require(sha(source[spec[0]]) == spec[1], "Changed frozen inner contract")
    return case, spec


def records(data):
    answer, cursor = [], 0
    while cursor < len(data):
        require(cursor + 8 <= len(data), "Truncated native record")
        size = struct.unpack_from("<I", data, cursor)[0]
        end = cursor + 4 + size
        require(
            end + 4 <= len(data) and data[cursor : cursor + 4] == data[end : end + 4],
            "Broken native framing",
        )
        answer.append(memoryview(data)[cursor + 4 : end])
        cursor = end + 4
    return answer


def numbers(data, dtype, count):
    require(len(data) == 8 * count, "Wrong native array size")
    result = np.frombuffer(data, dtype)
    require(np.isfinite(result).all(), "Nonfinite native array")
    return result


def symmetric(packed, dimension):
    require(packed.shape[-1] == dimension * (dimension + 1) // 2, "Wrong triangular length")
    result = np.zeros((*packed.shape[:-1], dimension, dimension))
    cursor = 0
    for i in range(dimension):
        for j in range(i + 1):
            result[..., i, j] = result[..., j, i] = packed[..., cursor]
            cursor += 1
    return result


def ci_export(data):
    rows = records(data)
    require(
        len(rows) >= 2
        and bytes(rows[0]) == b"********                CIDATA  "
        and len(rows[1]) >= 208
        and bytes(rows[1][16:136]).decode().strip() == "CO CI calculation - doublet B1",
        "Wrong scattering CIDATA",
    )
    centers, dimension, roots, _ = struct.unpack_from("<4q", rows[1], 136)
    require(
        (centers, dimension, roots) == (3, 87198, 2048)
        and struct.unpack_from("<ddq", rows[1], 168) == (0.5, 0.5, 15),
        "Wrong B1 CI header",
    )
    geometry = [struct.unpack("<4d", row[8:]) for row in rows[2:5]]
    np.testing.assert_allclose([r[3] for r in geometry], [6, 8, 0], atol=0, rtol=0)
    np.testing.assert_allclose(
        np.linalg.norm(np.asarray(geometry[0][:3]) - geometry[1][:3]), 2.5, atol=1e-12, rtol=0
    )
    np.testing.assert_allclose(geometry[2], [0, 0, 0, 0], atol=0, rtol=0)
    mixed = rows[5]
    require(len(mixed) == 8 * (2 * dimension + roots), "Wrong CI mixed record")
    phases = np.frombuffer(mixed[: 8 * dimension], "<i8")
    require(np.isin(phases, [-1, 1]).all(), "Wrong CI phase flags")
    nuclear = struct.unpack_from("<d", rows[1], 192)[0]
    energies = numbers(mixed[8 * dimension : 8 * (dimension + roots)], "<f8", roots) + nuclear
    require(len(rows) == 6 + roots and np.all(np.diff(energies) > 0), "Wrong CI root count/order")
    for index, row in enumerate(rows[6:], 1):
        require(struct.unpack_from("<q", row)[0] == index, "Wrong CI vector index")
        numbers(row[8:], "<f8", 2256)
    return energies


def native_arrays(files):
    ci = ci_export(files[PREFIX + "fort.25"])
    channel, boundary = records(files[PREFIX + "fort.10"]), records(files[PREFIX + "fort.21"])
    require(
        struct.unpack_from("<q", channel[0])[0] == 10
        and struct.unpack_from("<q", boundary[0])[0] == 11,
        "Wrong native file types",
    )
    require(
        bytes(channel[2]) == bytes(boundary[2])
        and struct.unpack("<4q", channel[2]) == (64, 0, 0, 400)
        and len(channel) == 468
        and bytes(channel[3]) == bytes(boundary[3]),
        "Wrong64/400 channel interface",
    )
    require(struct.unpack("<4q", channel[3][:32])[3] == 0, "Nonneutral interface")
    targets = np.asarray([struct.unpack_from("<d", r, 32)[0] for r in channel[4:68]])
    require(
        [struct.unpack_from("<q", r)[0] for r in channel[4:68]] == list(range(1, 65)),
        "Wrong target IDs",
    )
    labels = np.asarray([struct.unpack("<4q", r[:32]) for r in channel[68:]])
    thresholds = np.asarray([struct.unpack_from("<d", r, 32)[0] for r in channel[68:]])
    np.testing.assert_array_equal(labels[:, 0], np.arange(1, 401))
    require(
        labels[:, 1].min() == 1
        and labels[:, 1].max() <= 64
        and np.all(labels[:, 2] >= abs(labels[:, 3])),
        "Invalid channel labels",
    )
    np.testing.assert_allclose(
        thresholds, 2 * (targets[labels[:, 1] - 1] - targets[0]), atol=1e-10, rtol=0
    )
    props = [
        line.split()
        for line in files[PREFIX + "fort.24"].decode().splitlines()
        if line.split()[:1] == ["5"]
    ]
    require(len(props) == 64, "Wrong own properties")
    np.testing.assert_allclose(
        targets, [float(p[8].replace("D", "E")) for p in props], atol=1e-9, rtol=0
    )
    require(
        len(boundary) == 8 and struct.unpack("<4q", boundary[4][:32]) == (2, 2048, 0, 0),
        "Wrong boundary layout",
    )
    np.testing.assert_allclose(struct.unpack("<d", boundary[4][32:])[0], 18, atol=1e-12, rtol=0)
    poles = numbers(boundary[6], "<f8", 2048)
    np.testing.assert_allclose(poles, ci, atol=1e-10, rtol=0)
    amplitudes = numbers(boundary[7], "<f8", 2048 * 400).reshape(2048, 400)
    cf = symmetric(numbers(boundary[5], "<f8", 400 * 401).reshape(2, -1), 400)
    km = records(files[PREFIX + "fort.921"])
    require(
        len(km) == 104
        and struct.unpack_from("<2q", km[0]) == (11, 1)
        and bytes(km[2]) == bytes(channel[3])
        and struct.unpack("<5q", km[3]) == (64, 0, 0, 400, 1),
        "Wrong K metadata",
    )
    index, count, start, step = struct.unpack("<2q2d", km[4])
    require(index == 1 and count == 99, "Wrong K grid metadata")
    wave, packed = np.empty(99), np.empty((99, 21))
    for index, row in enumerate(km[5:]):
        require(struct.unpack_from("<3q", row) == (6, 0, 21), "Wrong open-channel K shape")
        values = numbers(row[24:], "<f8", 22)
        wave[index], packed[index] = values[0], values[1:]
    k = symmetric(packed, 6)
    np.testing.assert_allclose(wave, start + step * np.arange(99), atol=1e-12, rtol=0)
    np.testing.assert_allclose(wave / 0.0735, 0.1 + 0.05 * np.arange(99), atol=1e-11, rtol=0)
    np.testing.assert_allclose(thresholds[:6], 0, atol=0, rtol=0)
    require(np.all(thresholds[6:] > wave.max()), "Closed-channel threshold crosses diagnostic grid")
    phase = np.loadtxt(io.BytesIO(files[PREFIX + "fort.1021"]), ndmin=2)
    require(phase.shape == (99, 2) and np.isfinite(phase).all(), "Wrong phase export")
    np.testing.assert_allclose(phase[:, 0], 0.1 + 0.05 * np.arange(99), atol=5e-5, rtol=0)
    error = np.max(
        abs(
            (np.sum(np.arctan(np.linalg.eigvalsh(k)), axis=1) - phase[:, 1] + np.pi / 2) % np.pi
            - np.pi / 2
        )
    )
    require(error <= 1e-6, "K/native phase mismatch")
    arrays = dict(
        poles=poles,
        amplitudes=amplitudes,
        c2=cf[0] + np.diag(labels[:, 2] * (labels[:, 2] + 1)),
        c3=cf[1],
        thresholds=thresholds,
        saved=k,
    )
    return (
        arrays,
        dict(
            boundary_to_ci_energy_error_hartree=float(np.max(abs(poles - ci))),
            max_k_to_eigenphase_sum_difference_rad_mod_pi=float(error),
            first_closed_threshold_rydberg=float(thresholds[6:].min()),
            ground_hartree=float(targets[0]),
        ),
        labels,
    )


def disk_peak(periodic_peak, allocation):
    final = sum(row["allocated_bytes"] for row in allocation["files"].values())
    final -= allocation["files"]["resources.json"]["allocated_bytes"]
    return max(periodic_peak, final), final


def verify(path, collection, allocation):
    native = {
        PREFIX + n
        for n in ["fort.10", "fort.21", "fort.25", "fort.921", "fort.1021", "fort.24", "fort.26"]
    }
    files, digests, _ = stream_files(path, collection, capture_names=native)

    def load(name):
        return json.loads(files[name])

    source = files_in(files["source.tar.gz"])
    pins = load("source-manifest.json")
    case, spec = select_contract(source)
    require(set(source) == set(pins) and len(pins) == spec[2], "Wrong frozen source")
    for n, w in pins.items():
        require(sha(source[n]) == w, "Changed source: " + n)
    c = json.loads(source[spec[0]])
    references = collection["reference_files"]
    require(set(references) == set(c["input_sha256"]), "Wrong prerequisite set")
    for n, w in c["input_sha256"].items():
        require(digests[references[n]] == w, "Changed prerequisite: " + n)
    if case == "tz":
        prior = load(references[c["prior_inner_receipt"]])
        require(
            prior["case"] == "dz"
            and prior["inner_export_verdict"] is True
            and prior["public_reader_exit"] == 0
            and prior["original_exits"] == [0] * 6
            and prior["contract_sha256"] == CONTRACT_SHA
            and (prior["roots"], prior["target_count"], prior["boundary_channels"])
            == (2048, 64, 400),
            "Unqualified public DZ inner predecessor",
        )
    for n, w in c["source_prerequisite_sha256"].items():
        require(sha(source[n]) == w, "Changed source prerequisite")
    owner, result, resources = (
        load("execution.json"),
        load("outputs/inner/result.json"),
        load("outputs/inner/resources.json"),
    )
    exits = [
        result["scientific_exit"],
        owner["original_command_exit"],
        owner["original_container_exit"],
        owner["original_owner_exit"],
        load("service-return.json")["original_service_exit"],
        resources["exit_code"],
    ]
    require(
        exits == [0] * 6
        and result["inner_export_verdict"] is True
        and result["copied_native_inputs_unchanged"]
        and owner["source_integrity_after"]
        and not owner["timed_out"]
        and not owner["container_oom_killed"],
        "Failed original run",
    )
    require(
        digests["source.tar.gz"] == owner["source_archive_sha256"]
        and sha(source["projects/ukrmol_co/compressed_cas12_inner_owner.py"])
        == owner["owner_sha256"],
        "Wrong executing source",
    )
    container = load("container-inspect.json")
    require(
        container["Image"] == c["image_id"]
        and not container["State"]["Running"]
        and container["State"]["ExitCode"] == 0
        and not container["State"]["OOMKilled"],
        "Wrong container state",
    )
    require(
        owner["engine_verified_before"] == owner["engine_verified_after"]
        and owner["engine_verified_after"]["installed_files"] == 242
        and owner["engine_verified_after"]["frozen_counts"] == dict(source=55, inputs=192),
        "Changed qualified engine attestation",
    )
    require(
        result["runtime_gate"]["runtime_gate_passed"]
        and result["runtime_gate"]["numerical_library_sha256"] == c["numerical_library_sha256"]
        and result["runtime_gate"]["telemetry_sha256"] == c["telemetry_sha256"],
        "Runtime library/telemetry gate failed",
    )
    pilot = load("outputs/inner/native/pilot.json")
    row = pilot["sectors"]["B1"]
    require(
        list(pilot["sectors"]) == ["B1"]
        and pilot["all_stages_complete"]
        and pilot["pi_check_passed"] is None,
        "Unexpected counterpart",
    )
    require(
        [s["name"] for s in row["stages"]] == ["scatci", "rsolve", "eigenp", "reson"]
        and all(s["original_exit"] == 0 and not s["timed_out"] for s in row["stages"]),
        "Failed native stage",
    )
    require(
        all(
            result[k] is None
            for k in [
                "pi_pair_verdict",
                "pole_verdict",
                "electronic_model_verdict",
                "spectral_truncation_verdict",
            ]
        )
        and not result["full_sweep_released"]
        and result["route_holdout_s_matrices_analyzed"] == 0,
        "Unsupported scientific verdict",
    )
    for n, w in row["input_sha256"].items():
        original = c["preparation_directory"] + "/" + ("moints" if n == "fort.16" else n)
        require(digests[PREFIX + n] == w == c["input_sha256"][original], "Changed native copy")
    require(
        digests[PREFIX + "scatci.inp"] == sha(source[c["frozen_native_deck"]]),
        "Changed native deck",
    )
    own_target = json.loads(files[references[c["preparation_directory"] + "/target.json"]])
    energies, _ = native_energies(files[PREFIX + "fort.26"])
    own = inventory(files[PREFIX + "fort.24"].decode(), energies, own_target["ground_dipole_au"][2])
    props = [
        line.split()
        for line in files[PREFIX + "fort.24"].decode().splitlines()
        if line.split()[:1] == ["5"]
    ]
    mapping = [
        int(p[1]) for sector in SECTORS for p in props if p[-2].lower() + "." + p[-1] == sector
    ]
    require(mapping == pilot["target_map"] == c["own_target_map"], "Wrong contraction/property map")
    telemetry = [json.loads(line) for line in files[PREFIX + "eps-telemetry.jsonl"].splitlines()]
    require(len(telemetry) == 1, "Ambiguous telemetry")
    t = telemetry[0]
    residual = np.asarray([v["absolute_residual"] for v in t["eigenpairs"]])
    require(
        t["requested"] == 2048
        and t["converged"] >= 2048
        and t["reason"] > 0
        and t["dimension"] == 87198
        and t["matrix_type"] == "mpisbaij"
        and t["solver_type"] == "krylovschur"
        and len(residual) >= 2048
        and np.isfinite(residual).all()
        and residual.min() >= 0
        and residual.max() <= 1e-7,
        "Native residual/solver gate failed",
    )
    np.testing.assert_allclose(
        residual.max(), row["maximum_native_residual_hartree"], atol=1e-20, rtol=1e-12
    )
    arrays, metrics, labels = native_arrays(files)
    export = result["export_gate"]
    require(
        export["strict64_target_400_channel_assembly_verdict"]
        and export["open_labels"] == labels[:6].tolist(),
        "Reported export differs",
    )
    for n, array in arrays.items():
        require(
            list(array.shape) == export["arrays"][n]["shape"] and np.isfinite(array).all(),
            "Array export differs",
        )
    samples = list(csv.DictReader(io.StringIO(files["outputs/inner/resources.csv"].decode())))
    allocation_prefix = "outputs/inner/"
    regular = {allocation_prefix + n for n, v in allocation["files"].items() if "sha256" in v}
    require(
        regular == {n for n in digests if n.startswith(allocation_prefix)},
        "Wrong final allocation file set",
    )
    links = {
        allocation_prefix + n: v["symlink_target"]
        for n, v in allocation["files"].items()
        if "symlink_target" in v
    }
    require(
        links
        == {n: v for n, v in collection["symlinks"].items() if n.startswith(allocation_prefix)},
        "Wrong final allocation symlinks",
    )
    for n, v in allocation["files"].items():
        require(
            v["allocated_bytes"] >= 0 and v["allocated_bytes"] % 512 == 0,
            "Invalid filesystem block count",
        )
        if "sha256" in v:
            require(
                v["sha256"] == digests[allocation_prefix + n]
                and v["bytes"] == collection["files"][allocation_prefix + n]["bytes"],
                "Changed final allocation payload",
            )
    periodic_disk = max(int(r["run_disk_bytes"]) for r in samples)
    combined_disk, final_disk = disk_peak(periodic_disk, allocation)
    require(
        combined_disk == resources["sampled_peak_run_allocated_bytes"],
        "Periodic plus post-CSV final disk peak differs",
    )
    for column, key in [
        ("container_memory_bytes", "sampled_peak_container_memory_bytes"),
        ("scratch_bytes", "sampled_peak_scratch_allocated_bytes"),
    ]:
        require(
            samples and max(int(r[column]) for r in samples) == resources[key],
            "Resource peak mismatch",
        )
    require(
        resources["kernel_peak_container_memory_bytes"]
        >= resources["sampled_peak_container_memory_bytes"]
        and result["inputs_verified_before"] == result["inputs_verified_after"] == len(references)
        and load("outputs/controls.json")["original_exit"] == 0
        and f"Ran {spec[3]} tests".encode() in files["outputs/controls.log"],
        "Control/resource/integrity gate failed",
    )
    return dict(
        case=case,
        contract_sha256=spec[1],
        inner_export_verdict=True,
        original_exits=exits,
        roots=2048,
        contracted_dimension=87198,
        continuum_coefficients=2256,
        boundary_channels=400,
        target_count=64,
        source_commit=owner["source_commit"],
        source_entries_verified=len(pins),
        archived_files_verified=len(digests),
        reference_inputs_verified=len(references),
        resource_samples_verified=len(samples),
        periodic_csv_disk_peak_bytes=periodic_disk,
        retained_final_pre_summary_allocated_bytes=final_disk,
        allocation_witness_scope=(
            "Separate authenticated retained host st_blocks metadata; original profiler takes a "
            "final measurement after CSV closes and before resources.json is written. "
            "Peak is max(periodic CSV, final). Symlink storage is counted by the original "
            "lstat implementation."
        ),
        maximum_native_residual_hartree=float(residual.max()),
        native_k_energy_records=99,
        **metrics,
        array_sha256={n: sha(np.asarray(a, dtype="<f8").tobytes()) for n, a in arrays.items()},
        shapes={n: list(a.shape) for n, a in arrays.items()},
        own_properties=own,
        resources=resources,
        native_stages=row["stages"],
        finished_utc=owner["finished_utc"],
        owner_wall_seconds=owner["wall_seconds"],
        native_floating_point_note_count=files[PREFIX + "scatci.out"].count(b"IEEE_UNDERFLOW_FLAG"),
        new_native_calls=0,
        new_s_constructions=0,
        new_complex_python_maps=0,
        route_holdout_s_matrices_analyzed=0,
        pole_verdict=None,
        electronic_model_verdict=None,
        spectral_truncation_verdict=None,
        full_sweep_released=False,
        **(
            {
                "independent_dz_predecessor_authenticated": True,
                "dz_inner_receipt_sha256": digests[references[c["prior_inner_receipt"]]],
            }
            if case == "tz"
            else {}
        ),
        independence=(
            "Streaming full-byte authentication; independent CIDATA/channel/boundary/K decoding, "
            "six-array assembly, recorded residual and own-K/native phase arithmetic. "
            "No producer imports, native calls, S construction, Hamiltonian action, CIRMAT-phase "
            "audit, target identity, electronic-accuracy or omitted-spectrum claim."
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("collection", type=Path)
    parser.add_argument("allocation", type=Path)
    args = parser.parse_args()
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
