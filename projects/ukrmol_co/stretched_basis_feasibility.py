"""Reconstruct saved stretched target gates and retention inventory, never solving CI."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import struct
import tarfile
from pathlib import Path

import numpy as np

IRREPS = ("A1", "B1", "B2", "A2")
SECTORS = tuple(f"{spin}.{irrep}" for irrep in IRREPS for spin in ("singlet", "triplet"))
DZ = "runs/co-r2500-ccdz-sa12-slepc-roots8-dz80-residual-qualified-target"
TZ = "runs/co-r2500-cctz-sa12-slepc-roots8-dz80-residual-qualified-target-mpi12-package3"
DRIFT = "diagnostics/cas12-tz-seed-drift-20261010"
REPLAY = "diagnostics/cas12-tz-seed-drift-replay-20261010"
QUAL = "diagnostics/cas12-tz-reoptimized-target-qualification-20261010"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def records(data):
    """Require complete little-endian four-byte Fortran record framing."""
    cursor = 0
    while cursor < len(data):
        if cursor + 8 > len(data):
            raise ValueError("Truncated native record")
        size = struct.unpack_from("<I", data, cursor)[0]
        end = cursor + 4 + size
        if end + 4 > len(data) or struct.unpack_from("<I", data, end)[0] != size:
            raise ValueError("Invalid native record framing")
        yield memoryview(data)[cursor + 4 : end]
        cursor = end + 4


def read_cidata(data):
    """Read all eight eight-root sets, including phases and unmodified CI vectors."""
    blocks = iter(records(data))
    answer = {}
    for index, sector in enumerate(SECTORS, 1):
        if bytes(next(blocks)) != b"********                CIDATA  ":
            raise ValueError("Wrong CIDATA marker")
        header = next(blocks)
        if len(header) != 208:
            raise ValueError("Unknown CIDATA header")
        number, count = struct.unpack_from("<qq", header)
        centers, dimension, roots, symmetry = struct.unpack_from("<4q", header, 136)
        spin = index % 2 == 0
        if not (
            number == index
            and count == 11
            and centers == 2
            and roots == 8
            and dimension >= 8
            and symmetry == IRREPS.index(sector.split(".")[1])
            and bytes(header[16:136]).decode().strip()
            == "CO CI calculation - " + sector.replace(".", " ")
            and struct.unpack_from("<dd", header, 168) == (float(spin), float(spin))
            and struct.unpack_from("<q", header, 184)[0] == 14
        ):
            raise ValueError("Unexpected CIDATA set/root/electron inventory")
        nuclear = struct.unpack_from("<d", header, 192)[0]
        for charge in (6, 8):
            center = next(blocks)
            if len(center) != 40 or struct.unpack_from("<d", center, 32)[0] != charge:
                raise ValueError("Wrong physical center")
        mixed = next(blocks)
        if len(mixed) != 8 * (2 * dimension + roots):
            raise ValueError("Wrong mixed phase/energy record")
        phases = np.frombuffer(mixed, "<i8", count=dimension)
        if not np.isin(phases, [-1, 1]).all():
            raise ValueError("Invalid CSF phase")
        energies = np.frombuffer(mixed, "<f8", count=roots, offset=8 * dimension) + nuclear
        if not np.isfinite(energies).all() or np.any(np.diff(energies) < 0):
            raise ValueError("Invalid native energies")
        vectors = []
        for root in range(1, 9):
            record = next(blocks)
            if len(record) != 8 * (dimension + 1) or struct.unpack_from("<q", record)[0] != root:
                raise ValueError("Missing/reordered native vector")
            vectors.append(np.frombuffer(record, "<f8", offset=8))
        vectors = np.array(vectors)
        if not np.isfinite(vectors).all():
            raise ValueError("Nonfinite native CI vector")
        error = float(np.max(abs(vectors @ vectors.T - np.eye(8))))
        if error > 1e-9:
            raise ValueError("Native CI orthogonality failed")
        answer[sector] = dict(energies=energies, dimension=dimension, orthogonality_error=error)
    if next(blocks, None) is not None:
        raise ValueError("Trailing CIDATA set")
    return answer


def native_tokens(text):
    """Select only the newly solved set, excluding inherited CIDATA echoes."""
    matches = re.findall(
        r"CI data will be stored as set number\s+(\d+)\s*\n\s*EIGEN-ENERGIES\s*([\s\d.+-]+)", text
    )
    if len(matches) != 1:
        raise ValueError("No unique new native sector")
    number, block = matches[0]
    energies = [float(t) for t in re.findall(r"-\d+\.\d+", block)]
    if len(energies) != 8 or "Not all requested eigenpairs have converged" in text:
        raise ValueError("Incomplete native roots")
    return int(number), energies


def property_inventory(text, native):
    """Audit actual property IDs/counts against the same frame's CIDATA energies."""
    counts = dict.fromkeys(SECTORS, 0)
    ids, moments, header = [], set(), []
    for line in text.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "6":
            header.append(int(fields[4]))
        elif fields[0] == "5":
            sector = fields[-2].lower() + "." + fields[-1]
            counts[sector] += 1
            root = counts[sector] - 1
            if abs(float(fields[8].replace("D", "E")) - native[sector]["energies"][root]) > 1e-7:
                raise ValueError("Property/native energy mismatch")
            ids.append(int(fields[1]))
        elif fields[0] == "1":
            key = tuple(map(int, fields[1:8]))
            if key in moments or not np.isfinite(float(fields[8].replace("D", "E"))):
                raise ValueError("Duplicate/nonfinite property")
            moments.add(key)
    if header != [40] or ids != list(range(1, 41)) or sum(counts.values()) != 40:
        raise ValueError("Wrong existing property inventory")
    return dict(
        counts=counts,
        multipoles=len(moments),
        selected_all_eight_available=False,
        missing_states_for_all64=sum(8 - n for n in counts.values()),
    )


def reconstruct(archive, collection):
    """Public saved-record gate replay; no package import or Hamiltonian action."""
    if digest(archive) != collection["archive_sha256"]:
        raise ValueError("Archive hash differs")
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:gz") as stream:
        files = {m.name: stream.extractfile(m).read() for m in stream if m.isfile()}
    if set(files) != set(collection["files"]):
        raise ValueError("Archive inventory differs")
    for name, expected in collection["files"].items():
        if len(files[name]) != expected["bytes"] or digest(files[name]) != expected["sha256"]:
            raise ValueError("Saved file changed: " + name)

    def load(name):
        return json.loads(files[name])

    verified = set()
    for base in [
        "prepared/co-tz-reoptimized-qualification-source-20261010",
        "prepared/co-tz-seed-drift-replay-source-20261010",
    ]:
        freeze = load(base + "/frozen-inputs.json")
        for kind, mapping in freeze.items():
            for name, want in mapping.items():
                path = base + "/" + name if kind == "source" else name
                if digest(files[path]) != want:
                    raise ValueError("Original frozen prerequisite differs: " + path)
                verified.add(path)
    recovered = load(REPLAY + "/recovered-result.json")
    independent = load(REPLAY + "/independent-verification.json")
    qualification = load(QUAL + "/result.json")
    contract = load("prepared/co-tz-reoptimized-qualification-source-20261010/contract.json")
    for name, want in contract["replay_sha256"].items():
        assert digest(files[name]) == want
    assert (
        qualification["qualification_passed"]
        and not qualification["original_import_seed_preservation_verdict"]
    )
    assert qualification["original_import_controller_exit"] == 1
    assert (
        independent["probe_records_verified"] == 48
        and independent["widest_space_states_independently_reconstructed"] == 128
    )
    results = {}
    native_max = 0.0
    for label, base in [("dz", DZ), ("tz", TZ)]:
        geometry = base + "/output/CO/geom1/"
        target, diagnostics = load(base + "/target.json"), load(base + "/target-diagnostics.json")
        assert target["converged"] and diagnostics["orbital_converged"]
        assert all(diagnostics["ci_converged"]) and all(diagnostics["fresh_ci_converged"])
        assert len(target["states"]) == 40
        assert {(s["spin"], s["irrep"], s["root"]) for s in target["states"]} == {
            (spin, irrep, root)
            for spin in ("singlet", "triplet")
            for irrep in IRREPS
            for root in range(1, 6)
        }
        assert all(
            abs(s["spin_square"] - (0 if s["spin"] == "singlet" else 2)) <= 1e-6
            for s in target["states"]
        )
        assert target["mo_orthogonality_max_error"] <= 1e-9
        assert target["pi_energy_max_difference_hartree"] <= 1e-7
        assert all(abs(s["weight"] - 0.025) < 1e-15 for s in target["states"])
        assert (
            abs(
                sum(s["weight"] * s["energy_hartree"] for s in target["states"])
                - target["ensemble_energy_hartree"]
            )
            < 1e-10
        )
        last = diagnostics["iterations"][-1]
        assert last["orbital_gradient_norm"] <= 1e-7 and abs(last["energy_change_hartree"]) <= 1e-11
        assert digest(files[geometry + "co.molden"]) == target["molden_sha256"]
        native = read_cidata(files[geometry + "fort.26"])
        for index, sector in enumerate(SECTORS, 1):
            number, tokens = native_tokens(
                files[geometry + "outputs/target.scatci." + sector + ".out"].decode()
            )
            assert number == index
            np.testing.assert_allclose(tokens, native[sector]["energies"], atol=5.1e-11, rtol=0)
            if label == "tz":
                row = load(
                    DRIFT
                    + "/imported_orbitals/spin"
                    + ("0" if sector.startswith("singlet") else "2")
                    + "-"
                    + sector.split(".")[1]
                    + "-space240/result.json"
                )
                error = float(
                    np.max(
                        abs(
                            native[sector]["energies"]
                            - [s["energy_hartree"] for s in row["states"]]
                        )
                    )
                )
                assert error <= 1e-7
                native_max = max(native_max, error)
        stages = [line.split("\t") for line in files[base + "/stages.tsv"].decode().splitlines()]
        denprop = [float(row[-2]) for row in stages if row[1] == "denprop"]
        assert len(denprop) == 1 and all(row[-1] == "0" for row in stages)
        results[label] = dict(
            checkpoint_sha256=digest(files[geometry + "co.casscf.chk"]),
            native_roots=sum(len(s["energies"]) for s in native.values()),
            native_dimensions={s: row["dimension"] for s, row in native.items()},
            native_ci_orthogonality_error=max(
                row["orthogonality_error"] for row in native.values()
            ),
            properties=property_inventory(files[geometry + "fort.24"].decode(), native),
            original_denprop_seconds=denprop[0],
        )
    assert results["tz"]["checkpoint_sha256"] == contract["checkpoint_sha256"]
    assert (
        digest(
            files[
                "prepared/mpi12-import-successors-source-20261009/basis/projects/ukrmol_co/target.py"
            ]
        )
        == load(TZ + "/target.json")["source_sha256"]
    )
    probe_count = 0
    vectors_checked = 0
    max_saved_residual = 0.0
    for frame in recovered["frames"]:
        assert len(frame["probes"]) == 24 and frame["fixed_orbital_coverage_pass"]
        for row in frame["probes"]:
            path = (
                DRIFT
                + "/"
                + frame["label"]
                + f"/spin{row['spin']}-{row['irrep']}-space{row['space']}"
            )
            assert load(path + "/result.json") == row and all(row["ci_converged"])
            assert row["ci_orthogonality_error"] <= 1e-9
            for state in row["states"]:
                for key in [
                    "physical_residual_hartree",
                    "penalized_residual_hartree",
                    "spin_penalty_action_hartree",
                ]:
                    assert 0 <= state[key] <= 1e-9
                    max_saved_residual = max(max_saved_residual, state[key])
                assert abs(state["spin_square_error"]) <= 1e-6
            if row["space"] == 240:
                ci = np.load(io.BytesIO(files[path + "/ci.npy"]), allow_pickle=False).reshape(8, -1)
                np.testing.assert_allclose(ci @ ci.T, np.eye(8), atol=1e-9, rtol=0)
                vectors_checked += 8
            probe_count += 1
        assert (
            max(
                frame["maximum_space_energy_error_hartree"],
                frame["maximum_pi_energy_error_hartree"],
            )
            <= 1e-7
        )
    assert probe_count == 48 and vectors_checked == 128
    for key, value in independent["maxima"].items():
        assert value <= (1e-6 if key == "spin_square_error" else 1e-9)
    assert (
        min(v for group in independent["subspace_singular_values"].values() for v in group)
        >= 0.99999999
    )
    assert (
        independent["original_seed_gate_remains_rejected"]
        and recovered["original_controller_exit"] == 1
    )
    assert abs(independent["triplet_A1_root3_frame_difference_hartree"]) > 1e-7
    return dict(
        saved_proof_reconstruction_pass=True,
        archived_files_verified=len(files),
        original_frozen_entries_verified=len(verified),
        original_qualification_exit=load(
            "prepared/co-tz-reoptimized-qualification-source-20261010/execution.json"
        )["original_controller_exit"],
        original_import_seed_preservation_verdict=False,
        original_import_controller_exit=1,
        native_tz_max_saved_ci_energy_error_hartree=native_max,
        saved_probe_records=probe_count,
        saved_ci_vectors_norm_checked=vectors_checked,
        maximum_saved_probe_action_residual_hartree=max_saved_residual,
        models=results,
        selected_common_retention=dict(
            roots_per_sector=8,
            components=64,
            scope=(
                "Finite matched maximum-available native target model; "
                "not nine-A2 or production convergence"
            ),
        ),
        missing_property_exports=2,
        new_ci_calls=0,
        new_qc_calls=0,
        new_native_calls=0,
        new_maps=0,
        route_holdout_s_matrices_analyzed=0,
        fresh_hamiltonian_actions_recomputed=False,
        full_sweep_released=False,
        independence=(
            "Raw native energies/CIDATA vectors/properties, original frozen hashes, "
            "ordinary QC flags and saved128 CI-vector norms independently checked. "
            "Original physical/spin action metrics and subspace checks are authenticated "
            "saved independent evidence, not fresh Hamiltonian/spin contractions. "
            "No eigensolve or native engine call."
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("collection", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    report = reconstruct(args.archive.read_bytes(), json.loads(args.collection.read_text()))
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
