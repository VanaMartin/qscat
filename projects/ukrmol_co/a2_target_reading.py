"""Reconstruct completed fixed-orbital48-state target exports from checksum-fetched bytes."""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import FortranFile


def read_target(directory: Path) -> dict:
    """Read full CIDATA vectors and DENPROP independently of the producer inventory."""
    manifest = json.loads((directory / "manifest.json").read_text())
    for name, row in manifest.items():
        path = directory / name
        if (
            path.stat().st_size != row["bytes"]
            or hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]
        ):
            raise ValueError(f"Fetched input checksum failed: {name}")
    contract = json.loads((directory / "contract.json").read_text())
    pilot = json.loads((directory / "pilot.json").read_text())
    resources = json.loads((directory / "resources.json").read_text())
    owner = json.loads((directory / "execution.json").read_text())
    target = json.loads((directory / "target.json").read_text())
    assert owner["original_controller_exit"] == resources["exit_code"] == 0
    assert all(step["original_exit"] == 0 for step in owner["steps"])
    sectors = contract["native_set_order"]
    assert [s["name"] for s in pilot["stages"]] == [
        "target.scatci." + sector for sector in sectors
    ] + ["target.denprop"]
    assert all(s["original_exit"] == 0 and not s["timed_out"] for s in pilot["stages"])
    assert "has completed execution" in (directory / "target.denprop.out").read_text()
    irreps = ["A1", "B1", "B2", "A2"]
    energies, binary_rows, residuals = {}, [], []
    with FortranFile(directory / "fort.26", "r") as stream:
        for sector in sectors:
            assert stream.read_record("u1").tobytes() == b"********                CIDATA  "
            header = stream.read_record("u1").tobytes()
            assert header[16:136].decode().strip() == "CO CI calculation - " + sector.replace(
                ".", " "
            )
            centers, dimension, roots, irrep = np.frombuffer(header, "<i8", count=4, offset=136)
            spin, spin_z = np.frombuffer(header, "<f8", count=2, offset=168)
            electron_count = np.frombuffer(header, "<i8", count=1, offset=184)[0]
            assert centers == 3 and dimension > 0 and electron_count == 14
            assert roots == (9 if sector.endswith(".A2") else 5)
            assert irrep == irreps.index(sector.split(".")[1])
            assert spin == spin_z == (0 if sector.startswith("singlet") else 1)
            core = np.frombuffer(header, "<f8", count=1, offset=192)[0]
            charges = []
            for _ in range(centers):
                center = stream.read_record("u1").tobytes()
                assert len(center) == 40
                coordinates_and_charge = np.frombuffer(center, "<f8", offset=8)
                assert np.isfinite(coordinates_and_charge).all()
                charges.append(coordinates_and_charge[-1])
            # CIDATA includes the zero-charge scattering center as well as C/O.
            np.testing.assert_array_equal(charges, [6, 8, 0])
            record = stream.read_record("<f8")
            assert record.size == 2 * dimension + roots
            values = record[dimension : dimension + roots] + core
            assert np.isfinite(values).all() and np.all(np.diff(values) >= 0)
            independent = json.loads((directory / "independent" / (sector + ".json")).read_text())
            assert independent["gate_passed"] and all(independent["ci_converged"])
            for key in ("physical_residual_norms_hartree", "penalized_residual_norms_hartree"):
                assert np.isfinite(independent[key][:roots]).all()
                assert np.max(independent[key][:roots]) <= 1e-9
            np.testing.assert_allclose(
                independent["spin_square"][:roots], 0 if spin == 0 else 2, atol=1e-6, rtol=0
            )
            error = float(np.max(abs(values - np.asarray(independent["energy_hartree"][:roots]))))
            assert error <= 1e-7
            original = [
                s["energy_hartree"]
                for s in target["states"]
                if s["spin"] + "." + s["irrep"] == sector
            ]
            np.testing.assert_allclose(values[:5], original, atol=1e-7, rtol=0)
            vectors = []
            for index in range(roots):
                record = stream.read_record("u1").tobytes()
                assert np.frombuffer(record, "<i8", count=1)[0] == index + 1
                vector = np.frombuffer(record, "<f8", offset=8)
                assert vector.size == dimension and np.isfinite(vector).all()
                vectors.append(vector)
            gram = np.asarray(vectors) @ np.asarray(vectors).T
            orthogonality = float(np.max(abs(gram - np.eye(roots))))
            assert orthogonality <= 1e-9
            energies[sector] = values
            binary_rows.append(
                {
                    "sector": sector,
                    "dimension": int(dimension),
                    "exported_roots": int(roots),
                    "maximum_independent_energy_difference_hartree": error,
                    "vector_gram_maximum_error": orthogonality,
                }
            )
            telemetry = [
                json.loads(line)
                for line in (directory / ("target.scatci." + sector + ".telemetry.jsonl"))
                .read_text()
                .splitlines()
            ]
            assert len(telemetry) == 1
            t = telemetry[0]
            assert t["dimension"] == dimension and t["requested"] == roots
            assert t["converged"] >= roots and t["reason"] > 0
            assert t["which"] == 1 and t["solver_type"] == "krylovschur"
            assert len(t["eigenpairs"]) >= roots
            native_residuals = [row["absolute_residual"] for row in t["eigenpairs"]]
            assert np.isfinite(native_residuals).all() and max(native_residuals) <= 1e-9
            residuals.extend(native_residuals)
            np.testing.assert_allclose(
                [row["value"] + core for row in t["eigenpairs"][:roots]], values, atol=1e-10, rtol=0
            )
        assert stream._fp.tell() == (directory / "fort.26").stat().st_size
    states, moments, counts = [], {}, []
    for line in (directory / "fort.24").read_text().splitlines():
        t = line.split()
        if not t:
            continue
        if t[0] == "6":
            counts.append(int(t[4]))
        if t[0] == "5":
            states.append((int(t[1]), t[-2].lower() + "." + t[-1], float(t[8].replace("D", "E"))))
        if t[0] == "1":
            key, value = tuple(map(int, t[1:8])), float(t[8].replace("D", "E"))
            assert key not in moments and np.isfinite(value)
            moments[key] = value
    assert counts == [48] and [s[0] for s in states] == list(range(1, 49))
    assert np.isfinite([s[2] for s in states]).all()
    assert np.all(np.diff([s[2] for s in states]) >= -1e-10)
    mapping, labels, property_errors = [], {}, []
    for sector in sectors:
        group = [s for s in states if s[1] == sector]
        assert len(group) == len(energies[sector])
        property_errors.extend(abs(np.asarray([s[2] for s in group]) - energies[sector]))
        for root, (state_id, _, _) in enumerate(group, 1):
            mapping.append(state_id)
            labels[str(state_id)] = sector + "." + str(root)
    assert len(mapping) == len(set(mapping)) == 48 and max(property_errors) <= 1e-7
    referenced, ranks = set(), set()
    for key in moments:
        left, left_irrep, right, right_irrep, kind, rank, projection = key
        assert str(left) in labels and str(right) in labels and kind == 3
        assert rank in (1, 2) and abs(projection) <= rank
        assert left_irrep == irreps.index(labels[str(left)].split(".")[1])
        assert right_irrep == irreps.index(labels[str(right)].split(".")[1])
        referenced.update([left, right])
        ranks.add(rank)
    assert referenced == set(mapping) and ranks == {1, 2}
    assert mapping == pilot["property_inventory"]["target_map"]
    assert labels == pilot["property_inventory"]["state_labels_by_id"]
    assert len(moments) == pilot["property_inventory"]["multipole_record_count"]
    ground = next(int(i) for i, label in labels.items() if label == "singlet.A1.1")
    dipole = -moments[(ground, 0, ground, 0, 3, 1, 0)]
    dipole_difference = dipole - target["ground_dipole_au"][2]
    assert abs(dipole_difference) <= 1e-5
    pi_error = max(
        float(np.max(abs(energies[spin + ".B1"] - energies[spin + ".B2"])))
        for spin in ("singlet", "triplet")
    )
    assert pi_error <= 1e-7
    assert pilot["unchanged_orbital_hashes_passed"] is True
    assert pilot["channel_import_verdict"] is True and pilot["scattering_channel_verdict"] is None
    return {
        "verified_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "reader_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_manifest": manifest,
        "verified_file_count": len(manifest),
        "binary_sectors": binary_rows,
        "maximum_native_reference_energy_difference_hartree": max(
            row["maximum_independent_energy_difference_hartree"] for row in binary_rows
        ),
        "maximum_denprop_binary_energy_difference_hartree": float(max(property_errors)),
        "maximum_native_residual_hartree": float(max(residuals)),
        "maximum_pi_energy_difference_hartree": pi_error,
        "maximum_vector_gram_error": max(row["vector_gram_maximum_error"] for row in binary_rows),
        "target_component_count": 48,
        "multipole_record_count": len(moments),
        "target_map": mapping,
        "a2_root9_native_ids": {
            spin: next(int(i) for i, v in labels.items() if v == spin + ".A2.9")
            for spin in ("singlet", "triplet")
        },
        "ground_dipole_difference_au": dipole_difference,
        "resources": resources,
        "target_sectors_total_wall_seconds": sum(s["wall_seconds"] for s in pilot["stages"][:-1]),
        "denprop_wall_seconds": pilot["stages"][-1]["wall_seconds"],
        "channel_import_verdict": True,
        "scattering_channel_verdict": None,
        "scope": "Full native-vector/property reconstruction on one fixed orbital frame. "
        "Multipole inventory is not independent validation of every transition moment; "
        "full channel convergence, global state identity and scattering remain open.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = read_target(args.directory)
    with args.output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
