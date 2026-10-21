"""Qualify one fifty-component orbital frame against covered native target roots."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import FortranFile

from projects.ukrmol_co.analyze import qc_properties, target_properties

IRREPS = ("A1", "B1", "B2", "A2")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def covered_sectors(scan, checkpoint_sha256):
    """Require the complete two-count/three-space coverage, not an aggregate flag."""
    expected = {
        (spin, irrep, roots, space)
        for spin in (0, 2)
        for irrep in IRREPS
        for roots in ((10, 13) if irrep == "A1" else (5, 8))
        for space in (40, 80, 160)
    }
    probes = scan["probes"]
    keys = [(p["spin"], p["irrep"], p["roots"], p["space"]) for p in probes]
    if len(keys) != len(expected) or set(keys) != expected or not scan["all_roots_pass"]:
        raise ValueError("Incomplete fifty-component coverage")
    for probe in probes:
        states = probe["states"]
        if not (
            probe["checkpoint_sha256"] == checkpoint_sha256
            and probe["cycles"] == 600
            and probe["all_roots_pass"]
            and len(states) == probe["roots"]
            and len(probe["ci_converged"]) == probe["roots"]
            and all(probe["ci_converged"])
        ):
            raise ValueError("Unqualified fixed-orbital coverage probe")
        for state in states:
            residuals = [
                state[k]
                for k in (
                    "physical_residual_norm_hartree",
                    "penalized_residual_norm_hartree",
                    "spin_penalty_action_norm_hartree",
                )
            ]
            if not np.isfinite(residuals).all() or min(residuals) < 0 or max(residuals) > 1e-9:
                raise ValueError("Independent coverage residual fails")
            if (
                not np.isfinite(state["spin_square"])
                or abs(state["spin_square"] - probe["spin"]) > 1e-6
            ):
                raise ValueError("Independent coverage spin fails")
    selected = {}
    for spin in (0, 2):
        for irrep in IRREPS:
            roots = 13 if irrep == "A1" else 8
            group = [p for p in probes if (p["spin"], p["irrep"]) == (spin, irrep)]
            reference = next(p for p in group if p["roots"] == roots and p["space"] == 160)
            energies = np.array([s["energy_hartree"] for s in reference["states"]])
            for probe in group:
                np.testing.assert_allclose(
                    [s["energy_hartree"] for s in probe["states"]],
                    energies[: probe["roots"]],
                    atol=1e-7,
                    rtol=0,
                )
            selected[(spin, irrep)] = energies
    return selected


def verify(root, contract):
    from pyscf import lib

    for name, want in contract["input_sha256"].items():
        if digest(root / name) != want:
            raise ValueError(f"Original input changed: {name}")
    run = root / "runs" / contract["run"]
    reference = root / "runs" / contract["reference"]

    def load(p):
        return json.loads(p.read_text())

    config, seed_config = [load(p / "config.json") for p in (run, reference)]
    qc_properties(reference, seed_config)
    for key in (
        "basis",
        "bond_length",
        "active_orbitals",
        "frozen_orbitals",
        "target_energy_tolerance",
        "target_gradient_tolerance",
        "target_ci_tolerance",
        "target_ci_residual_tolerance",
        "target_ci_lindep",
        "target_ah_lindep",
        "target_ah_tolerance",
        "target_ah_start_tolerance",
        "target_max_cycles",
    ):
        if config[key] != seed_config[key]:
            raise ValueError(f"Changed physical/QC setting: {key}")
    if load(run / "resources.json")["exit_code"] != 0:
        raise ValueError("Original native command did not pass")
    for spin in ("singlet", "triplet"):
        if config[f"sa_{spin}_roots"] != [10, 5, 5, 5] or config[f"target_{spin}_roots"] != [
            13,
            8,
            8,
            8,
        ]:
            raise ValueError("Requested orbital/native roots differ")
    scan = next(
        s for s in load(root / contract["coverage"])["full_scans"] if s["run"] == reference.name
    )
    independent = covered_sectors(scan, digest(reference / "co.casscf.chk"))
    native = target_properties(run, config)
    seed, fresh = [load(p / "target.json") for p in (reference, run)]
    optimization = fresh["optimization"]
    for key in ("ci_converged", "fresh_ci_converged"):
        assert len(optimization[key]) == 8 and all(optimization[key])
    assert fresh["converged"] and optimization["orbital_converged"]
    last = optimization["iterations"][-1]
    gradient = last.get("orbital_gradient_norm", last.get("orbital_ci_gradient_norm"))
    assert gradient is not None and np.isfinite(gradient) and gradient <= 1e-7
    assert len(fresh["states"]) == 50
    assert [(s["spin"], s["irrep"], s["root"]) for s in fresh["states"]] == [
        (s["spin"], s["irrep"], s["root"]) for s in seed["states"]
    ]
    np.testing.assert_allclose(
        optimization["fresh_state_energies_hartree"],
        [s["energy_hartree"] for s in fresh["states"]],
        atol=1e-7,
        rtol=0,
    )
    np.testing.assert_allclose(
        [s["energy_hartree"] for s in fresh["states"]],
        [s["energy_hartree"] for s in seed["states"]],
        atol=1e-7,
        rtol=0,
    )
    np.testing.assert_allclose(
        fresh["ground_dipole_au"], seed["ground_dipole_au"], atol=1e-5, rtol=0
    )
    np.testing.assert_allclose(
        native["ukrmol_ground_state_dipole_z_au"], seed["ground_dipole_au"][2], atol=1e-5, rtol=0
    )
    path = run / "output/CO/geom1/fort.26"
    sectors = []
    with FortranFile(path, "r") as stream:
        for irrep in IRREPS:
            for spin, label in ((0, "singlet"), (2, "triplet")):
                assert stream.read_record("u1").tobytes() == b"********                CIDATA  "
                header = stream.read_record("u1").tobytes()
                assert len(header) == 208
                assert header[200:] == bytes(8)
                assert header[16:136].decode().strip() == f"CO CI calculation - {label} {irrep}"
                centers, dimension, roots, symmetry = np.frombuffer(
                    header, "<i8", count=4, offset=136
                )
                assert centers == 2 and dimension > 0 and symmetry == IRREPS.index(irrep)
                assert roots == len(independent[(spin, irrep)])
                np.testing.assert_array_equal(
                    np.frombuffer(header, "<f8", count=2, offset=168), [spin / 2] * 2
                )
                assert np.frombuffer(header, "<i8", count=1, offset=184)[0] == 14
                core = np.frombuffer(header, "<f8", count=1, offset=192)[0]
                charges = []
                for _ in range(centers):
                    record = stream.read_record("u1").tobytes()
                    assert len(record) == 40
                    values = np.frombuffer(record, "<f8", offset=8)
                    assert np.isfinite(values).all()
                    charges.append(values[-1])
                np.testing.assert_array_equal(charges, [6, 8])
                record = stream.read_record("<f8")
                assert len(record) == 2 * dimension + roots
                energies = record[dimension : dimension + roots] + core
                assert np.isfinite(energies).all() and np.all(np.diff(energies) >= 0)
                error = float(np.max(abs(energies - independent[(spin, irrep)])))
                assert error <= 1e-7
                vectors = []
                for index in range(roots):
                    record = stream.read_record("u1").tobytes()
                    assert np.frombuffer(record, "<i8", count=1)[0] == index + 1
                    vector = np.frombuffer(record, "<f8", offset=8)
                    assert len(vector) == dimension and np.isfinite(vector).all()
                    vectors.append(vector)
                gram_error = float(
                    np.max(abs(np.asarray(vectors) @ np.asarray(vectors).T - np.eye(roots)))
                )
                assert gram_error <= 1e-9
                sectors.append(
                    dict(
                        sector=f"{label}.{irrep}",
                        roots=int(roots),
                        dimension=int(dimension),
                        energy_error_hartree=error,
                        vector_gram_error=gram_error,
                    )
                )
        assert stream._fp.tell() == path.stat().st_size
    old, new = reference / "co.casscf.chk", run / "output/CO/geom1/co.casscf.chk"
    lib.num_threads(1)
    mol = lib.chkfile.load_mol(str(old))
    np.testing.assert_allclose(
        mol.atom_coords(), lib.chkfile.load_mol(str(new)).atom_coords(), atol=1e-12, rtol=0
    )
    coeffs = [lib.chkfile.load(str(p), "mcscf/mo_coeff") for p in (old, new)]
    metric = mol.intor_symmetric("int1e_ovlp")
    overlaps = {
        name: np.linalg.svd(coeffs[0][:, s].T @ metric @ coeffs[1][:, s], compute_uv=False).tolist()
        for name, s in (("core", slice(0, 2)), ("active", slice(2, 13)))
    }
    assert min(overlaps["core"] + overlaps["active"]) >= 1 - 1e-8
    return dict(
        import_verdict=True,
        roots_verified=sum(s["roots"] for s in sectors),
        sectors=sectors,
        subspace_singular_values=overlaps,
        native_properties=native,
        input_sha256=contract["input_sha256"],
        native_cidata_sha256=digest(path),
        scattering_ensemble_verdict=None,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root, json.loads(args.contract.read_text()))
    with args.output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
