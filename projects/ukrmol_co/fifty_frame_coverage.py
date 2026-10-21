"""Independently qualify the completed fifty-component frame without reoptimization."""

from __future__ import annotations

import argparse
import contextlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

from projects.ukrmol_co.analyze import target_properties
from projects.ukrmol_co.coverage_reference import benchmark, solve
from projects.ukrmol_co.fifty_import_verification import (
    IRREPS,
    covered_sectors,
    digest,
    native_binary_sectors,
)


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def checkpoint_export_overlap(checkpoint, exported):
    """Compare the new fixed checkpoint with its imported Molden, not the old seed."""
    from pyscf import gto, lib
    from pyscf.tools import molden

    lib.num_threads(4)
    mol = lib.chkfile.load_mol(str(checkpoint))
    coeff = lib.chkfile.load(str(checkpoint), "mcscf/mo_coeff")
    assert int(lib.chkfile.load(str(checkpoint), "mcscf/ncore")) == 2
    assert int(lib.chkfile.load(str(checkpoint), "mcscf/ncas")) == 11
    other, _, imported, _, _, _ = molden.load(str(exported))
    np.testing.assert_allclose(mol.atom_coords(), other.atom_coords(), atol=1e-12, rtol=0)
    np.testing.assert_array_equal(mol.atom_charges(), other.atom_charges())
    metric = gto.intor_cross("int1e_ovlp", mol, other)
    values = {
        label: np.linalg.svd(coeff[:, s].T @ metric @ imported[:, s], compute_uv=False).tolist()
        for label, s in (("core", slice(0, 2)), ("active", slice(2, 13)))
    }
    if (
        not np.isfinite(values["core"] + values["active"]).all()
        or min(values["core"] + values["active"]) < 1 - 1e-8
    ):
        raise ValueError("Completed checkpoint and imported Molden frame differ")
    return values


def verify_inputs(root, contract):
    for name, want in contract["input_sha256"].items():
        if digest(root / name) != want:
            raise ValueError(f"Pinned original file changed: {name}")
    previous = json.loads((root / contract["original_execution"]).read_text())
    assert previous["original_owner_exit"] == 1
    assert [s["original_command_exit"] for s in previous["steps"]] == [0, 0, 0, 1]


def run(root, output, contract):
    verify_inputs(root, contract)
    output.mkdir(exist_ok=False)
    native_run = root / "runs" / contract["run"]
    checkpoint = native_run / "output/CO/geom1/co.casscf.chk"
    checksum = digest(checkpoint)
    target = json.loads((native_run / "target.json").read_text())
    config = json.loads((native_run / "config.json").read_text())
    props = target_properties(native_run, config)
    property_counts = dict(
        Counter(
            t[-2].lower() + "." + t[-1]
            for line in (native_run / "output/CO/geom1/fort.24").read_text().splitlines()
            if (t := line.split()) and t[0] == "5"
        )
    )
    assert property_counts == contract["saved_property_counts"]
    overlaps = checkpoint_export_overlap(checkpoint, native_run / "output/CO/geom1/co.molden")
    report = dict(
        residual_analytic_benchmark=benchmark(),
        checkpoint_sha256=checksum,
        checkpoint_export_subspace_singular_values=overlaps,
        original_seed_preservation_passed=False,
        original_seed_verifier_and_owner_exits=[1, 1, 1],
        new_frame_coverage_verdict=None,
        new_frame_import_verdict=None,
        scattering_ensemble_verdict=None,
        saved_property_counts=property_counts,
        matched_five_A2_property_manifold_complete=False,
        probes=[],
    )
    save(output / "progress.json", report)
    for spin in (0, 2):
        for irrep in IRREPS:
            count = 10 if irrep == "A1" else 5
            expected = [
                s["energy_hartree"]
                for s in target["states"]
                if s["spin"] == ("singlet" if spin == 0 else "triplet") and s["irrep"] == irrep
            ]
            assert len(expected) == count
            for roots in (count, count + 3):
                for space in (40, 80, 160):
                    name = f"spin{spin}-{irrep}-roots{roots}-space{space}"
                    with (
                        (output / (name + ".log")).open("x") as log,
                        contextlib.redirect_stdout(log),
                    ):
                        row = solve(checkpoint, spin, irrep, roots, space, 600)
                    assert digest(checkpoint) == checksum
                    energies = np.array([s["energy_hartree"] for s in row["states"][:count]])
                    error = float(np.max(abs(energies - expected)))
                    row.update(
                        checkpoint_sha256=checksum, maximum_ensemble_energy_difference_hartree=error
                    )
                    row["all_roots_pass"] = bool(
                        row["all_roots_pass"]
                        and error <= 1e-7
                        and all(
                            s["spin_penalty_action_norm_hartree"] <= 1e-9 for s in row["states"]
                        )
                    )
                    report["probes"].append(row)
                    save(output / (name + ".json"), row)
                    save(output / "progress.json", report)
    report["all_roots_pass"] = all(p["all_roots_pass"] for p in report["probes"])
    save(output / "coverage.json", report)
    independent = covered_sectors(report, checksum)
    report["new_frame_coverage_verdict"] = True
    report["binary_sectors"] = native_binary_sectors(
        native_run / "output/CO/geom1/fort.26", independent
    )
    assert sum(s["roots"] for s in report["binary_sectors"]) == 74
    report["native_properties"] = props
    report["new_frame_import_verdict"] = True
    report["scope"] = (
        "Independent qualification of the completed new fixed orbital frame only; "
        "original seed-preservation rejection unchanged. No orbitals optimized, "
        "no native engine called, no matched scattering result."
    )
    verify_inputs(root, contract)
    save(output / "result.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    args = parser.parse_args()
    run(args.root, args.output, json.loads(args.contract.read_text()))


if __name__ == "__main__":
    main()
