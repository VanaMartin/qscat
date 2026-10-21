"""Own eight fixed-orbital CI/physical-overlap comparisons, with no QC rerun."""

from __future__ import annotations

import argparse
import contextlib
import json
import time
from math import comb
from pathlib import Path

import numpy as np

from projects.ukrmol_co.ci_overlap import frozen_core_overlaps
from projects.ukrmol_co.coverage_reference import benchmark, state_residual
from projects.ukrmol_co.ensemble_frame_overlaps import (
    decision_controls,
    diagnose_overlap,
    full_minor_overlaps,
)
from projects.ukrmol_co.fifty_import_verification import IRREPS, covered_sectors, digest


def save(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def load(path):
    return json.loads(path.read_text())


def verify_inputs(root, contract):
    for name, want in contract["input_sha256"].items():
        if digest(root / name) != want:
            raise ValueError(f"Changed pinned input: {name}")


def validate_cached(row, vectors, checkpoint_sha256, native, spin, irrep):
    """Reconstruct saved CI gates; residual actions themselves are not rerun."""
    label = f"{'singlet' if spin == 0 else 'triplet'}.{irrep}"
    assert row["sector"] == label and row["checkpoint_sha256"] == checkpoint_sha256
    assert row["gate_passed"] is True and row["ci_converged"] == [True] * 12
    assert row["ci_file"] == label + ".ci.npy" and row["native_comparison_root_count"] == 5
    electrons = ((10 + spin) // 2, (10 - spin) // 2)
    assert vectors.shape == (12, *(comb(11, n) for n in electrons))
    assert np.isfinite(vectors).all()
    energies = np.asarray(row["energy_hartree"])
    assert energies.shape == (12,) and np.isfinite(energies).all()
    assert np.all(np.diff(energies) >= -1e-10)
    for key in ("physical_residual_norms_hartree", "penalized_residual_norms_hartree"):
        values = np.asarray(row[key])
        assert values.shape == (12,) and np.isfinite(values).all()
        assert min(values) >= 0 and max(values) <= 1e-9
    np.testing.assert_allclose(row["spin_square"], np.full(12, spin), atol=1e-6, rtol=0)
    np.testing.assert_allclose(
        row["original_eight_energy_errors_hartree"], np.zeros(8), atol=1e-7, rtol=0
    )
    weights = np.asarray(row["original_eight_energy_matched_projection_weights"])
    assert weights.shape == (8,) and np.isfinite(weights).all()
    assert min(weights) >= 0.99999999 and max(weights) <= 1 + 1e-9
    errors = energies[:5] - [native[f"{label}.{i + 1}"] for i in range(5)]
    np.testing.assert_allclose(errors, np.zeros(5), atol=1e-7, rtol=0)
    np.testing.assert_allclose(errors, row["native_root_errors_hartree"], atol=1e-12, rtol=0)
    flat = vectors.reshape(12, -1)
    gram_error = float(np.max(abs(flat.conj() @ flat.T - np.eye(12))))
    assert gram_error <= 1e-9
    return dict(
        cached_ci_gates_reconstructed=True,
        ci_gram_error=gram_error,
        maximum_native_energy_error_hartree=float(max(abs(errors))),
        residual_actions_recomputed=False,
    )


def qualified_inputs(root, contract):
    old_checkpoint = root / contract["old_checkpoint"]
    new_checkpoint = root / contract["new_checkpoint"]
    old_owner, old_result = (load(root / contract[k]) for k in ("old_execution", "old_result"))
    assert old_owner["original_controller_exit"] == 0
    assert len(old_owner["steps"]) == 1 and old_owner["steps"][0]["original_exit"] == 0
    assert old_result["all_sectors_and_pairs_complete"] is True
    old_anchor = next(a for a in old_result["anchors"] if a["run"] == contract["old_run"])
    assert len(old_anchor["sectors"]) == 8
    for path, flag in (
        ("frame_execution", "new_frame_import_verdict"),
        ("property_execution", "new_property_export_verdict"),
    ):
        owner = load(root / contract[path])
        assert owner[flag] is True and owner["original_owner_exit"] == 0
        assert [s["original_command_exit"] for s in owner["steps"]] == [0]
        assert [s["original_container_exit"] for s in owner["steps"]] == [0]
    frame = load(root / contract["frame_result"])
    assert frame["new_frame_import_verdict"] is True
    assert frame["original_seed_preservation_passed"] is False
    reference = covered_sectors(frame, digest(new_checkpoint))
    properties = load(root / contract["property_result"])
    assert properties["new_property_export_verdict"] is True
    assert properties["frame_checkpoint_sha256"] == digest(new_checkpoint)
    assert properties["property_inventory"]["target_count"] == 40
    assert properties["selected_property_sha256"] == digest(root / contract["selected_properties"])
    rejected = load(root / contract["original_import_execution"])
    assert rejected["original_owner_exit"] == 1
    configs = [
        load(root / "runs" / contract[key] / "config.json") for key in ("old_run", "new_run")
    ]
    for c, counts in zip(configs, ([5, 5, 5, 5], [10, 5, 5, 5]), strict=True):
        assert c["basis"] == "cc-pVDZ" and c["bond_length"] == 2.1323
        assert c["active_orbitals"] == [5, 3, 3, 0] and c["frozen_orbitals"] == 2
        assert c["sa_singlet_roots"] == c["sa_triplet_roots"] == counts
    table = root / contract["old_native_energies"]
    header = next(s for s in table.read_text().splitlines() if "R_bohr" in s).split()[2:]
    data = np.loadtxt(table, ndmin=2)
    assert data.shape == (1, len(header) + 1) and abs(data[0, 0] - 2.1323) <= 1e-10
    native = dict(zip(header, data[0, 1:], strict=True))
    assert np.isfinite(list(native.values())).all()
    return old_checkpoint, new_checkpoint, old_anchor, native, reference


def physical_metrics(old_checkpoint, new_checkpoint, output):
    from pyscf import gto, lib

    molecules, coefficients, self_metrics = [], [], []
    for name, checkpoint in (("old", old_checkpoint), ("new", new_checkpoint)):
        mol = lib.chkfile.load_mol(str(checkpoint))
        assert int(lib.chkfile.load(str(checkpoint), "mcscf/ncore")) == 2
        assert int(lib.chkfile.load(str(checkpoint), "mcscf/ncas")) == 11
        c = np.asarray(lib.chkfile.load(str(checkpoint), "mcscf/mo_coeff"))[:, :13]
        metric = mol.intor_symmetric("int1e_ovlp")
        own = c.conj().T @ metric @ c
        np.testing.assert_allclose(own, np.eye(13), atol=1e-9, rtol=0)
        for suffix, value in (("mo", c), ("ao", metric), ("self", own)):
            np.save(output / f"{name}.{suffix}.npy", value)
        molecules.append(mol)
        coefficients.append(c)
        self_metrics.append(own)
    left, right = molecules
    assert left.ao_labels() == right.ao_labels()
    np.testing.assert_array_equal(left.atom_charges(), right.atom_charges())
    np.testing.assert_allclose(left.atom_coords(), right.atom_coords(), atol=1e-12, rtol=0)
    cross = gto.intor_cross("int1e_ovlp", left, right)
    reverse = gto.intor_cross("int1e_ovlp", right, left)
    np.testing.assert_allclose(cross, reverse.conj().T, atol=1e-12, rtol=0)
    a, b = coefficients
    forward, backward = a.conj().T @ cross @ b, b.conj().T @ reverse @ a
    np.testing.assert_allclose(forward, backward.conj().T, atol=1e-12, rtol=0)
    for name, value in (
        ("cross.ao", cross),
        ("reverse.ao", reverse),
        ("cross.mo", forward),
        ("reverse.mo", backward),
    ):
        np.save(output / f"{name}.npy", value)
    assert np.linalg.svd(forward, compute_uv=False)[0] <= 1 + 1e-9
    return forward, backward, self_metrics


def new_wavefunctions(checkpoint, spin, irrep, expected, output):
    from pyscf import fci, lib, mcscf, scf

    roots = 13 if irrep == "A1" else 8
    label = f"{'singlet' if spin == 0 else 'triplet'}.{irrep}"
    started = time.monotonic()
    with (output / f"{label}.log").open("x") as log, contextlib.redirect_stdout(log):
        mol = lib.chkfile.load_mol(str(checkpoint))
        mol.verbose = 4
        c = lib.chkfile.load(str(checkpoint), "mcscf/mo_coeff")
        electrons = ((10 + spin) // 2, (10 - spin) // 2)
        cas = mcscf.CASCI(scf.RHF(mol), 11, electrons, ncore=2)
        cas.canonicalization = False
        solver = fci.direct_spin1_symm.FCI(mol)
        solver.spin, solver.wfnsym = spin, irrep
        solver.nroots, solver.max_space, solver.max_cycle = roots, 160, 600
        solver.conv_tol, solver.conv_tol_residual, solver.lindep = 1e-12, 1e-9, 1e-22
        cas.fcisolver = fci.addons.fix_spin_(solver, shift=1.0, ss=spin)
        cas.kernel(c)
        vectors, energies = np.asarray(cas.ci), np.atleast_1d(cas.e_tot)
        np.save(output / f"{label}.ci.npy", vectors)
        assert vectors.shape == (roots, *(comb(11, n) for n in electrons))
        h1, ecore = cas.get_h1eff(c)
        h2 = cas.get_h2eff(c)
        effective = cas.fcisolver.absorb_h1e(h1, h2, 11, electrons, 0.5)
        states = [
            state_residual(cas.fcisolver, effective, ci, energy, ecore, 11, electrons)
            for ci, energy in zip(vectors, energies, strict=True)
        ]
    flat = vectors.reshape(roots, -1)
    gram_error = float(np.max(abs(flat.conj() @ flat.T - np.eye(roots))))
    flags = np.atleast_1d(cas.fcisolver.converged).tolist()
    errors = energies - expected
    passed = bool(
        len(flags) == roots
        and all(flags)
        and gram_error <= 1e-9
        and np.isfinite(vectors).all()
        and np.isfinite(errors).all()
        and max(abs(errors)) <= 1e-7
        and all(
            np.isfinite(list(s.values())).all()
            and abs(s["spin_square"] - spin) <= 1e-6
            and all(
                0 <= s[k] <= 1e-9
                for k in (
                    "physical_residual_norm_hartree",
                    "penalized_residual_norm_hartree",
                    "spin_penalty_action_norm_hartree",
                )
            )
            for s in states
        )
    )
    row = dict(
        sector=label,
        roots=roots,
        space=160,
        cycles=600,
        states=states,
        ci_converged=flags,
        ci_vector_gram_error=gram_error,
        reference_energy_errors_hartree=errors.tolist(),
        numerical_gates_passed=passed,
        checkpoint_sha256=digest(checkpoint),
        ci_sha256=digest(output / f"{label}.ci.npy"),
        wall_seconds=time.monotonic() - started,
    )
    save(output / f"{label}.json", row)
    if not passed:
        raise ValueError(f"New fixed-orbital wavefunction gate failed: {label}")
    return vectors, row


def compare_sector(left, right, forward, backward, self_metrics, electrons, output, label):
    """Persist independent routes before asserting their numerical agreement."""
    started = time.monotonic()
    own_errors = []
    for ci, metric in zip((left, right), self_metrics, strict=True):
        own = frozen_core_overlaps(ci, ci, metric, 2, electrons)
        own_errors.append(float(np.max(abs(own - np.eye(len(ci))))))
    overlap = frozen_core_overlaps(left, right, forward, 2, electrons)
    explicit = full_minor_overlaps(left, right, forward, 2, electrons)
    reverse = frozen_core_overlaps(right, left, backward, 2, electrons)
    for suffix, value in (("overlap", overlap), ("full-minors", explicit), ("reverse", reverse)):
        np.save(output / f"{label}.{suffix}.npy", value)
    row = dict(
        maximum_full_minor_error=float(np.max(abs(overlap - explicit))),
        maximum_reciprocity_error=float(np.max(abs(overlap - reverse.conj().T))),
        own_frame_physical_gram_errors=own_errors,
    )
    save(output / f"{label}.numerics.json", row)
    np.testing.assert_allclose(overlap, explicit, atol=1e-10, rtol=0)
    np.testing.assert_allclose(overlap, reverse.conj().T, atol=1e-10, rtol=0)
    assert max(own_errors) <= 1e-9
    row.update(
        diagnose_overlap(overlap),
        numerical_gates_passed=True,
        overlap_sha256=digest(output / f"{label}.overlap.npy"),
        wall_seconds=time.monotonic() - started,
    )
    save(output / f"{label}.comparison.json", row)
    return row


def run(root, output, contract):
    from pyscf import lib

    lib.num_threads(4)
    verify_inputs(root, contract)
    output.mkdir(exist_ok=False)
    report = dict(
        original_seed_preservation_passed=False,
        original_seed_exits=[1, 1, 1],
        retained_subspace_dominance_verdict=None,
        scattering_ensemble_verdict=None,
        global_state_identity_verdict=None,
        sectors=[],
    )
    try:
        report["decision_controls"] = decision_controls()
        report["residual_analytic_benchmark"] = benchmark()
        old_checkpoint, new_checkpoint, anchor, native, reference = qualified_inputs(root, contract)
        forward, backward, self_metrics = physical_metrics(old_checkpoint, new_checkpoint, output)
        save(output / "progress.json", report)
        for spin in (0, 2):
            for irrep in IRREPS:
                label = f"{'singlet' if spin == 0 else 'triplet'}.{irrep}"
                cached = root / contract["old_cache"]
                metadata = load(cached / f"{label}.json")
                assert metadata == next(s for s in anchor["sectors"] if s["sector"] == label)
                left = np.load(cached / f"{label}.ci.npy", allow_pickle=False)
                old_checks = validate_cached(
                    metadata, left, digest(old_checkpoint), native, spin, irrep
                )
                right, new_checks = new_wavefunctions(
                    new_checkpoint, spin, irrep, reference[(spin, irrep)], output
                )
                electrons = ((10 + spin) // 2, (10 - spin) // 2)
                comparison = compare_sector(
                    left, right, forward, backward, self_metrics, electrons, output, label
                )
                report["sectors"].append(
                    dict(sector=label, cached=old_checks, new_ci=new_checks, comparison=comparison)
                )
                save(output / "progress.json", report)
        report["all_eight_sectors_complete"] = True
        report["retained_subspace_dominance_verdict"] = all(
            s["comparison"]["retained_subspace_dominance_passed"] for s in report["sectors"]
        )
        report["scope"] = (
            "Finite retained-five subspace dominance only; no global identity, "
            "model convergence or scattering-error bound. No QC/native/scattering call."
        )
        verify_inputs(root, contract)
        save(output / "result.json", report)
    except BaseException as error:
        report["error"] = repr(error)
        save(output / "failure.json", report)
        raise
    finally:
        verify_inputs(root, contract)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.root, args.output, load(args.contract))
    raise SystemExit(0 if result["retained_subspace_dominance_verdict"] else 1)


if __name__ == "__main__":
    main()
