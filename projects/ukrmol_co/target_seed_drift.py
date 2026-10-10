"""Separate fixed-orbital CI error from seed-to-import orbital-frame drift.

This diagnostic preserves a failed import gate. It never resumes an import queue.
The contract freezes both checkpoints, solver controls and output location.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import itertools
import json
import os
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.run import run_profiled


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def _actions(fci, solver, effective, vectors, energies, ecore, ncas, electrons):
    rows = []
    expected_spin = solver.spin * (solver.spin + 2) / 4
    for vector, energy in zip(vectors, energies, strict=True):
        norm = np.linalg.norm(vector)
        physical = fci.direct_spin1.contract_2e(effective, vector, ncas, electrons)
        penalized = solver.contract_2e(effective, vector, ncas, electrons)
        rayleigh = float(np.vdot(vector, physical).real / norm**2 + ecore)
        rows.append(
            {
                "energy_hartree": float(energy),
                "physical_rayleigh_energy_hartree": rayleigh,
                "physical_residual_hartree": float(
                    np.linalg.norm(physical - (energy - ecore) * vector) / norm
                ),
                "penalized_residual_hartree": float(
                    np.linalg.norm(penalized - (energy - ecore) * vector) / norm
                ),
                "spin_penalty_action_hartree": float(np.linalg.norm(penalized - physical) / norm),
                "spin_square_error": float(
                    abs(fci.spin_op.spin_square(vector, ncas, electrons)[0] - expected_spin)
                ),
            }
        )
    return rows


def analytic_controls():
    """Gate the action/residual calculation on an offset Hubbard dimer."""
    from pyscf import fci

    h1 = np.array([[0.0, -1.0], [-1.0, 0.0]])
    h2 = np.zeros((2, 2, 2, 2))
    h2[0, 0, 0, 0] = h2[1, 1, 1, 1] = 4.0
    solver = fci.direct_spin1.FCI()
    solver.spin = 0
    solver = fci.addons.fix_spin_(solver, shift=1.0, ss=0.0)
    solver.conv_tol, solver.conv_tol_residual = 1e-14, 1e-12
    energy, ci = solver.kernel(h1, h2, 2, (1, 1), ecore=3.0)
    expected = 3.0 + (4.0 - np.sqrt(32.0)) / 2
    np.testing.assert_allclose(energy, expected, atol=1e-12, rtol=0)
    effective = fci.direct_spin1.absorb_h1e(h1, h2, 2, (1, 1), 0.5)
    good = _actions(fci, solver, effective, [ci], [energy], 3.0, 2, (1, 1))[0]
    assert good["physical_residual_hartree"] <= 1e-12
    assert good["penalized_residual_hartree"] <= 1e-12
    assert good["spin_penalty_action_hartree"] <= 1e-12
    wrong = np.array(ci, copy=True)
    wrong[0, 0] += 0.01
    bad = _actions(fci, solver, effective, [wrong], [energy], 3.0, 2, (1, 1))[0]
    assert bad["physical_residual_hartree"] >= 1e-3
    return {"exact_ground_energy_hartree": float(expected), "exact": good, "perturbed": bad}


def _probe(mol, coefficients, contract, spin, irrep, space, directory):
    from pyscf import fci, mcscf, scf

    electrons = ((10 + spin) // 2, (10 - spin) // 2)
    mc = mcscf.CASCI(scf.RHF(mol), 12, electrons, ncore=2)
    mc.canonicalization = False
    solver = fci.direct_spin1_symm.FCI(mol)
    solver.spin, solver.wfnsym, solver.nroots = spin, irrep, 8
    solver.max_space, solver.max_cycle = space, contract["ci_cycles"]
    solver.conv_tol = contract["ci_energy_tolerance_hartree"]
    solver.conv_tol_residual = contract["ci_residual_tolerance_hartree"]
    solver.lindep = contract["ci_lindep"]
    mc.fcisolver = fci.addons.fix_spin_(
        solver, shift=contract["spin_penalty_shift_hartree"], ss=spin * (spin + 2) / 4
    )
    started = time.monotonic()
    with (directory / "casci.log").open("x") as stream, contextlib.redirect_stdout(stream):
        mc.kernel(coefficients)
        h1, ecore = mc.get_h1eff(coefficients)
        h2 = mc.get_h2eff(coefficients)
        effective = fci.direct_spin1.absorb_h1e(h1, h2, 12, electrons, 0.5)
        energies, vectors = np.asarray(mc.e_tot), np.asarray(mc.ci)
        assert len(energies) == len(vectors) == 8
        actions = _actions(fci, mc.fcisolver, effective, vectors, energies, ecore, 12, electrons)
    orthogonality = float(
        np.max(abs(vectors.reshape(8, -1) @ vectors.reshape(8, -1).T - np.eye(8)))
    )
    flags = np.atleast_1d(mc.fcisolver.converged).tolist()
    gates = contract["gates"]
    passed = bool(
        all(flags)
        and orthogonality <= gates["ci_and_mo_orthogonality_error"]
        and all(
            np.isfinite(list(row.values())).all()
            and row["spin_square_error"] <= gates["spin_square_error"]
            and max(
                row[k]
                for k in (
                    "physical_residual_hartree",
                    "penalized_residual_hartree",
                    "spin_penalty_action_hartree",
                )
            )
            <= gates["physical_penalized_and_action_residual_hartree"]
            for row in actions
        )
    )
    if space == max(contract["trial_spaces"]):
        np.save(directory / "ci.npy", vectors)
    return {
        "spin": spin,
        "irrep": irrep,
        "space": space,
        "roots": 8,
        "cycles": contract["ci_cycles"],
        "states": actions,
        "ci_converged": flags,
        "ci_orthogonality_error": orthogonality,
        "physical_gates_pass": passed,
        "wall_seconds": time.monotonic() - started,
    }


def compute(root: Path, work: Path, contract_path: Path) -> dict:
    """Run both frozen frames and retain seed/native comparisons separately."""
    from pyscf import lib

    contract = json.loads(contract_path.read_text())
    lib.num_threads(contract["resources"]["library_threads"])
    assert contract["ncore"] == 2 and contract["ncas"] == 12
    assert contract["root_counts"] == [8] and contract["probes"] == 48
    report = {
        "analytic_controls": analytic_controls(),
        "frames": [],
        "contract_sha256": _sha(contract_path),
    }
    frames = []
    for frame in contract["frames"]:
        run = root / "runs" / frame["run"]
        checkpoint = run / frame["checkpoint"]
        assert _sha(checkpoint) == frame["checkpoint_sha256"]
        config = json.loads((run / "config.json").read_text())
        assert config["basis"] == contract["basis"]
        assert config["bond_length"] == contract["geometry_bohr"]
        assert sum(config["active_orbitals"]) == 12
        assert json.loads((run / "resources.json").read_text())["exit_code"] == 0
        target = json.loads((run / "target.json").read_text())
        mol = lib.chkfile.load_mol(str(checkpoint))
        mol.verbose = 0
        coefficients = lib.chkfile.load(str(checkpoint), "mcscf/mo_coeff")
        assert int(lib.chkfile.load(str(checkpoint), "mcscf/ncore")) == 2
        assert int(lib.chkfile.load(str(checkpoint), "mcscf/ncas")) == 12
        metric = mol.intor_symmetric("int1e_ovlp")
        np.testing.assert_allclose(
            coefficients.T @ metric @ coefficients, np.eye(coefficients.shape[1]), atol=1e-9, rtol=0
        )
        result = {
            "label": frame["label"],
            "checkpoint_sha256": frame["checkpoint_sha256"],
            "probes": [],
        }
        report["frames"].append(result)
        for spin, irrep, space in itertools.product(
            (0, 2), ("A1", "B1", "B2", "A2"), contract["trial_spaces"]
        ):
            path = work / frame["label"] / f"spin{spin}-{irrep}-space{space}"
            path.mkdir(parents=True)
            row = _probe(mol, coefficients, contract, spin, irrep, space, path)
            saved = [
                s["energy_hartree"]
                for s in target["states"]
                if s["spin"] == ("singlet" if spin == 0 else "triplet") and s["irrep"] == irrep
            ]
            row["maximum_saved_target_error_hartree"] = float(
                np.max(abs(np.asarray([s["energy_hartree"] for s in row["states"][:5]]) - saved))
            )
            _save(path / "result.json", row)
            result["probes"].append(row)
            _save(work / "progress.json", report)
        energies = {
            (row["spin"], row["irrep"], row["space"]): np.asarray(
                [s["energy_hartree"] for s in row["states"]]
            )
            for row in result["probes"]
        }
        widest = max(contract["trial_spaces"])
        result["maximum_space_energy_error_hartree"] = max(
            float(np.max(abs(values - energies[spin, irrep, widest])))
            for (spin, irrep, space), values in energies.items()
        )
        result["maximum_pi_energy_error_hartree"] = max(
            float(np.max(abs(energies[spin, "B1", space] - energies[spin, "B2", space])))
            for spin, space in itertools.product((0, 2), contract["trial_spaces"])
        )
        gate = contract["gates"]["saved_target_space_and_pi_energy_error_hartree"]
        result["fixed_orbital_coverage_pass"] = bool(
            all(
                row["physical_gates_pass"] and row["maximum_saved_target_error_hartree"] <= gate
                for row in result["probes"]
            )
            and max(
                result["maximum_space_energy_error_hartree"],
                result["maximum_pi_energy_error_hartree"],
            )
            <= gate
        )
        if (run / "output/CO/target.energies").exists():
            header = next(
                line
                for line in (run / "output/CO/target.energies").read_text().splitlines()
                if "R_bohr" in line
            ).split()[2:]
            native = dict(
                zip(
                    header,
                    np.loadtxt(run / "output/CO/target.energies", ndmin=2)[0, 1:],
                    strict=True,
                )
            )
            result["maximum_native_root_error_hartree"] = float(
                max(
                    abs(
                        float(e)
                        - native[f"{'singlet' if spin == 0 else 'triplet'}.{irrep}.{i + 1}"]
                    )
                    for (spin, irrep, space), values in energies.items()
                    if space == widest
                    for i, e in enumerate(values)
                )
            )
            result["native_root_gate_pass"] = (
                result["maximum_native_root_error_hartree"]
                <= contract["gates"]["native_root_error_hartree"]
            )
        assert _sha(checkpoint) == frame["checkpoint_sha256"]
        frames.append((mol, coefficients, energies))
    mol, old, previous = frames[0]
    new_mol, new, current = frames[1]
    np.testing.assert_allclose(mol.atom_coords(), new_mol.atom_coords(), atol=1e-12, rtol=0)
    assert mol.ao_labels() == new_mol.ao_labels()
    metric = mol.intor_symmetric("int1e_ovlp")
    report["subspace_singular_values"] = {
        label: np.linalg.svd(
            old[:, selection].T @ metric @ new[:, selection], compute_uv=False
        ).tolist()
        for label, selection in (("core", slice(0, 2)), ("active", slice(2, 14)))
    }
    space = max(contract["trial_spaces"])
    report["fixed_orbital_frame_energy_differences_hartree"] = {
        f"{'singlet' if spin == 0 else 'triplet'}.{irrep}": (
            current[spin, irrep, space] - previous[spin, irrep, space]
        ).tolist()
        for spin, irrep in itertools.product((0, 2), ("A1", "B1", "B2", "A2"))
    }
    report["original_seed_gate_remains_rejected"] = True
    report["diagnostic_gates_pass"] = all(
        f["fixed_orbital_coverage_pass"] and f.get("native_root_gate_pass", True)
        for f in report["frames"]
    )
    report["scope"] = contract["release"]
    _save(work / "result.json", report)
    return report


def main() -> None:
    """Run a bounded, profiled fixed-checkpoint diagnostic."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("/work"))
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--profiled", action="store_true")
    args = parser.parse_args()
    if not args.profiled:
        args.workdir.mkdir()
        raise SystemExit(
            run_profiled(
                args.workdir,
                os.environ.copy(),
                [
                    "python3",
                    "-m",
                    "projects.ukrmol_co.target_seed_drift",
                    "--root",
                    str(args.root),
                    "--workdir",
                    str(args.workdir),
                    "--contract",
                    str(args.contract),
                    "--profiled",
                ],
            )
        )
    result = compute(args.root, args.workdir, args.contract)
    print(
        json.dumps(
            {
                "diagnostic_gates_pass": result["diagnostic_gates_pass"],
                "original_seed_gate_remains_rejected": True,
            }
        )
    )
    raise SystemExit(0 if result["diagnostic_gates_pass"] else 1)


if __name__ == "__main__":
    main()
