"""Build common CO orbitals by multi-spin/multi-irrep state-averaged CASSCF.

PySCF is an optional dependency of the CO Docker layers. The upstream scripts
invoke this module through their quantum-chemistry adapter; its JSON and Molden
outputs are independently checked against UKRmol's target diagonalizations.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import time
from pathlib import Path

import numpy as np

IRREPS = ("A1", "B1", "B2", "A2")


def ci_trial_space(roots: int, override: int | None) -> int:
    """Size the CI search space independently of the orbital ensemble."""
    space = max(40, 8 * roots) if override is None else override
    if roots < 1 or space <= roots:
        raise ValueError("CI trial space must exceed the requested ensemble roots")
    return space


def _fresh_ci_kernel(kernel):
    def solve(h1e, eri, norb, nelec, ci0=None, **kwargs):
        return kernel(h1e, eri, norb, nelec, ci0=None, **kwargs)

    return solve


def _sector_ci_driver(fci, spin: int, irrep: str, config: dict):
    """Select the experimental singlet-A1 driver without altering other sectors."""
    choice = config.get("target_singlet_a1_driver", "spin1")
    if choice not in ("spin1", "spin0"):
        raise ValueError("Singlet-A1 CI driver must be spin1 or spin0")
    if choice == "spin0" and spin == 0 and irrep == "A1":
        return fci.direct_spin0_symm.FCI
    return fci.direct_spin1_symm.FCI


def project_initial_orbitals(mc, previous_mo: np.ndarray, previous_mol) -> np.ndarray:
    """Project a checkpoint's core/active spaces, completing destination virtuals."""
    from pyscf import mcscf

    # A full larger-basis MO matrix exceeds PySCF's destination-column limit.
    # Only inactive core + active orbitals define this CAS; PySCF constructs
    # the orthogonal virtual complement in the destination basis itself.
    if previous_mo.shape[1] > mc._scf.mo_coeff.shape[1]:
        previous_mo = previous_mo[:, : mc.ncore + mc.ncas]
    return mcscf.project_init_guess(mc, previous_mo, prev_mol=previous_mol)


def build_target(config: dict, molden_path: Path) -> dict:
    """Optimize an equally weighted ensemble and export core/active/virtual MOs."""
    import pyscf
    from pyscf import fci, gto, lib, mcscf, scf, symm
    from pyscf.tools import molden

    started = time.monotonic()
    lib.num_threads(config["ranks"])
    # Match correct_cm=1 in the pinned upstream geometry generator (its masses
    # are 12.0110 and 15.9994), so the continuum sphere has the same origin.
    radius = config["bond_length"]
    carbon_z = -radius * 15.9994 / (12.0110 + 15.9994)
    oxygen_z = carbon_z + radius
    mol = gto.M(
        atom=[("C", (0, 0, carbon_z)), ("O", (0, 0, oxygen_z))],
        unit="Bohr",
        basis=config["basis"],
        symmetry="C2v",
        cart=False,
        verbose=config.get("target_verbosity", 4),
        max_memory=config["target_memory_mb"],
    )
    hf = scf.RHF(mol)
    hf.conv_tol = 1e-10
    hf.conv_tol_grad = 1e-8
    hf.chkfile = str(molden_path.with_suffix(".rhf.chk"))
    hf.kernel()
    if not hf.converged or not np.isfinite(hf.e_tot):
        raise ValueError("PySCF RHF did not converge")
    active = dict(zip(IRREPS, config["active_orbitals"], strict=True))
    ncas = sum(active.values())
    mc = mcscf.CASSCF(hf, ncas, 10, ncore=2)
    if config.get("target_optimizer", "one-step") == "newton":
        mc = mc.newton()
    mc.conv_tol = config["target_energy_tolerance"]
    mc.conv_tol_grad = config["target_gradient_tolerance"]
    mc.max_cycle_macro = config["target_max_cycles"]
    mc.ah_lindep = config["target_ah_lindep"]
    mc.ah_conv_tol = config.get("target_ah_tolerance", 1e-12)
    mc.ah_start_tol = config["target_ah_start_tolerance"]
    mc.chkfile = str(molden_path.with_suffix(".casscf.chk"))
    # A single canonicalization preserves the three subspaces. Active orbitals
    # are not reordered globally by orbital energy when imported by UKRmol.
    mc.natorb = False
    initial = mcscf.sort_mo_by_irrep(mc, hf.mo_coeff, active, {"A1": 2})
    initial_provenance = None
    if config.get("target_initial_checkpoint"):
        checkpoint = Path(config["target_initial_checkpoint"])
        previous_mol = lib.chkfile.load_mol(str(checkpoint))
        previous_mo = lib.chkfile.load(str(checkpoint), "mcscf/mo_coeff")
        previous_ncas = lib.chkfile.load(str(checkpoint), "mcscf/ncas")
        previous_ncore = lib.chkfile.load(str(checkpoint), "mcscf/ncore")
        if int(previous_ncas) != ncas or int(previous_ncore) != 2:
            raise ValueError("Initial checkpoint must have the same core/active-space sizes")
        initial = project_initial_orbitals(mc, previous_mo, previous_mol)
        initial_symmetries = symm.label_orb_symm(
            mol, mol.irrep_name, mol.symm_orb, initial[:, 2 : 2 + ncas]
        )
        if any(
            np.count_nonzero(initial_symmetries == irrep) != count
            for irrep, count in active.items()
        ):
            raise ValueError("Initial checkpoint must have the requested active irreps")
        initial_provenance = {
            "checkpoint_sha256": hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
            "source": config.get("target_initial_source_checkpoint", str(checkpoint)),
            "source_mo_columns": previous_mo.shape[1],
            "projected_mo_columns": (
                2 + ncas if previous_mo.shape[1] > hf.mo_coeff.shape[1] else previous_mo.shape[1]
            ),
        }
    solvers = []
    sectors = []
    labels = []
    for spin, counts in ((0, config["sa_singlet_roots"]), (2, config["sa_triplet_roots"])):
        for irrep, roots in zip(IRREPS, counts, strict=True):
            if not roots:
                continue
            solver = _sector_ci_driver(fci, spin, irrep, config)(mol)
            solver.spin = spin
            solver.wfnsym = irrep
            solver.nroots = roots
            solver.conv_tol = config["target_ci_tolerance"]
            solver.lindep = config.get("target_ci_lindep", 1e-14)
            if config.get("target_ci_residual_tolerance") is not None:
                solver.conv_tol_residual = config["target_ci_residual_tolerance"]
            solver.max_cycle = 200
            solver.max_space = ci_trial_space(roots, config.get("target_ci_max_space"))
            # M_s alone does not exclude higher-spin eigenstates. Penalize them
            # and verify every resulting root's S^2 before exporting anything.
            solver = fci.addons.fix_spin_(solver, shift=1.0, ss=spin * (spin + 2) / 4)
            solvers.append(solver)
            sectors.append(f"{'singlet' if spin == 0 else 'triplet'}.{irrep}")
            labels.extend((spin, irrep, root + 1) for root in range(roots))
    weights = np.full(len(labels), 1 / len(labels))
    # The mixed solver always returns a CI list; PySCF's Newton code expects
    # an array for one component. Its ordinary single-state path is the same
    # unit-weight variational objective and supplies that representation.
    if len(labels) == 1:
        mc.fcisolver = solvers[0]
    else:
        mcscf.state_average_mix_(mc, solvers, weights)
    # Construct the mixer first: it copies the first solver's instance fields.
    # An earlier instance-level kernel override would replace the mixer itself.
    if config.get("target_ci_fresh_start", False):
        for solver in solvers:
            solver.kernel = _fresh_ci_kernel(solver.kernel)
    history = {}

    def record_iteration(env: dict) -> None:
        # PySCF also calls this during microiterations, before rotation metrics
        # exist. Later microiterations are overwritten by the completed macro.
        if "de" in env and "norm_gall" in env:
            history[env["imacro"]] = {
                "macro_iteration": int(env["imacro"]),
                "energy_hartree": float(env["e_tot"]),
                "energy_change_hartree": float(env["de"]),
                "orbital_ci_gradient_norm": float(env["norm_gall"]),
                "max_rotation": float(np.max(np.abs(np.triu(env["u"], 1)))),
            }
        elif "de" in env and "max_offdiag_u" in env:
            history[env["imacro"]] = {
                "macro_iteration": int(env["imacro"]),
                "energy_hartree": float(env["e_tot"]),
                "energy_change_hartree": float(env["de"]),
                "orbital_gradient_norm": float(env["norm_gorb0"]),
                "max_rotation": float(env["max_offdiag_u"]),
            }

    mc.callback = record_iteration
    mc.kernel(initial)
    state_energies = np.atleast_1d(mc.e_tot if len(labels) == 1 else mc.e_states)
    ci_vectors = [mc.ci] if len(labels) == 1 else mc.ci
    diagnostics = {
        "orbital_converged": bool(mc.converged),
        "ci_converged": [bool(np.all(s.converged)) for s in solvers],
        "ci_max_space_by_sector": {
            sector: solver.max_space for sector, solver in zip(sectors, solvers, strict=True)
        },
        "ci_driver_by_sector": {
            sector: _sector_ci_driver(fci, solver.spin, sector.split(".")[1], config).__module__
            for sector, solver in zip(sectors, solvers, strict=True)
        },
        "state_energies_hartree": [float(e) for e in state_energies],
        "iterations": list(history.values()),
    }
    run_dir = Path(config["workdir"])
    (run_dir / "target-diagnostics.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    if (
        not mc.converged
        or not all(np.all(s.converged) for s in solvers)
        or not np.all(np.isfinite(state_energies))
    ):
        raise ValueError(
            "State-averaged CASSCF or a target CI solver did not converge; "
            f"orbital={diagnostics['orbital_converged']}, CI={diagnostics['ci_converged']}"
        )
    # Warm CI vectors can follow a higher eigenpair through an orbital change.
    # Audit lowest-root coverage before paying for UKRmol diagonalization.
    h1eff, core_energy = mc.get_h1eff()
    h2eff = mc.get_h2eff()
    fresh_energies = []
    fresh_spins = []
    for solver in solvers:
        nelec = ((10 + solver.spin) // 2, (10 - solver.spin) // 2)
        energy, vectors = solver.kernel(h1eff, h2eff, ncas, nelec, ci0=None, ecore=core_energy)
        fresh_energies.extend(float(e) for e in np.atleast_1d(energy))
        if solver.nroots == 1:
            vectors = [vectors]
        fresh_spins.extend(float(solver.spin_square(ci, ncas, nelec)[0]) for ci in vectors)
    differences = np.array(fresh_energies) - state_energies
    diagnostics.update(
        fresh_ci_converged=[bool(np.all(s.converged)) for s in solvers],
        fresh_state_energies_hartree=fresh_energies,
        fresh_spin_square=fresh_spins,
        lowest_root_max_energy_difference_hartree=float(np.max(np.abs(differences))),
    )
    (run_dir / "target-diagnostics.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    if not all(diagnostics["fresh_ci_converged"]) or not np.all(np.isfinite(fresh_energies)):
        raise ValueError("Fresh fixed-orbital CI audit did not converge")
    expected_spins = [spin * (spin + 2) / 4 for spin, _, _ in labels]
    np.testing.assert_allclose(fresh_spins, expected_spins, atol=1e-6, rtol=0)
    np.testing.assert_allclose(differences, 0, atol=1e-7, rtol=0)
    states = []
    offset = 0
    ground_density = None
    for solver in solvers:
        spin = solver.spin
        nelec = ((10 + spin) // 2, (10 - spin) // 2)
        for _ in range(solver.nroots):
            ci = ci_vectors[offset]
            ss = float(solver.spin_square(ci, ncas, nelec)[0])
            expected_ss = spin * (spin + 2) / 4
            if not np.isfinite(ss) or abs(ss - expected_ss) > 1e-6:
                raise ValueError(f"Wrong spin for {labels[offset]}: S^2={ss}")
            density = solver.make_rdm1(ci, ncas, nelec)
            if abs(np.trace(density) - 10) > 1e-8:
                raise ValueError("Active density does not contain ten electrons")
            states.append(
                {
                    "spin": "singlet" if spin == 0 else "triplet",
                    "irrep": labels[offset][1],
                    "root": labels[offset][2],
                    "energy_hartree": float(state_energies[offset]),
                    "spin_square": ss,
                    "weight": float(weights[offset]),
                }
            )
            if labels[offset] == (0, "A1", 1):
                ground_density = density
            offset += 1
    if ground_density is None:
        raise ValueError("The state average must include the ground A1 singlet")
    ground = next(
        s["energy_hartree"]
        for s in states
        if (s["spin"], s["irrep"], s["root"]) == ("singlet", "A1", 1)
    )
    if min(s["energy_hartree"] for s in states) < ground - 1e-8:
        raise ValueError("The A1 singlet is not the lowest target state")
    for state in states:
        state["excitation_hartree"] = state["energy_hartree"] - ground
    pi_error = 0.0
    for spin in ("singlet", "triplet"):
        b1 = [s["energy_hartree"] for s in states if s["spin"] == spin and s["irrep"] == "B1"]
        b2 = [s["energy_hartree"] for s in states if s["spin"] == spin and s["irrep"] == "B2"]
        np.testing.assert_allclose(b1, b2, atol=1e-7, rtol=0)
        if b1:
            pi_error = max(pi_error, float(np.max(np.abs(np.array(b1) - b2))))
    coeff = mc.mo_coeff
    overlap = hf.get_ovlp()
    active_overlap = np.linalg.svd(
        initial[:, 2 : 2 + ncas].T @ overlap @ coeff[:, 2 : 2 + ncas], compute_uv=False
    )
    orthogonality = float(np.max(np.abs(coeff.T @ overlap @ coeff - np.eye(coeff.shape[1]))))
    if orthogonality > 1e-9:
        raise ValueError("Exported MOs are not orthonormal in the AO metric")
    core = coeff[:, :2]
    cas = coeff[:, 2 : 2 + ncas]
    dm = 2 * core @ core.T + cas @ ground_density @ cas.T
    dipole = hf.dip_moment(mol, dm, unit="AU")
    orbsym = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, coeff)
    # The block order is the import contract, and is recorded explicitly.
    occupations = np.zeros(coeff.shape[1])
    occupations[:2] = 2
    occupations[2 : 2 + ncas] = np.diag(mc.fcisolver.make_rdm1(mc.ci, ncas, 10))
    molden.from_mo(
        mol,
        str(molden_path),
        coeff,
        symm=orbsym,
        ene=mc.mo_energy,
        occ=occupations,
    )
    # GBTOlib's shell scan treats lowercase d/f/g in cards as shell headers.
    # Its spherical-basis cards use uppercase, unlike PySCF's Molden writer.
    exported = molden_path.read_text()
    for lower, upper in (("[5d]", "[5D]"), ("[7f]", "[7F]"), ("[9g]", "[9G]")):
        exported = exported.replace(lower, upper)
    molden_path.write_text(exported)
    inventory = []
    counts = dict.fromkeys(IRREPS, 0)
    for index, (irrep, energy) in enumerate(zip(orbsym, mc.mo_energy, strict=True)):
        counts[irrep] += 1
        inventory.append(
            {
                "irrep": str(irrep),
                "index_in_irrep": counts[irrep],
                "energy_hartree": float(energy),
                "space": "core" if index < 2 else "active" if index < 2 + ncas else "virtual",
            }
        )
    return {
        "pyscf_version": pyscf.__version__,
        "atoms_bohr": [["C", 0, 0, carbon_z], ["O", 0, 0, oxygen_z]],
        "method": "equal-weight multi-spin/multi-irrep SA-CASSCF",
        "initial_orbitals": initial_provenance,
        "optimization": diagnostics,
        "converged": True,
        "rhf_energy_hartree": float(hf.e_tot),
        "ensemble_energy_hartree": float(mc.e_tot),
        "ground_energy_hartree": ground,
        "ground_dipole_au": [float(x) for x in dipole],
        "pi_energy_max_difference_hartree": pi_error,
        "mo_orthogonality_max_error": orthogonality,
        "initial_to_final_active_overlap_singular_values": [float(x) for x in active_overlap],
        "states": states,
        "orbitals": inventory,
        "molden_sha256": hashlib.sha256(molden_path.read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "wall_seconds": time.monotonic() - started,
    }


def main() -> None:
    """Read the run configuration and retain both successful and failed QC logs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    # PySCF writes through its own stdout object as well as ordinary print.
    with args.output.open("w", buffering=1) as log, contextlib.redirect_stdout(log):
        result = build_target(json.loads(args.config.read_text()), Path("co.molden"))
        args.config.with_name("target.json").write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
