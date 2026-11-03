"""One four-atom ROHF/UCCSD(T) endpoint pilot; no molecular or scattering successor."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path

CONTRACT = "projects/ukrmol_co/atomic-asymptote-release-20261103.json"
P_IRREPS = ("B1u", "B2u", "B3u")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    temporary = path.with_name(path.name + ".new")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def occupations(state):
    """Freeze one orientation of the high-spin 2s2/2p^n atomic configuration."""
    if (
        4 + sum(state["p_occupations"]) != state["electrons"]
        or sum(n == 1 for n in state["p_occupations"]) != state["spin_2s"]
        or any(n not in (0, 1, 2) for n in state["p_occupations"])
    ):
        raise ValueError("Inconsistent atomic valence/electron/spin specification")
    result = {name: (0, 0) for name in ("Ag", "B1g", "B2g", "B3g", "Au", *P_IRREPS)}
    result["Ag"] = (2, 2)
    result.update(
        {
            ir: (int(n > 0), int(n == 2))
            for ir, n in zip(P_IRREPS, state["p_occupations"], strict=True)
        }
    )
    return result


def authenticate(source, contract):
    """Reject altered or expanded inputs before importing a quantum implementation."""
    for name, pin in contract["source_sha256"].items():
        if digest(source / name) != pin:
            raise ValueError("Changed atomic prerequisite " + name)
    inputs = json.loads((source / contract["inputs"]).read_text())
    basis = contract.get("basis", "aug-cc-pVQZ")
    if (
        contract.get("execution_authorized") is not True
        or contract["maximum_atoms"] != 4
        or basis not in ("aug-cc-pVQZ", "aug-cc-pV5Z")
        or inputs["basis"] != basis
        or [s["id"] for s in inputs["states"]] != ["C", "O", "O-", "C-"]
        or inputs["frozen_spatial_orbitals_per_atom"] != 1
        or contract["automatic_second_basis"]
        or contract["full_sweep_released"]
    ):
        raise ValueError("Only the reviewed four-atom QZ or separately released5Z is supported")
    if basis == "aug-cc-pV5Z":
        previous = json.loads((source / contract["predecessor_result"]).read_text())
        original = json.loads((source / contract["predecessor_inputs"]).read_text())
        proof = json.loads((source / contract["predecessor_public_proof"]).read_text())
        if not (
            previous.get("finite_method_endpoint_reconstruction_verdict") is True
            and previous.get("original_exits") == [0] * 8
            and previous.get("original_atom_profile_exits") == [0] * 4
            and proof.get("finite_method_endpoint_reconstruction_verdict") is True
            and proof.get("fresh_public_payloads_verified") == 5
            and proof.get("archive_sha256") == previous["archive_sha256"]
            and inputs["states"] == original["states"]
            and original["basis"] == "aug-cc-pVQZ"
        ):
            raise ValueError("Independently qualified QZ predecessor and identical states required")
    for state in inputs["states"]:
        z = {"C": 6, "O": 8}[state["element"]]
        if z - state["charge"] != state["electrons"]:
            raise ValueError("Atomic charge/electron inventory differs")
        occupations(state)
    return inputs


def endpoint_sums(totals):
    """Exact ledger identities; external EA agreement is a separate physical check."""
    if set(totals) != {"C", "O", "O-", "C-"} or not all(math.isfinite(x) for x in totals.values()):
        raise ValueError("All four finite atomic totals are required")
    return dict(
        fragment_sum_hartree={
            "C+O": totals["C"] + totals["O"],
            "C+O-": totals["C"] + totals["O-"],
            "C-+O": totals["C-"] + totals["O"],
        },
        electron_affinity_hartree={
            "O": totals["O"] - totals["O-"],
            "C": totals["C"] - totals["C-"],
        },
    )


def calculate(state, inputs, resources, work):
    """Use explicit ROHF orbitals without a subsequent unrestricted SCF relaxation."""
    import numpy as np
    import pyscf
    from pyscf import cc, gto, lib, scf, symm

    lib.num_threads(resources["physical_cores"])
    report = dict(
        id=state["id"],
        term_target=state["term"],
        basis=inputs["basis"],
        pyscf_version=pyscf.__version__,
        method="ROHF-reference conventional UCCSD(T)",
        frozen_spatial_orbitals=1,
        external_rohf_stability=None,
        external_stability_reason="Not available in pinned PySCF2.11",
        stages=[],
    )

    def stage(name, fn):
        started = time.monotonic()
        item = dict(name=name, original_exit=1)
        report["stages"].append(item)
        try:
            value = fn()
            item["original_exit"] = 0
            return value
        finally:
            item["wall_seconds"] = time.monotonic() - started
            save(work / "diagnostics.json", report)

    try:
        mol = gto.M(
            atom=[(state["element"], (0, 0, 0))],
            unit="Bohr",
            basis=inputs["basis"],
            charge=state["charge"],
            spin=state["spin_2s"],
            symmetry="D2h",
            cart=False,
            verbose=4,
            max_memory=resources["pyscf_memory_mb"],
        )
        if pyscf.__version__ != "2.11.0" or mol.nelectron != state["electrons"]:
            raise ValueError("Unreviewed PySCF runtime or electron count")
        hf = scf.ROHF(mol)
        hf.irrep_nelec = occupations(state)
        hf.conv_tol = 1e-11
        hf.conv_tol_grad = 1e-8
        hf.max_cycle = 100
        hf.chkfile = str(work / "rohf.chk")
        stage("rohf", hf.kernel)
        report.update(
            rohf_converged=bool(hf.converged),
            rohf_energy_hartree=float(hf.e_tot),
            reference_s2=float(hf.spin_square()[0]),
            nao=mol.nao_nr(),
        )
        if not hf.converged:
            raise ValueError("Atomic ROHF did not converge")
        spin = state["spin_2s"] / 2
        np.testing.assert_allclose(report["reference_s2"], spin * (spin + 1), atol=1e-8, rtol=0)
        labels = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, hf.mo_coeff)
        measured = {
            ir: (
                int(np.count_nonzero((labels == ir) & (hf.mo_occ > 0))),
                int(np.count_nonzero((labels == ir) & (hf.mo_occ == 2))),
            )
            for ir in mol.irrep_name
        }
        if measured != occupations(state):
            raise ValueError("Final atomic occupation changed")
        report["irrep_electrons_alpha_beta"] = measured
        overlap = hf.get_ovlp()
        vals, vec = np.linalg.eigh(overlap)
        if vals.min() <= 1e-9:
            raise ValueError("Atomic AO overlap rank is unresolved")
        root = (vec * np.sqrt(vals)) @ vec.T
        ao_l = [label[2][-1] for label in mol.ao_labels(fmt=False)]
        angular = []
        for i in np.flatnonzero(hf.mo_occ > 0):
            weights = (root @ hf.mo_coeff[:, i]) ** 2
            target = "s" if labels[i] == "Ag" else "p"
            fraction = float(
                sum(w for w, kind in zip(weights, ao_l, strict=True) if kind == target)
            )
            angular.append(
                dict(
                    orbital=int(i),
                    irrep=str(labels[i]),
                    occupation=float(hf.mo_occ[i]),
                    dominant_angular_target=target,
                    target_weight=fraction,
                )
            )
            if fraction <= 0.5:
                raise ValueError("Occupied orbital is not majority intended s/p character")
        report.update(
            occupied_angular_character=angular, minimum_overlap_eigenvalue=float(vals.min())
        )
        if labels[0] != "Ag" or hf.mo_occ[0] != 2 or hf.mo_energy[0] >= -3:
            raise ValueError("Frozen orbital is not the deep occupied s core")
        _, _, internal, _ = stage(
            "rohf_internal_stability",
            lambda: hf.stability(internal=True, external=False, return_status=True),
        )
        report["rohf_internally_stable"] = bool(internal)
        if not internal:
            raise ValueError("Atomic ROHF internally unstable;no orbital retry released")
        uhf = hf.to_uhf()
        np.testing.assert_allclose(uhf.e_tot, hf.e_tot, atol=1e-12, rtol=0)
        np.testing.assert_allclose(uhf.mo_coeff, np.array([hf.mo_coeff] * 2), atol=1e-12, rtol=0)
        coupled = cc.UCCSD(uhf, frozen=1)
        coupled.conv_tol = 1e-10
        coupled.conv_tol_normt = 1e-8
        coupled.max_cycle = 100
        coupled.chkfile = str(work / "uccsd.chk")
        stage("uccsd", coupled.kernel)
        report.update(
            uccsd_converged=bool(coupled.converged),
            uccsd_energy_hartree=float(coupled.e_tot),
            t1_norms=[float(np.linalg.norm(t)) for t in coupled.t1],
            reference_conversion_class=uhf.__class__.__module__ + "." + uhf.__class__.__name__,
        )
        if not coupled.converged:
            raise ValueError("Atomic UCCSD did not converge")
        coupled.dump_chk()
        triples = float(stage("triples", coupled.ccsd_t))
        total = float(coupled.e_tot + triples)
        if not all(math.isfinite(x) for x in [total, triples, *report["t1_norms"]]):
            raise ValueError("Nonfinite atomic energy/amplitude")
        report.update(
            triples_correction_hartree=triples,
            total_energy_hartree=total,
            finite_method_energy_complete=True,
            correlated_L2_purity_proved=False,
            absolute_accuracy_qualified=False,
            molecular_dissociation_qualified=False,
        )
        save(work / "atom.json", report)
        return report
    finally:
        save(work / "diagnostics.json", report)


def run_queue(source, output, contract, contract_path=CONTRACT):
    """Profile sequential atoms; first failed original child exit stops without fallback."""
    from projects.ukrmol_co.run import run_profiled

    inputs = authenticate(source, contract)
    started = time.monotonic()
    completed = {}
    report = dict(
        scientific_exit=2,
        atoms=[],
        new_scf_calls_reserved=0,
        automatic_successor=False,
        full_sweep_released=False,
        route_holdout_s_matrices_analyzed=0,
    )
    try:
        for state in inputs["states"]:
            remaining = contract["resources"]["total_science_timeout_seconds"] - (
                time.monotonic() - started
            )
            if remaining <= 0:
                raise TimeoutError("Total atomic science cap exhausted")
            work = output / state["id"]
            work.mkdir()
            (work / "scratch").mkdir()
            save(work / "state.json", state)
            report["atoms"].append(
                dict(id=state["id"], state="running", original_profile_exit=None)
            )
            report["new_scf_calls_reserved"] += 1
            save(output / "result.json", report)
            env = os.environ | {"TMPDIR": str(work / "scratch")}
            command = [
                "timeout",
                "--signal=TERM",
                "--kill-after=30s",
                str(min(contract["resources"]["per_atom_timeout_seconds"], int(remaining))) + "s",
                sys.executable,
                "-m",
                "projects.ukrmol_co.atomic_references",
                "--contract",
                str(source / contract_path),
                "--output",
                str(output),
                "--atom",
                state["id"],
            ]
            exit_code = run_profiled(work, env, command)
            report["atoms"][-1].update(state="finished", original_profile_exit=exit_code)
            save(output / "result.json", report)
            if exit_code:
                raise ValueError("Original atomic worker failed;no peer/retry released")
            completed[state["id"]] = json.loads((work / "atom.json").read_text())[
                "total_energy_hartree"
            ]
        from qscat.units import HARTREE_TO_EV

        report.update(
            endpoint_sums(completed),
            total_energy_hartree=completed,
            scientific_exit=0,
            finite_method_energies_complete=True,
            atomic_absolute_accuracy_qualified=False,
            molecular_dissociation_qualified=False,
        )
        report["electron_affinity_ev"] = {
            k: v * HARTREE_TO_EV for k, v in report["electron_affinity_hartree"].items()
        }
    except Exception as error:
        report["error"] = repr(error)
    finally:
        try:
            authenticate(source, contract)
        except Exception as error:
            report.update(scientific_exit=2, integrity_error=repr(error))
        report["wall_seconds"] = time.monotonic() - started
        save(output / "result.json", report)
    return report["scientific_exit"]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--contract", type=Path, default=Path(CONTRACT))
    p.add_argument("--output", type=Path, default=Path("/out"))
    p.add_argument("--atom")
    args = p.parse_args()
    source = Path(__file__).resolve().parents[2]
    c = json.loads(args.contract.read_text())
    inputs = authenticate(source, c)
    if args.atom is None:
        raise SystemExit(run_queue(source, args.output, c, args.contract))
    state = next(s for s in inputs["states"] if s["id"] == args.atom)
    calculate(state, inputs, c["resources"], args.output / state["id"])
