"""Saved-state oxygen triples control: same determinant and CCSD state, new MO gauge."""

from __future__ import annotations

import argparse
import gc
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.atomic_references import digest, save

CONTRACT = "projects/ukrmol_co/oxygen-affinity-orbital-release-20261103.json"


def off_diagonal(matrix):
    return float(np.max(np.abs(matrix - np.diag(np.diag(matrix))), initial=0))


def block_rotation(fock, nocc, labels):
    """Diagonalize real occupied/virtual Fock blocks without mixing spatial irreps."""
    fock = np.asarray(fock)
    labels = np.asarray(labels)
    if fock.shape != (len(labels), len(labels)) or not 0 < nocc < len(labels):
        raise ValueError("Fock/occupation/symmetry dimensions disagree")
    np.testing.assert_allclose(fock, fock.T, atol=1e-10, rtol=0)
    rotation = np.eye(len(labels))
    for lo, hi in [(0, nocc), (nocc, len(labels))]:
        for label in np.unique(labels[lo:hi]):
            indices = np.flatnonzero(labels[lo:hi] == label) + lo
            _, vectors = np.linalg.eigh(fock[np.ix_(indices, indices)])
            rotation[np.ix_(indices, indices)] = vectors
    return rotation


def rotate_amplitudes(t1, t2, rotations, nocc):
    """Transform real UCCSD amplitudes for C_new=C_old U in the active MO spaces."""
    oa, ob = [u[:n, :n] for u, n in zip(rotations, nocc, strict=True)]
    va, vb = [u[n:, n:] for u, n in zip(rotations, nocc, strict=True)]
    singles = (oa.T @ t1[0] @ va, ob.T @ t1[1] @ vb)
    doubles = tuple(
        np.ascontiguousarray(np.einsum("ijab,iI,jJ,aA,bB->IJAB", t, oi, oj, vi, vj, optimize=True))
        for t, oi, oj, vi, vj in [
            (t2[0], oa, oa, va, va),
            (t2[1], oa, ob, va, vb),
            (t2[2], ob, ob, vb, vb),
        ]
    )
    return singles, doubles


def authenticate(root, previous, contract):
    for name, pin in contract["source_sha256"].items():
        if digest(root / name) != pin:
            raise ValueError("Changed diagnostic source " + name)
    for name, pin in contract["input_sha256"].items():
        if digest(previous / name) != pin:
            raise ValueError("Changed original oxygen evidence " + name)
    if not (
        contract["execution_authorized"]
        and contract["states"] == ["O", "O-"]
        and contract["new_scf_calls"] == contract["new_ccsd_iterations"] == 0
    ):
        raise ValueError("Only two saved-state oxygen controls are released")


def diagnose(previous, output, species):
    """Replay saved CCSD/triples, then semicanonicalize without reference relaxation."""
    import pyscf
    from pyscf import cc, lib, scf, symm

    if pyscf.__version__ != "2.11.0":
        raise ValueError("Pinned PySCF2.11 required")
    lib.num_threads(4)
    old = previous / "outputs" / species
    original = json.loads((old / "atom.json").read_text())
    mol, state = scf.chkfile.load_scf(str(old / "rohf.chk"))
    mol.max_memory = 12000
    mol.verbose = 4
    hf = scf.ROHF(mol)
    hf.__dict__.update(state)
    hf.converged = True
    np.testing.assert_allclose(hf.e_tot, original["rohf_energy_hartree"], atol=1e-11, rtol=0)
    if mol.spin != (2 if species == "O" else 1) or mol.nelectron != (8 if species == "O" else 9):
        raise ValueError("Original oxygen state changed")
    mf = hf.to_uhf()
    coupled = cc.UCCSD(mf, frozen=1)
    stored = lib.chkfile.load(str(old / "uccsd.chk"), "ccsd")
    coupled.t1 = tuple(stored["t1"])
    coupled.t2 = tuple(stored["t2"])
    coupled.e_hf = float(hf.e_tot)
    eris = coupled.ao2mo()
    e_corr = float(coupled.energy(coupled.t1, coupled.t2, eris))
    np.testing.assert_allclose(e_corr, stored["e_corr"], atol=1e-8, rtol=0)
    np.testing.assert_allclose(
        hf.e_tot + e_corr, original["uccsd_energy_hartree"], atol=1e-8, rtol=0
    )
    triples = float(coupled.ccsd_t(eris=eris))
    np.testing.assert_allclose(triples, original["triples_correction_hartree"], atol=1e-9, rtol=0)
    labels = symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, hf.mo_coeff)
    masks = coupled.get_frozen_mask()
    active_labels = [labels[m] for m in masks]
    fock = (eris.focka.copy(), eris.fockb.copy())
    nocc = coupled.nocc
    rotations = tuple(
        block_rotation(f, n, label) for f, n, label in zip(fock, nocc, active_labels, strict=True)
    )
    before = [
        max(off_diagonal(f[:n, :n]), off_diagonal(f[n:, n:]))
        for f, n in zip(fock, nocc, strict=True)
    ]
    transformed = [u.T @ f @ u for u, f in zip(rotations, fock, strict=True)]
    after = [
        max(off_diagonal(f[:n, :n]), off_diagonal(f[n:, n:]))
        for f, n in zip(transformed, nocc, strict=True)
    ]
    if max(after) > 1e-9:
        raise ValueError("Semicanonical Fock blocks do not diagonalize")
    full_rotations = []
    for mask, u in zip(masks, rotations, strict=True):
        full = np.eye(len(mask))
        indices = np.flatnonzero(mask)
        full[np.ix_(indices, indices)] = u
        np.testing.assert_allclose(full.T @ full, np.eye(len(mask)), atol=1e-10, rtol=0)
        np.testing.assert_allclose(full[0], np.eye(len(mask))[0], atol=1e-10, rtol=0)
        full_rotations.append(full)
    new_mf = mf.copy()
    new_mf.mo_coeff = np.asarray([c @ u for c, u in zip(mf.mo_coeff, full_rotations, strict=True)])
    np.testing.assert_allclose(new_mf.make_rdm1(), mf.make_rdm1(), atol=1e-10, rtol=0)
    new_t1, new_t2 = rotate_amplitudes(coupled.t1, coupled.t2, rotations, nocc)
    np.savez_compressed(
        output / "orbital-witnesses.npz",
        focka=fock[0],
        fockb=fock[1],
        rotationa=rotations[0],
        rotationb=rotations[1],
        nocc=nocc,
        density_old=mf.make_rdm1(),
        density_new=new_mf.make_rdm1(),
        t1a_old=coupled.t1[0],
        t1b_old=coupled.t1[1],
        t1a_new=new_t1[0],
        t1b_new=new_t1[1],
    )
    del eris, coupled
    gc.collect()
    new_cc = cc.UCCSD(new_mf, frozen=1)
    new_cc.e_hf = float(hf.e_tot)
    new_cc.t1, new_cc.t2 = new_t1, new_t2
    new_eris = new_cc.ao2mo()
    for actual, expected in zip((new_eris.focka, new_eris.fockb), transformed, strict=True):
        np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=0)
    rotated_e_corr = float(new_cc.energy(new_t1, new_t2, new_eris))
    np.testing.assert_allclose(rotated_e_corr, e_corr, atol=1e-8, rtol=0)
    new_triples = float(new_cc.ccsd_t(eris=new_eris))
    report = dict(
        species=species,
        basis=original["basis"],
        new_scf_calls=0,
        new_ccsd_iterations=0,
        reference_energy_hartree=float(hf.e_tot),
        original_ccsd_correlation_hartree=e_corr,
        rotated_ccsd_correlation_hartree=rotated_e_corr,
        original_triples_hartree=triples,
        semicanonical_triples_hartree=new_triples,
        original_total_hartree=float(hf.e_tot + e_corr + triples),
        semicanonical_total_hartree=float(hf.e_tot + rotated_e_corr + new_triples),
        fock_off_diagonal_before_hartree=before,
        fock_off_diagonal_after_hartree=after,
        density_and_core_invariance_passed=True,
        ccsd_energy_invariance_passed=True,
        absolute_accuracy_qualified=False,
    )
    save(output / "result.json", report)
    print(json.dumps(report, indent=2), flush=True)


def queue(previous, output, contract):
    from qscat.units import HARTREE_TO_EV

    from projects.ukrmol_co.run import run_profiled

    result = dict(
        original_science_exit=2,
        states=[],
        new_scf_calls=0,
        new_ccsd_iterations=0,
        automatic_successor=False,
    )
    start = time.monotonic()
    try:
        for species in contract["states"]:
            work = output / species
            work.mkdir()
            (work / "scratch").mkdir()
            command = [
                "timeout",
                "--signal=TERM",
                "--kill-after=30s",
                "900s",
                sys.executable,
                "-m",
                "projects.ukrmol_co.oxygen_affinity_orbitals",
                "--previous",
                str(previous),
                "--output",
                str(work),
                "--species",
                species,
            ]
            code = run_profiled(work, os.environ | {"TMPDIR": str(work / "scratch")}, command)
            result["states"].append(dict(species=species, original_profile_exit=code))
            save(output / "result.json", result)
            if code:
                raise ValueError("Original saved-state diagnostic failed;no peer/retry")
            if time.monotonic() - start > 1800:
                raise TimeoutError("Total science cap exceeded")
        atoms = {
            s: json.loads((output / s / "result.json").read_text()) for s in contract["states"]
        }
        for kind in ["original", "semicanonical"]:
            result[kind + "_electron_affinity_ev"] = (
                atoms["O"][kind + "_total_hartree"] - atoms["O-"][kind + "_total_hartree"]
            ) * HARTREE_TO_EV
        result["orbital_convention_ea_shift_ev"] = (
            result["semicanonical_electron_affinity_ev"] - result["original_electron_affinity_ev"]
        )
        result["original_science_exit"] = 0
    except Exception as error:
        result["error"] = repr(error)
    finally:
        result["wall_seconds"] = time.monotonic() - start
        save(output / "result.json", result)
    return result["original_science_exit"]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--previous", type=Path, default=Path("/previous"))
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--species", choices=["O", "O-"])
    p.add_argument("--queue", action="store_true")
    a = p.parse_args()
    root = Path(__file__).resolve().parents[2]
    contract = json.loads((root / CONTRACT).read_text())
    authenticate(root, a.previous, contract)
    try:
        if a.queue:
            code = queue(a.previous, a.output, contract)
        elif a.species is not None:
            diagnose(a.previous, a.output, a.species)
            code = 0
        else:
            p.error("Specify one species or the frozen queue")
    finally:
        authenticate(root, a.previous, contract)
    raise SystemExit(code)
