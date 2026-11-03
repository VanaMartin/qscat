"""One fixed-QZ diffuse-function control, reusing the original atomic algorithm."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.atomic_references import calculate, digest, occupations, save

CONTRACT = "projects/ukrmol_co/oxygen-affinity-diffuse-release-20261103.json"
INPUTS = "projects/ukrmol_co/oxygen-diffuse-inputs-20261103.json"


def authenticate(root, previous, contract):
    for name, pin in contract["source_sha256"].items():
        if digest(root / name) != pin:
            raise ValueError("Changed frozen source " + name)
    for name, pin in contract["input_sha256"].items():
        if digest(previous / name) != pin:
            raise ValueError("Changed original QZ input " + name)
    inputs = json.loads((root / INPUTS).read_text())
    if (
        not contract["execution_authorized"]
        or contract["states"] != ["O", "O-"]
        or inputs["parent_basis"] != "aug-cc-pVQZ"
        or contract["maximum_new_scf_calls"] != 2
        or contract["automatic_successor"]
    ):
        raise ValueError("Only one two-state fixed-QZ diffuse control is released")
    for state in inputs["states"]:
        occupations(state)
    return inputs


def basis_control(inputs):
    """Prove the parent AO space is unchanged and exactly 25 functions are added."""
    from pyscf import gto

    base = gto.basis.load(inputs["parent_basis"], "O")
    basis = base + inputs["second_diffuse_shells"]
    parent = gto.M(atom="O 0 0 0", basis={"O": base}, spin=2, verbose=0)
    extended = gto.M(atom="O 0 0 0", basis={"O": basis}, spin=2, verbose=0)
    if parent.nao_nr() != 80 or extended.nao_nr() != inputs["expected_spherical_ao_count"]:
        raise ValueError("Unexpected AO count")
    # PySCF sorts shells by angular momentum: map original shells, not AO positions.
    indices = []
    offsets = extended.ao_loc_nr()
    used = set()
    for shell in range(parent.nbas):
        for other in range(extended.nbas):
            if other not in used and parent.bas_angular(shell) == extended.bas_angular(other):
                if np.array_equal(
                    parent.bas_exp(shell), extended.bas_exp(other)
                ) and np.array_equal(parent.bas_ctr_coeff(shell), extended.bas_ctr_coeff(other)):
                    indices.extend(range(offsets[other], offsets[other + 1]))
                    used.add(other)
                    break
        else:
            raise ValueError("An original parent shell was changed")
    overlap = extended.intor("int1e_ovlp")
    np.testing.assert_allclose(
        overlap[np.ix_(indices, indices)], parent.intor("int1e_ovlp"), atol=1e-13, rtol=0
    )
    if np.linalg.eigvalsh(overlap).min() <= 1e-9:
        raise ValueError("Diffuse overlap rank unresolved")
    return basis


def queue(root, previous, output, contract, inputs):
    from qscat.units import HARTREE_TO_EV

    from projects.ukrmol_co.run import run_profiled

    result = dict(original_science_exit=2, states=[], automatic_successor=False)
    start = time.monotonic()
    try:
        for state in inputs["states"]:
            species = state["id"]
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
                "projects.ukrmol_co.oxygen_affinity_diffuse",
                "--output",
                str(work),
                "--species",
                species,
            ]
            code = run_profiled(work, os.environ | {"TMPDIR": str(work / "scratch")}, command)
            result["states"].append(dict(species=species, original_profile_exit=code))
            save(output / "result.json", result)
            if code or time.monotonic() - start > 1800:
                raise ValueError("Original worker failed/time cap; no peer/retry")
        totals = {
            s: json.loads((output / s / "atom.json").read_text())["total_energy_hartree"]
            for s in contract["states"]
        }
        old = {
            s: json.loads((previous / "outputs" / s / "atom.json").read_text())[
                "total_energy_hartree"
            ]
            for s in contract["states"]
        }
        result.update(
            total_energy_hartree=totals,
            extended_electron_affinity_ev=(totals["O"] - totals["O-"]) * HARTREE_TO_EV,
            original_qz_electron_affinity_ev=(old["O"] - old["O-"]) * HARTREE_TO_EV,
        )
        result["diffuse_shift_ev"] = (
            result["extended_electron_affinity_ev"] - result["original_qz_electron_affinity_ev"]
        )
        result["original_science_exit"] = 0
    except Exception as error:
        result["error"] = repr(error)
    finally:
        authenticate(root, previous, contract)
        result["wall_seconds"] = time.monotonic() - start
        save(output / "result.json", result)
    return result["original_science_exit"]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--species", choices=["O", "O-"])
    p.add_argument("--queue", action="store_true")
    a = p.parse_args()
    root = Path(__file__).resolve().parents[2]
    previous = Path("/previous")
    contract = json.loads((root / CONTRACT).read_text())
    inputs = authenticate(root, previous, contract)
    if a.queue:
        raise SystemExit(queue(root, previous, a.output, contract, inputs))
    if a.species is None:
        p.error("Choose the queue or one frozen species")
    basis = basis_control(inputs)
    try:
        state = next(s for s in inputs["states"] if s["id"] == a.species)
        calculate(state, {"basis": {"O": basis}}, contract["resources"], a.output)
    finally:
        authenticate(root, previous, contract)
