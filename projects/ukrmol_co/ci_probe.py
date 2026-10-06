"""Check lowest-root coverage at fixed CASSCF orbitals with PySCF CASCI.

Converged CI eigenpairs can omit a lower root. Vary requested roots and Davidson
space without changing the orbital subspaces, then compare with independent
UKRmol energies. PySCF is supplied by the optional CO Docker layer.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import shutil
import time
from pathlib import Path

import numpy as np


def probe(checkpoint: Path, spin: int, irrep: str, roots: int, space: int, ranks: int) -> dict:
    """Solve one fixed-orbital CAS sector; report every eigenpair and its spin."""
    from pyscf import fci, lib, mcscf, scf

    started = time.monotonic()
    lib.num_threads(ranks)
    mol = lib.chkfile.load_mol(str(checkpoint))
    mol.verbose = 4
    coeff = lib.chkfile.load(str(checkpoint), "mcscf/mo_coeff")
    ncas = int(lib.chkfile.load(str(checkpoint), "mcscf/ncas"))
    ncore = int(lib.chkfile.load(str(checkpoint), "mcscf/ncore"))
    if ncore != 2:
        raise ValueError("CO CI probes require two core orbitals and ten active electrons")
    nelec = ((10 + spin) // 2, (10 - spin) // 2)
    casci = mcscf.CASCI(scf.RHF(mol), ncas, nelec, ncore=ncore)
    casci.canonicalization = False
    solver = fci.direct_spin1_symm.FCI(mol)
    solver.spin = spin
    solver.wfnsym = irrep
    solver.nroots = roots
    solver.max_space = space
    solver.max_cycle = 200
    solver.conv_tol = 1e-12
    solver.conv_tol_residual = 1e-9
    solver.lindep = 1e-22
    expected_spin = spin * (spin + 2) / 4
    casci.fcisolver = fci.addons.fix_spin_(solver, shift=1.0, ss=expected_spin)
    casci.kernel(coeff)
    vectors = [casci.ci] if roots == 1 else casci.ci
    energies = np.atleast_1d(casci.e_tot)
    spins = [float(casci.fcisolver.spin_square(ci, ncas, nelec)[0]) for ci in vectors]
    return {
        "roots": roots,
        "max_space": space,
        "ci_converged": np.atleast_1d(casci.fcisolver.converged).tolist(),
        "finite_energies": bool(np.all(np.isfinite(energies))),
        "energy_hartree": [float(e) for e in energies],
        "spin_square": spins,
        "correct_spin": bool(np.all(np.abs(np.array(spins) - expected_spin) <= 1e-6)),
        "wall_seconds": time.monotonic() - started,
    }


def main() -> None:
    """Preserve a fresh checkpoint-based CI coverage scan, including failed probes."""
    import pyscf

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--spin", type=int, choices=[0, 2], default=2)
    parser.add_argument("--irrep", choices=["A1", "B1", "B2", "A2"], default="A1")
    parser.add_argument("--roots", type=int, nargs="+", default=[5, 8])
    parser.add_argument("--spaces", type=int, nargs="+", default=[40, 80, 160])
    parser.add_argument("--ranks", type=int, default=4)
    args = parser.parse_args()
    if min(args.roots + args.spaces + [args.ranks]) < 1 or min(args.spaces) <= max(args.roots):
        parser.error("Roots/ranks must be positive; every trial space must exceed all root counts")
    args.output.mkdir(parents=True, exist_ok=False)
    checkpoint = args.output / "checkpoint.chk"
    shutil.copyfile(args.checkpoint, checkpoint)
    records = []
    report = {
        "scope": "Fixed-orbital CI coverage diagnostic; no orbital or electronic qualification",
        "pyscf_version": pyscf.__version__,
        "checkpoint_source": str(args.checkpoint),
        "checkpoint_sha256": hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "spin": args.spin,
        "irrep": args.irrep,
        "ranks": args.ranks,
        "requested_roots": args.roots,
        "requested_spaces": args.spaces,
        "ci_max_cycles": 200,
        "ci_energy_tolerance": 1e-12,
        "ci_residual_tolerance": 1e-9,
        "ci_lindep": 1e-22,
        "spin_penalty_shift": 1.0,
        "probes": records,
    }
    for roots in args.roots:
        for space in args.spaces:
            with (args.output / f"roots{roots}-space{space}.log").open("w") as log:
                with contextlib.redirect_stdout(log):
                    try:
                        record = probe(checkpoint, args.spin, args.irrep, roots, space, args.ranks)
                    except Exception as error:
                        record = {"roots": roots, "max_space": space, "error": repr(error)}
                records.append(record)
            (args.output / "result.json").write_text(json.dumps(report, indent=2) + "\n")
            print(json.dumps(record), flush=True)
    raise SystemExit(
        int(
            any(
                "error" in record
                or not all(record["ci_converged"])
                or not record["finite_energies"]
                or not record["correct_spin"]
                for record in records
            )
        )
    )


if __name__ == "__main__":
    main()
