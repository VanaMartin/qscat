"""Profile an independently correlated CO neutral-curve pilot in the CO image.

Energies are frozen-core conventional RHF/CCSD(T); the density dipole is CCSD.
PySCF and the optional independent Psi4 control are supplied by Docker.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.run import run_profiled


def calculate(config: dict) -> dict:
    """Retain energies, reference stability and lambda-density diagnostics."""
    import pyscf
    from pyscf import cc, gto, lib, scf

    workdir = Path(config["workdir"])
    lib.num_threads(config["ranks"])
    radius = config["bond_length"]
    carbon_z = -radius * 15.9994 / (12.0110 + 15.9994)
    atoms = [["C", 0, 0, carbon_z], ["O", 0, 0, carbon_z + radius]]
    report = {
        "method": "frozen-core conventional RHF/CCSD(T)",
        "pyscf_version": pyscf.__version__,
        "basis": config["basis"],
        "bond_length_bohr": radius,
        "atoms_bohr": atoms,
        "frozen_spatial_orbitals": 2,
        "correlated_electrons": 10,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }

    def timed(program, operation):
        started = time.monotonic()
        status = 1
        try:
            result = operation()
            status = 0
            return result
        finally:
            with (workdir / "stages.tsv").open("a") as log:
                log.write(f"neutral\t{program}\t\t\t{time.monotonic() - started}\t{status}\n")

    try:
        mol = gto.M(
            atom=[(symbol, xyz) for symbol, *xyz in atoms],
            unit="Bohr",
            basis=config["basis"],
            symmetry="C2v",
            cart=False,
            verbose=4,
            max_memory=config["memory_mb"],
        )
        hf = scf.RHF(mol)
        hf.conv_tol = 1e-11
        hf.conv_tol_grad = 1e-8
        hf.chkfile = str(workdir / "neutral.rhf.chk")
        timed("rhf", hf.kernel)
        report.update(rhf_converged=bool(hf.converged), rhf_energy_hartree=float(hf.e_tot))
        if not hf.converged:
            raise ValueError("Neutral RHF did not converge")
        _, _, internal, external = timed(
            "rhf_stability", lambda: hf.stability(external=True, return_status=True)
        )
        report.update(rhf_internally_stable=bool(internal), rhf_externally_stable=bool(external))
        if not internal or not external:
            raise ValueError("Neutral RHF reference is unstable")
        coupled = cc.CCSD(hf, frozen=2)
        coupled.conv_tol = 1e-10
        coupled.conv_tol_normt = 1e-8
        coupled.max_cycle = 100
        coupled.chkfile = str(workdir / "neutral.ccsd.chk")
        timed("ccsd", coupled.kernel)
        report.update(
            ccsd_converged=bool(coupled.converged),
            ccsd_energy_hartree=float(coupled.e_tot),
            t1_frobenius_over_sqrt_correlated_electrons=float(
                np.linalg.norm(coupled.t1) / np.sqrt(10)
            ),
            t1_largest_singular_value=float(np.linalg.svd(coupled.t1, compute_uv=False)[0]),
        )
        if not coupled.converged:
            raise ValueError("Neutral CCSD did not converge")
        coupled.dump_chk()
        triples = float(timed("triples", coupled.ccsd_t))
        report.update(
            triples_correction_hartree=triples,
            neutral_energy_hartree=float(coupled.e_tot + triples),
        )
        timed("lambda", coupled.solve_lambda)
        report["lambda_converged"] = bool(coupled.converged_lambda)
        if not coupled.converged_lambda:
            raise ValueError("Neutral CCSD lambda equations did not converge")
        density = coupled.make_rdm1(ao_repr=True)
        electrons = float(np.einsum("ij,ji->", density, hf.get_ovlp()))
        np.testing.assert_allclose(electrons, 14, atol=1e-7, rtol=0)
        report.update(
            density_electrons=electrons,
            ccsd_dipole_au=[float(x) for x in hf.dip_moment(mol, density, unit="AU")],
            dipole_method="CCSD lambda density; not a CCSD(T) energy derivative",
        )
        if config["reference_check"]:
            report["independent_reference"] = timed(
                "psi4_reference", lambda: psi4_reference(config, atoms)
            )
            for key in ("rhf_energy_hartree", "neutral_energy_hartree"):
                np.testing.assert_allclose(
                    report[key], report["independent_reference"][key], atol=1e-7, rtol=0
                )
        return report
    finally:
        (workdir / "neutral-diagnostics.json").write_text(json.dumps(report, indent=2) + "\n")


def psi4_reference(config: dict, atoms: list) -> dict:
    """Compute a separate conventional, frozen-core CCSD(T) energy control."""
    workdir = Path(config["workdir"])
    geometry = "0 1\n" + "\n".join(f"{symbol} {x} {y} {z}" for symbol, x, y, z in atoms)
    geometry += "\nunits bohr\nno_com\nno_reorient\n"
    options = {
        "basis": config["basis"],
        "puream": True,
        "scf_type": "pk",
        "cc_type": "conv",
        "freeze_core": True,
        "e_convergence": 1e-10,
        "d_convergence": 1e-10,
        "r_convergence": 1e-8,
    }
    # The pinned Psi4 executable has its own Python environment. Importing it
    # from the PySCF interpreter would bypass that deployment contract.
    input_file = workdir / "psi4-reference.inp"
    input_file.write_text(
        "import json\nimport psi4\n"
        f"psi4.set_memory({str(config['memory_mb']) + ' MB'!r})\n"
        f"psi4.set_num_threads({config['ranks']})\n"
        f"mol = psi4.geometry({geometry!r})\n"
        f"psi4.set_options({options!r})\n"
        "value = psi4.energy('ccsd(t)', molecule=mol)\n"
        "with open('psi4-reference.json', 'w') as output:\n"
        "    json.dump({'psi4_version': psi4.__version__, "
        "'rhf_energy_hartree': float(psi4.variable('SCF TOTAL ENERGY')), "
        "'neutral_energy_hartree': float(value)}, output, indent=2)\n"
    )
    subprocess.run(
        ["psi4", "--input", str(input_file), "--output", str(workdir / "psi4-reference.out")],
        cwd=workdir,
        check=True,
    )
    return json.loads((workdir / "psi4-reference.json").read_text()) | {
        "method": "independent Psi4 conventional RHF/CCSD(T), frozen core",
        "comparison_tolerance_hartree": 1e-7,
    }


def main() -> None:
    """Run one fresh, profiled neutral calculation, preserving failed diagnostics."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--bond-length", type=float, default=2.1323)
    parser.add_argument(
        "--basis",
        choices=["cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ", "aug-cc-pV5Z"],
        default="aug-cc-pVTZ",
    )
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--memory-mb", type=int, default=12000)
    parser.add_argument("--reference-check", action="store_true")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not all(math.isfinite(x) and x > 0 for x in (args.bond_length, args.ranks, args.memory_mb)):
        parser.error("Geometry, threads and memory must be finite and positive")
    workdir = args.workdir.resolve()
    if args.worker:
        config = json.loads((workdir / "config.json").read_text())
        report = calculate(config)
        (workdir / "neutral.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2), flush=True)
        return
    workdir.mkdir(parents=True, exist_ok=False)
    (workdir / "scratch").mkdir()
    config = {k: v for k, v in vars(args).items() if k != "worker"} | {
        "workdir": str(workdir),
        "model": "neutral-ccsd(t)",
        "neutral_only": True,
    }
    (workdir / "config.json").write_text(json.dumps(config, default=str, indent=2) + "\n")
    env = os.environ | {
        "TMPDIR": str(workdir / "scratch"),
        "PSI_SCRATCH": str(workdir / "scratch"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "PYTHONPATH": str(Path(__file__).resolve().parents[2]),
    }
    raise SystemExit(
        run_profiled(
            workdir,
            env,
            [
                "python3",
                "-m",
                "projects.ukrmol_co.neutral",
                "--worker",
                "--workdir",
                str(workdir),
            ],
        )
    )


if __name__ == "__main__":
    main()
