"""Exact non-A1 controls, then at most two diagnostic fixed-orbital A2 solves."""

from __future__ import annotations

import argparse
import datetime
import json
import math
import os
import shutil
import subprocess
import sys
import time
import traceback
from itertools import combinations
from pathlib import Path

from projects.ukrmol_co.compressed_cas12_properties import digest

CONTRACT = "projects/ukrmol_co/singlet-a2-pilot-contract-20261026.json"
NAME = "co-singlet-a2-pilot-20261026"


def save(path, value):
    partial = path.with_suffix(".partial")
    partial.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    partial.replace(path)


def authenticate(root, source, contract):
    for base, key in [(root, "input_sha256"), (source, "source_sha256")]:
        for name, want in contract[key].items():
            if digest(base / name) != want:
                raise ValueError(f"Changed frozen prerequisite: {name}")


def authenticate_runtime(contract):
    import pyscf
    from pyscf.fci import (
        addons,
        cistring,
        direct_spin0,
        direct_spin0_symm,
        direct_spin1,
        direct_spin1_symm,
        spin_op,
    )

    if pyscf.__version__ != "2.11.0":
        raise ValueError("Pinned PySCF2.11.0 required")
    for module in [
        addons,
        cistring,
        direct_spin0,
        direct_spin0_symm,
        direct_spin1,
        direct_spin1_symm,
        spin_op,
    ]:
        name = module.__name__.rsplit(".", 1)[-1]
        if digest(Path(module.__file__)) != contract["package_sha256"][name]:
            raise ValueError(f"Changed numerical package source: {name}")


def state_passes(row):
    """Singlet support is required even for a small full-Hamiltonian residual."""
    return (
        all(math.isfinite(value) for value in row.values())
        and abs(row["spin_square"]) <= 1e-6
        and max(
            row[key]
            for key in [
                "full_residual_hartree",
                "penalized_residual_hartree",
                "spin_penalty_action_hartree",
                "restricted_full_action_difference_hartree",
                "spin_action_norm",
                "alpha_beta_exchange_error",
                "spatial_leakage",
            ]
        )
        <= 1e-9
    )


def fixed_a2(checkpoint, roots, output):
    """One fresh fixed-orbital request; full spin1 action has no spatial compression."""
    import numpy as np
    from pyscf import ao2mo, fci, lib, mcscf, scf, symm

    started = time.monotonic()
    lib.num_threads(4)
    mol = lib.chkfile.load_mol(str(checkpoint))
    mol.verbose = 4
    coefficients = lib.chkfile.load(str(checkpoint), "mcscf/mo_coeff")
    if list(mol.atom_charges()) != [6, 8]:
        raise ValueError("C/O diagnostic geometry required")
    metric = scf.RHF(mol).get_ovlp()
    mo_error = float(
        np.max(abs(coefficients.T @ metric @ coefficients - np.eye(len(coefficients))))
    )
    if mo_error > 1e-9:
        raise ValueError("Frozen MOs are not orthonormal")
    if (
        lib.chkfile.load(str(checkpoint), "mcscf/ncore") != 2
        or lib.chkfile.load(str(checkpoint), "mcscf/ncas") != 12
    ):
        raise ValueError("The frozen CAS12 checkpoint is required")
    np.testing.assert_allclose(
        np.linalg.norm(mol.atom_coords()[1] - mol.atom_coords()[0]), 1.9, atol=1e-12, rtol=0
    )
    if str(mol.basis).lower() != "cc-pvtz" or mol.cart:
        raise ValueError("Wrong basis/convention in diagnostic checkpoint")
    labels = symm.label_orb_symm(mol, mol.irrep_name, mol.symm_orb, coefficients[:, 2:14])
    if {ir: int(sum(labels == ir)) for ir in ["A1", "B1", "B2", "A2"]} != {
        "A1": 6,
        "B1": 3,
        "B2": 3,
        "A2": 0,
    }:
        raise ValueError("Wrong fixed active space")
    cas = mcscf.CASCI(scf.RHF(mol), 12, (5, 5), ncore=2)
    cas.canonicalization = False
    solver = fci.direct_spin0_symm.FCI(mol)
    solver.spin, solver.wfnsym = 0, "A2"
    solver.nroots, solver.max_space, solver.max_cycle = roots, 160, 200
    solver.conv_tol, solver.conv_tol_residual, solver.lindep = 1e-12, 1e-9, 1e-22
    cas.fcisolver = fci.addons.fix_spin_(solver, shift=1.0, ss=0.0)
    cas.kernel(coefficients)
    h1, ecore = cas.get_h1eff(coefficients)
    eri = ao2mo.restore(1, cas.get_h2eff(coefficients), 12)
    effective = fci.direct_spin1.absorb_h1e(h1, eri, 12, (5, 5), 0.5)
    ids = np.array([{"A1": 0, "A2": 1, "B1": 2, "B2": 3}[str(label)] for label in labels])
    strings = sorted(sum(1 << i for i in occupied) for occupied in combinations(range(12), 5))
    np.testing.assert_array_equal(strings, fci.cistring.gen_strings4orblist(range(12), 5))
    determinants = np.array(
        [np.bitwise_xor.reduce([ids[i] for i in range(12) if state >> i & 1]) for state in strings]
    )
    forbidden = (determinants[:, None] ^ determinants) != 1
    np.savez_compressed(
        output / f"inputs-{roots}.npz",
        h1=h1,
        eri=eri,
        ecore=ecore,
        orbsym=ids,
        mo_coefficients=coefficients,
        mo_overlap=metric,
    )
    rows = []
    vectors = []
    flags = np.atleast_1d(cas.fcisolver.converged).tolist()
    for index, (energy, ci) in enumerate(zip(np.atleast_1d(cas.e_tot), cas.ci, strict=True)):
        ci = np.asarray(ci)
        norm = float(np.linalg.norm(ci))
        full = fci.direct_spin1.contract_2e(effective, ci, 12, (5, 5))
        spin_action = fci.spin_op.contract_ss(ci, 12, (5, 5))
        penalized = cas.fcisolver.contract_2e(effective, ci, 12, (5, 5))
        restricted = cas.fcisolver.undo_fix_spin().contract_2e(effective, ci, 12, (5, 5))
        row = dict(
            root=index + 1,
            energy_hartree=float(energy),
            spin_square=float(np.vdot(ci, spin_action).real / norm**2),
            full_residual_hartree=float(np.linalg.norm(full - (energy - ecore) * ci) / norm),
            penalized_residual_hartree=float(
                np.linalg.norm(penalized - (energy - ecore) * ci) / norm
            ),
            spin_penalty_action_hartree=float(np.linalg.norm(penalized - full) / norm),
            restricted_full_action_difference_hartree=float(
                np.linalg.norm(restricted - full) / norm
            ),
            spin_action_norm=float(np.linalg.norm(spin_action) / norm),
            alpha_beta_exchange_error=float(np.linalg.norm(ci - ci.T) / norm),
            spatial_leakage=float(np.linalg.norm(ci[forbidden]) / norm),
        )
        np.savez_compressed(
            output / f"root-{roots}-{index + 1}.npz",
            ci=ci,
            full_action=full,
            spin_action=spin_action,
            penalized_action=penalized,
            restricted_action=restricted,
        )
        rows.append(row)
        vectors.append(ci.ravel() / norm)
    orthogonality = float(np.max(abs(np.array(vectors) @ np.array(vectors).T - np.eye(roots))))
    passed = (
        all(flags)
        and len(rows) == roots
        and orthogonality <= 1e-9
        and all(state_passes(row) for row in rows)
    )
    result = dict(
        roots=roots,
        space=160,
        cycles=200,
        ci_converged=flags,
        states=rows,
        orthogonality_error=orthogonality,
        mo_orthogonality_error=mo_error,
        passed=bool(passed),
        wall_seconds=time.monotonic() - started,
    )
    save(output / f"case-{roots}.json", result)
    return result


def compute(root, source, work, contract):
    from projects.ukrmol_co.singlet_a2_reference import exact_controls, mixed_optimizer_control

    report = dict(
        scientific_exit=2,
        pilot_verdict=None,
        molecular_ci_calls_reserved=0,
        new_qc_calls=0,
        qc_count_scope="CO only; H4 mixed-optimizer validation is separate",
        native_calls=0,
        new_maps=0,
        failed_optimization_qualified=False,
        new_qc_released=False,
        full_sweep_released=False,
        route_holdout_s_matrices_analyzed=0,
        cases=[],
    )
    started = time.monotonic()
    try:
        authenticate(root, source, contract)
        authenticate_runtime(contract)
        failed = json.loads((root / contract["failed_result"]).read_text())
        if not (
            failed["scientific_exit"] == 2
            and failed["original_runner_exit"] == 1
            and failed["qc_candidate_verdict"] is None
        ):
            raise ValueError("The original rejected QC must remain rejected")
        report["exact_controls"] = exact_controls()
        save(work / "result.json", report)
        report["mixed_optimizer_control"] = mixed_optimizer_control()
        save(work / "result.json", report)
        checkpoint = root / contract["checkpoint"]
        for roots in [5, 8]:
            report["molecular_ci_calls_reserved"] += 1
            save(work / "result.json", report)
            case = fixed_a2(checkpoint, roots, work)
            report["cases"].append(case)
            save(work / "result.json", report)
            if not case["passed"]:
                raise ValueError("Diagnostic A2 gate failed; no retry or QC release")
        left, right = report["cases"]
        error = max(
            abs(a["energy_hartree"] - b["energy_hartree"])
            for a, b in zip(left["states"], right["states"][:5], strict=True)
        )
        report["common_first_five_max_difference_hartree"] = error
        if error > 1e-7:
            raise ValueError("Requested-count agreement failed")
        old = json.loads((root / contract["failed_diagnostics"]).read_text())[
            "state_energies_hartree"
        ][15:20]
        report["old_rejected_a2_reported_energies_hartree"] = old
        report["diagnostic_objective_substitution_shift_hartree"] = (
            sum(a["energy_hartree"] - b for a, b in zip(left["states"], old, strict=True)) / 40
        )
        report["diagnostic_substitution_is_not_a_reoptimized_objective"] = True
        report.update(pilot_verdict=True, scientific_exit=0)
    except Exception as error:
        report.update(error=repr(error), traceback=traceback.format_exc())
    finally:
        try:
            authenticate(root, source, contract)
        except Exception as error:
            report.update(integrity_error=repr(error), pilot_verdict=None, scientific_exit=2)
        report["wall_seconds"] = time.monotonic() - started
        save(work / "result.json", report)
    return report["scientific_exit"]


def container_command(root, source, output, contract):
    return [
        "docker",
        "run",
        "--name",
        NAME,
        "--cpuset-cpus=0-3",
        "--memory=8g",
        "--memory-swap=8g",
        "--network=none",
        "--entrypoint",
        "/opt/ukrmolp/entrypoint.sh",
        "-e",
        "OMP_NUM_THREADS=1",
        "-e",
        "OPENBLAS_NUM_THREADS=1",
        "-e",
        "MKL_NUM_THREADS=1",
        "-e",
        "PYTHONDONTWRITEBYTECODE=1",
        "-e",
        "PYTHONPATH=/opt/qmodeling",
        "-v",
        str(root) + ":/work:ro",
        "-v",
        str(source) + ":/opt/qmodeling:ro",
        "-v",
        str(output) + ":/out",
        "-w",
        "/opt/qmodeling",
        contract["image"],
        "python3",
        "-m",
        "projects.ukrmol_co.singlet_a2_pilot",
        "--inside",
    ]


def own(capsule, archive_sha256, commit):
    root, source = capsule.parent.parent, capsule / "source"
    contract = json.loads((source / CONTRACT).read_text())
    pins = json.loads((capsule / "source-manifest.json").read_text())
    started = time.monotonic()
    report = dict(
        source_commit=commit,
        source_archive_sha256=archive_sha256,
        owner_sha256=digest(Path(__file__)),
        started_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        original_owner_exit=2,
    )

    def source_ok():
        return digest(capsule / "source.tar.gz") == archive_sha256 and all(
            digest(source / n) == w for n, w in pins.items()
        )

    try:
        if not source_ok():
            raise ValueError("Frozen source integrity failure")
        report["topology"] = subprocess.check_output(
            ["lscpu", "-e=CPU,CORE,SOCKET,NODE"], text=True
        )
        rows = [
            s.split() for s in report["topology"].splitlines()[1:] if int(s.split()[0]) in range(4)
        ]
        if len(rows) != 4 or len({(s[1], s[2]) for s in rows}) != 4:
            raise ValueError("Not four distinct physical cores")
        ram = (
            int(
                next(
                    s.split()[1]
                    for s in Path("/proc/meminfo").read_text().splitlines()
                    if s.startswith("MemAvailable:")
                )
            )
            * 1024
        )
        report.update(
            available_ram_bytes_before=ram,
            free_disk_bytes_before=shutil.disk_usage(root).free,
            prior_workloads=subprocess.check_output(
                ["docker", "ps", "--format", "{{.Names}} {{.Status}}"], text=True
            ),
        )
        if ram < 16 * 2**30 or report["free_disk_bytes_before"] < 2 * 2**30:
            raise ValueError("RAM/disk launch gate failed")
        output = capsule / "outputs"
        output.mkdir(mode=0o777)
        output.chmod(0o777)
        command = container_command(root, source, output, contract)
        report.update(state="running", command=command)
        save(capsule / "execution.json", report)
        with (capsule / "owner.log").open("x") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            try:
                process.wait(timeout=contract["container_timeout_seconds"])
                report["timed_out"] = False
            except subprocess.TimeoutExpired:
                report["timed_out"] = True
                subprocess.run(["docker", "kill", NAME], check=False, timeout=30)
                process.wait(timeout=30)
        report["original_command_exit"] = process.returncode
        state = json.loads(subprocess.check_output(["docker", "inspect", NAME], text=True))[0]
        save(capsule / "container-inspect.json", state)
        report.update(
            original_container_exit=state["State"]["ExitCode"],
            container_oom_killed=state["State"]["OOMKilled"],
        )
        report["original_owner_exit"] = (
            0
            if process.returncode == 0
            and not report["timed_out"]
            and not report["container_oom_killed"]
            else 2
        )
    except Exception as error:
        report["error"] = repr(error)
    finally:
        report.update(
            state="finished",
            source_integrity_after=source_ok(),
            finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),
            wall_seconds=time.monotonic() - started,
        )
        if not report["source_integrity_after"]:
            report["original_owner_exit"] = 2
        save(capsule / "execution.json", report)
        print(json.dumps(report, indent=2))
    return report["original_owner_exit"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inside", action="store_true")
    parser.add_argument("--compute", action="store_true")
    parser.add_argument("--capsule", type=Path)
    parser.add_argument("--archive-sha256")
    parser.add_argument("--commit")
    args = parser.parse_args()
    if args.inside or args.compute:
        root, source, work = Path("/work"), Path("/opt/qmodeling"), Path("/out/pilot")
        contract = json.loads((source / CONTRACT).read_text())
        if args.compute:
            raise SystemExit(compute(root, source, work, contract))
        from projects.ukrmol_co.run import run_profiled

        work.mkdir()
        (work / "scratch").mkdir()
        raise SystemExit(
            run_profiled(
                work,
                os.environ | {"TMPDIR": str(work / "scratch")},
                [sys.executable, "-m", "projects.ukrmol_co.singlet_a2_pilot", "--compute"],
            )
        )
    if not all([args.capsule, args.archive_sha256, args.commit]):
        parser.error("Owner needs capsule/archive-sha256/commit")
    raise SystemExit(own(args.capsule, args.archive_sha256, args.commit))
