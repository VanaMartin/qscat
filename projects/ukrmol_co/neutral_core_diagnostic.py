"""Matched FC/AE neutral diagnostic; never emits accepted neutral-curve rows.

One stable real RHF reference is shared by conventional CCSD(T) with two frozen
occupied orbitals and with none. See docs/physics/co-neutral-core-diagnostic.md.
PySCF imports are lazy; generated backends can qualify orchestration without QC.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from projects.ukrmol_co.run import allocated_bytes, cgroup_anon_bytes, cgroup_bytes

BASIS = "aug-cc-pwCVQZ"
BASIS_SHA256 = "8d53c8167cfa30f6c7bbaebb37d42c3516b0cc5fdedfae84c3e2ed90926305b5"
RADII_BOHR = (2.45, 2.50)
PHASE_CAPS = {"reference": 300, "fc": 1500, "ae": 1500}


def _save(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")
    temporary.replace(path)


def _backend():
    import pyscf
    from pyscf import cc, gto, lib, scf

    return SimpleNamespace(__version__=pyscf.__version__, cc=cc, gto=gto, lib=lib, scf=scf)


def _validate(config):
    radius = config["bond_length"]
    if isinstance(radius, bool) or not isinstance(radius, (int, float)) or radius not in RADII_BOHR:
        raise ValueError("Diagnostic geometry must be exactly 2.45 or 2.50 bohr")
    if config.get("basis", BASIS) != BASIS:
        raise ValueError("Diagnostic basis is fixed")
    if type(config["ranks"]) is not int or config["ranks"] != 4:
        raise ValueError("Need the declared four threads")
    if type(config["memory_mb"]) is not int or not 0 < config["memory_mb"] <= 20000:
        raise ValueError("Memory must be an integer in (0, 20000] MB")
    if "source_identity" in config and config["source_identity"] != source_identity():
        raise ValueError("Parent/child source identity mismatch")


def source_identity():
    """Record actual module paths/digests, catching bundled-source shadowing."""
    from projects.ukrmol_co import run

    return {
        str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (Path(__file__), Path(run.__file__))
    }


def _basis_identity(backend):
    path = Path(backend.gto.basis.__file__).with_name("aug-cc-pwcvqz.dat")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != BASIS_SHA256:
        raise ValueError("Stock core-capable basis digest mismatch")
    return digest


def _prepare(config, backend):
    _validate(config)
    if backend.__version__ != "2.11.0":
        raise ValueError("This diagnostic requires PySCF 2.11.0")
    digest = _basis_identity(backend)
    radius = config["bond_length"]
    z = -radius * 15.9994 / (12.0110 + 15.9994)
    atoms = [["C", 0, 0, z], ["O", 0, 0, z + radius]]
    backend.lib.num_threads(config["ranks"])
    mol = backend.gto.M(
        atom=atoms,
        unit="Bohr",
        charge=0,
        spin=0,
        basis=BASIS,
        symmetry="C2v",
        cart=False,
        verbose=4,
        max_memory=config["memory_mb"],
    )
    if (mol.nelectron, mol.spin, mol.charge, mol.cart, mol.nao_nr()) != (14, 0, 0, False, 218):
        raise ValueError("Unexpected molecule, spherical AO count or electron count")
    metadata = {
        "schema": "co-neutral-core-pair-diagnostic-v1",
        "diagnostic_only": True,
        "eligible_for_delivered_curve": False,
        "basis": BASIS,
        "basis_data_sha256": digest,
        "pyscf_version": backend.__version__,
        "bond_length_bohr": radius,
        "atoms_bohr": atoms,
        "spherical_AOs": 218,
        "source_identity": source_identity(),
        "hamiltonian": "nonrelativistic, spherical, all-electron molecule, real RHF",
    }
    return mol, metadata


def preflight(config, *, backend=None):
    """Construct/check the molecule and basis only; never call SCF or CC."""
    _, metadata = _prepare(config, _backend() if backend is None else backend)
    return metadata


def _finite(value, name):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"Nonfinite {name}")
    return value


def _reference_identity(hf):
    digest = hashlib.sha256()
    for name, shape in [("mo_coeff", (218, 218)), ("mo_occ", (218,)), ("mo_energy", (218,))]:
        array = np.asarray(getattr(hf, name))
        if array.shape != shape or not np.isrealobj(array) or not np.isfinite(array).all():
            raise ValueError(f"Invalid RHF {name}")
        digest.update(array.dtype.str.encode())
        digest.update(np.ascontiguousarray(array).tobytes())
    if not np.array_equal(hf.mo_occ, np.r_[np.full(7, 2.0), np.zeros(211)]):
        raise ValueError("Expected seven doubly occupied reference orbitals")
    if set(np.argsort(hf.mo_energy[:7])[:2]) != {0, 1}:
        raise ValueError("First two orbitals must be the lowest occupied core pair")
    if getattr(hf, "with_df", None) is not None:
        raise ValueError("Density-fitted reference is not this conventional diagnostic")
    digest.update(_finite(hf.e_tot, "RHF energy").hex().encode())
    return digest.hexdigest()


def _check_mask(coupled, frozen):
    expected = np.ones(218, dtype=bool)
    expected[: frozen or 0] = False
    if (
        coupled.frozen != frozen
        or coupled._nmo is not None
        or coupled._nocc is not None
        or not np.array_equal(coupled.get_frozen_mask(), expected)
        or coupled.nocc != 7 - (frozen or 0)
        or coupled.nmo != 218 - (frozen or 0)
    ):
        raise ValueError("Core mask/occupied count/full virtual space mismatch")


def calculate_pair(config, *, backend=None):
    """Run one fresh matched pair, preserving failures and successful partials."""
    _validate(config)
    backend = _backend() if backend is None else backend
    work = Path(config["workdir"])
    report_path = work / "diagnostic-pair.json"
    # No resume or silent overwrite of a prior success/failure or accepted row.
    for name in [
        "diagnostic-pair.json",
        "neutral.json",
        "neutral.rhf.chk",
        "fc.ccsd.chk",
        "ae.ccsd.chk",
    ]:
        if (work / name).exists():
            raise ValueError("Pair work directory already contains scientific output")
    report = {
        "status": "running_diagnostic",
        "diagnostic_only": True,
        "eligible_for_delivered_curve": False,
        "core_treatments": [],
        "stages": [],
    }
    _save(report_path, report)

    def phase(name):
        _save(
            work / "active-stage.json",
            {"phase": name, "started_monotonic": time.monotonic(), "cap_seconds": PHASE_CAPS[name]},
        )

    def timed(name, operation):
        start = time.monotonic()
        status = 1
        try:
            result = operation()
            status = 0
            return result
        finally:
            report["stages"].append(
                {"name": name, "original_exit": status, "wall_seconds": time.monotonic() - start}
            )
            _save(report_path, report)

    try:
        mol, metadata = _prepare(config, backend)
        report.update(metadata)
        phase("reference")
        hf = backend.scf.RHF(mol)
        hf.conv_tol, hf.conv_tol_grad = 1e-11, 1e-8
        hf.chkfile = str(work / "neutral.rhf.chk")
        timed("rhf", hf.kernel)
        report["rhf_converged"] = bool(hf.converged)
        if not hf.converged:
            raise ValueError("Diagnostic RHF did not converge")
        _, _, internal, external = timed(
            "rhf_stability", lambda: hf.stability(external=True, return_status=True)
        )
        report.update(rhf_internally_stable=bool(internal), rhf_externally_stable=bool(external))
        if not internal or not external:
            raise ValueError("Diagnostic RHF reference is unstable")
        reference = _reference_identity(hf)
        report.update(rhf_reference_sha256=reference, rhf_energy_hartree=float(hf.e_tot))
        for label, frozen in [("fc", 2), ("ae", None)]:
            phase(label)
            if _reference_identity(hf) != reference:
                raise ValueError("Shared RHF reference changed")
            row = {
                "core_treatment": label,
                "status": "running_diagnostic",
                "diagnostic_only": True,
                "eligible_for_delivered_curve": False,
                "frozen_spatial_orbitals": frozen or 0,
                "correlated_electrons": 14 - 2 * (frozen or 0),
                "rhf_reference_sha256": reference,
            }
            report["core_treatments"].append(row)
            coupled = backend.cc.CCSD(hf, frozen=frozen)
            _check_mask(coupled, frozen)
            if _reference_identity(coupled._scf) != reference:
                raise ValueError("CCSD factory changed the shared reference")
            coupled.conv_tol, coupled.conv_tol_normt, coupled.max_cycle = 1e-10, 1e-8, 100
            coupled.chkfile = str(work / f"{label}.ccsd.chk")
            timed(f"{label}.ccsd", coupled.kernel)
            row["ccsd_converged"] = bool(coupled.converged)
            if not coupled.converged:
                raise ValueError("Diagnostic CCSD did not converge")
            _check_mask(coupled, frozen)
            if _reference_identity(hf) != reference:
                raise ValueError("Correlation operation changed shared RHF reference")
            nocc = 7 - (frozen or 0)
            if (
                coupled.t1.shape != (nocc, 211)
                or coupled.t2.shape != (nocc, nocc, 211, 211)
                or not np.isrealobj(coupled.t1)
                or not np.isrealobj(coupled.t2)
                or not np.isfinite(coupled.t1).all()
                or not np.isfinite(coupled.t2).all()
            ):
                raise ValueError("Invalid diagnostic amplitudes")
            row.update(
                ccsd_energy_hartree=_finite(coupled.e_tot, "CCSD energy"),
                correlated_occupied_spatial_orbitals=nocc,
                retained_virtual_orbitals=211,
                t1_frobenius_over_sqrt_correlated_electrons=_finite(
                    np.linalg.norm(coupled.t1) / math.sqrt(row["correlated_electrons"]), "T1 norm"
                ),
                t1_largest_singular_value=_finite(
                    np.linalg.svd(coupled.t1, compute_uv=False)[0], "T1 singular value"
                ),
            )
            coupled.dump_chk()
            triples = _finite(timed(f"{label}.triples", coupled.ccsd_t), "triples")
            row.update(
                triples_correction_hartree=triples,
                total_energy_hartree=_finite(coupled.e_tot + triples, "CCSD(T) energy"),
            )
            timed(f"{label}.lambda", coupled.solve_lambda)
            row["lambda_converged"] = bool(coupled.converged_lambda)
            if not coupled.converged_lambda:
                raise ValueError("Diagnostic lambda did not converge")
            density = coupled.make_rdm1(ao_repr=True)
            if (
                density.shape != (218, 218)
                or not np.isrealobj(density)
                or not np.isfinite(density).all()
            ):
                raise ValueError("Invalid diagnostic AO density")
            electrons = _finite(np.einsum("ij,ji->", density, hf.get_ovlp()), "density trace")
            if abs(electrons - 14) > 1e-7:
                raise ValueError("Diagnostic density must contain all fourteen electrons")
            if _reference_identity(hf) != reference:
                raise ValueError("Shared RHF reference changed during triples/lambda/density")
            row.update(density_electrons=electrons, status="numerically_complete_diagnostic")
            _save(report_path, report)
        report["status"] = "numerically_complete_diagnostic"
        return report
    except Exception as error:
        report.update(status="failed_diagnostic", error=repr(error))
        if (
            report["core_treatments"]
            and report["core_treatments"][-1]["status"] == "running_diagnostic"
        ):
            report["core_treatments"][-1].update(status="failed_diagnostic", error=repr(error))
        raise
    finally:
        _save(report_path, report)


def _pair_complete(path, identity):
    try:
        value = json.loads(path.read_text())
        rows = value["core_treatments"]
        for row, frozen, nocc in zip(rows, [2, 0], [5, 7], strict=True):
            if not (
                row["frozen_spatial_orbitals"] == frozen
                and row["correlated_occupied_spatial_orbitals"] == nocc
                and row["retained_virtual_orbitals"] == 211
                and row["ccsd_converged"] is True
                and row["lambda_converged"] is True
                and math.isfinite(row["density_electrons"])
                and abs(row["density_electrons"] - 14) <= 1e-7
                and all(
                    math.isfinite(row[k])
                    for k in [
                        "total_energy_hartree",
                        "ccsd_energy_hartree",
                        "triples_correction_hartree",
                    ]
                )
                and abs(
                    row["total_energy_hartree"]
                    - row["ccsd_energy_hartree"]
                    - row["triples_correction_hartree"]
                )
                <= 1e-12
            ):
                return False
        names = [
            "rhf",
            "rhf_stability",
            "fc.ccsd",
            "fc.triples",
            "fc.lambda",
            "ae.ccsd",
            "ae.triples",
            "ae.lambda",
        ]
        return (
            value["schema"] == "co-neutral-core-pair-diagnostic-v1"
            and value["status"] == "numerically_complete_diagnostic"
            and value["eligible_for_delivered_curve"] is False
            and value["diagnostic_only"] is True
            and value["source_identity"] == identity
            and value["rhf_internally_stable"] is True
            and value["rhf_externally_stable"] is True
            and [r["core_treatment"] for r in rows] == ["fc", "ae"]
            and [r["correlated_electrons"] for r in rows] == [10, 14]
            and all(
                r["status"] == "numerically_complete_diagnostic"
                and r["eligible_for_delivered_curve"] is False
                and r["rhf_reference_sha256"] == value["rhf_reference_sha256"]
                for r in rows
            )
            and [s["name"] for s in value["stages"]] == names
            and all(s["original_exit"] == 0 for s in value["stages"])
        )
    except (OSError, ValueError, KeyError, TypeError):
        return False


def _stop_child(process):
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass
    # Also remove surviving descendants if their process-group leader exited.
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait(timeout=5)


def _terminate(signum, frame):
    raise SystemExit(128 + signum)


def launch_pair(
    workdir,
    config,
    *,
    command=None,
    wall_seconds=3600,
    poll_seconds=0.2,
    disk_limit_bytes=20 * 2**30,
):
    """Profile and cap only this child group; outer container/host admission is separate.

    ``command`` permits inert subprocess tests, never an alternative scientific
    method through the CLI. Original child exits are retained independently of
    watchdog/schema/source rejection. No retry or second geometry is launched.
    """
    _validate(config)
    if (
        not 0 < wall_seconds <= 3600
        or not 0 < poll_seconds <= 1
        or not 0 < disk_limit_bytes <= 20 * 2**30
    ):
        raise ValueError("Invalid or enlarged local resource caps")
    work = Path(workdir).resolve()
    work.mkdir(parents=True, exist_ok=False)
    (work / "scratch").mkdir()
    identity = source_identity()
    config = config | {"workdir": str(work), "source_identity": identity, "execute_reviewed": True}
    _save(work / "config.json", config)
    command = (
        command
        if command is not None
        else [
            sys.executable,
            "-m",
            "projects.ukrmol_co.neutral_core_diagnostic",
            "--worker",
            "--workdir",
            str(work),
        ]
    )
    source_root = Path(__file__).resolve().parents[2]
    env = os.environ | {
        "PYTHONPATH": str(source_root),
        "TMPDIR": str(work / "scratch"),
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    start = time.monotonic()
    process = None
    report = {
        "command": command,
        "stop_reason": None,
        "source_identity": identity,
        "eligible_for_delivered_curve": False,
        "diagnostic_only": True,
        "peak_container_memory_bytes": 0,
        "peak_anonymous_bytes": 0,
        "peak_allocated_disk_bytes": 0,
    }
    previous_handlers = {
        sig: signal.signal(sig, _terminate) for sig in [signal.SIGTERM, signal.SIGHUP]
    }
    try:
        with (work / "run.log").open("x") as log, (work / "resources.csv").open("x") as samples:
            writer = csv.writer(samples)
            writer.writerow(
                ["elapsed_s", "container_memory_bytes", "anonymous_bytes", "disk_bytes"]
            )
            process = subprocess.Popen(
                command,
                cwd=work,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
            while True:
                now = time.monotonic()
                charge, anon, disk = (
                    cgroup_bytes("memory.current"),
                    cgroup_anon_bytes(),
                    allocated_bytes(work),
                )
                report["peak_container_memory_bytes"] = max(
                    report["peak_container_memory_bytes"], charge
                )
                report["peak_anonymous_bytes"] = max(report["peak_anonymous_bytes"], anon)
                report["peak_allocated_disk_bytes"] = max(report["peak_allocated_disk_bytes"], disk)
                writer.writerow([now - start, charge, anon, disk])
                samples.flush()
                reason = None
                if now - start > wall_seconds:
                    reason = "overall_wall_cap"
                elif disk > disk_limit_bytes:
                    reason = "allocated_disk_cap"
                elif charge > 24 * 2**30:
                    reason = "container_memory_cap"
                elif (work / "active-stage.json").exists():
                    try:
                        stage = json.loads((work / "active-stage.json").read_text())
                        phase = stage["phase"]
                        began = float(stage["started_monotonic"])
                        if phase not in PHASE_CAPS or not math.isfinite(began) or began > now + 1:
                            reason = "invalid_stage_metadata"
                        elif now - began > PHASE_CAPS[phase]:
                            reason = f"{phase}_wall_cap"
                    except (OSError, ValueError, KeyError, TypeError):
                        reason = "invalid_stage_metadata"
                elif now - start > 300:
                    reason = "preflight_wall_cap"
                if reason:
                    report["stop_reason"] = reason
                    _stop_child(process)
                    break
                if process.poll() is not None:
                    break
                time.sleep(poll_seconds)
    except BaseException as error:
        report["stop_reason"] = "launcher_exception"
        report["error"] = repr(error)
        raise
    finally:
        try:
            if process is not None:
                _stop_child(process)
            report.update(
                original_child_exit=None if process is None else process.returncode,
                wall_seconds=time.monotonic() - start,
                source_integrity_after=source_identity() == identity,
            )
            if not report["source_integrity_after"] and report["stop_reason"] is None:
                report["stop_reason"] = "source_changed"
            if (
                report["original_child_exit"] == 0
                and report["stop_reason"] is None
                and not _pair_complete(work / "diagnostic-pair.json", identity)
            ):
                report["stop_reason"] = "missing_or_invalid_pair_output"
            _save(work / "diagnostic-resources.json", report)
        finally:
            for sig, previous in previous_handlers.items():
                signal.signal(sig, previous)
    if report["stop_reason"]:
        return (
            124
            if report["stop_reason"].endswith("_cap")
            or report["stop_reason"] == "invalid_stage_metadata"
            else 2
        )
    return process.returncode


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path)
    parser.add_argument("--bond-length", type=float, choices=RADII_BOHR)
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--memory-mb", type=int, default=20000)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--execute-reviewed", action="store_true")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        if args.workdir is None:
            parser.error("Worker requires its parent workdir")
        config = json.loads((args.workdir / "config.json").read_text())
        if (
            not config.get("execute_reviewed")
            or "source_identity" not in config
            or Path(config["workdir"]).resolve() != args.workdir.resolve()
        ):
            parser.error("Worker requires the exact reviewed parent configuration")
        calculate_pair(config)
        return
    if args.bond_length is None:
        parser.error("An explicit diagnostic geometry is required")
    config = {"bond_length": args.bond_length, "ranks": args.ranks, "memory_mb": args.memory_mb}
    if args.preflight:
        print(json.dumps(preflight(config), indent=2))
        return
    if not args.execute_reviewed or args.workdir is None:
        parser.error("Explicit reviewed execution and a fresh workdir are required")
    raise SystemExit(launch_pair(args.workdir, config))


if __name__ == "__main__":
    main()
