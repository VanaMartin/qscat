"""Finite selected-root scattering pilot from retained native integral/CONGEN inputs.

Unlike the dense-reference replay, this experiment has no same-model dense oracle.
Its export, residual and Pi checks are necessary; truncation verdicts remain unset.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import time
from pathlib import Path

import numpy as np

from projects.ukrmol_co.analyze import native_resonances
from projects.ukrmol_co.run import run_profiled
from projects.ukrmol_co.sparse_scattering import (
    boundary_spectrum,
    ci_spectrum,
    mpi_boundary_deck,
    selected_deck,
)


def target_map(prop: str, states: list[dict], roots: int = 5) -> list[int]:
    """Map spin/irrep contraction order onto DENPROP's energy-ordered state IDs."""
    rows = []
    for line in prop.splitlines():
        if line.split()[:1] != ["5"]:
            continue
        tokens = line.split()
        rows.append(
            (int(tokens[1]), tokens[-2].lower(), tokens[-1], float(tokens[8].replace("D", "E")))
        )
    if len(rows) != 8 * roots or [r[0] for r in rows] != list(range(1, 8 * roots + 1)):
        raise ValueError("Incomplete or reordered DENPROP state table")
    if not np.isfinite([r[3] for r in rows]).all() or np.any(
        np.diff([r[3] for r in rows]) < -1e-10
    ):
        raise ValueError("DENPROP states are not finite and energy ordered")
    mapping = []
    for irrep in ("A1", "B1", "B2", "A2"):
        for spin in ("singlet", "triplet"):
            group = [r for r in rows if r[1:3] == (spin, irrep)]
            expected = [
                s["energy_hartree"] for s in states if (s["spin"], s["irrep"]) == (spin, irrep)
            ]
            if len(group) != roots or len(expected) != roots:
                raise ValueError("Incomplete retained target sector")
            np.testing.assert_allclose([r[3] for r in group], expected, atol=1e-7, rtol=0)
            mapping.extend(r[0] for r in group)
    return mapping


def _save(path: Path, value: dict) -> None:
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def audit_sector(
    directory: Path, roots: int, dimension: int, continuum: int, boundary_channels: int = 250
) -> dict:
    """Independently read emitted CI/boundary records and check the native residual log."""
    telemetry = [
        json.loads(line) for line in (directory / "eps-telemetry.jsonl").read_text().splitlines()
    ]
    if len(telemetry) != 1:
        raise ValueError("No unique solver telemetry")
    t = telemetry[0]
    if not (
        t["requested"] == roots
        and t["converged"] >= roots
        and t["reason"] > 0
        and t["dimension"] == dimension
        and t["matrix_type"] == "mpisbaij"
        and t["solver_type"] == "krylovschur"
    ):
        raise ValueError("Selected-root solver contract failed")
    residuals = np.asarray([r["absolute_residual"] for r in t["eigenpairs"]])
    if len(residuals) < roots or not np.isfinite(residuals).all() or np.max(residuals) > 1e-7:
        raise ValueError("Native eigenpair residual gate failed")
    energies, vectors = ci_spectrum(directory / "fort.25")
    poles, amplitudes = boundary_spectrum(directory / "fort.21")
    if len(energies) != roots or len(poles) != roots or vectors.shape != (roots, continuum):
        raise ValueError("Incomplete selected-spectrum export")
    if amplitudes.shape != (roots, boundary_channels) or not all(
        np.isfinite(v).all() for v in (energies, vectors, poles, amplitudes)
    ):
        raise ValueError("Nonfinite or inconsistent boundary export")
    if np.any(np.diff(energies) <= 0):
        raise ValueError("Selected eigenvalues are not distinct and energy ordered")
    np.testing.assert_allclose(poles, energies, atol=1e-10, rtol=0)
    unit = int(re.search(r"luphso\s*=\s*(\d+)", (directory / "eigenp.inp").read_text())[1])
    phases = np.loadtxt(directory / f"fort.{unit}", ndmin=2)
    if phases.shape != (99, 2) or not np.isfinite(phases).all():
        raise ValueError("Incomplete or nonfinite pilot phase grid")
    np.testing.assert_allclose(phases[:, 0], 0.1 + 0.05 * np.arange(99), atol=5e-5, rtol=0)
    candidates = native_resonances(directory / "reson.out")
    return {
        "export_checks_passed": True,
        "maximum_native_residual_hartree": float(np.max(residuals)),
        "energy_range_hartree": [float(energies[0]), float(energies[-1])],
        "boundary_energy_error_hartree": float(np.max(abs(poles - energies))),
        "native_resonances": candidates,
        "native_maxfit_diagnostic": "MAXFIT" in (directory / "reson.out").read_text().upper(),
        "spectral_truncation_verdict": None,
        "observable_qualification_verdict": None,
    }


def compute(args: argparse.Namespace, *, sectors=("B1", "B2")) -> None:
    if sectors not in (("B1",), ("B1", "B2")):
        raise ValueError("Only the declared B1 pilot or complete B1/B2 pair is supported")
    geometry = args.preflight / "output/CO/geom1"
    preflight = json.loads((args.preflight / "preflight-result.json").read_text())
    if preflight["original_exit"] != 0:
        raise ValueError("Native preflight failed")
    target = json.loads((args.preflight / "target.json").read_text())
    mapping = target_map((geometry / "prop.out").read_text(), target["states"])
    report = {
        "roots": args.roots,
        "target_map": mapping,
        "sectors": {},
        "spectral_truncation_verdict": None,
    }
    for sector in sectors:
        directory = args.workdir / sector
        directory.mkdir()
        text = (geometry / f"inputs/scattering.scatci.doublet.{sector}.inp").read_text()
        unit = int(re.search(r"megul\s*=\s*(\d+)", text)[1])
        row = {"stages": [], "input_sha256": {}}
        report["sectors"][sector] = row
        for name, source in {
            "fort.16": geometry / f"doublet.{sector}" / "fort.16",
            "fort.26": geometry / "fort.26",
            "fort.24": geometry / "fort.24",
            f"fort.{unit}": geometry / f"doublet.{sector}" / f"fort.{unit}",
        }.items():
            shutil.copy2(source, directory / name)
            row["input_sha256"][name] = hashlib.sha256(source.read_bytes()).hexdigest()
        text = selected_deck(text, args.roots, 1e-12, 1000)
        boundary = f"ntarg = 40, idtarg = {','.join(map(str, mapping))}, rmatr = 18,"
        (directory / "scatci.inp").write_text(mpi_boundary_deck(text, boundary))
        for stage in ("rsolve", "eigenp", "reson"):
            shutil.copy2(
                args.outer_inputs / f"scattering.{stage}.doublet.{sector}.inp",
                directory / f"{stage}.inp",
            )
        (directory / "reson_message").write_text("")
        for stage in ("scatci", "rsolve", "eigenp", "reson"):
            command = (
                [
                    "mpiexec",
                    "-n",
                    str(args.ranks),
                    "/opt/ukrmolp/bin.double/mpi-scatci",
                    "scatci.inp",
                ]
                if stage == "scatci"
                else [f"/opt/ukrmolp/bin.double/{stage}"]
            )
            env = os.environ.copy()
            if stage == "scatci":
                env.update(
                    LD_PRELOAD=str(args.telemetry),
                    UKRMOL_EPS_AUDIT=str(directory / "eps-telemetry.jsonl"),
                )
            started = time.monotonic()
            with (
                (directory / f"{stage}.inp").open() as inp,
                (directory / f"{stage}.out").open("x") as out,
            ):
                process = subprocess.Popen(
                    command,
                    cwd=directory,
                    env=env,
                    stdin=inp,
                    stdout=out,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )
                timed_out = False
                try:
                    process.wait(
                        timeout=args.stage_timeout
                        if stage == "scatci"
                        else getattr(args, "outer_stage_timeout", args.stage_timeout)
                    )
                except subprocess.TimeoutExpired:
                    timed_out = True
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=20)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
            row["stages"].append(
                {
                    "name": stage,
                    "original_exit": process.returncode,
                    "timed_out": timed_out,
                    "wall_seconds": time.monotonic() - started,
                }
            )
            _save(args.workdir / "pilot.json", report)
            output = (directory / f"{stage}.out").read_text(errors="replace")
            if (
                process.returncode
                or timed_out
                or any(
                    s in output
                    for s in ("MPI_ABORT", "***DUE TO AN ERROR", "NOT IMPLEMENTED FOR PARTITIONED")
                )
            ):
                raise ValueError(f"Native {sector} {stage} failed")
        row.update(
            audit_sector(
                directory,
                args.roots,
                preflight["sectors"][sector]["contracted_dimension"],
                preflight["sectors"][sector]["contracted_continuum_dimension"],
            )
        )
        _save(args.workdir / "pilot.json", report)
    if sectors == ("B1",):
        report.update(
            all_stages_complete=True,
            pi_phase_difference_rad_mod_pi=None,
            pi_check_passed=None,
            scope="B1 only; no counterpart, pole, electronic-model or omitted-spectrum verdict",
        )
        _save(args.workdir / "pilot.json", report)
        return
    phases = []
    for sector in ("B1", "B2"):
        p = args.workdir / sector
        unit = int(re.search(r"luphso\s*=\s*(\d+)", (p / "eigenp.inp").read_text())[1])
        phases.append(np.loadtxt(p / f"fort.{unit}", ndmin=2))
    np.testing.assert_allclose(phases[0][:, 0], phases[1][:, 0], atol=5e-5, rtol=0)
    delta = float(
        np.max(abs((phases[0][:, 1:] - phases[1][:, 1:] + np.pi / 2) % np.pi - np.pi / 2))
    )
    report.update(
        all_stages_complete=True,
        pi_phase_difference_rad_mod_pi=delta,
        pi_check_passed=delta <= 0.002,
    )
    _save(args.workdir / "pilot.json", report)
    if delta > 0.002:
        raise ValueError("Pi phase symmetry failed")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("preflight", "outer-inputs", "workdir", "telemetry"):
        parser.add_argument(f"--{field}", type=Path, required=True)
    parser.add_argument("--roots", type=int, default=2048)
    parser.add_argument("--ranks", type=int, default=8)
    parser.add_argument("--stage-timeout", type=float, default=86400)
    parser.add_argument("--profiled", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for field in ("preflight", "outer_inputs", "workdir", "telemetry"):
        setattr(args, field, getattr(args, field).resolve())
    if args.roots <= 0 or args.ranks <= 0 or not 0 < args.stage_timeout < np.inf:
        parser.error("Positive finite counts and timeout required")
    if args.profiled:
        compute(args)
    else:
        args.workdir.mkdir()
        (args.workdir / "scratch").mkdir()
        _save(
            args.workdir / "config.json",
            {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        )
        env = os.environ | {
            "TMPDIR": str(args.workdir / "scratch"),
            "PSI_SCRATCH": str(args.workdir / "scratch"),
        }
        import sys

        raise SystemExit(
            run_profiled(
                args.workdir,
                env,
                ["python3", "-m", "projects.ukrmol_co.sparse_pilot", *sys.argv[1:], "--profiled"],
            )
        )


if __name__ == "__main__":
    main()
