"""Replay native selected-root scattering against an immutable dense reference.

This experimental audit reuses integrals, CONGEN and target CI data. It measures
the scattering solve/export/outer-region stages, not a fresh electronic pipeline.
Run in the pinned SLEPc image with a separately built telemetry interposer.
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
from scipy.io import FortranEOFError, FortranFile

from projects.ukrmol_co.analyze import native_resonances, target_properties
from projects.ukrmol_co.run import run_profiled


def ci_spectrum(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read one double/int64 CIDATA set, returning total energies and continuum vectors."""
    with FortranFile(path, "r") as stream:
        if stream.read_record("S8").tobytes() != b"********                CIDATA  ":
            raise ValueError("Not native CIDATA")
        header = stream.read_record("u1").tobytes()
        # SolutionHandler.write_header_sol: two int64s, a 120-character name,
        # four int64s, spin/spin-z, electron count, then frozen-core energy.
        nnuc, dimension, roots, _ = np.frombuffer(header, "<i8", count=4, offset=136)
        core_energy = np.frombuffer(header, "<f8", count=1, offset=192)[0]
        for _ in range(nnuc):
            stream.read_record("u1")
        record = stream.read_record("<f8")
        if record.size != roots + dimension * 2:
            raise ValueError("Unexpected CIDATA eigenvalue/diagonal record")
        energies = record[dimension : dimension + roots] + core_energy
        vectors = []
        for i in range(roots):
            record = stream.read_record("u1").tobytes()
            index = np.frombuffer(record, "<i8", count=1)[0]
            if index != i + 1:
                raise ValueError("Incomplete CIDATA vector sequence")
            vectors.append(np.frombuffer(record, "<f8", offset=8).copy())
    vectors = np.asarray(vectors)
    if not np.all(np.isfinite(energies)) or not np.all(np.isfinite(vectors)):
        raise ValueError("Nonfinite native eigenpair data")
    return energies, vectors


def selected_deck(text: str, eigenpairs: int, tolerance: float, max_cycles: int) -> str:
    """Replace the unique native solver slot while retaining the contraction contract."""
    if eigenpairs <= 0 or tolerance <= 0 or not np.isfinite(tolerance) or max_cycles <= 0:
        raise ValueError("Positive finite selected-root controls are required")
    uncommented = "\n".join(line.split("!", 1)[0] for line in text.splitlines())
    if re.search(r"(?:^|,)\s*nstat\s*=", uncommented, re.M | re.I):
        raise ValueError("Reference already selects scattering eigenpairs")
    text, count = re.subn(
        r"^!\s*igh\s*=\s*1,\s*$",
        f"  igh = -1, nstat = {eigenpairs}, crite = {tolerance:.17g}, maxiter = {max_cycles},",
        text,
        flags=re.M,
    )
    if count != 1:
        raise ValueError("No unique native scattering diagonalizer slot")
    return text


def dense_root_indices(energies: np.ndarray, dense_energies: np.ndarray) -> np.ndarray:
    """Match energy-ordered oracle roots without allocating a roots-by-dimension matrix."""
    if np.any(np.diff(dense_energies) < 0):
        raise ValueError("Dense oracle spectrum is not energy ordered")
    right = np.clip(np.searchsorted(dense_energies, energies), 0, len(dense_energies) - 1)
    left = np.maximum(right - 1, 0)
    matches = np.where(
        abs(energies - dense_energies[left]) <= abs(energies - dense_energies[right]), left, right
    )
    if len(set(matches.tolist())) != len(matches):
        raise ValueError("Selected roots do not match distinct dense states")
    return matches


def save(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n")


def boundary_spectrum(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read the pinned non-partitioned unformatted R-matrix pole/amplitude layout."""
    with FortranFile(path, "r") as stream:
        header = stream.read_record("u1").tobytes()
        if len(header) != 40 or np.frombuffer(header, "<i8", count=1)[0] != 11:
            raise ValueError("Not native RMTAMP")
        stream.read_record("u1")  # Dataset title.
        _, _, _, channels = stream.read_record("<i8")
        stream.read_record("u1")  # Symmetry, charge and nuclear data.
        record = stream.read_record("u1").tobytes()
        multipoles, roots, _, buttle = np.frombuffer(record, "<i8", count=4)
        if buttle != 0:
            raise ValueError("Partitioned/Buttle boundary layout is unqualified")
        if multipoles > 0:
            stream.read_record("u1")
        energies = stream.read_record("<f8")
        amplitudes = stream.read_record("<f8")
        if energies.size != roots or amplitudes.size != roots * channels:
            raise ValueError("Incomplete boundary pole/amplitude export")
        try:
            stream.read_record("u1")
        except FortranEOFError:
            pass
        else:
            raise ValueError("Unexpected trailing boundary records")
    if not np.all(np.isfinite(energies)) or not np.all(np.isfinite(amplitudes)):
        raise ValueError("Nonfinite native boundary data")
    return energies, amplitudes.reshape(roots, channels)


def interface_digest(path: Path, *, boundary: bool = False) -> str:
    """Hash native channel/threshold/multipole records independently of title/poles."""
    digest = hashlib.sha256()
    pole_record = 6
    with FortranFile(path, "r") as stream:
        index = 0
        while True:
            try:
                record = stream.read_record("u1").tobytes()
            except FortranEOFError:
                break
            if boundary and index == 4:
                multipoles = np.frombuffer(record, "<i8", count=1)[0]
                pole_record = 5 + int(multipoles > 0)
            if index != 1 and not (boundary and index >= pole_record):
                if boundary and index == 4:
                    record = record[:8] + bytes(8) + record[16:]
                digest.update(len(record).to_bytes(8, "little"))
                digest.update(record)
            index += 1
    return digest.hexdigest()


def truncate_boundary(source: Path, destination: Path, roots: int) -> None:
    """Write a dense-oracle low-energy pole subset without a missing-spectrum correction."""
    energies, amplitudes = boundary_spectrum(source)
    if not 0 < roots <= energies.size or np.any(np.diff(energies) < 0):
        raise ValueError("An energy-ordered boundary oracle and valid root count are required")
    with FortranFile(source, "r") as stream:
        metadata = [stream.read_record("u1") for _ in range(5)]
        if np.frombuffer(metadata[4], "<i8", count=1)[0] > 0:
            metadata.append(stream.read_record("u1"))
    # Header record five stores ismax/nstat/npole/ibut/radius (int64/double).
    np.frombuffer(metadata[4], "<i8", count=4)[1] = roots
    with destination.open("xb") as raw, FortranFile(raw, "w") as stream:
        for record in metadata:
            stream.write_record(record)
        stream.write_record(energies[:roots])
        stream.write_record(amplitudes[:roots].ravel())


def mpi_boundary_deck(text: str, swinterf: str) -> str:
    """Use native in-memory boundary export with an explicit uncorrected pole sum."""
    text, count = re.subn(r"\bvecstore\s*=\s*1,", "vecstore = 3,", text, flags=re.I)
    if count != 1 or "&outer_interface" in text.lower():
        raise ValueError("No unique continuum/in-memory export slot")
    target_ids = re.search(r"\bidtarg\s*=\s*([\d,]+)", swinterf, re.I)[1].rstrip(",")
    targets = int(re.search(r"\bntarg\s*=\s*(\d+)", swinterf, re.I)[1])
    radius = float(re.search(r"\brmatr\s*=\s*([\d.]+)", swinterf, re.I)[1])
    if len(target_ids.split(",")) != targets:
        raise ValueError("Incomplete retained-state map for native boundary export")
    return text + (
        "\n&outer_interface\n"
        " write_amp = .true., write_dip = .false., write_rmt = .false.,\n"
        " auto_nvo = .false., lutarg = 24, luchan = 10, lurmt = 21,\n"
        f" rmatr = {radius:.17g}, ntarg = {targets}, idtarg = {target_ids},\n"
        " cform = 'U', rform = 'U',\n/\n"
    )


def replay(args: argparse.Namespace) -> None:
    """Execute both Pi sectors, preserving stage exits before applying differential gates."""
    reference = args.reference.resolve()
    baseline = json.loads((reference / "result.json").read_text())
    if baseline.get("resources", {}).get("exit_code") != 0:
        raise ValueError("Dense reference is not qualified")
    target_properties(reference, json.loads((reference / "config.json").read_text()))
    geometry = reference / "output/CO/geom1"
    report = {
        "reference": str(reference),
        "eigenpairs": args.eigenpairs,
        "boundary_export": args.boundary_export,
        "omitted_spectrum": "uncorrected truncated pole sum"
        if args.boundary_export == "mpi"
        else "native SWINTERF partitioned path (unqualified)",
        "sectors": {},
    }
    for irrep in args.sectors:
        directory = args.workdir / irrep
        directory.mkdir()
        retained = geometry / f"doublet.{irrep}"
        original = geometry / f"inputs/scattering.scatci.doublet.{irrep}.inp"
        text = original.read_text()
        unit = int(re.search(r"megul\s*=\s*(\d+)", text)[1])
        inputs = {
            "fort.16": retained / "fort.16",
            "fort.26": retained / "fort.26",
            f"fort.{unit}": retained / f"fort.{unit}",
            "fort.24": geometry / "fort.24",
        }
        row = {"stages": [], "input_sha256": {}, "differential_passed": False}
        report["sectors"][irrep] = row
        for name, source in inputs.items():
            if not source.is_file():
                raise FileNotFoundError(source)
            # Copy: engine routines may open or position their unformatted inputs.
            shutil.copy2(source, directory / name)
            row["input_sha256"][name] = hashlib.sha256(source.read_bytes()).hexdigest()
        text = selected_deck(text, args.eigenpairs, args.tolerance, args.max_cycles)
        if args.boundary_export == "mpi":
            text = mpi_boundary_deck(
                text, (geometry / f"inputs/scattering.swinterf.doublet.{irrep}.inp").read_text()
            )
        (directory / "scatci.inp").write_text(text)
        env = os.environ | {
            "LD_PRELOAD": str(args.telemetry.resolve()),
            "UKRMOL_EPS_AUDIT": str(directory.resolve() / "eps-telemetry.jsonl"),
        }
        commands = [
            (
                "scatci",
                [
                    "mpiexec",
                    "-n",
                    str(args.ranks),
                    "/opt/ukrmolp/bin.double/mpi-scatci",
                    "scatci.inp",
                ],
            ),
            ("swinterf", ["/opt/ukrmolp/bin.double/swinterf"]),
            ("rsolve", ["/opt/ukrmolp/bin.double/rsolve"]),
            ("eigenp", ["/opt/ukrmolp/bin.double/eigenp"]),
            ("reson", ["/opt/ukrmolp/bin.double/reson"]),
        ]
        if args.boundary_export == "mpi":
            commands = [(stage, command) for stage, command in commands if stage != "swinterf"]
        (directory / "reson_message").write_text("")
        for stage, command in commands:
            if stage != "scatci":
                shutil.copy2(
                    geometry / f"inputs/scattering.{stage}.doublet.{irrep}.inp",
                    directory / f"{stage}.inp",
                )
            started = time.monotonic()
            with (
                (directory / f"{stage}.inp").open() as inp,
                (directory / f"{stage}.out").open("x") as out,
            ):
                process = subprocess.Popen(
                    command,
                    cwd=directory,
                    env=env if stage == "scatci" else os.environ,
                    stdin=inp,
                    stdout=out,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )
                timed_out = False
                try:
                    process.wait(timeout=args.stage_timeout)
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
            save(args.workdir / "comparison.json", report)
            output = (directory / f"{stage}.out").read_text(errors="replace")
            if (
                process.returncode
                or timed_out
                or "MPI_ABORT" in output
                or "***DUE TO AN ERROR" in output
                or "NOT IMPLEMENTED FOR PARTITIONED" in output
            ):
                raise ValueError(f"Native {irrep} {stage} failed")
            if stage == "scatci":
                telemetry = [
                    json.loads(line)
                    for line in (directory / "eps-telemetry.jsonl").read_text().splitlines()
                ]
                if len(telemetry) != 1:
                    raise ValueError("No unique scattering telemetry record")
                audit = telemetry[0]
                row["solver"] = audit
                if (
                    audit["matrix_type"] not in ("seqsbaij", "mpisbaij")
                    or audit["solver_type"] != "krylovschur"
                    or audit["requested"] != args.eigenpairs
                    or audit["converged"] < args.eigenpairs
                    or audit["reason"] <= 0
                ):
                    raise ValueError("Sparse scattering solver contract failed")
                energies, vectors = ci_spectrum(directory / "fort.25")
                dense_energies, dense_vectors = ci_spectrum(retained / "fort.25")
                if len(energies) != args.eigenpairs or vectors.shape[1] != dense_vectors.shape[1]:
                    raise ValueError("Incomplete scattering eigenpair/continuum export")
                # Match each emitted state to the dense spectrum; selection need not be lowest-real.
                matches = dense_root_indices(energies, dense_energies)
                row["dense_root_indices_one_based"] = (matches + 1).tolist()
                row["energy_max_difference_hartree"] = float(
                    np.max(abs(energies - dense_energies[matches]))
                )
                signs = np.sign(np.sum(vectors * dense_vectors[matches], axis=1))
                row["continuum_vector_max_difference"] = float(
                    np.max(abs(vectors * signs[:, None] - dense_vectors[matches]))
                )
                row["energy_range_hartree"] = [float(min(energies)), float(max(energies))]
                row["max_absolute_residual_hartree"] = max(
                    pair["absolute_residual"] for pair in audit["eigenpairs"]
                )
                gaps = np.diff(dense_energies)
                if np.any(gaps < 0):
                    raise ValueError("Dense oracle spectrum is not energy ordered")
                nearest_gap = np.minimum(np.r_[np.inf, gaps], np.r_[gaps, np.inf])
                isolated = nearest_gap[matches] > 1e-6
                row["nonisolated_selected_roots"] = int(np.count_nonzero(~isolated))
                row["isolated_continuum_vector_max_difference"] = (
                    float(
                        np.max(
                            abs(
                                vectors[isolated] * signs[isolated, None]
                                - dense_vectors[matches[isolated]]
                            )
                        )
                    )
                    if np.any(isolated)
                    else None
                )
                save(args.workdir / "comparison.json", report)
                if row["energy_max_difference_hartree"] > 1e-7:
                    raise ValueError("Selected scattering eigenvalues fail dense differential")
                if row["max_absolute_residual_hartree"] > 1e-7:
                    raise ValueError("Selected scattering eigenpairs fail residual gate")
                if (
                    row["isolated_continuum_vector_max_difference"] is not None
                    and row["isolated_continuum_vector_max_difference"] > 1e-5
                ):
                    raise ValueError("Isolated continuum coefficients fail dense differential")
                if args.boundary_export == "mpi":
                    data = geometry.parent / "collected_scattering_data"
                    for name, current, oracle, boundary in (
                        (
                            "channel",
                            directory / "fort.10",
                            data / f"channels/channels.geom1.doublet.{irrep}",
                            False,
                        ),
                        (
                            "boundary_static",
                            directory / "fort.21",
                            data / f"rmat_amplitudes/ramps.geom1.doublet.{irrep}",
                            True,
                        ),
                    ):
                        row[f"{name}_sha256"] = interface_digest(current, boundary=boundary)
                        if row[f"{name}_sha256"] != interface_digest(oracle, boundary=boundary):
                            raise ValueError(
                                "Native channels/thresholds/multipoles differ from dense reference"
                            )
                    pole_e, amp = boundary_spectrum(directory / "fort.21")
                    dense_e, dense_amp = boundary_spectrum(
                        geometry.parent
                        / "collected_scattering_data/rmat_amplitudes"
                        / f"ramps.geom1.doublet.{irrep}"
                    )
                    if len(pole_e) != args.eigenpairs or amp.shape[1] != dense_amp.shape[1]:
                        raise ValueError("Incomplete native boundary export")
                    row["boundary_energy_max_difference_hartree"] = float(
                        np.max(abs(pole_e - dense_e[matches]))
                    )
                    row["boundary_amplitude_max_difference"] = float(
                        np.max(abs(amp * signs[:, None] - dense_amp[matches]))
                    )
                    row["isolated_boundary_amplitude_max_difference"] = (
                        float(
                            np.max(
                                abs(
                                    amp[isolated] * signs[isolated, None]
                                    - dense_amp[matches[isolated]]
                                )
                            )
                        )
                        if np.any(isolated)
                        else None
                    )
                    save(args.workdir / "comparison.json", report)
                    if row["boundary_energy_max_difference_hartree"] > 1e-7:
                        raise ValueError("Native boundary energies fail dense differential")
                    if (
                        row["isolated_boundary_amplitude_max_difference"] is not None
                        and row["isolated_boundary_amplitude_max_difference"] > 1e-5
                    ):
                        raise ValueError("Isolated boundary amplitudes fail dense differential")
        eigen_input = (directory / "eigenp.inp").read_text()
        phase_unit = int(re.search(r"luphso\s*=\s*(\d+)", eigen_input)[1])
        phases = np.loadtxt(directory / f"fort.{phase_unit}", ndmin=2)
        dense_phases = np.loadtxt(geometry / f"eigenph.doublet.{irrep}", ndmin=2)
        if phases.shape != dense_phases.shape or not np.all(np.isfinite(phases)):
            raise ValueError("Incomplete or nonfinite native phase grid")
        np.testing.assert_allclose(phases[:, 0], dense_phases[:, 0], atol=5e-5, rtol=0)
        delta = phases[:, 1:] - dense_phases[:, 1:]
        row["phase_max_difference_rad_mod_pi"] = float(
            np.max(abs((delta + np.pi / 2) % np.pi - np.pi / 2))
        )
        save(args.workdir / "comparison.json", report)
        row["native_resonances"] = native_resonances(directory / "reson.out")
        # Position/width qualification uses independent fixed-window fits later.
        row["differential_passed"] = row["phase_max_difference_rad_mod_pi"] <= 0.05
        save(args.workdir / "comparison.json", report)
    report["all_stages_complete"] = True
    report["phase_gate_passed"] = all(
        row["differential_passed"] for row in report["sectors"].values()
    )
    report["resonance_fit_verdict"] = None
    save(args.workdir / "comparison.json", report)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--telemetry", type=Path, required=True)
    parser.add_argument("--eigenpairs", type=int, required=True)
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--tolerance", type=float, default=1e-12)
    parser.add_argument("--max-cycles", type=int, default=1000)
    parser.add_argument("--stage-timeout", type=float, default=21600)
    parser.add_argument("--boundary-export", choices=["swinterf", "mpi"], default="swinterf")
    parser.add_argument("--sectors", choices=["B1", "B2"], nargs="+", default=["B1", "B2"])
    parser.add_argument("--profiled", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    for field in ("reference", "workdir", "telemetry"):
        setattr(args, field, getattr(args, field).resolve())
    if (
        args.eigenpairs <= 0
        or args.ranks <= 0
        or args.max_cycles <= 0
        or not 0 < args.tolerance < np.inf
        or not 0 < args.stage_timeout < np.inf
    ):
        parser.error("Eigenpairs, ranks, cycles and finite tolerance/timeout must be positive")
    if args.profiled:
        replay(args)
    else:
        args.workdir.mkdir(parents=True)
        (args.workdir / "scratch").mkdir()
        save(
            args.workdir / "config.json",
            {
                key: str(value) if isinstance(value, Path) else value
                for key, value in vars(args).items()
            },
        )
        command = [
            "python3",
            "-m",
            "projects.ukrmol_co.sparse_scattering",
            "--reference",
            str(args.reference),
            "--workdir",
            str(args.workdir),
            "--telemetry",
            str(args.telemetry),
            "--eigenpairs",
            str(args.eigenpairs),
            "--ranks",
            str(args.ranks),
            "--tolerance",
            str(args.tolerance),
            "--max-cycles",
            str(args.max_cycles),
            "--stage-timeout",
            str(args.stage_timeout),
            "--boundary-export",
            args.boundary_export,
            "--sectors",
            *args.sectors,
            "--profiled",
        ]
        env = os.environ | {
            "TMPDIR": str(args.workdir / "scratch"),
            "PSI_SCRATCH": str(args.workdir / "scratch"),
        }
        raise SystemExit(run_profiled(args.workdir, env, command))


if __name__ == "__main__":
    main()
