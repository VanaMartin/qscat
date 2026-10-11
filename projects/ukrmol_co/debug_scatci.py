"""One profiled native B1 SCATCI diagnostic from the preserved failed deck and inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

from projects.ukrmol_co.run import run_profiled


def compute(args: argparse.Namespace) -> int:
    directory = args.workdir / "B1"
    directory.mkdir()
    copied = {}
    for name in ("fort.16", "fort.24", "fort.26", "fort.721", "scatci.inp"):
        shutil.copy2(args.reference / name, directory / name)
        copied[name] = hashlib.sha256((directory / name).read_bytes()).hexdigest()
    (directory / "debugger").mkdir()
    env = os.environ | {"CO_DEBUG_DIRECTORY": str(directory / "debugger")}
    command = [
        "mpiexec",
        "-n",
        str(args.ranks),
        "python3",
        "-m",
        "projects.ukrmol_co.rank_debugger",
        "/opt/ukrmolp/bin.double/mpi-scatci",
        "scatci.inp",
    ]
    started = time.monotonic()
    timed_out = False
    with (directory / "scatci.out").open("x") as stream:
        process = subprocess.Popen(
            command,
            cwd=directory,
            env=env,
            stdout=stream,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        try:
            process.wait(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=20)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
    ranks = [
        json.loads(p.read_text()) for p in sorted((directory / "debugger").glob("rank-*.json"))
    ]
    exports_present = all((directory / name).is_file() for name in ("fort.25", "fort.21"))
    report = {
        "command": command,
        "mpi_debugger_original_exit": process.returncode,
        "timed_out": timed_out,
        "wall_seconds": time.monotonic() - started,
        "input_sha256": copied,
        "ranks": ranks,
        "native_export_files_present": exports_present,
        "diagnostic_native_completion": bool(
            process.returncode == 0
            and not timed_out
            and len(ranks) == args.ranks
            and exports_present
            and all(
                r.get("inferior_exit_code") == 0 and not r.get("observed_stop_signal")
                for r in ranks
            )
        ),
        "root_pair_stability_verdict": None,
        "spectral_truncation_verdict": None,
        "scope": "One instrumented B1 SCATCI; no B2,outer stages or convergence qualification",
    }
    (args.workdir / "diagnostic.json").write_text(json.dumps(report, indent=2) + "\n")
    return 0 if report["diagnostic_native_completion"] else 1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("reference", "workdir"):
        parser.add_argument(f"--{field}", type=Path, required=True)
    parser.add_argument("--ranks", type=int, default=8)
    parser.add_argument("--timeout", type=float, default=86400)
    parser.add_argument("--profiled", action="store_true")
    args = parser.parse_args()
    if args.ranks != 8 or not 0 < args.timeout <= 86400:
        parser.error("Declared eight ranks and a positive bounded timeout required")
    args.reference, args.workdir = args.reference.resolve(), args.workdir.resolve()
    if args.profiled:
        raise SystemExit(compute(args))
    args.workdir.mkdir()
    (args.workdir / "scratch").mkdir()
    env = os.environ | {
        "TMPDIR": str(args.workdir / "scratch"),
        "PSI_SCRATCH": str(args.workdir / "scratch"),
    }
    raise SystemExit(
        run_profiled(
            args.workdir,
            env,
            ["python3", "-m", "projects.ukrmol_co.debug_scatci", *sys.argv[1:], "--profiled"],
        )
    )


if __name__ == "__main__":
    main()
