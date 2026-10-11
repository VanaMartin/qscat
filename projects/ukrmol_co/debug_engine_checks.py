"""Gate the diagnostic engine and prove normal-exit, first-signal and bounds capture."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    args = parser.parse_args()
    root = args.build
    results = {}
    gate = [
        "ctest",
        "--test-dir",
        str(root / "build"),
        "--output-on-failure",
        "-j1",
        "-R",
        "C2v_target_HF.*_(serial|parallel)$",
        "--no-tests=error",
    ]
    with (root / "upstream-tests.log").open("x") as stream:
        completed = subprocess.run(gate, stdout=stream, stderr=subprocess.STDOUT, timeout=2700)
    results["upstream_original_exit"] = completed.returncode
    (root / "engine-checks.json").write_text(json.dumps(results, indent=2) + "\n")
    if completed.returncode != 0:
        raise RuntimeError("Diagnostic engine upstream gate failed")
    text = (root / "upstream-tests.log").read_text()
    if "100% tests passed" not in text or "out of 12" not in text:
        raise RuntimeError("Upstream gate did not run and pass all twelve cases")
    smoke = root / "debugger-smoke"
    smoke.mkdir()
    sources = {
        "normal": ("c", "int main(void) {return 0;}\n"),
        "signal": ("c", "#include <signal.h>\nint main(void) {raise(SIGSEGV);return 0;}\n"),
        "bounds": ("f90", "program bad\ninteger :: a(1),i\ni=2\na(i)=7\nprint *,a\nend program\n"),
    }
    for name, (extension, source) in sources.items():
        directory = smoke / name
        directory.mkdir()
        path = directory / f"smoke.{extension}"
        path.write_text(source)
        compiler = "gcc" if extension == "c" else "gfortran"
        command = [compiler, "-g", "-O0"]
        if extension == "f90":
            command.append("-fcheck=bounds")
        command.extend([str(path), "-o", str(directory / "smoke")])
        subprocess.run(command, check=True, timeout=60)
        env = os.environ | {"CO_DEBUG_DIRECTORY": str(directory), "OMPI_COMM_WORLD_RANK": "0"}
        result = subprocess.run(
            ["python3", "-m", "projects.ukrmol_co.rank_debugger", str(directory / "smoke")],
            env=env,
            check=False,
            timeout=60,
        )
        record = json.loads((directory / "rank-0.json").read_text())
        log = (directory / "rank-0.gdb.log").read_text()
        results[name] = record
        (root / "engine-checks.json").write_text(json.dumps(results, indent=2) + "\n")
        if name == "normal" and not (result.returncode == 0 and record["inferior_exit_code"] == 0):
            raise RuntimeError("Normal native exit was not preserved")
        if name == "signal" and not (
            result.returncode != 0
            and record.get("observed_stop_signal") == "SIGSEGV"
            and record.get("inferior_exit_code") is None
            and "smoke.c" in log
        ):
            raise RuntimeError("First-signal stack was not captured")
        if name == "bounds" and not (
            result.returncode != 0
            and record.get("inferior_exit_code") == 2
            and "above upper bound" in log
            and "smoke.f90" in log
        ):
            raise RuntimeError("Bounds failure/source backtrace was not captured")
    results["all_checks_passed"] = True
    (root / "engine-checks.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
