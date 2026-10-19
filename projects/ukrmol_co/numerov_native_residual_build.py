"""Bounded locked native build, Rust checks, Python differential tests and benchmarks."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import time
import zipfile
from pathlib import Path


def main():
    output = Path("/work/build")
    output.mkdir()
    crate = "native/qscat-extended/Cargo.toml"
    commands = [
        ["rustc", "--version", "--verbose"],
        ["cargo", "--version"],
        ["cc", "--version"],
        ["cargo", "fmt", "--manifest-path", crate, "--check"],
        ["cargo", "test", "--locked", "--manifest-path", crate],
        [
            "cargo",
            "clippy",
            "--locked",
            "--manifest-path",
            crate,
            "--all-targets",
            "--",
            "-D",
            "warnings",
        ],
        [
            "/app/.venv/bin/maturin",
            "build",
            "--release",
            "--locked",
            "--manifest-path",
            crate,
            "--out",
            "/work/wheels",
            "-i",
            "/app/.venv/bin/python",
        ],
        ["cargo", "bench", "--locked", "--manifest-path", crate, "--bench", "residual"],
    ]
    report = dict(
        commands=[],
        build_passed=False,
        python=sys.version,
        compiler_environment={
            key: os.environ.get(key)
            for key in ["CC", "CFLAGS", "RUSTFLAGS", "CARGO_BUILD_JOBS", "PYO3_PYTHON"]
        },
    )

    def execute(command):
        name = f"command-{len(report['commands']):02d}.log"
        begin = time.monotonic()
        with (output / name).open("x") as stream:
            exit_code = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT).returncode
        report["commands"].append(
            dict(command=command, log=name, exit=exit_code, wall_seconds=time.monotonic() - begin)
        )
        (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        if exit_code:
            raise SystemExit(exit_code)

    for command in commands:
        execute(command)
    wheels = list(Path("/work/wheels").glob("*.whl"))
    if len(wheels) != 1:
        raise ValueError("Expected exactly one native wheel")
    report["wheel"] = dict(
        file=wheels[0].name, sha256=hashlib.sha256(wheels[0].read_bytes()).hexdigest()
    )
    with zipfile.ZipFile(wheels[0]) as archive:
        archive.extractall("/work/extension")
    os.environ["PYTHONPATH"] = "/work/extension:/work/source"
    execute(
        [
            "/app/.venv/bin/python",
            "-m",
            "pytest",
            "-q",
            "projects/ukrmol_co/test_numerov_native_residual.py",
        ]
    )
    report["build_passed"] = True
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
