"""One missing compressed CAS12/TZ QC target, using unchanged target numerics."""

from __future__ import annotations

import argparse
import datetime
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

from projects.ukrmol_co.compressed_cas12_properties import digest

CONTRACT = "projects/ukrmol_co/compressed-cas12-tz-qc-contract-20261026.json"
NAME = "co-compressed-cas12-tz-qc-20261026"


def save(path, value):
    temporary = path.with_suffix(".partial")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def qc_command(root, output, contract):
    """Use the existing QC-only CLI, with no native or import continuation."""
    args = contract["runner_args"]
    if (
        args.count("--qc-only") != 1
        or any(
            flag in args
            for flag in [
                "--target-only",
                "--workdir",
                "--scripts-cache",
                "--target-initial-checkpoint",
            ]
        )
        or args[args.index("--basis") + 1] != "cc-pVTZ"
    ):
        raise ValueError("Only the frozen missing TZ QC calculation is allowed")
    return [
        sys.executable,
        "-m",
        "projects.ukrmol_co.run",
        "--workdir",
        str(output / "run"),
        "--scripts-cache",
        str(output / "cache"),
        *args,
        "--target-initial-checkpoint",
        str(root / contract["seed_checkpoint"]),
    ]


def execute_commands(command, output, report):
    """A failed runner stops before analysis; preserve each original exit separately."""
    for name, argv in [
        ("runner", command),
        ("analyzer", [sys.executable, "-m", "projects.ukrmol_co.analyze", str(output / "run")]),
    ]:
        with (output / (name + ".log")).open("x") as log:
            status = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT).returncode
        report["original_" + name + "_exit"] = status
        save(output / "result.json", report)
        if status:
            raise ValueError(f"Original {name} failed; no retry or successor")


def inside(root, source, output):
    """Authenticate the new physical input and dispatch exactly one existing QC run."""
    contract = json.loads((source / CONTRACT).read_text())
    report = dict(
        qc_candidate_verdict=None,
        scientific_exit=2,
        qc_calls=0,
        original_runner_exit=None,
        original_analyzer_exit=None,
        native_calls=0,
        new_pole_maps=0,
        route_holdout_s_matrices_analyzed=0,
        electronic_model_verdict=None,
        full_sweep_released=False,
    )
    started = time.monotonic()

    def authenticate():
        for base, key in [(root, "input_sha256"), (source, "source_sha256")]:
            for name, want in contract[key].items():
                if digest(base / name) != want:
                    raise ValueError(f"Changed frozen prerequisite: {name}")

    try:
        authenticate()
        import pyscf

        if pyscf.__version__ != contract["pyscf_version"]:
            raise ValueError("Pinned PySCF runtime required; no dependency installation")
        report["pyscf_version"] = pyscf.__version__
        native = json.loads((root / contract["native_dz_result"]).read_text())
        model = json.loads((root / contract["completed_model_result"]).read_text())
        seed = json.loads((root / contract["seed_target"]).read_text())
        if not (
            native["inner_export_verdict"] is True
            and native["scientific_exit"] == 0
            and model["disposition"] == "reject"
            and model["scientific_exit"] == 1
            and model["original_pole_scientific_exit"] == 0
            and seed["converged"]
        ):
            raise ValueError("Qualified DZ input and resolved model-review prerequisites required")
        report["inputs_verified_before"] = len(contract["input_sha256"])
        cache = output / "cache"
        cache.mkdir()
        shutil.copy2(
            root / contract["cached_scripts_archive"], cache / "UKRmol-scripts-release-1.0.zip"
        )
        report["qc_calls"] = 1
        save(output / "result.json", report)
        execute_commands(qc_command(root, output, contract), output, report)
        result = json.loads((output / "run/result.json").read_text())
        if result.get("qc_only") is not True or result["basis"] != "cc-pVTZ":
            raise ValueError("Missing existing QC analyzer verdict")
        report.update(
            qc_candidate_verdict=True,
            scientific_exit=0,
            checkpoint_sha256=digest(output / "run/co.casscf.chk"),
            molden_sha256=digest(output / "run/co.molden"),
            validation_scope=result["validation_scope"],
        )
    except Exception as error:
        report["error"] = repr(error)
    finally:
        try:
            authenticate()
            report["inputs_verified_after"] = len(contract["input_sha256"])
        except Exception as error:
            report.update(integrity_error=repr(error), qc_candidate_verdict=None, scientific_exit=2)
        report["wall_seconds"] = time.monotonic() - started
        save(output / "result.json", report)
    return report["scientific_exit"]


def container_command(root, source, output, contract):
    inside_code = """from pathlib import Path
import json,subprocess,sys
with Path('/out/controls.log').open('x') as log:
 control=subprocess.run([sys.executable,'-m','projects.ukrmol_co.test_compressed_cas12_tz_qc'],stdout=log,stderr=subprocess.STDOUT,timeout=120)
Path('/out/controls.json').write_text(json.dumps({'original_exit':control.returncode})+'\\n')
if control.returncode:raise SystemExit(2)
command=[sys.executable,'-m','projects.ukrmol_co.compressed_cas12_tz_qc','--inside','--root','/work','--source','/opt/qmodeling','--output','/out']
raise SystemExit(subprocess.call(command))
"""
    return [
        "docker",
        "run",
        "--name",
        NAME,
        "--cpuset-cpus=0-3",
        "--memory=24g",
        "--memory-swap=24g",
        "--shm-size=1g",
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
        "QSCAT_REQUIRE_PYSCF_CONTROLS=1",
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
        "-c",
        inside_code,
    ]


def own(capsule, archive_sha256, commit):
    """Durable owner: four physical cores, bounded container, no queued successors."""
    root, source = capsule.parent.parent, capsule / "source"
    contract = json.loads((source / CONTRACT).read_text())
    pins = json.loads((capsule / "source-manifest.json").read_text())
    started = time.monotonic()
    report = dict(
        source_commit=commit,
        source_archive_sha256=archive_sha256,
        owner_sha256=digest(Path(__file__)),
        original_owner_exit=2,
        started_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        native_calls=0,
    )

    def source_ok():
        return digest(capsule / "source.tar.gz") == archive_sha256 and all(
            digest(source / n) == w for n, w in pins.items()
        )

    try:
        if not source_ok():
            raise ValueError("Frozen source capsule integrity failure")
        report["topology"] = subprocess.check_output(
            ["lscpu", "-e=CPU,CORE,SOCKET,NODE"], text=True
        )
        selected = [
            line.split()
            for line in report["topology"].splitlines()[1:]
            if int(line.split()[0]) in range(4)
        ]
        if len(selected) != 4 or len({(s[1], s[2]) for s in selected}) != 4:
            raise ValueError("Not four distinct physical cores")
        report["prior_workloads"] = subprocess.check_output(
            ["docker", "ps", "--format", "{{.Names}} {{.Status}}"], text=True
        )
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
            available_ram_bytes_before=ram, free_disk_bytes_before=shutil.disk_usage(root).free
        )
        if ram < 32 * 2**30 or report["free_disk_bytes_before"] < 8 * 2**30:
            raise ValueError("Available RAM/disk launch gate failed")
        output = capsule / "outputs"
        output.mkdir(mode=0o777)
        output.chmod(0o777)
        command = container_command(root, source, output, contract)
        report.update(state="running", command=command, image=contract["image"])
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
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inside", action="store_true")
    for name in ["root", "source", "output", "capsule"]:
        p.add_argument("--" + name, type=Path)
    p.add_argument("--archive-sha256")
    p.add_argument("--commit")
    a = p.parse_args()
    if a.inside:
        if not all([a.root, a.source, a.output]):
            p.error("Inside mode requires root/source/output")
        raise SystemExit(inside(a.root, a.source, a.output))
    if not all([a.capsule, a.archive_sha256, a.commit]):
        p.error("Owner mode requires capsule/archive-sha256/commit")
    raise SystemExit(own(a.capsule, a.archive_sha256, a.commit))
