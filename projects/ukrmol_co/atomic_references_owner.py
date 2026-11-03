"""Durable, capped owner and ten-minute status checks for the four-atom pilot."""

from __future__ import annotations

import argparse
import datetime
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from projects.ukrmol_co.atomic_references import CONTRACT, authenticate, digest, save


def container_command(source, output, contract):
    controls = "projects.ukrmol_co.test_atomic_references"
    inside = (
        "from pathlib import Path\nimport json,subprocess\n"
        "with Path('/out/controls.log').open('x') as f:\n"
        f" p=subprocess.run(['python3','-m',{controls!r}],"
        "stdout=f,stderr=subprocess.STDOUT,timeout=120)\n"
        "Path('/out/controls.json').write_text(json.dumps({'original_exit':p.returncode})+'\\n')\n"
        "if p.returncode:raise SystemExit(2)\n"
        "raise SystemExit(subprocess.call(['python3','-m',"
        "'projects.ukrmol_co.atomic_references']))\n"
    )
    return [
        "docker",
        "run",
        "--pull=never",
        "--name",
        contract["container_name"],
        "--cpuset-cpus=0-3",
        "--memory=16g",
        "--memory-swap=16g",
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
        "PYTHONPATH=/opt/qmodeling:/opt/qmodeling/libs/qscat",
        "-v",
        str(source) + ":/opt/qmodeling:ro",
        "-v",
        str(output) + ":/out",
        "-w",
        "/opt/qmodeling",
        contract["image"],
        "python3",
        "-c",
        inside,
    ]


def own(capsule, archive_sha256, commit):
    source = capsule / "source"
    c = json.loads((source / CONTRACT).read_text())
    pins = json.loads((capsule / "source-manifest.json").read_text())
    start = time.monotonic()
    report = dict(
        source_commit=commit,
        source_archive_sha256=archive_sha256,
        owner_sha256=digest(Path(__file__)),
        started_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        original_owner_exit=2,
        automatic_successor=False,
        owner_pid=os.getpid(),
    )
    process = None

    def source_ok():
        return digest(capsule / "source.tar.gz") == archive_sha256 and all(
            digest(source / n) == h for n, h in pins.items()
        )

    try:
        if not source_ok():
            raise ValueError("Changed frozen atomic source")
        authenticate(source, c)
        topology = subprocess.check_output(
            ["lscpu", "-e=CPU,CORE,SOCKET,NODE"], text=True, timeout=15
        )
        rows = [r.split() for r in topology.splitlines()[1:] if int(r.split()[0]) in range(4)]
        if len(rows) != 4 or len({(r[1], r[2]) for r in rows}) != 4:
            raise ValueError("Need four distinct physical cores0-3")
        ram = (
            int(
                next(
                    r.split()[1]
                    for r in Path("/proc/meminfo").read_text().splitlines()
                    if r.startswith("MemAvailable:")
                )
            )
            * 1024
        )
        disk = shutil.disk_usage(capsule).free
        report.update(
            topology=topology,
            available_ram_bytes_before=ram,
            free_disk_bytes_before=disk,
            prior_workloads=subprocess.check_output(
                ["docker", "ps", "--format", "{{.Names}} {{.Status}}"], text=True, timeout=15
            ),
        )
        if ram < 24 * 2**30 or disk < 10 * 2**30:
            raise ValueError("Atomic RAM/disk floor failed")
        output = capsule / "outputs"
        output.mkdir(mode=0o777)
        output.chmod(0o777)
        command = container_command(source, output, c)
        report.update(state="running", command=command)
        save(capsule / "execution.json", report)
        with (
            (capsule / "owner.log").open("x") as log,
            (capsule / "monitor.jsonl").open("x") as monitor,
        ):
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            deadline = time.monotonic() + c["resources"]["container_timeout_seconds"]
            report["timed_out"] = False
            while True:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    report["timed_out"] = True
                    subprocess.run(["docker", "kill", c["container_name"]], check=False, timeout=30)
                    process.wait(timeout=30)
                else:
                    try:
                        process.wait(
                            timeout=min(c["resources"]["observation_interval_seconds"], remaining)
                        )
                    except subprocess.TimeoutExpired:
                        pass
                row = dict(
                    observed_utc=datetime.datetime.now(datetime.UTC).isoformat(),
                    owner_elapsed_seconds=time.monotonic() - start,
                    container_command_exit=process.poll(),
                )
                if (output / "result.json").exists():
                    row["science"] = json.loads((output / "result.json").read_text())
                monitor.write(json.dumps(row) + "\n")
                monitor.flush()
                print(json.dumps(row), flush=True)
                if process.poll() is not None:
                    break
        state = json.loads(
            subprocess.check_output(
                ["docker", "inspect", c["container_name"]], text=True, timeout=30
            )
        )[0]
        save(capsule / "container-inspect.json", state)
        report.update(
            original_command_exit=process.returncode,
            original_container_exit=state["State"]["ExitCode"],
            container_oom_killed=state["State"]["OOMKilled"],
        )
        if (
            process.returncode == state["State"]["ExitCode"] == 0
            and not report["timed_out"]
            and not state["State"]["OOMKilled"]
        ):
            report["original_owner_exit"] = 0
    except Exception as error:
        report["error"] = repr(error)
    finally:
        if process is not None and process.poll() is None:
            subprocess.run(["docker", "kill", c["container_name"]], check=False, timeout=30)
            process.wait(timeout=30)
        try:
            valid = source_ok()
        except Exception:
            valid = False
        report.update(
            state="finished",
            source_integrity_after=valid,
            wall_seconds=time.monotonic() - start,
            finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        )
        if not valid:
            report["original_owner_exit"] = 2
        save(capsule / "execution.json", report)
    return report["original_owner_exit"]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--capsule", type=Path, required=True)
    p.add_argument("--archive-sha256", required=True)
    p.add_argument("--commit", required=True)
    p.add_argument("--execute-reviewed", action="store_true")
    p.add_argument("--under-service", action="store_true")
    a = p.parse_args()
    a.capsule = a.capsule.resolve()
    if not a.execute_reviewed:
        p.error("Explicit reviewed execution release required")
    if a.under_service:
        raise SystemExit(own(a.capsule, a.archive_sha256, a.commit))
    c = json.loads((a.capsule / "source" / CONTRACT).read_text())
    command = [
        "systemd-run",
        "--user",
        "--unit=" + c["container_name"],
        "--property=Type=exec",
        "--property=RuntimeMaxSec=7800",
        "--property=TimeoutStopSec=60",
        "--property=WorkingDirectory=" + str(a.capsule / "source"),
        "--wait",
        sys.executable,
        "-m",
        "projects.ukrmol_co.atomic_references_owner",
        *sys.argv[1:],
        "--under-service",
    ]
    save(a.capsule / "service-launch.json", dict(command=command, source_commit=a.commit))
    result = subprocess.run(command, capture_output=True, text=True)
    save(
        a.capsule / "service-return.json",
        dict(
            original_service_exit=result.returncode,
            stdout=result.stdout,
            stderr=result.stderr,
            finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        ),
    )
    raise SystemExit(result.returncode)
