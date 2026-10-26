"""Durably own the one frozen compressed CAS12 integral/two-CONGEN preparation."""

from __future__ import annotations

import argparse
import datetime
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

from projects.ukrmol_co.compressed_cas12_properties import digest

CONTRACT = "projects/ukrmol_co/compressed-cas12-preparation-contract-20261026.json"
NAME = "co-compressed-cas12-preparation-20261026"


def container_command(root, source, output, contract):
    """Exact mounts/resources; only the new output tree is writable."""
    inside = """from pathlib import Path
import json,subprocess,sys
with Path('/out/controls.log').open('x') as log:
 control=subprocess.run([sys.executable,'-m','projects.ukrmol_co.test_compressed_cas12_preparation'],stdout=log,stderr=subprocess.STDOUT,timeout=120)
Path('/out/controls.json').write_text(json.dumps({'original_exit':control.returncode})+'\\n')
if control.returncode:raise SystemExit(2)
command=[sys.executable,'-m','projects.ukrmol_co.compressed_cas12_preparation','--root','/work','--source','/opt/qmodeling','--work','/out/preparation','--contract','/opt/qmodeling/projects/ukrmol_co/compressed-cas12-preparation-contract-20261026.json']
raise SystemExit(subprocess.call(command))
"""
    return [
        "docker",
        "run",
        "--name",
        NAME,
        "--cpuset-cpus=0-3",
        "--memory=16g",
        "--memory-swap=16g",
        "--shm-size=1g",
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
        inside,
    ]


def run(capsule, archive_sha256, commit):
    root, source = capsule.parent.parent, capsule / "source"
    contract = json.loads((source / CONTRACT).read_text())
    limits = contract["resources"]
    pins = json.loads((capsule / "source-manifest.json").read_text())
    report = dict(
        source_commit=commit,
        source_archive_sha256=archive_sha256,
        owner_sha256=digest(Path(__file__)),
        owner_pid=os.getpid(),
        started_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        original_owner_exit=2,
        integral_calls_max=1,
        congen_calls_max=2,
        qc_calls=0,
        target_diagonalizations=0,
        property_exports=0,
        scattering_inner_solves=0,
    )
    started = time.monotonic()

    def save():
        p = capsule / "execution.partial"
        p.write_text(json.dumps(report, indent=2) + "\n")
        p.replace(capsule / "execution.json")

    def verify():
        return digest(capsule / "source.tar.gz") == archive_sha256 and all(
            digest(source / n) == w for n, w in pins.items()
        )

    try:
        if not verify():
            raise ValueError("Changed frozen preparation source")
        report["topology"] = subprocess.check_output(
            ["lscpu", "-e=CPU,CORE,SOCKET,NODE"], text=True
        )
        rows = [line.split() for line in report["topology"].splitlines()[1:]]
        chosen = [row for row in rows if int(row[0]) in range(4)]
        if len(chosen) != 4 or len({(row[1], row[2]) for row in chosen}) != 4:
            raise ValueError("CPU0-3 are not four physical cores")
        ram = (
            int(
                next(
                    line.split()[1]
                    for line in Path("/proc/meminfo").read_text().splitlines()
                    if line.startswith("MemAvailable:")
                )
            )
            * 1024
        )
        disk = shutil.disk_usage(root).free
        report.update(
            available_ram_bytes_before=ram,
            free_disk_bytes_before=disk,
            prior_workloads=subprocess.check_output(
                ["docker", "ps", "--format", "{{.Names}} {{.Status}}"], text=True
            ),
        )
        if (
            ram < limits["minimum_available_ram_gib"] * 2**30
            or disk < limits["minimum_free_disk_gib"] * 2**30
        ):
            raise ValueError("Insufficient available RAM/disk")
        output = capsule / "outputs"
        output.mkdir(mode=0o777)
        output.chmod(0o777)
        command = container_command(root, source, output, contract)
        report.update(state="launching", command=command, image=contract["image"])
        save()
        with (capsule / "check.log").open("x") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            report.update(state="running", docker_client_pid=process.pid)
            save()
            try:
                process.wait(timeout=limits["owner_timeout_seconds"] - (time.monotonic() - started))
                report["timed_out"] = False
            except subprocess.TimeoutExpired:
                report["timed_out"] = True
                subprocess.run(["docker", "kill", NAME], check=False, timeout=30)
                process.wait(timeout=30)
        report["original_command_exit"] = process.returncode
        state = json.loads(subprocess.check_output(["docker", "inspect", NAME], text=True))[0]
        (capsule / "container-inspect.json").write_text(json.dumps(state, indent=2) + "\n")
        report.update(
            original_container_exit=state["State"]["ExitCode"],
            container_oom_killed=state["State"]["OOMKilled"],
        )
        report["original_owner_exit"] = (
            process.returncode if process.returncode in [0, 1, 2] and not report["timed_out"] else 2
        )
    except Exception as error:
        report["error"] = repr(error)
    finally:
        report.update(
            state="finished",
            wall_seconds=time.monotonic() - started,
            finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),
            source_integrity_after=verify(),
        )
        if not report["source_integrity_after"]:
            report["original_owner_exit"] = 2
        save()
        print(json.dumps(report, indent=2))
    return report["original_owner_exit"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capsule", type=Path, required=True)
    parser.add_argument("--archive-sha256", required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    raise SystemExit(run(args.capsule, args.archive_sha256, args.commit))
