"""Durable container owner for the necessary compressed CAS12 property inventory."""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import time
from pathlib import Path


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def run(capsule, archive_sha256, commit):
    root, source = capsule.parent.parent, capsule / "source"
    contract = json.loads(
        (
            source / "projects/ukrmol_co/compressed-cas12-matched-properties-contract-20261025.json"
        ).read_text()
    )
    limits = contract["resources"]
    pins = json.loads((capsule / "source-manifest.json").read_text())
    if digest(capsule / "source.tar.gz") != archive_sha256 or not all(
        digest(source / n) == w for n, w in pins.items()
    ):
        raise ValueError("Source capsule integrity failure")
    report = dict(
        source_commit=commit,
        source_archive_sha256=archive_sha256,
        owner_sha256=digest(Path(__file__)),
        owner_pid=os.getpid(),
        started_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        original_owner_exit=2,
        native_denprop_calls_max=1,
        qc_calls=0,
        target_diagonalizations=0,
        scattering_calls=0,
    )

    def save():
        temporary = capsule / "execution.partial"
        temporary.write_text(json.dumps(report, indent=2) + "\n")
        temporary.replace(capsule / "execution.json")

    report["topology"] = subprocess.check_output(["lscpu", "-e=CPU,CORE,SOCKET,NODE"], text=True)
    rows = [s.split() for s in report["topology"].splitlines()[1:]]
    chosen = [row for row in rows if 4 <= int(row[0]) <= 11]
    if len(chosen) != 8 or len({(row[1], row[2]) for row in chosen}) != 8:
        raise ValueError("Declared CPU IDs are not eight physical cores")
    available = (
        int(
            next(
                line
                for line in Path("/proc/meminfo").read_text().splitlines()
                if line.startswith("MemAvailable:")
            ).split()[1]
        )
        * 1024
    )
    report.update(
        available_ram_bytes_before=available,
        free_disk_bytes_before=shutil.disk_usage(root).free,
        prior_workloads=subprocess.check_output(
            ["docker", "ps", "--format", "{{.Names}} {{.Status}}"], text=True
        ),
    )
    if (
        available < limits["minimum_available_ram_gib"] * 2**30
        or report["free_disk_bytes_before"] < limits["minimum_free_disk_gib"] * 2**30
    ):
        raise ValueError("Insufficient available RAM/disk")
    output = capsule / "outputs"
    output.mkdir(mode=0o777)
    output.chmod(0o777)
    inside = """from pathlib import Path
import json,subprocess,sys
with Path('/out/controls.log').open('x') as log:
 control=subprocess.run([sys.executable,'-m','projects.ukrmol_co.test_compressed_cas12_properties'],stdout=log,stderr=subprocess.STDOUT,timeout=120)
Path('/out/controls.json').write_text(json.dumps({'original_exit':control.returncode})+'\\n')
if control.returncode:raise SystemExit(2)
command=[sys.executable,'-m','projects.ukrmol_co.compressed_cas12_properties','--root','/work','--source','/opt/qmodeling','--work','/out/properties','--contract','/opt/qmodeling/projects/ukrmol_co/compressed-cas12-matched-properties-contract-20261025.json']
raise SystemExit(subprocess.call(command))
"""
    name = "co-compressed-cas12-properties-20261025"
    command = [
        "docker",
        "run",
        "--name",
        name,
        "--cpuset-cpus=" + limits["physical_cpu_group"],
        "--memory=64g",
        "--memory-swap=64g",
        "--shm-size=1g",
        "--entrypoint",
        "/opt/ukrmolp/entrypoint.sh",
    ]
    for value in [
        "OMP_NUM_THREADS=1",
        "OPENBLAS_NUM_THREADS=1",
        "MKL_NUM_THREADS=1",
        "PYTHONDONTWRITEBYTECODE=1",
        "PYTHONPATH=/opt/qmodeling:/opt/qmodeling/libs/qscat",
    ]:
        command += ["-e", value]
    command += [
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
    report.update(command=command, image=contract["image"], state="launching")
    save()
    started = time.monotonic()
    try:
        with (capsule / "check.log").open("x") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            report.update(state="running", docker_client_pid=process.pid)
            save()
            try:
                process.wait(timeout=limits["owner_timeout_seconds"])
                report["timed_out"] = False
            except subprocess.TimeoutExpired:
                report["timed_out"] = True
                subprocess.run(["docker", "kill", name], check=False)
                process.wait(timeout=20)
        report["original_command_exit"] = process.returncode
        state = json.loads(subprocess.check_output(["docker", "inspect", name], text=True))[0]
        (capsule / "container-inspect.json").write_text(json.dumps(state, indent=2) + "\n")
        report.update(
            original_container_exit=state["State"]["ExitCode"],
            container_oom_killed=state["State"]["OOMKilled"],
        )
        report["original_owner_exit"] = (
            process.returncode if process.returncode in [0, 1, 2] and not report["timed_out"] else 2
        )
    finally:
        report.update(
            state="finished",
            wall_seconds=time.monotonic() - started,
            finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),
            source_integrity_after=all(digest(source / n) == w for n, w in pins.items()),
        )
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
