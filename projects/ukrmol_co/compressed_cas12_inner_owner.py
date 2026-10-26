"""Durable one-B1 owner over the qualified immutable CAS12 scattering inputs."""

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
from projects.ukrmol_co.engine_baseline_owner import verify_engine

CONTRACT = "projects/ukrmol_co/compressed-cas12-inner-contract-20261026.json"
NAME = "co-compressed-cas12-b1-2048-20261026"


def container_command(root, source, output, contract):
    engine = root / contract["engine_source_directory"] / "build/install"
    inside = """from pathlib import Path
import json,subprocess,sys
with Path('/out/controls.log').open('x') as log:
 control=subprocess.run([sys.executable,'-m','projects.ukrmol_co.test_compressed_cas12_inner'],stdout=log,stderr=subprocess.STDOUT,timeout=120)
Path('/out/controls.json').write_text(json.dumps({'original_exit':control.returncode})+'\\n')
if control.returncode:raise SystemExit(2)
command=[sys.executable,'-m','projects.ukrmol_co.compressed_cas12_inner','--root','/work','--source','/opt/qmodeling','--work','/out/inner','--contract','/opt/qmodeling/projects/ukrmol_co/compressed-cas12-inner-contract-20261026.json']
raise SystemExit(subprocess.call(command))
"""
    return [
        "docker",
        "run",
        "--name",
        NAME,
        "--cpuset-cpus=0-11",
        "--memory=64g",
        "--memory-swap=64g",
        "--shm-size=8g",
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
        "-v",
        str(engine / "bin") + ":/opt/ukrmolp/bin.double:ro",
        "-v",
        str(engine / "lib") + ":/opt/ukrmolp/lib.double:ro",
        "-w",
        "/opt/qmodeling",
        contract["image_id"],
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
        sectors=["B1"],
        roots=2048,
        mpi_ranks=12,
        pi_pair_verdict=None,
        pole_verdict=None,
        full_sweep_released=False,
    )
    started = time.monotonic()

    def save():
        p = capsule / "execution.partial"
        p.write_text(json.dumps(report, indent=2) + "\n")
        p.replace(capsule / "execution.json")

    def source_ok():
        return digest(capsule / "source.tar.gz") == archive_sha256 and all(
            digest(source / n) == w for n, w in pins.items()
        )

    try:
        if not source_ok():
            raise ValueError("Changed frozen inner source capsule")
        report["engine_verified_before"] = verify_engine(root, contract)
        report["topology"] = subprocess.check_output(
            ["lscpu", "-e=CPU,CORE,SOCKET,NODE"], text=True
        )
        rows = [line.split() for line in report["topology"].splitlines()[1:]]
        chosen = [row for row in rows if int(row[0]) in range(12)]
        if len(chosen) != 12 or len({(row[1], row[2]) for row in chosen}) != 12:
            raise ValueError("CPU0-11 are not twelve physical cores")
        for cid in subprocess.check_output(["docker", "ps", "-q"], text=True).split():
            other = json.loads(subprocess.check_output(["docker", "inspect", cid], text=True))[0]
            if "comfyui" not in other["Name"]:
                raise ValueError("Another scientific container is active")
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
            raise ValueError("RAM/disk launch gate failed")
        output = capsule / "outputs"
        output.mkdir(mode=0o777)
        output.chmod(0o777)
        command = container_command(root, source, output, contract)
        report.update(state="launching", command=command, image=contract["image_id"])
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
        report["engine_verified_after"] = verify_engine(root, contract)
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
            source_integrity_after=source_ok(),
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
