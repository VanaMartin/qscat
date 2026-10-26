"""Persistent four-core owner for one frozen saved CAS12 B1 pole/model contrast."""

from __future__ import annotations

import argparse
import datetime
import json
import shutil
import subprocess
import time
from pathlib import Path

from projects.ukrmol_co.compressed_cas12_properties import digest

CONTRACT = "projects/ukrmol_co/compressed-cas12-pole-contract-20261026.json"
NAME = "co-compressed-cas12-pole-20261026"


def container_command(root, source, output, contract):
    inside = """from pathlib import Path
import json,subprocess,sys
with Path('/out/controls.log').open('x') as log:
 control=subprocess.run([sys.executable,'-m','projects.ukrmol_co.test_compressed_cas12_pole_contrast'],stdout=log,stderr=subprocess.STDOUT,timeout=120)
Path('/out/controls.json').write_text(json.dumps({'original_exit':control.returncode})+'\\n')
if control.returncode:raise SystemExit(2)
command=[sys.executable,'-m','projects.ukrmol_co.compressed_cas12_pole_contrast','--root','/work','--source','/opt/qmodeling','--output','/out/check','--contract','/opt/qmodeling/projects/ukrmol_co/compressed-cas12-pole-contract-20261026.json']
raise SystemExit(subprocess.call(command))
"""
    return [
        "docker",
        "run",
        "--name",
        NAME,
        "--cpuset-cpus=0-3",
        "--memory=8g",
        "--memory-swap=8g",
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
    pins = json.loads((capsule / "source-manifest.json").read_text())
    started = time.monotonic()
    record = dict(
        source_commit=commit,
        source_archive_sha256=archive_sha256,
        owner_sha256=digest(Path(__file__)),
        started_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        original_owner_exit=2,
        native_calls=0,
    )

    def save():
        p = capsule / "execution.partial"
        p.write_text(json.dumps(record, indent=2) + "\n")
        p.replace(capsule / "execution.json")

    def source_ok():
        return digest(capsule / "source.tar.gz") == archive_sha256 and all(
            digest(source / n) == w for n, w in pins.items()
        )

    try:
        if not source_ok():
            raise ValueError("Source capsule integrity failure")
        record["topology"] = subprocess.check_output(
            ["lscpu", "-e=CPU,CORE,SOCKET,NODE"], text=True
        )
        rows = [line.split() for line in record["topology"].splitlines()[1:]]
        chosen = [row for row in rows if int(row[0]) in range(4)]
        if len(chosen) != 4 or len({(row[1], row[2]) for row in chosen}) != 4:
            raise ValueError("Not four distinct physical cores")
        record["prior_workloads"] = subprocess.check_output(
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
        record.update(
            available_ram_bytes_before=ram, free_disk_bytes_before=shutil.disk_usage(root).free
        )
        if ram < 16 * 2**30 or record["free_disk_bytes_before"] < 2 * 2**30:
            raise ValueError("RAM/disk launch gate failed")
        output = capsule / "outputs"
        output.mkdir(mode=0o777)
        output.chmod(0o777)
        command = container_command(root, source, output, contract)
        record.update(state="running", command=command, image=contract["image"])
        save()
        with (capsule / "check.log").open("x") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            try:
                process.wait(
                    timeout=contract["owner_timeout_seconds"] - (time.monotonic() - started)
                )
                record["timed_out"] = False
            except subprocess.TimeoutExpired:
                record["timed_out"] = True
                subprocess.run(["docker", "kill", NAME], check=False, timeout=30)
                process.wait(timeout=30)
        record["original_command_exit"] = process.returncode
        state = json.loads(subprocess.check_output(["docker", "inspect", NAME], text=True))[0]
        (capsule / "container-inspect.json").write_text(json.dumps(state, indent=2) + "\n")
        record.update(
            original_container_exit=state["State"]["ExitCode"],
            container_oom_killed=state["State"]["OOMKilled"],
        )
        record["original_owner_exit"] = (
            process.returncode if process.returncode in [0, 1, 2] and not record["timed_out"] else 2
        )
    except Exception as error:
        record["error"] = repr(error)
    finally:
        record.update(
            state="finished",
            wall_seconds=time.monotonic() - started,
            finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),
            source_integrity_after=source_ok(),
        )
        if not record["source_integrity_after"]:
            record["original_owner_exit"] = 2
        save()
        print(json.dumps(record, indent=2))
    return record["original_owner_exit"]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--capsule", type=Path, required=True)
    p.add_argument("--archive-sha256", required=True)
    p.add_argument("--commit", required=True)
    args = p.parse_args()
    raise SystemExit(run(args.capsule, args.archive_sha256, args.commit))
