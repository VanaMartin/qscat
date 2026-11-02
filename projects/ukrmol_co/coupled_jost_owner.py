"""One isolated frozen coupled-Jost allocation; no CO inputs or successor."""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path

CONTRACT = "projects/ukrmol_co/coupled-jost-contract-20261102.json"
NAME = "co-coupled-jost-control-20261102"
SOURCE_FILES = {
    "projects/__init__.py",
    "projects/ukrmol_co/__init__.py",
    CONTRACT,
    "docs/physics/coupled-jost-boundary-control.md",
    *[
        "projects/ukrmol_co/" + name
        for name in ["coupled_jost.py", "coupled_jost_owner.py", "test_coupled_jost.py"]
    ],
}


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def command(source, output, contract):
    code = """from pathlib import Path
import json,subprocess,sys
with Path('/out/controls.log').open('x') as log:
 r=subprocess.run([sys.executable,'-m','projects.ukrmol_co.test_coupled_jost'],
                  stdout=log,stderr=subprocess.STDOUT,timeout=60)
Path('/out/controls.json').write_text(json.dumps({'original_exit':r.returncode})+'\\n')
if r.returncode:raise SystemExit(2)
raise SystemExit(subprocess.call([sys.executable,'-m','projects.ukrmol_co.coupled_jost',
 '/source/projects/ukrmol_co/coupled-jost-contract-20261102.json','/out/control']))
"""
    return [
        "docker",
        "run",
        "--name",
        NAME,
        "--network=none",
        "--cpuset-cpus=0-3",
        "--memory=2g",
        "--memory-swap=2g",
        "--entrypoint=/opt/ukrmolp/entrypoint.sh",
        "-e",
        "OMP_NUM_THREADS=1",
        "-e",
        "OPENBLAS_NUM_THREADS=1",
        "-e",
        "MKL_NUM_THREADS=1",
        "-e",
        "PYTHONDONTWRITEBYTECODE=1",
        "-e",
        "PYTHONPATH=/source",
        "-v",
        str(source) + ":/source:ro",
        "-v",
        str(output) + ":/out",
        "-w",
        "/source",
        contract["image"],
        "python3",
        "-c",
        code,
    ]


def own(capsule, archive_sha, commit):
    source = capsule / "source"
    pins = json.loads((capsule / "source-manifest.json").read_text())
    contract = json.loads((source / CONTRACT).read_text())
    started = time.monotonic()
    record = dict(
        source_commit=commit,
        source_archive_sha256=archive_sha,
        original_owner_exit=2,
        co_maps=0,
        native_calls=0,
        automatic_successor=False,
        started_utc=datetime.datetime.now(datetime.UTC).isoformat(),
    )

    def integrity():
        actual = {p.relative_to(source).as_posix() for p in source.rglob("*") if p.is_file()}
        return (
            actual == set(pins) == SOURCE_FILES
            and not any(p.is_symlink() for p in source.rglob("*"))
            and digest(capsule / "source.tar.gz") == archive_sha
            and all(digest(source / n) == pin for n, pin in pins.items())
            and digest(Path(__file__)) == pins["projects/ukrmol_co/coupled_jost_owner.py"]
        )

    def save():
        partial = capsule / "execution.partial"
        partial.write_text(json.dumps(record, indent=2) + "\n")
        partial.replace(capsule / "execution.json")

    try:
        if not integrity():
            raise ValueError("Frozen source mismatch")
        topology = subprocess.check_output(["lscpu", "-e=CPU,CORE,SOCKET,NODE"], text=True)
        selected = [row.split() for row in topology.splitlines()[1:] if int(row.split()[0]) < 4]
        if len(selected) != 4 or len({(r[1], r[2]) for r in selected}) != 4:
            raise ValueError("Not four physical cores")
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
        if ram < 4 * 2**30 or shutil.disk_usage(capsule).free < 2**30:
            raise ValueError("RAM/disk floor")
        record.update(
            topology=topology,
            available_ram_bytes=ram,
            prior_workloads=subprocess.check_output(
                ["docker", "ps", "--format", "{{.Names}} {{.Status}}"], text=True
            ),
        )
        output = capsule / "outputs"
        output.mkdir(mode=0o777)
        output.chmod(0o777)
        record.update(state="running", command=command(source, output, contract))
        save()
        with (capsule / "command.log").open("x") as log:
            process = subprocess.Popen(record["command"], stdout=log, stderr=subprocess.STDOUT)
            try:
                process.wait(timeout=contract["owner_seconds"])
                record["timed_out"] = False
            except subprocess.TimeoutExpired:
                record["timed_out"] = True
                subprocess.run(["docker", "kill", NAME], check=False, timeout=30)
                process.wait(timeout=30)
        inspected = json.loads(subprocess.check_output(["docker", "inspect", NAME], text=True))[0]
        (capsule / "container-inspect.json").write_text(json.dumps(inspected, indent=2) + "\n")
        record.update(
            original_command_exit=process.returncode,
            original_container_exit=inspected["State"]["ExitCode"],
            oom_killed=inspected["State"]["OOMKilled"],
        )
        if (
            process.returncode == inspected["State"]["ExitCode"] == 0
            and not record["timed_out"]
            and not record["oom_killed"]
        ):
            record["original_owner_exit"] = 0
    except Exception as error:
        record["error"] = repr(error)
    finally:
        record.update(
            state="finished",
            source_integrity_after=integrity(),
            owner_wall_seconds=time.monotonic() - started,
            finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        )
        if not record["source_integrity_after"]:
            record["original_owner_exit"] = 2
        save()
    return record["original_owner_exit"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capsule", type=Path, required=True)
    parser.add_argument("--archive-sha256", required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    raise SystemExit(own(args.capsule, args.archive_sha256, args.commit))
