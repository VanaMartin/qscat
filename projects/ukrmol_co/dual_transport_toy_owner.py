"""Own one frozen mathematical toy with no CO/campaign input mount or successor."""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path

CONTRACT = "projects/ukrmol_co/dual-transport-toy-contract-20261102.json"
NAME = "co-dual-transport-toy-20261102"
SOURCE_FILES = {
    "projects/__init__.py",
    "projects/ukrmol_co/__init__.py",
    CONTRACT,
    *[
        "projects/ukrmol_co/" + n
        for n in [
            "dual_transport_toy.py",
            "dual_transport_toy_owner.py",
            "test_dual_transport_toy.py",
        ]
    ],
    "docs/physics/inverse-jost-dual-transport-control.md",
}


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def container_command(source, output, contract):
    code = """from pathlib import Path
import json,subprocess,sys
with Path('/out/controls.log').open('x') as log:
 r=subprocess.run([sys.executable,'-m','projects.ukrmol_co.test_dual_transport_toy'],
                  stdout=log,stderr=subprocess.STDOUT,timeout=60)
Path('/out/controls.json').write_text(json.dumps({'original_exit':r.returncode})+'\\n')
if r.returncode:raise SystemExit(2)
raise SystemExit(subprocess.call([
 sys.executable,'-m','projects.ukrmol_co.dual_transport_toy',
 '/source/projects/ukrmol_co/dual-transport-toy-contract-20261102.json','/out/toy']))
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


def own(capsule, archive_sha256, commit):
    source = capsule / "source"
    pins = json.loads((capsule / "source-manifest.json").read_text())
    contract = json.loads((source / CONTRACT).read_text())
    started = time.monotonic()
    record = dict(
        source_commit=commit,
        source_archive_sha256=archive_sha256,
        owner_sha256=digest(Path(__file__)),
        started_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        original_owner_exit=2,
        campaign_input_mounted=False,
        co_maps=0,
        native_calls=0,
        automatic_successor=False,
    )

    def integrity():
        actual = {p.relative_to(source).as_posix() for p in source.rglob("*") if p.is_file()}
        return (
            set(pins) == actual == SOURCE_FILES
            and not any(p.is_symlink() for p in source.rglob("*"))
            and digest(capsule / "source.tar.gz") == archive_sha256
            and all(digest(source / n) == w for n, w in pins.items())
            and digest(Path(__file__)) == pins["projects/ukrmol_co/dual_transport_toy_owner.py"]
        )

    def save():
        q = capsule / "execution.partial"
        q.write_text(json.dumps(record, indent=2) + "\n")
        q.replace(capsule / "execution.json")

    try:
        if CONTRACT not in pins or not integrity():
            raise ValueError("Changed frozen toy source")
        topology = subprocess.check_output(["lscpu", "-e=CPU,CORE,SOCKET,NODE"], text=True)
        selected = [
            row.split() for row in topology.splitlines()[1:] if int(row.split()[0]) in range(4)
        ]
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
            raise ValueError("Toy RAM/disk floor failed")
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
        command = container_command(source, output, contract)
        record.update(state="running", command=command)
        save()
        with (capsule / "command.log").open("x") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            try:
                process.wait(timeout=360)
                record["timed_out"] = False
            except subprocess.TimeoutExpired:
                record["timed_out"] = True
                subprocess.run(["docker", "kill", NAME], check=False, timeout=30)
                process.wait(timeout=30)
        state = json.loads(subprocess.check_output(["docker", "inspect", NAME], text=True))[0]
        (capsule / "container-inspect.json").write_text(json.dumps(state, indent=2) + "\n")
        record.update(
            original_command_exit=process.returncode,
            original_container_exit=state["State"]["ExitCode"],
            container_oom_killed=state["State"]["OOMKilled"],
        )
        if (
            process.returncode == state["State"]["ExitCode"] == 0
            and not state["State"]["OOMKilled"]
            and not record["timed_out"]
        ):
            record["original_owner_exit"] = 0
    except Exception as error:
        record["error"] = repr(error)
    finally:
        record.update(
            state="finished",
            wall_seconds=time.monotonic() - started,
            finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),
            source_integrity_after=integrity(),
        )
        if not record["source_integrity_after"]:
            record["original_owner_exit"] = 2
        save()
    return record["original_owner_exit"]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--capsule", type=Path, required=True)
    p.add_argument("--archive-sha256", required=True)
    p.add_argument("--commit", required=True)
    a = p.parse_args()
    raise SystemExit(own(a.capsule, a.archive_sha256, a.commit))
