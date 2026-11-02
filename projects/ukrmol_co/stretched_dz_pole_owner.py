"""Immutable eight-physical-core owner for one bounded own-DZ pole attempt."""

from __future__ import annotations

import argparse
import datetime
import json
import shutil
import subprocess
import time
from pathlib import Path

from projects.ukrmol_co.compressed_cas12_pole_owner import container_command as base_command
from projects.ukrmol_co.compressed_cas12_properties import digest

CONTRACT = "projects/ukrmol_co/stretched-dz-pole-execution-20261102.json"


def command(root, source, output, c):
    result = base_command(root, source, output, c, CONTRACT)
    result[result.index("--cpuset-cpus=0-3")] = "--cpuset-cpus=0-7"
    return result


def run(capsule, archive_sha, commit):
    root, source = capsule.parent.parent, capsule / "source"
    c = json.loads((source / CONTRACT).read_text())
    pins = json.loads((capsule / "source-manifest.json").read_text())
    started = time.monotonic()
    record = dict(
        source_commit=commit,
        source_archive_sha256=archive_sha,
        owner_sha256=digest(Path(__file__)),
        started_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        original_owner_exit=2,
        native_calls=0,
    )

    def save():
        path = capsule / "execution.partial"
        path.write_text(json.dumps(record, indent=2) + "\n")
        path.replace(capsule / "execution.json")

    def source_ok():
        actual = {p.relative_to(source).as_posix() for p in source.rglob("*") if p.is_file()}
        return (
            digest(capsule / "source.tar.gz") == archive_sha
            and actual == set(pins)
            and not any(p.is_symlink() for p in source.rglob("*"))
            and all(digest(source / n) == want for n, want in pins.items())
            and set(c["source_sha256"]) <= set(pins)
        )

    try:
        if CONTRACT not in pins or not source_ok():
            raise ValueError("Incomplete frozen source capsule")
        topology = subprocess.check_output(["lscpu", "-e=CPU,CORE,SOCKET,NODE"], text=True)
        rows = [line.split() for line in topology.splitlines()[1:]]
        chosen = [r for r in rows if int(r[0]) in range(8)]
        if len(chosen) != 8 or len({(r[1], r[2]) for r in chosen}) != 8:
            raise ValueError("CPU0-7 are not eight distinct physical cores")
        record["topology"] = topology
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
            raise ValueError("RAM/disk launch floor")
        output = capsule / "outputs"
        output.mkdir(mode=0o777)
        output.chmod(0o777)
        cmd = command(root, source, output, c)
        record.update(state="running", command=cmd, image=c["image"])
        save()
        with (capsule / "check.log").open("x") as log:
            process = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
            try:
                process.wait(timeout=c["owner_timeout_seconds"] - (time.monotonic() - started))
                record["timed_out"] = False
            except subprocess.TimeoutExpired:
                record["timed_out"] = True
                subprocess.run(["docker", "kill", c["container_name"]], check=False, timeout=30)
                process.wait(timeout=30)
        state = json.loads(
            subprocess.check_output(["docker", "inspect", c["container_name"]], text=True)
        )[0]
        (capsule / "container-inspect.json").write_text(json.dumps(state, indent=2) + "\n")
        record.update(
            original_command_exit=process.returncode,
            original_container_exit=state["State"]["ExitCode"],
            container_oom_killed=state["State"]["OOMKilled"],
        )
        if (
            process.returncode in [0, 1, 2]
            and process.returncode == state["State"]["ExitCode"]
            and not record["timed_out"]
            and not record["container_oom_killed"]
        ):
            record["original_owner_exit"] = process.returncode
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
    return record["original_owner_exit"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capsule", type=Path, required=True)
    parser.add_argument("--archive-sha256", required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    raise SystemExit(run(args.capsule, args.archive_sha256, args.commit))
