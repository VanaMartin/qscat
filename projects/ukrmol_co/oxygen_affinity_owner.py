"""Durable bounded owner for two saved-state oxygen orbital-convention diagnostics."""

from __future__ import annotations

import argparse
import datetime
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

from projects.ukrmol_co.atomic_references import digest, save

CONTRACT = "projects/ukrmol_co/oxygen-affinity-orbital-release-20261103.json"


def own(capsule, commit, archive_hash):
    source = capsule / "source"
    c = json.loads((source / CONTRACT).read_text())
    pins = json.loads((capsule / "source-manifest.json").read_text())
    previous = Path(c["previous_capsule"])
    start = time.monotonic()
    result = dict(
        source_commit=commit,
        source_archive_sha256=archive_hash,
        original_owner_exit=2,
        owner_pid=os.getpid(),
        automatic_successor=False,
        started_utc=datetime.datetime.now(datetime.UTC).isoformat(),
    )

    def integrity():
        return (
            digest(capsule / "source.tar.gz") == archive_hash
            and all(digest(source / n) == h for n, h in pins.items())
            and all(digest(previous / n) == h for n, h in c["input_sha256"].items())
        )

    process = None
    try:
        if not integrity():
            raise ValueError("Frozen source/previous evidence changed")
        topology = subprocess.check_output(
            ["lscpu", "-e=CPU,CORE,SOCKET,NODE"], text=True, timeout=15
        )
        rows = [r.split() for r in topology.splitlines()[1:] if int(r.split()[0]) in range(4)]
        if len(rows) != 4 or len({(r[1], r[2]) for r in rows}) != 4:
            raise ValueError("Four physical cores0-3 required")
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
        result.update(
            topology=topology,
            available_ram_bytes_before=ram,
            free_disk_bytes_before=disk,
            prior_workloads=subprocess.check_output(
                ["docker", "ps", "--format", "{{.Names}} {{.Status}}"], text=True, timeout=15
            ),
        )
        if ram < 24 * 2**30 or disk < 10 * 2**30:
            raise ValueError("RAM/disk floors failed")
        output = capsule / "outputs"
        output.mkdir(mode=0o777)
        output.chmod(0o777)
        inside = (
            "from pathlib import Path\nimport subprocess,json\n"
            "with Path('/out/controls.log').open('x') as f:\n"
            " p=subprocess.run(['python3','-m','projects.ukrmol_co.test_oxygen_affinity_orbitals'],"
            "stdout=f,stderr=subprocess.STDOUT,timeout=120)\n"
            "Path('/out/controls.json').write_text(json.dumps({'original_exit':p.returncode})+'\\n')\n"
            "if p.returncode:raise SystemExit(2)\n"
            "raise SystemExit(subprocess.call(['python3','-m',"
            "'projects.ukrmol_co.oxygen_affinity_orbitals',"
            "'--queue','--output','/out']))\n"
        )
        command = [
            "docker",
            "run",
            "--pull=never",
            "--name",
            c["container_name"],
            "--cpuset-cpus=0-3",
            "--memory=16g",
            "--memory-swap=16g",
            "--shm-size=1g",
            "--network=none",
            "--entrypoint",
            "/opt/ukrmolp/entrypoint.sh",
        ]
        for key, value in {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONPATH": "/opt/qmodeling:/opt/qmodeling/libs/qscat",
        }.items():
            command += ["-e", key + "=" + value]
        command += [
            "-v",
            str(source) + ":/opt/qmodeling:ro",
            "-v",
            str(output) + ":/out",
            "-v",
            str(previous) + ":/previous:ro",
            "-w",
            "/opt/qmodeling",
            c["image"],
            "python3",
            "-c",
            inside,
        ]
        result.update(state="running", command=command)
        save(capsule / "execution.json", result)
        with (
            (capsule / "owner.log").open("x") as log,
            (capsule / "monitor.jsonl").open("x") as monitor,
        ):
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            deadline = time.monotonic() + 1950
            result["timed_out"] = False
            while True:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    result["timed_out"] = True
                    subprocess.run(["docker", "kill", c["container_name"]], check=False, timeout=30)
                    process.wait(timeout=30)
                else:
                    try:
                        process.wait(timeout=min(600, remaining))
                    except subprocess.TimeoutExpired:
                        pass
                row = dict(
                    observed_utc=datetime.datetime.now(datetime.UTC).isoformat(),
                    elapsed_seconds=time.monotonic() - start,
                    original_command_exit=process.poll(),
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
        result.update(
            original_command_exit=process.returncode,
            original_container_exit=state["State"]["ExitCode"],
            container_oom_killed=state["State"]["OOMKilled"],
        )
        if process.returncode == state["State"]["ExitCode"] == 0 and not (
            result["timed_out"] or state["State"]["OOMKilled"]
        ):
            result["original_owner_exit"] = 0
    except Exception as error:
        result["error"] = repr(error)
    finally:
        if process is not None and process.poll() is None:
            subprocess.run(["docker", "kill", c["container_name"]], check=False, timeout=30)
            process.wait(timeout=30)
        try:
            valid = integrity()
        except Exception:
            valid = False
        result.update(
            state="finished",
            source_integrity_after=valid,
            finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),
            wall_seconds=time.monotonic() - start,
        )
        if not valid:
            result["original_owner_exit"] = 2
        save(capsule / "execution.json", result)
    return result["original_owner_exit"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capsule", type=Path, required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--archive-sha256", required=True)
    parser.add_argument("--execute-reviewed", action="store_true")
    a = parser.parse_args()
    if not a.execute_reviewed:
        parser.error("Exact reviewed release required")
    raise SystemExit(own(a.capsule.resolve(), a.commit, a.archive_sha256))
