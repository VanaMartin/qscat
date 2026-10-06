"""Run a JSON list of CO jobs in CPU-pinned Docker slots on the remote host.

Each entry has a unique ``name`` and an ``args`` list passed to ``run.py``.
Every calculation uses a fresh container and a fresh bind-mounted run directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import queue
import re
import shlex
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def main() -> None:
    """Execute a finite batch and preserve commands, logs and batch wall time."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("jobs", type=Path)
    parser.add_argument("--root", type=Path, required=True, help="Host artifact and scratch root")
    parser.add_argument("--cpu-groups", nargs="+", default=["0-3", "4-7", "8-11"])
    parser.add_argument("--memory", default="16g")
    parser.add_argument("--image", default="qmodeling/ukrmol-co:source")
    args = parser.parse_args()
    assigned = set()
    for group in args.cpu_groups:
        cpus = set()
        try:
            for part in group.split(","):
                ends = [int(value) for value in part.split("-")]
                if len(ends) == 1:
                    cpus.add(ends[0])
                elif len(ends) == 2 and 0 <= ends[0] <= ends[1]:
                    cpus.update(range(ends[0], ends[1] + 1))
                else:
                    raise ValueError
        except ValueError:
            parser.error("CPU groups must be comma-separated nonnegative IDs or ascending ranges")
        if not cpus or min(cpus) < 0 or assigned.intersection(cpus):
            parser.error("CPU groups must be nonempty and disjoint")
        assigned.update(cpus)
    jobs = json.loads(args.jobs.read_text())
    names = [job["name"] for job in jobs]
    if len(set(names)) != len(names) or any(
        n in (".", "..") or not re.fullmatch(r"[\w.-]+", n) for n in names
    ):
        parser.error("Job names must be unique and contain only letters, digits, _, . or -")
    root = args.root.resolve()
    if any((root / "runs" / name).exists() for name in names):
        parser.error("Every job needs a fresh run directory; choose new names for repeat runs")
    logs = root / "batches" / args.jobs.stem
    logs.mkdir(parents=True, exist_ok=False)
    (logs / "jobs.json").write_text(json.dumps(jobs, indent=2) + "\n")
    image_id = subprocess.check_output(
        ["docker", "image", "inspect", "--format", "{{.Id}}", args.image], text=True
    ).strip()
    snapshot = logs / "source"
    source = snapshot / "projects/ukrmol_co"
    shutil.copytree(
        Path(__file__).parent, source, ignore=shutil.ignore_patterns("__pycache__", "._*")
    )
    shutil.copy(Path(__file__).parents[1] / "__init__.py", snapshot / "projects/__init__.py")
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source.glob("*.py")}
    hashes["co.pl"] = hashlib.sha256((source / "co.pl").read_bytes()).hexdigest()
    hashes["state_average.pm"] = hashlib.sha256(
        (source / "state_average.pm").read_bytes()
    ).hexdigest()
    pending = queue.Queue()
    for job in jobs:
        pending.put(job)

    def worker(cpus: str) -> list[dict]:
        results = []
        while True:
            try:
                job = pending.get_nowait()
            except queue.Empty:
                return results
            workdir = f"/work/runs/{job['name']}"
            calculation = [
                "python3",
                "-m",
                "projects.ukrmol_co.run",
                "--workdir",
                workdir,
                *job["args"],
            ]
            analysis = ["python3", "-m", "projects.ukrmol_co.analyze", workdir]
            command = [
                "docker",
                "run",
                "--rm",
                f"--cpuset-cpus={cpus}",
                f"--memory={args.memory}",
                f"--memory-swap={args.memory}",
                "--shm-size=1g",
                "-v",
                f"{root}:/work",
                "-v",
                f"{snapshot}:/opt/qmodeling:ro",
                "-w",
                "/opt/qmodeling",
                "--entrypoint",
                "/opt/ukrmolp/entrypoint.sh",
                image_id,
                "bash",
                "-c",
                shlex.join(calculation) + " && " + shlex.join(analysis),
            ]
            started = time.monotonic()
            with (logs / f"{job['name']}.log").open("w") as log:
                status = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
            result = {
                "name": job["name"],
                "cpus": cpus,
                "exit_code": status,
                "container_wall_seconds": time.monotonic() - started,
                "command": command,
            }
            results.append(result)
            print(
                f"{job['name']}: exit={status}, cpus={cpus}, "
                f"container wall={result['container_wall_seconds']:.2f}s",
                flush=True,
            )

    started = time.monotonic()
    with ThreadPoolExecutor(max_workers=len(args.cpu_groups)) as pool:
        results = [result for slot in pool.map(worker, args.cpu_groups) for result in slot]
    summary = {
        "wall_seconds": time.monotonic() - started,
        "image_id": image_id,
        "source_sha256": hashes,
        "cpu_groups": args.cpu_groups,
        "memory_limit_per_job": args.memory,
        "jobs": results,
    }
    (logs / "result.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Batch wall={summary['wall_seconds']:.2f}s; records: {logs / 'result.json'}", flush=True)
    raise SystemExit(int(any(job["exit_code"] for job in results)))


if __name__ == "__main__":
    main()
