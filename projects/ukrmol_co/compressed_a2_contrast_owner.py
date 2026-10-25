"""Persistent Sadaharu container owner for the single frozen compressed A2 contrast.

Invoked as a module from the checksum-verified source capsule, not a mutable
checkout. No native executable, image rebuild, repair or automatic successor.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def run(capsule, archive_sha256, commit):
    root = capsule.parent.parent
    source = capsule / "source"
    execution = json.loads(
        (source / "projects/ukrmol_co/compressed-a2-contrast-execution-20261025.json").read_text()
    )
    resources = execution["resources"]
    pins = json.loads((capsule / "source-manifest.json").read_text())
    if digest(capsule / "source.tar.gz") != archive_sha256 or not all(
        digest(source / name) == want for name, want in pins.items()
    ):
        raise ValueError("Source capsule integrity failure before execution")
    report = dict(
        source_commit=commit,
        source_archive_sha256=archive_sha256,
        owner_sha256=digest(Path(__file__)),
        started_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        original_owner_exit=2,
        native_calls=0,
    )
    report["topology"] = subprocess.check_output(["lscpu", "-e=CPU,CORE,SOCKET,NODE"], text=True)
    topology = [line.split() for line in report["topology"].splitlines()[1:]]
    selected = [row for row in topology if int(row[0]) in range(8)]
    if len(selected) != 8 or len({(row[1], row[2]) for row in selected}) != 8:
        raise ValueError("CPU0-7 are not eight distinct physical cores")
    report["prior_workloads"] = subprocess.check_output(
        ["docker", "ps", "--format", "{{.Names}} {{.Status}}"], text=True
    )
    available = int(
        next(
            line
            for line in Path("/proc/meminfo").read_text().splitlines()
            if line.startswith("MemAvailable:")
        ).split()[1]
    )
    report.update(
        mem_available_kib_before=available, disk_free_bytes_before=shutil.disk_usage(root).free
    )
    if available < 16 * 2**20 or report["disk_free_bytes_before"] < 2 * 2**30:
        raise ValueError("Insufficient available RAM/disk for bounded contrast")
    output = capsule / "outputs"
    output.mkdir(mode=0o777)
    output.chmod(0o777)
    image = execution["image"]
    name = "co-compressed-a2-contrast-20261025"
    inside = """from pathlib import Path
import json,subprocess,sys
controls={}
modules=['test_compressed_a2_contrast_run','test_compressed_a2_pole_contrast','test_pole_character']
for module in modules:
 with Path('/out/'+module+'.log').open('x') as log:
  process=subprocess.run([sys.executable,'-m','projects.ukrmol_co.'+module],stdout=log,stderr=subprocess.STDOUT,timeout=120)
 controls[module]=process.returncode
 Path('/out/controls.json').write_text(json.dumps(controls)+'\\n')
 if process.returncode:raise SystemExit(2)
command=[sys.executable,'-m','projects.ukrmol_co.compressed_a2_contrast_run','--root','/work','--source','/opt/qmodeling','--output','/out/check','--contract','/opt/qmodeling/projects/ukrmol_co/compressed-a2-contrast-execution-20261025.json']
raise SystemExit(subprocess.call(command))
"""
    command = [
        "docker",
        "run",
        "--name",
        name,
        "--cpuset-cpus=" + resources["cpuset_cpus"],
        "--memory=8g",
        "--memory-swap=8g",
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
        image,
        "python3",
        "-c",
        inside,
    ]
    report.update(command=command, image=image)
    started = time.monotonic()
    try:
        with (capsule / "check.log").open("x") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            try:
                process.wait(timeout=resources["combined_hard_cap_seconds"])
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
            wall_seconds=time.monotonic() - started,
            finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),
            source_integrity_after=all(digest(source / n) == w for n, w in pins.items()),
        )
        (capsule / "execution.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
    return report["original_owner_exit"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capsule", type=Path, required=True)
    parser.add_argument("--archive-sha256", required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    raise SystemExit(run(args.capsule, args.archive_sha256, args.commit))
