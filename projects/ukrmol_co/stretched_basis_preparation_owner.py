"""Durably own one reviewed stretched B1 preparation; never launch the peer/inner job."""

from __future__ import annotations

import argparse
import datetime
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from projects.ukrmol_co.stretched_basis_preparation import (
    CONTRACT,
    authenticate,
    reuse_prerequisites,
)
from projects.ukrmol_co.stretched_basis_properties import digest, save


def name(case):
    return "co-stretched-" + case + "-preparation-20261031"


def container_command(
    root, source, output, case, contract, receipt=None, receipt_sha256=None, reuse_dz_csf=False
):
    science = ["python3", "-m", "projects.ukrmol_co.stretched_basis_preparation", "--case", case]
    if receipt is not None:
        science += ["--previous-receipt", receipt, "--previous-receipt-sha256", receipt_sha256]
    if reuse_dz_csf:
        science += ["--reuse-dz-csf"]
    inside = (
        "from pathlib import Path\nimport json,subprocess\n"
        "with Path('/out/controls.log').open('x') as f:\n"
        " p=subprocess.run(['python3','-m','projects.ukrmol_co.test_stretched_basis_preparation'],"
        "stdout=f,stderr=subprocess.STDOUT,timeout=120)\n"
        "Path('/out/controls.json').write_text(json.dumps({'original_exit':p.returncode})+'\\n')\n"
        "if p.returncode:raise SystemExit(2)\n"
        f"raise SystemExit(subprocess.call({science!r}))\n"
    )
    command = [
        "docker",
        "run",
        "--name",
        name(case),
        "--cpuset-cpus=0-3",
        "--memory=32g",
        "--memory-swap=32g",
        "--shm-size=1g",
        "--network=none",
        "--entrypoint",
        "/opt/ukrmolp/entrypoint.sh",
    ]
    for item in [
        "OMP_NUM_THREADS=1",
        "OPENBLAS_NUM_THREADS=1",
        "MKL_NUM_THREADS=1",
        "PYTHONDONTWRITEBYTECODE=1",
        "PYTHONPATH=/opt/qmodeling:/opt/qmodeling/libs/qscat",
    ]:
        command += ["-e", item]
    return [
        *command,
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


def own(
    capsule, archive_sha256, commit, case, receipt=None, receipt_sha256=None, reuse_dz_csf=False
):
    root, source = capsule.parent.parent, capsule / "source"
    contract = json.loads((source / CONTRACT).read_text())
    pins = json.loads((capsule / "source-manifest.json").read_text())
    started = time.monotonic()
    report = dict(
        case=case,
        source_commit=commit,
        source_archive_sha256=archive_sha256,
        owner_sha256=digest(Path(__file__)),
        owner_pid=os.getpid(),
        original_owner_exit=2,
        started_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        max_integral_calls=1,
        max_congen_calls=1,
        automatic_successor=False,
    )

    def source_ok():
        return digest(capsule / "source.tar.gz") == archive_sha256 and all(
            digest(source / n) == want for n, want in pins.items()
        )

    try:
        if not source_ok():
            raise ValueError("Changed frozen preparation source")
        _, _, previous = authenticate(root, source, case, contract, receipt, receipt_sha256)
        if reuse_dz_csf:
            reuse_prerequisites(root, source, case, previous)
        topology = subprocess.check_output(
            ["lscpu", "-e=CPU,CORE,SOCKET,NODE"], text=True, timeout=15
        )
        rows = [r.split() for r in topology.splitlines()[1:] if int(r.split()[0]) in range(4)]
        if len(rows) != 4 or len({(r[1], r[2]) for r in rows}) != 4:
            raise ValueError("CPU0-3 must be four distinct physical cores")
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
            topology=topology,
            available_ram_bytes_before=ram,
            free_disk_bytes_before=disk,
            prior_workloads=subprocess.check_output(
                ["docker", "ps", "--format", "{{.Names}} {{.Status}}"], text=True, timeout=15
            ),
        )
        if ram < 40 * 2**30 or disk < 20 * 2**30:
            raise ValueError("Preparation RAM/disk floor failed")
        output = capsule / "outputs"
        output.mkdir(mode=0o777)
        output.chmod(0o777)
        command = container_command(
            root, source, output, case, contract, receipt, receipt_sha256, reuse_dz_csf
        )
        report.update(state="running", command=command)
        save(capsule / "execution.json", report)
        with (capsule / "owner.log").open("x") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            try:
                process.wait(timeout=contract["resources"]["container_timeout_seconds"])
                report["timed_out"] = False
            except subprocess.TimeoutExpired:
                report["timed_out"] = True
                subprocess.run(["docker", "kill", name(case)], check=False, timeout=30)
                process.wait(timeout=30)
        state = json.loads(
            subprocess.check_output(["docker", "inspect", name(case)], text=True, timeout=30)
        )[0]
        save(capsule / "container-inspect.json", state)
        report.update(
            original_command_exit=process.returncode,
            original_container_exit=state["State"]["ExitCode"],
            container_oom_killed=state["State"]["OOMKilled"],
        )
        if (
            process.returncode == state["State"]["ExitCode"] == 0
            and not report["timed_out"]
            and not report["container_oom_killed"]
        ):
            report["original_owner_exit"] = 0
    except Exception as error:
        report["error"] = repr(error)
    finally:
        report.update(
            state="finished",
            wall_seconds=time.monotonic() - started,
            source_integrity_after=source_ok(),
            finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        )
        if not report["source_integrity_after"]:
            report["original_owner_exit"] = 2
        save(capsule / "execution.json", report)
    return report["original_owner_exit"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=["dz", "tz"], required=True)
    parser.add_argument("--capsule", type=Path, required=True)
    parser.add_argument("--archive-sha256", required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--previous-receipt")
    parser.add_argument("--previous-receipt-sha256")
    parser.add_argument("--reuse-dz-csf", action="store_true")
    parser.add_argument("--execute-reviewed", action="store_true")
    parser.add_argument("--under-service", action="store_true")
    args = parser.parse_args()
    args.capsule = args.capsule.resolve()
    if not args.execute_reviewed or (args.previous_receipt is None) != (
        args.previous_receipt_sha256 is None
    ):
        parser.error("Separate reviewed release and paired predecessor path/hash required")
    if args.under_service:
        raise SystemExit(
            own(
                args.capsule,
                args.archive_sha256,
                args.commit,
                args.case,
                args.previous_receipt,
                args.previous_receipt_sha256,
                args.reuse_dz_csf,
            )
        )
    command = [
        "systemd-run",
        "--user",
        "--unit=" + name(args.case),
        "--property=Type=exec",
        "--property=RuntimeMaxSec=3600",
        "--property=TimeoutStopSec=60",
        "--property=WorkingDirectory=" + str(args.capsule / "source"),
        "--wait",
        sys.executable,
        "-m",
        "projects.ukrmol_co.stretched_basis_preparation_owner",
        *sys.argv[1:],
        "--under-service",
    ]
    launch = dict(
        command=command,
        started_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        source_commit=args.commit,
    )
    with (args.capsule / "service-launch.json").open("x") as stream:
        stream.write(json.dumps(launch, indent=2) + "\n")
    result = subprocess.run(command, text=True, capture_output=True)
    save(
        args.capsule / "service-return.json",
        launch
        | dict(
            original_service_exit=result.returncode,
            stdout=result.stdout,
            stderr=result.stderr,
            finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        ),
    )
    raise SystemExit(result.returncode)
