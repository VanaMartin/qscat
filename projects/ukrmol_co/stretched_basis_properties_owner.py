"""Own one explicitly released stretched property export; never launch its peer."""

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

from projects.ukrmol_co.stretched_basis_properties import CONTRACT, authenticate, digest, save


def service_command(capsule, limits, arguments):
    return [
        "systemd-run",
        "--user",
        "--unit=co-stretched-basis-properties-" + arguments[1] + "-20261026",
        "--property=Type=exec",
        "--property=RuntimeMaxSec=" + str(limits["owner_timeout_seconds"]),
        "--property=TimeoutStopSec=60",
        "--property=WorkingDirectory=" + str(capsule / "source"),
        "--wait",
        sys.executable,
        "-m",
        "projects.ukrmol_co.stretched_basis_properties_owner",
        *arguments,
        "--under-service",
    ]


def container_command(root, source, output, case, contract, receipt=None, receipt_sha256=None):
    name = "co-stretched-basis-properties-" + case + "-20261026"
    controls = ["python3", "-m", "projects.ukrmol_co.test_stretched_basis_properties"]
    science = ["python3", "-m", "projects.ukrmol_co.stretched_basis_properties", "--case", case]
    if receipt is not None:
        science += ["--dz-receipt", receipt, "--dz-receipt-sha256", receipt_sha256]
    inside = (
        "from pathlib import Path\nimport json,subprocess\n"
        "with Path('/out/controls.log').open('x') as log:\n"
        f" p=subprocess.run({controls!r},stdout=log,stderr=subprocess.STDOUT,timeout=120)\n"
        "Path('/out/controls.json').write_text(json.dumps({'original_exit':p.returncode})+'\\n')\n"
        "if p.returncode:raise SystemExit(2)\n"
        f"raise SystemExit(subprocess.call({science!r}))\n"
    )
    command = [
        "docker",
        "run",
        "--name",
        name,
        "--cpuset-cpus=4-11",
        "--memory=64g",
        "--memory-swap=64g",
        "--shm-size=1g",
        "--network=none",
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


def own(capsule, archive_sha256, commit, case, receipt=None, receipt_sha256=None):
    root, source = capsule.parent.parent, capsule / "source"
    contract = json.loads((source / CONTRACT).read_text())
    pins = json.loads((capsule / "source-manifest.json").read_text())
    limits = contract["resources"]
    started = time.monotonic()
    report = dict(
        case=case,
        source_commit=commit,
        source_archive_sha256=archive_sha256,
        owner_sha256=digest(Path(__file__)),
        owner_pid=os.getpid(),
        original_owner_exit=2,
        native_denprop_calls_max=1,
        peer_job_automatically_released=False,
        started_utc=datetime.datetime.now(datetime.UTC).isoformat(),
    )

    def source_ok():
        return digest(capsule / "source.tar.gz") == archive_sha256 and all(
            digest(source / name) == want for name, want in pins.items()
        )

    try:
        if not source_ok():
            raise ValueError("Frozen source capsule differs")
        authenticate(root, source, case, contract, receipt, receipt_sha256)
        report["topology"] = subprocess.check_output(
            ["lscpu", "-e=CPU,CORE,SOCKET,NODE"], text=True, timeout=15
        )
        rows = [
            row.split()
            for row in report["topology"].splitlines()[1:]
            if 4 <= int(row.split()[0]) <= 11
        ]
        if len(rows) != 8 or len({(row[1], row[2]) for row in rows}) != 8:
            raise ValueError("CPU IDs4-11 must be eight distinct physical cores")
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
        report.update(
            available_ram_bytes_before=ram,
            free_disk_bytes_before=shutil.disk_usage(root).free,
            prior_workloads=subprocess.check_output(
                ["docker", "ps", "--format", "{{.Names}} {{.Status}}"], text=True, timeout=15
            ),
        )
        if (
            ram < limits["minimum_available_ram_gib"] * 2**30
            or report["free_disk_bytes_before"] < limits["minimum_free_disk_gib"] * 2**30
        ):
            raise ValueError("RAM/disk eligibility floor failed")
        output = capsule / "outputs"
        output.mkdir(mode=0o777)
        output.chmod(0o777)
        command = container_command(root, source, output, case, contract, receipt, receipt_sha256)
        report.update(state="running", command=command)
        save(capsule / "execution.json", report)
        with (capsule / "owner.log").open("x") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            try:
                process.wait(timeout=limits["container_timeout_seconds"])
                report["timed_out"] = False
            except subprocess.TimeoutExpired:
                report["timed_out"] = True
                subprocess.run(["docker", "kill", command[3]], check=False, timeout=30)
                process.wait(timeout=30)
        state = json.loads(
            subprocess.check_output(["docker", "inspect", command[3]], text=True, timeout=30)
        )[0]
        save(capsule / "container-inspect.json", state)
        report.update(
            original_command_exit=process.returncode,
            original_container_exit=state["State"]["ExitCode"],
            container_oom_killed=state["State"]["OOMKilled"],
        )
        report["original_owner_exit"] = (
            0
            if process.returncode == 0
            and not report["timed_out"]
            and not report["container_oom_killed"]
            else 2
        )
    except Exception as error:
        report["error"] = repr(error)
    finally:
        report.update(
            state="finished",
            source_integrity_after=source_ok(),
            wall_seconds=time.monotonic() - started,
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
    parser.add_argument("--dz-receipt")
    parser.add_argument("--dz-receipt-sha256")
    parser.add_argument("--execute-reviewed", action="store_true")
    parser.add_argument("--under-service", action="store_true")
    args = parser.parse_args()
    args.capsule = args.capsule.resolve()
    if (args.dz_receipt is None) != (args.dz_receipt_sha256 is None):
        parser.error("DZ receipt path and SHA256 must be supplied together")
    if not args.execute_reviewed:
        parser.error("Execution requires a separate reviewed release; this freeze does not launch")
    if args.under_service:
        raise SystemExit(
            own(
                args.capsule,
                args.archive_sha256,
                args.commit,
                args.case,
                args.dz_receipt,
                args.dz_receipt_sha256,
            )
        )
    contract = json.loads((args.capsule / "source" / CONTRACT).read_text())
    # Normalize the case position; no second service or model is constructed.
    arguments = [
        "--case",
        args.case,
        "--capsule",
        str(args.capsule),
        "--archive-sha256",
        args.archive_sha256,
        "--commit",
        args.commit,
        "--execute-reviewed",
    ]
    if args.dz_receipt is not None:
        arguments += [
            "--dz-receipt",
            args.dz_receipt,
            "--dz-receipt-sha256",
            args.dz_receipt_sha256,
        ]
    command = service_command(args.capsule, contract["resources"], arguments)
    launch = dict(
        command=command,
        source_commit=args.commit,
        started_utc=datetime.datetime.now(datetime.UTC).isoformat(),
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
    print(result.stdout, result.stderr)
    raise SystemExit(result.returncode)
