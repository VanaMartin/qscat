"""Own one finite CAS12 root comparator and its predeclared independent reading.

Run from the fresh frozen source directory on the campaign host. The immutable
contract and frozen-inputs manifest must exist before the owner starts.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import shutil
import signal
import subprocess
import time
import traceback
from pathlib import Path


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save(source: Path, record: dict) -> None:
    temporary = source / "execution.json.tmp"
    temporary.write_text(json.dumps(record, indent=2) + "\n")
    temporary.replace(source / "execution.json")


def check_freeze(source: Path, campaign: Path) -> dict:
    frozen = json.loads((source / "frozen-inputs.json").read_text())
    for base, key in ((source, "source"), (campaign, "inputs")):
        for name, digest in frozen[key].items():
            if sha(base / name) != digest:
                raise ValueError(f"Changed frozen {key}: {name}")
    return {key: len(files) for key, files in frozen.items()}


def command(source: Path, campaign: Path, contract: dict, *, reading: bool) -> list[str]:
    resources = contract["resources"]
    name = contract["run"] + ("-reading" if reading else "")
    memory = resources["independent_reader_memory_gib" if reading else "container_memory_gib"]
    cpus = resources["independent_reader_cpu_group" if reading else "physical_cpu_group"]
    return [
        "docker",
        "run",
        "--rm",
        "--name",
        name,
        f"--cpuset-cpus={cpus}",
        f"--memory={memory}g",
        f"--memory-swap={memory}g",
        f"--shm-size={resources['shared_memory_gib']}g",
        "--entrypoint",
        "/opt/ukrmolp/entrypoint.sh",
        "-e",
        "PYTHONPATH=/opt/qmodeling",
        "-e",
        "OMP_NUM_THREADS=1",
        "-e",
        "OPENBLAS_NUM_THREADS=1",
        "-e",
        "MKL_NUM_THREADS=1",
        "-v",
        f"{campaign}:/work:ro",
        "-v",
        f"{source}:/opt/qmodeling:ro",
        "-v",
        f"{source}/native:/outputs",
        "-w",
        "/opt/qmodeling",
        contract["image_id"],
    ]


def step(source: Path, record: dict, label: str, cmd: list[str], timeout: float) -> None:
    row = {"label": label, "command": cmd, "started_unix": time.time()}
    record["steps"].append(row)
    save(source, record)
    with (source / f"{label}.launch.log").open("x") as stream:
        process = subprocess.Popen(cmd, stdout=stream, stderr=subprocess.STDOUT)
        try:
            while process.poll() is None:
                elapsed = time.time() - row["started_unix"]
                if elapsed > timeout or (source / "STOP").exists():
                    row.update(
                        timed_out=elapsed > timeout, administrative_stop=(source / "STOP").exists()
                    )
                    raise RuntimeError("Finite step stopped or timed out")
                record["heartbeat_unix"] = time.time()
                save(source, record)
                try:
                    process.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    pass
        finally:
            if process.poll() is None:
                subprocess.run(
                    ["docker", "stop", "-t", "20", cmd[cmd.index("--name") + 1]],
                    check=False,
                    timeout=45,
                )
            row.update(original_exit=process.wait(timeout=60), finished_unix=time.time())
            save(source, record)
    if row["original_exit"] != 0:
        raise RuntimeError(f"{label} failed: {row}")


def preflight(source: Path, campaign: Path, contract: dict, record: dict) -> None:
    baseline = campaign / contract["baseline_directory"]
    owner = json.loads((baseline.parent.parent / "execution.json").read_text())
    reader_name = contract.get("baseline_reader_filename", "independent-verification.json")
    verdict = contract.get("baseline_reader_verdict", "independent_raw_reading_passed")
    reader = json.loads((baseline / reader_name).read_text())
    if owner["original_controller_exit"] != 0 or reader[verdict] is not True:
        raise ValueError("Baseline owner/reader gate failed")
    record["frozen_file_counts_verified"] = check_freeze(source, campaign)
    record["cpu_topology"] = subprocess.check_output(
        ["lscpu", "-e=CPU,CORE,SOCKET,NODE"], text=True
    )
    for cpu in range(4, 12):
        siblings = (
            Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")
            .read_text()
            .strip()
        )
        if siblings != f"{cpu},{cpu + 16}":
            raise ValueError("Physical core topology changed")
    for cid in subprocess.check_output(["docker", "ps", "-q"], text=True).split():
        current = json.loads(subprocess.check_output(["docker", "inspect", cid], text=True))[0]
        if "comfyui" in current["Name"]:
            continue
        groups = current["HostConfig"]["CpusetCpus"]
        if not groups:
            raise ValueError("Another workload has unbounded CPUs")
        cpus = set()
        for group in groups.split(","):
            ends = list(map(int, group.split("-")))
            cpus.update(range(ends[0], ends[-1] + 1))
        if cpus.intersection(set(range(4, 12)) | set(range(20, 28))):
            raise ValueError("Physical comparator CPUs already allocated")
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
    disk = shutil.disk_usage(campaign).free
    record.update(available_ram_bytes=ram, free_disk_bytes=disk)
    limits = contract["resources"]
    if (
        ram < limits["minimum_available_ram_gib"] * 2**30
        or disk < limits["minimum_free_disk_gib"] * 2**30
    ):
        raise ValueError("RAM/disk launch gate failed")
    image = subprocess.check_output(
        ["docker", "image", "inspect", "--format", "{{.Id}}", contract["image_id"]], text=True
    ).strip()
    if image != contract["image_id"]:
        raise ValueError("Image pin changed")
    record["image_id"] = image
    if (source / "STOP").exists():
        raise ValueError("Administrative STOP exists")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    args = parser.parse_args()
    source, campaign = Path.cwd(), args.campaign.resolve()
    contract = json.loads((source / "contract.json").read_text())
    if (source / "execution.json").exists():
        raise ValueError("Owner output already exists")
    record = {
        "supervisor_pid": os.getpid(),
        "started_unix": time.time(),
        "steps": [],
        "scope": contract["scientific_question"],
    }
    save(source, record)

    def terminate(signum, frame):
        raise RuntimeError(f"Owner received signal {signum}")

    signal.signal(signal.SIGTERM, terminate)
    try:
        preflight(source, campaign, contract, record)
        (source / "native").mkdir()
        alias = campaign / "diagnostics" / contract["run"]
        if alias.exists() or alias.is_symlink():
            raise ValueError("Diagnostic output already exists")
        alias.symlink_to(source / "native" / contract["run"])
        limits = contract["resources"]
        native = [
            *command(source, campaign, contract, reading=False),
            "python3",
            "-m",
            "projects.ukrmol_co.sparse_pilot",
            "--preflight",
            "/work/" + contract["native_preflight"],
            "--outer-inputs",
            "/work/" + contract["outer_deck_reference"],
            "--workdir",
            "/outputs/" + contract["run"],
            "--telemetry",
            "/work/diagnostics/sparse-scattering-telemetry-v6-20261007/scattering_telemetry.so",
            "--roots",
            str(contract["candidate_roots"]),
            "--ranks",
            str(limits["mpi_ranks"]),
            "--stage-timeout",
            str(limits["native_stage_timeout_seconds"]),
        ]
        step(source, record, "native", native, limits["whole_native_timeout_seconds"])
        reading = [
            *command(source, campaign, contract, reading=True),
            "python3",
            "-m",
            "projects.ukrmol_co.sparse_comparison",
            "--reference",
            "/work/" + contract["baseline_directory"],
            "--candidate",
            "/outputs/" + contract["run"],
            "--contract",
            "/opt/qmodeling/contract.json",
            "--output",
            "/outputs/" + contract["run"] + "/independent-comparison.json",
        ]
        step(
            source,
            record,
            "independent-reading",
            reading,
            limits["independent_reader_timeout_seconds"],
        )
        record["final_frozen_file_counts_verified"] = check_freeze(source, campaign)
        record["original_controller_exit"] = 0
    except BaseException as error:
        record.update(
            original_controller_exit=1, error=repr(error), traceback=traceback.format_exc()
        )
        raise
    finally:
        record.update(
            finished_unix=time.time(),
            heartbeat_unix=time.time(),
            finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        )
        save(source, record)


if __name__ == "__main__":
    main()
