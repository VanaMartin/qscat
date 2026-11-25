"""Capped sequential owner for two diagnostic-only neutral core pairs.

Standard-library only on the host. The frozen external manifest and a fresh
explicit allocation authorize execution, not this module or its CLI flag.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import signal
import subprocess
import time
from pathlib import Path

IMAGE = "sha256:48b9af05d30ce3c289e519433d9fc72a399b8caff6c109ea6c361e8a6b3599b6"
RADII = (2.45, 2.50)
OWNER_PATH = "projects/ukrmol_co/neutral_core_owner.py"
REQUIRED_SOURCE = (
    "projects/__init__.py",
    "projects/ukrmol_co/__init__.py",
    "projects/ukrmol_co/run.py",
    "projects/ukrmol_co/neutral_core_diagnostic.py",
    OWNER_PATH,
)


def _save(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
    temporary.replace(path)


def _output(command):
    return subprocess.check_output(command, text=True, timeout=20)


def _allocated(root):
    return int(_output(["du", "-s", "-B1", str(root)]).split()[0])


def _source_matches(root, manifest):
    try:
        if not set(REQUIRED_SOURCE).issubset(manifest["source_files"]):
            return False
        for name, digest in manifest["source_files"].items():
            path = Path(name)
            if path.is_absolute() or ".." in path.parts:
                return False
            target = root / "source" / path
            if target.is_symlink() or hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                return False
        return (
            hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
            == manifest["source_files"][OWNER_PATH]
        )
    except (OSError, KeyError, TypeError):
        return False


def host_ready(root):
    """Check physical topology, image, workload metadata and RAM/disk floors."""
    topology = _output(["lscpu", "-e=CPU,CORE,SOCKET,NODE"])
    rows = [r.split() for r in topology.splitlines()[1:] if int(r.split()[0]) in [8, 9, 10, 11]]
    if len(rows) != 4 or len({(r[1], r[2]) for r in rows}) != 4:
        raise ValueError("Four distinct physical cores8-11 required")
    memory = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
    available = int(memory["MemAvailable"].split()[0]) * 1024
    free_disk = shutil.disk_usage(root).free
    if available < 32 * 2**30 or free_disk < 40 * 2**30:
        raise ValueError("Host RAM/disk floor failed")
    image = json.loads(_output(["docker", "image", "inspect", IMAGE]))[0]
    if image["Id"] != IMAGE:
        raise ValueError("Image identity mismatch")

    def sample():
        return {
            int(r.split()[0][3:]): list(map(int, r.split()[1:]))
            for r in Path("/proc/stat").read_text().splitlines()
            if r.startswith("cpu") and r.split()[0][3:].isdigit()
        }

    first = sample()
    time.sleep(0.5)
    second = sample()
    selected_cores = {(r[1], r[2]) for r in rows}
    logical = [
        int(r.split()[0])
        for r in topology.splitlines()[1:]
        if tuple(r.split()[1:3]) in selected_cores
    ]
    loads = {}
    for cpu in logical:
        diff = [b - a for a, b in zip(first[cpu], second[cpu], strict=True)]
        loads[str(cpu)] = 1 - (diff[3] + diff[4]) / max(sum(diff[:8]), 1)
    if max(loads.values()) > 0.3:
        raise ValueError(f"Selected cores/siblings busy: {loads}")
    return {
        "topology": topology,
        "available_ram_bytes": available,
        "free_disk_bytes": free_disk,
        "sampled_logical_cpu_loads": loads,
        "image": image,
        "workloads": _output(["docker", "ps", "--format", "{{.Names}} {{.Status}} {{.Image}}"]),
    }


def container_command(root, name, radius, budget):
    """Fixed image/core/basis route; no install or parameter ladder."""
    if isinstance(radius, bool) or radius not in RADII or not 0 < budget <= 3600:
        raise ValueError("Unadmitted geometry or budget")
    if not re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9_.-]+", name):
        raise ValueError("Invalid owned container name")
    output = root / "outputs" / name
    return [
        "docker",
        "run",
        "--name",
        name,
        "--cidfile",
        str(root / f"{name}.cid"),
        "--pull=never",
        "--network=none",
        "--memory=24g",
        "--memory-swap=24g",
        "--cpuset-cpus=8-11",
        "--cpus=4",
        "--pids-limit=256",
        "--read-only",
        "--shm-size=1g",
        "--tmpfs",
        "/tmp:rw,size=64m",
        "--user",
        f"{os.getuid()}:{os.getgid()}",
        "-v",
        str(root / "source") + ":/work/source:ro",
        "-v",
        str(output) + ":/out",
        "-e",
        "PYTHONPATH=/work/source",
        "-e",
        "PYTHONDONTWRITEBYTECODE=1",
        "-e",
        "TMPDIR=/out",
        "-e",
        "OMP_NUM_THREADS=1",
        "-e",
        "OPENBLAS_NUM_THREADS=1",
        "-e",
        "MKL_NUM_THREADS=1",
        "--workdir=/work/source",
        "--entrypoint=/opt/ukrmolp/entrypoint.sh",
        IMAGE,
        "/usr/bin/timeout",
        "--kill-after=5s",
        f"{int(budget)}s",
        "python3",
        "-m",
        "projects.ukrmol_co.neutral_core_diagnostic",
        "--execute-reviewed",
        "--bond-length",
        str(radius),
        "--workdir",
        "/out/pair",
        "--ranks",
        "4",
        "--memory-mb",
        "20000",
    ]


def _kill_owned(root, name):
    path = root / f"{name}.cid"
    if path.exists():
        cid = path.read_text().strip()
        if not re.fullmatch(r"[0-9a-f]{64}", cid):
            raise ValueError("Invalid owned container ID")
        subprocess.run(["docker", "kill", cid], capture_output=True, timeout=20, check=False)


def _stop_pair(root, name, process):
    _kill_owned(root, name)
    try:
        process.wait(timeout=20)
    except subprocess.TimeoutExpired:
        process.terminate()
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait(timeout=5)


def _completed(root, name, radius):
    try:
        path = root / "outputs" / name / "pair"
        value = json.loads((path / "diagnostic-pair.json").read_text())
        resources = json.loads((path / "diagnostic-resources.json").read_text())
        if not (
            value["status"] == "numerically_complete_diagnostic"
            and value["bond_length_bohr"] == radius
            and value["basis"] == "aug-cc-pwCVQZ"
            and value["eligible_for_delivered_curve"] is False
            and all(
                value[k] is True
                for k in ["rhf_converged", "rhf_internally_stable", "rhf_externally_stable"]
            )
            and resources["original_child_exit"] == 0
            and resources["stop_reason"] is None
            and resources["source_integrity_after"] is True
        ):
            return False
        expected = {
            "/work/source/" + p: hashlib.sha256((root / "source" / p).read_bytes()).hexdigest()
            for p in ["projects/ukrmol_co/run.py", "projects/ukrmol_co/neutral_core_diagnostic.py"]
        }
        if value["source_identity"] != expected:
            return False
        for row, label, electrons in zip(
            value["core_treatments"], ["fc", "ae"], [10, 14], strict=True
        ):
            if not (
                row["core_treatment"] == label
                and row["status"] == "numerically_complete_diagnostic"
                and row["eligible_for_delivered_curve"] is False
                and row["correlated_electrons"] == electrons
                and row["rhf_reference_sha256"] == value["rhf_reference_sha256"]
                and row["ccsd_converged"] is True
                and row["lambda_converged"] is True
                and math.isfinite(row["density_electrons"])
                and abs(row["density_electrons"] - 14) < 1e-7
                and all(
                    math.isfinite(row[k])
                    for k in [
                        "total_energy_hartree",
                        "ccsd_energy_hartree",
                        "triples_correction_hartree",
                    ]
                )
                and abs(
                    row["total_energy_hartree"]
                    - row["ccsd_energy_hartree"]
                    - row["triples_correction_hartree"]
                )
                < 1e-12
            ):
                return False
        return True
    except (OSError, KeyError, ValueError, TypeError, OverflowError):
        return False


def run_pair(root, radius, name, budget, deadline):
    """Watch one container to terminal state, with no automatic numerical retry."""
    command = container_command(root, name, radius, budget)
    existing = subprocess.run(
        ["docker", "container", "inspect", name], capture_output=True, timeout=20
    )
    if existing.returncode == 0:
        raise ValueError("Owned container name already exists; no reuse")
    output = root / "outputs" / name
    output.mkdir(parents=True, exist_ok=False)
    record = {
        "radius_bohr": radius,
        "container_name": name,
        "command": command,
        "budget_seconds": budget,
        "observations": [],
    }
    start = time.monotonic()
    process = None
    try:
        with (
            (output / "container.stdout").open("x") as stdout,
            (output / "container.stderr").open("x") as stderr,
        ):
            process = subprocess.Popen(command, stdout=stdout, stderr=stderr)
            next_observation = 0.0
            while process.poll() is None:
                if time.time() >= deadline - 5 or time.monotonic() - start > budget + 10:
                    record["guard_stop"] = "wall_deadline"
                    _stop_pair(root, name, process)
                    break
                if _allocated(root) > 20 * 2**30:
                    record["guard_stop"] = "allocated_disk_cap"
                    _stop_pair(root, name, process)
                    break
                if time.monotonic() >= next_observation:
                    row = {"elapsed_seconds": time.monotonic() - start}
                    for file in ["active-stage.json", "diagnostic-pair.json"]:
                        path = output / "pair" / file
                        if path.exists():
                            value = json.loads(path.read_text())
                            row[file] = {
                                k: value[k]
                                for k in [
                                    "phase",
                                    "status",
                                    "rhf_converged",
                                    "rhf_externally_stable",
                                ]
                                if k in value
                            }
                    record["observations"].append(row)
                    _save(output / "owner-observations.json", record)
                    print(json.dumps({"radius_bohr": radius, **row}), flush=True)
                    next_observation = time.monotonic() + 120
                try:
                    process.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    pass
    except BaseException as error:
        record.update(guard_stop="pair_owner_exception", error=repr(error))
    finally:
        if process is not None:
            _stop_pair(root, name, process)
    cid_path = root / f"{name}.cid"
    record.update(
        original_command_exit=None if process is None else process.returncode,
        wall_seconds=time.monotonic() - start,
        diagnostic_complete=False,
        container_running=False,
        container_exit=None,
        OOMKilled=None,
    )
    if cid_path.exists():
        state = json.loads(_output(["docker", "inspect", cid_path.read_text().strip()]))[0]
        _save(output / "container-inspect.json", state)
        record.update(
            container_exit=state["State"]["ExitCode"],
            OOMKilled=state["State"]["OOMKilled"],
            container_running=state["State"]["Running"],
        )
        limits = state["HostConfig"]
        if (
            limits["Memory"] != 24 * 2**30
            or limits["MemorySwap"] != 24 * 2**30
            or limits["CpusetCpus"] != "8-11"
            or state["Image"] != IMAGE
        ):
            record["guard_stop"] = "container_limits_or_image_changed"
        record["diagnostic_complete"] = _completed(root, name, radius)
    _save(output / "owner-result.json", record)
    return record


def owned_running(names):
    return "".join(
        _output(["docker", "ps", "--filter", "name=^/" + name + "$", "--format", "{{.Names}}"])
        for name in names
    )


def _stop_signal(signum, frame):
    raise SystemExit(128 + signum)


def own(root, manifest_sha256):
    """Run only the two fixed pairs under a frozen manifest and absolute deadline."""
    root = Path(root).resolve()
    path = root / "execution.json"
    with path.open("x") as stream:
        json.dump({"status": "starting"}, stream)
    report = {
        "status": "stopped",
        "pairs": [],
        "automatic_successor": False,
        "started_unix": time.time(),
    }
    manifest = None
    names = []
    previous = {sig: signal.signal(sig, _stop_signal) for sig in [signal.SIGTERM, signal.SIGHUP]}
    try:
        data = (root / "source-manifest.json").read_bytes()
        if hashlib.sha256(data).hexdigest() != manifest_sha256:
            raise ValueError("Manifest digest mismatch")
        manifest = json.loads(data)
        if (
            manifest["image"] != IMAGE
            or not re.fullmatch(r"[0-9a-f]{40}", manifest["git_sha"])
            or not 0 < manifest["deadline_unix"] - manifest["received_unix"] <= 7200
            or time.time() >= manifest["deadline_unix"] - 20
            or not _source_matches(root, manifest)
        ):
            raise ValueError("Source/image/allocation admission failed")
        report.update(
            source_commit=manifest["git_sha"],
            source_manifest_sha256=manifest_sha256,
            deadline_unix=manifest["deadline_unix"],
        )
        for radius in RADII:
            if not _source_matches(root, manifest):
                raise ValueError("Frozen source changed")
            remaining = int(manifest["deadline_unix"] - time.time()) - 20
            if remaining < 60:
                raise ValueError("Insufficient remaining bounded allocation")
            readiness = host_ready(root)
            report.setdefault("host_checks", []).append(readiness)
            name = root.name + ("-r2450" if radius == 2.45 else "-r2500")
            names.append(name)
            pair = run_pair(root, radius, name, min(3600, remaining), manifest["deadline_unix"])
            report["pairs"].append(pair)
            _save(path, report)
            if (
                pair["original_command_exit"] != 0
                or pair["container_exit"] != 0
                or pair["OOMKilled"]
                or pair["container_running"]
                or not pair["diagnostic_complete"]
                or pair.get("guard_stop")
            ):
                raise ValueError("Pair failed; stop before any successor")
        report["status"] = "numerically_complete_diagnostics"
    except BaseException as error:
        report["error"] = repr(error)
    finally:
        try:
            report["source_integrity_after"] = manifest is not None and _source_matches(
                root, manifest
            )
            report["own_running_after"] = owned_running(names) if names else ""
            report["finished_unix"] = time.time()
            if not report["source_integrity_after"] or report["own_running_after"]:
                report["status"] = "stopped"
            _save(path, report)
        finally:
            for sig, handler in previous.items():
                signal.signal(sig, handler)
    return 0 if report["status"] == "numerically_complete_diagnostics" else 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capsule", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--execute-reviewed", action="store_true")
    args = parser.parse_args()
    if not args.execute_reviewed:
        parser.error("Fresh explicit reviewed data/resource admission required")
    raise SystemExit(own(args.capsule, args.manifest_sha256))


if __name__ == "__main__":
    main()
