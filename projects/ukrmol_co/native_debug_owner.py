"""Finite build → engine/debugger gates → one B1 diagnostic, without follow-on models."""

from __future__ import annotations

import argparse
import datetime
import json
import os
import signal
import time
import traceback
from pathlib import Path

from projects.ukrmol_co.sparse_comparison_owner import (
    check_freeze,
    command,
    preflight,
    save,
    sha,
    step,
)


def container(
    source: Path, campaign: Path, contract: dict, suffix: str, memory: str, cpus: str
) -> list[str]:
    result = command(source, campaign, contract, reading=False)
    result[result.index("--name") + 1] += suffix
    for prefix, value in (
        ("--memory=", memory),
        ("--memory-swap=", memory),
        ("--cpuset-cpus=", cpus),
    ):
        index = next(i for i, word in enumerate(result) if word.startswith(prefix))
        result[index] = prefix + value
    result[-1:-1] = ["-v", f"{source}/build:/home/ukrmol/debug-source"]
    return result


def runtime(
    source: Path, campaign: Path, contract: dict, suffix: str, memory: str, cpus: str
) -> list[str]:
    result = container(source, campaign, contract, suffix, memory, cpus)
    result[-1:-1] = [
        "--cap-add=SYS_PTRACE",
        "--security-opt=seccomp=unconfined",
        "-v",
        f"{source}/build/install/bin:/opt/ukrmolp/bin.double:ro",
        "-v",
        f"{source}/build/install/lib:/opt/ukrmolp/lib.double:ro",
    ]
    return result


def verify_build(source: Path) -> dict:
    provenance = json.loads((source / "build/build-provenance.json").read_text())
    for name, digest in provenance["installed_sha256"].items():
        if sha(source / "build/install" / name) != digest:
            raise ValueError(f"Changed diagnostic engine artifact: {name}")
    return {
        "provenance_sha256": sha(source / "build/build-provenance.json"),
        "installed_files": len(provenance["installed_sha256"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    args = parser.parse_args()
    source, campaign = Path.cwd(), args.campaign.resolve()
    if (source / "execution.json").exists():
        raise ValueError("Owner output already exists")
    declaration = json.loads((source / "contract.json").read_text())
    contract = declaration | {
        "image_id": json.loads((source / "tools-image.json").read_text())["image_id"]
    }
    record = {
        "supervisor_pid": os.getpid(),
        "started_unix": time.time(),
        "steps": [],
        "scope": declaration["scientific_question"],
    }
    save(source, record)

    def terminate(signum, frame):
        raise RuntimeError(f"Diagnostic owner received signal {signum}")

    signal.signal(signal.SIGTERM, terminate)
    try:
        preflight(source, campaign, contract, record)
        (source / "native").mkdir()
        (source / "build").mkdir()
        alias = campaign / "diagnostics" / contract["run"]
        if alias.exists() or alias.is_symlink():
            raise ValueError("Diagnostic alias already exists")
        alias.symlink_to(source / "native" / contract["run"])
        build = [
            *container(source, campaign, contract, "-build", "32g", "4-7"),
            "bash",
            "/opt/qmodeling/projects/ukrmol_co/build_debug_native.sh",
        ]
        step(source, record, "diagnostic-build", build, 7200)
        record["diagnostic_build"] = verify_build(source)
        checks = [
            *runtime(source, campaign, contract, "-checks", "16g", "0-3"),
            "python3",
            "-m",
            "projects.ukrmol_co.debug_engine_checks",
            "--build",
            "/home/ukrmol/debug-source",
        ]
        step(source, record, "engine-debugger-checks", checks, 3600)
        record["checks_sha256"] = sha(source / "build/engine-checks.json")
        if not json.loads((source / "build/engine-checks.json").read_text())["all_checks_passed"]:
            raise ValueError("Engine/debugger gates failed")
        record["frozen_file_counts_reverified_before_native"] = check_freeze(source, campaign)
        preflight(source, campaign, contract, record)
        record["diagnostic_build"] = verify_build(source)
        native = [
            *runtime(source, campaign, contract, "", "80g", "4-11"),
            "python3",
            "-m",
            "projects.ukrmol_co.debug_scatci",
            "--reference",
            "/work/" + contract["failed_B1_directory"],
            "--workdir",
            "/outputs/" + contract["run"],
            "--ranks",
            "8",
            "--timeout",
            "86400",
        ]
        step(source, record, "native-B1-diagnostic", native, 90000)
        record["original_controller_exit"] = 0
    except BaseException as error:
        record.update(
            original_controller_exit=1, error=repr(error), traceback=traceback.format_exc()
        )
        raise
    finally:
        record["final_frozen_file_counts_verified"] = check_freeze(source, campaign)
        record.update(
            finished_unix=time.time(),
            heartbeat_unix=time.time(),
            finished_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        )
        save(source, record)


if __name__ == "__main__":
    main()
