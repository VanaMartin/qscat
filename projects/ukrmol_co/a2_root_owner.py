"""Bound one expanded-A2 root check on the exact completed audited engine."""

from __future__ import annotations

import argparse
import datetime
import json
import os
import signal
import time
import traceback
from pathlib import Path

from projects.ukrmol_co.engine_baseline_owner import engine_command, verify_engine
from projects.ukrmol_co.sparse_comparison_owner import check_freeze, preflight, save, step


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
        raise RuntimeError(f"A2 root owner received signal {signum}")

    signal.signal(signal.SIGTERM, terminate)
    output = source / "native" / contract["run"]
    try:
        preflight(source, campaign, contract, record)
        record["engine_verified"] = verify_engine(campaign, contract)
        (source / "native").mkdir()
        runtime = [
            *engine_command(source, campaign, contract, reading=True),
            "python3",
            "-m",
            "projects.ukrmol_co.engine_baseline",
            "runtime-gate",
            "--contract",
            "/opt/qmodeling/contract.json",
            "--output",
            "/outputs/runtime-gate.json",
        ]
        step(source, record, "runtime-gate", runtime, 300)
        for action, reading, timeout in (
            ("check-inputs", True, 1800),
            ("compute", False, contract["resources"]["whole_native_timeout_seconds"]),
            ("read", True, contract["resources"]["independent_reader_timeout_seconds"]),
        ):
            preflight(source, campaign, contract, record)
            command = [
                *engine_command(source, campaign, contract, reading=reading),
                "python3",
                "-m",
                "projects.ukrmol_co.a2_root_control",
                action,
                "--root",
                "/work",
                "--work",
                "/outputs/" + (contract["run"] if action != "check-inputs" else ""),
                "--contract",
                "/opt/qmodeling/contract.json",
            ]
            step(source, record, action, command, timeout)
        record["engine_verified_after_native"] = verify_engine(campaign, contract)
        record["original_controller_exit"] = 0
    except BaseException as error:
        record.update(
            original_controller_exit=1, error=repr(error), traceback=traceback.format_exc()
        )
        raise
    finally:
        result = output / "independent-root-decision.json"
        if result.exists():
            record["resolution_decision_verdict"] = json.loads(result.read_text())[
                "resolution_decision_verdict"
            ]
        record["final_frozen_file_counts_verified"] = check_freeze(source, campaign)
        record.update(
            finished_unix=time.time(), finished_utc=datetime.datetime.now(datetime.UTC).isoformat()
        )
        save(source, record)


if __name__ == "__main__":
    main()
