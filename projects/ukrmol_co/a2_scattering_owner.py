"""Bound the matched A2 control and gate expanded scattering on a dense baseline."""

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
        raise RuntimeError(f"A2 scattering owner received signal {signum}")

    signal.signal(signal.SIGTERM, terminate)
    try:
        preflight(source, campaign, contract, record)
        record["engine_verified"] = verify_engine(campaign, contract)
        (source / "native").mkdir()
        run = source / "native" / contract["run"]
        run.mkdir()
        limits = contract["resources"]
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
        actions = [
            ("check-inputs", None, True, 600),
            ("model", "dense", False, 3600),
            ("model", "five", False, limits["whole_model_timeout_seconds"]),
            ("baseline-gate", None, True, limits["independent_reader_timeout_seconds"]),
            ("model", "nine", False, limits["whole_model_timeout_seconds"]),
            ("compare", None, True, limits["independent_reader_timeout_seconds"]),
        ]
        for action, model, reading, timeout in actions:
            preflight(source, campaign, contract, record)
            command = [
                *engine_command(source, campaign, contract, reading=reading),
                "python3",
                "-m",
                "projects.ukrmol_co.a2_scattering",
                action,
                "--root",
                "/work",
                "--work",
                "/outputs/" + contract["run"],
                "--contract",
                "/opt/qmodeling/contract.json",
            ]
            if model:
                command.extend(["--model", model])
            step(source, record, action + ("-" + model if model else ""), command, timeout)
            if action in ("baseline-gate", "compare"):
                result = json.loads(
                    (
                        run
                        / (
                            "baseline-qualification.json"
                            if action == "baseline-gate"
                            else "independent-comparison.json"
                        )
                    ).read_text()
                )
                key = (
                    "dense_baseline_verdict"
                    if action == "baseline-gate"
                    else "a2_channel_sensitivity_verdict"
                )
                record[key] = result[key]
                if result[key] is not True:
                    raise ValueError(f"{key} failed or remains unset")
        record["engine_verified_after_native"] = verify_engine(campaign, contract)
        record["original_controller_exit"] = 0
    except BaseException as error:
        record.update(
            original_controller_exit=1, error=repr(error), traceback=traceback.format_exc()
        )
        raise
    finally:
        run = source / "native" / contract["run"]
        for name, key in (
            ("baseline-qualification.json", "dense_baseline_verdict"),
            ("independent-comparison.json", "a2_channel_sensitivity_verdict"),
        ):
            if (run / name).exists():
                record[key] = json.loads((run / name).read_text())[key]
        record["final_frozen_file_counts_verified"] = check_freeze(source, campaign)
        record.update(
            finished_unix=time.time(), finished_utc=datetime.datetime.now(datetime.UTC).isoformat()
        )
        save(source, record)


if __name__ == "__main__":
    main()
