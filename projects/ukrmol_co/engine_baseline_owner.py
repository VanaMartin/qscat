"""Own one bounded baseline or root-pair control on the pinned completed audit engine."""

from __future__ import annotations

import argparse
import datetime
import json
import os
import signal
import time
import traceback
from pathlib import Path

from projects.ukrmol_co.native_debug_owner import verify_build
from projects.ukrmol_co.sparse_comparison_owner import (
    check_freeze,
    command,
    preflight,
    save,
    sha,
    step,
)


def verify_engine(campaign: Path, contract: dict) -> dict:
    """Recheck the completed diagnostic, installed binaries and immutable evidence."""
    engine = campaign / contract["engine_source_directory"]
    for name, digest in contract["engine_evidence_sha256"].items():
        if sha(engine / name) != digest:
            raise ValueError(f"Changed engine evidence: {name}")
    owner = json.loads((engine / "execution.json").read_text())
    checks = json.loads((engine / "build/engine-checks.json").read_text())
    if owner["original_controller_exit"] != 0 or not checks["all_checks_passed"]:
        raise ValueError("Completed engine owner/check gate failed")
    return verify_build(engine) | {"frozen_counts": check_freeze(engine, campaign)}


def engine_command(source: Path, campaign: Path, contract: dict, *, reading: bool) -> list[str]:
    """Mount the completed audit engine read-only; retain the scientific MPI driver."""
    result = command(source, campaign, contract, reading=reading)
    engine = campaign / contract["engine_source_directory"] / "build/install"
    result[-1:-1] = [
        "-v",
        f"{engine}/bin:/opt/ukrmolp/bin.double:ro",
        "-v",
        f"{engine}/lib:/opt/ukrmolp/lib.double:ro",
    ]
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    args = parser.parse_args()
    source, campaign = Path.cwd(), args.campaign.resolve()
    if (source / "execution.json").exists():
        raise ValueError("Owner output already exists")
    contract = json.loads((source / "contract.json").read_text())
    kind = contract.get("control_kind", "matched-baseline")
    if kind not in ("matched-baseline", "root-comparison"):
        raise ValueError("Unknown audited-engine control kind")
    root_comparison = kind == "root-comparison"
    native_label = "native-root-comparison" if root_comparison else "native-matched-baseline"
    reader_label = "independent-root-reading" if root_comparison else "independent-engine-reading"
    result_name = (
        "independent-comparison.json" if root_comparison else "independent-engine-baseline.json"
    )
    verdict_key = (
        "root_pair_stability_verdict" if root_comparison else "matched_engine_baseline_verdict"
    )
    record = {
        "supervisor_pid": os.getpid(),
        "started_unix": time.time(),
        "steps": [],
        "scope": contract["scientific_question"],
    }
    save(source, record)

    def terminate(signum, frame):
        raise RuntimeError(f"Engine baseline owner received signal {signum}")

    signal.signal(signal.SIGTERM, terminate)
    try:
        preflight(source, campaign, contract, record)
        record["engine_verified"] = verify_engine(campaign, contract)
        (source / "native").mkdir()
        alias = campaign / "diagnostics" / contract["run"]
        if alias.exists() or alias.is_symlink():
            raise ValueError("Baseline alias already exists")
        alias.symlink_to(source / "native" / contract["run"])
        limits = contract["resources"]
        gate = [
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
        step(source, record, "runtime-gate", gate, 300)
        self_read = [
            *engine_command(source, campaign, contract, reading=True),
            "python3",
            "-m",
            "projects.ukrmol_co.engine_baseline",
            "compare",
            "--reference",
            "/work/" + contract["baseline_directory"],
            "--candidate",
            "/work/" + contract["baseline_directory"],
            "--contract",
            (
                "/work/" + contract["baseline_self_read_contract"]
                if root_comparison
                else "/opt/qmodeling/contract.json"
            ),
            "--output",
            "/outputs/baseline-self-reading.json",
        ]
        step(source, record, "baseline-self-reading", self_read, 3600)
        preflight(source, campaign, contract, record)
        record["engine_verified_before_native"] = verify_engine(campaign, contract)
        native = [
            *engine_command(source, campaign, contract, reading=False),
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
            "/work/" + contract["telemetry_directory"] + "/scattering_telemetry.so",
            "--roots",
            str(contract["candidate_roots"]),
            "--ranks",
            str(limits["mpi_ranks"]),
            "--stage-timeout",
            str(limits["native_stage_timeout_seconds"]),
        ]
        step(
            source,
            record,
            native_label,
            native,
            limits["whole_native_timeout_seconds"],
        )
        reading = [
            *engine_command(source, campaign, contract, reading=True),
            "python3",
            "-m",
            "projects.ukrmol_co.engine_baseline",
            "root-compare" if root_comparison else "compare",
            "--reference",
            "/work/" + contract["baseline_directory"],
            "--candidate",
            "/outputs/" + contract["run"],
            "--contract",
            "/opt/qmodeling/contract.json",
            "--output",
            "/outputs/" + contract["run"] + "/" + result_name,
        ]
        step(
            source,
            record,
            reader_label,
            reading,
            limits["independent_reader_timeout_seconds"],
        )
        result = source / "native" / contract["run"] / result_name
        record[verdict_key] = json.loads(result.read_text())[verdict_key]
        if record[verdict_key] is not True:
            raise ValueError("Audited-engine control verdict failed or remains unset")
        record["engine_verified_after_native"] = verify_engine(campaign, contract)
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
