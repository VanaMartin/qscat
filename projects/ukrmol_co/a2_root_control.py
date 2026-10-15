"""One model-specific root check to test persistence of the measured A2 phase rejection."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np

from projects.ukrmol_co.a2_channel_preflight import stage
from projects.ukrmol_co.a2_scattering import (
    input_gate,
    phases,
    save,
    verify_model_inputs,
)
from projects.ukrmol_co.engine_baseline import audit_root_pair
from projects.ukrmol_co.run import run_profiled
from projects.ukrmol_co.sparse_comparison import compare_phases, read_run
from projects.ukrmol_co.sparse_pilot import audit_sector


def refined_deck(text: str) -> str:
    """Only the actual2048-root count may change in the completed nine-A2 deck."""
    counts = re.findall(r"\bnstat\s*=\s*(\d+)", text)
    if counts != ["2048"]:
        raise ValueError("Expected exactly one2048-root selection")
    return re.sub(r"(\bnstat\s*=\s*)2048\b", r"\g<1>4096", text)


def check_reference(root: Path, contract: dict) -> dict:
    """Accept completed native evidence while preserving its scientific rejection."""
    source = root / contract["a2_source_directory"]
    original = json.loads((source / "contract.json").read_text())
    work = source / "native" / original["run"]
    owner = json.loads((source / "execution.json").read_text())
    required = [
        "runtime-gate",
        "check-inputs",
        "model-dense",
        "model-five",
        "baseline-gate",
        "model-nine",
        "compare",
    ]
    if (
        not owner.get("finished_utc")
        or owner["original_controller_exit"] != 1
        or [(s["label"], s["original_exit"]) for s in owner["steps"]]
        != [(name, int(name == "compare")) for name in required]
        or any(s.get("timed_out") or s.get("administrative_stop") for s in owner["steps"])
        or owner["dense_baseline_verdict"] is not True
        or owner["a2_channel_sensitivity_verdict"] is not False
    ):
        raise ValueError("Expected completed native runs with preserved phase rejection")
    result = json.loads((work / "independent-comparison.json").read_text())
    if (
        result["a2_channel_sensitivity_verdict"] is not False
        or result["independent_raw_reading_passed"] is not True
    ):
        raise ValueError("Original channel rejection prerequisite failed")
    inventory = input_gate(root, original)
    verify_model_inputs(root, work, original, inventory, "nine")
    verify_model_inputs(root, work, original, inventory, "five")
    _, nine = read_run(work / "nine", 2048, contract)
    _, five = read_run(work / "five", 2048, original | original["models"]["five"])
    for sector in contract["sectors"]:
        comparison = compare_phases(five[sector], nine[sector], original)
        if comparison["phase_budget_passed"] or not all(
            row["passed"] is True for row in comparison["fixed_window_fits"]
        ):
            raise ValueError("Named phase-only discrepancy no longer reproduces")
    return {"original_phase_rejection_preserved": True, "original_native_runs_passed": True}


def compute(root: Path, work: Path, contract: dict) -> None:
    check_reference(root, contract)
    reference = root / contract["root_reference_directory"]
    original = json.loads((reference / "pilot.json").read_text())
    report = {"target_map": original["target_map"], "roots": 4096, "sectors": {}}
    for sector in contract["sectors"]:
        directory = work / sector
        directory.mkdir()
        left = reference / sector
        text = (left / "scatci.inp").read_text()
        unit = int(re.search(r"megul\s*=\s*(\d+)", text)[1])
        for name in (
            "fort.16",
            "fort.24",
            "fort.26",
            f"fort.{unit}",
            "rsolve.inp",
            "eigenp.inp",
            "reson.inp",
        ):
            shutil.copy2(left / name, directory / name)
        (directory / "scatci.inp").write_text(refined_deck(text))
        (directory / "reson_message").write_text("")
        native = {"stages": []}
        for name in ("scatci", "rsolve", "eigenp", "reson"):
            env = os.environ.copy()
            command = [f"/opt/ukrmolp/bin.double/{name}"]
            if name == "scatci":
                command = [
                    "mpiexec",
                    "-n",
                    str(contract["resources"]["mpi_ranks"]),
                    "/opt/ukrmolp/bin.double/mpi-scatci",
                    "scatci.inp",
                ]
                env.update(
                    LD_PRELOAD=str(
                        root / contract["telemetry_directory"] / "scattering_telemetry.so"
                    ),
                    UKRMOL_EPS_AUDIT=str(directory / "eps-telemetry.jsonl"),
                )
            stage(
                directory,
                native,
                name,
                command,
                directory / f"{name}.inp",
                contract["resources"]["native_stage_timeout_seconds"],
                env,
            )
        native.update(audit_sector(directory, 4096, 27818, 1682, 298))
        report["sectors"][sector] = native
        save(work / "pilot.json", report)
    pi = float(
        np.max(
            abs(
                (phases(work / "B1")[:, 1] - phases(work / "B2")[:, 1] + np.pi / 2) % np.pi
                - np.pi / 2
            )
        )
    )
    report.update(all_stages_complete=True, pi_phase_difference_rad_mod_pi=pi)
    save(work / "pilot.json", report)
    if pi > contract["raw_checks"]["pi_phase_rad_mod_pi"]:
        raise ValueError("Native Pi symmetry failed")


def decision(root: Path, work: Path, contract: dict) -> dict:
    check_reference(root, contract)
    pair = audit_root_pair(root / contract["root_reference_directory"], work, contract)
    five = root / contract["five_reference_directory"]
    original = json.loads((root / contract["a2_source_directory"] / "contract.json").read_text())
    _, five_phases = read_run(five, 2048, original | original["models"]["five"])
    cross, refinement = {}, {}
    for sector in contract["sectors"]:
        cross[sector] = compare_phases(five_phases[sector], phases(work / sector), original)
        cross[sector]["observable_budget_verdict"] = cross[sector].pop(
            "root_pair_stability_verdict"
        )
        refinement[sector] = pair["sectors"][sector]["maximum_phase_difference_rad_mod_pi"]
    return assess_decision(pair, cross, refinement, contract["attribution_phase_budget_rad"])


def assess_decision(pair: dict, cross: dict, refinement: dict, attribution_budget: float) -> dict:
    """A stable refined rejection answers the named question; other outcomes stop unset/false."""
    valid = pair["root_pair_stability_verdict"] is not None and all(
        row["observable_budget_verdict"] is not None for row in cross.values()
    )
    established = bool(
        pair["root_pair_stability_verdict"] is True
        and all(value <= attribution_budget for value in refinement.values())
        and all(row["phase_budget_passed"] is False for row in cross.values())
    )
    return {
        "resolution_decision_verdict": established if valid else None,
        "original_a2_channel_sensitivity_verdict": False,
        "root_pair_reading": pair,
        "refined_channel_comparison": cross,
        "attribution_phase_budget_rad": attribution_budget,
        "spectral_truncation_verdict": None,
        "electronic_model_verdict": None,
        "global_continuity_verdict": None,
        "scope": (
            "Persistence of channel phase rejection under nine-A2 root refinement; "
            "no absolute omission bound or full channel convergence"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("check-inputs", "compute", "read"))
    for name in ("root", "work", "contract"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--profiled", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    root, work = args.root.resolve(), args.work.resolve()
    contract = json.loads(args.contract.read_text())
    if args.action == "check-inputs":
        save(work / "reference-gate.json", check_reference(root, contract))
    elif args.action == "compute":
        if args.profiled:
            compute(root, work, contract)
        else:
            work.mkdir()
            (work / "scratch").mkdir()
            env = os.environ | {
                "TMPDIR": str(work / "scratch"),
                "PSI_SCRATCH": str(work / "scratch"),
            }
            raise SystemExit(
                run_profiled(
                    work,
                    env,
                    [
                        "python3",
                        "-m",
                        "projects.ukrmol_co.a2_root_control",
                        *sys.argv[1:],
                        "--profiled",
                    ],
                )
            )
    else:
        result = decision(root, work, contract)
        save(work / "independent-root-decision.json", result)
        if result["resolution_decision_verdict"] is not True:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
