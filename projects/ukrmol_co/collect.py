"""Collect completed calibration batches, including failed jobs and phase differences.

Run locally with the QSCAT workspace environment after copying lightweight
artifacts from the host. Positions and full widths remain in Hartree in the
run records; comparisons are additionally labelled in eV.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from qscat.units import HARTREE_TO_EV


def collect(root: Path, pairs: list[list[str]], provenance: dict | None = None) -> dict:
    """Gather provenance and compare phases modulo pi at shared native energies."""
    batches = []
    records = {}
    for path in sorted((root / "batches").glob("*/result.json")):
        batch = json.loads(path.read_text())
        batch["name"] = path.parent.name
        for job in batch["jobs"]:
            job.pop("command", None)
        batches.append(batch)
        for job in batch["jobs"]:
            workdir = root / "runs" / job["name"]
            if not workdir.exists() and job["exit_code"] != 0:
                requested = json.loads(path.with_name("jobs.json").read_text())
                records[job["name"]] = {
                    "batch": batch["name"],
                    "initial_batch_exit_code": job["exit_code"],
                    "requested_args": next(
                        j["args"] for j in requested if j["name"] == job["name"]
                    ),
                    "status": "setup_failed",
                    "failure": "Run directory absent; see preserved batch log",
                    "refits": [],
                }
                continue
            config = json.loads((workdir / "config.json").read_text())
            resources = json.loads((workdir / "resources.json").read_text())
            record = {"config": config, "resources": resources, "batch": batch["name"]}
            target_path = workdir / "target.json"
            if target_path.exists():
                record["target"] = json.loads(target_path.read_text())
            diagnostics_path = workdir / "target-diagnostics.json"
            if diagnostics_path.exists():
                record["target_diagnostics"] = json.loads(diagnostics_path.read_text())
            record["initial_batch_exit_code"] = job["exit_code"]
            result_path = workdir / "result.json"
            if resources["exit_code"] == 0 and result_path.exists():
                record["result"] = json.loads(result_path.read_text())
                record["result"].pop("resources", None)
                record["status"] = "validated"
            else:
                record["failure"] = "See preserved run.log, geometry/log_file.0 and batch log"
                record["status"] = "engine_failed" if resources["exit_code"] else "analysis_failed"
            with (workdir / "stages.tsv").open() as stages:
                totals = {}
                for task, program, _spin, _ir, seconds, _status in csv.reader(
                    stages, delimiter="\t"
                ):
                    key = f"{task}.{program}"
                    totals[key] = totals.get(key, 0.0) + float(seconds)
                record["stage_totals_seconds"] = totals
            records[job["name"]] = record
            record["refits"] = [
                json.loads(p.read_text())
                for p in sorted((workdir / "refits").glob("*/result.json"))
            ]
    comparisons = []
    for left, right in pairs:
        lp = root / "runs" / left / "output/CO/geom1/eigenph.doublet.B1"
        rp = root / "runs" / right / "output/CO/geom1/eigenph.doublet.B1"
        a, b = np.loadtxt(lp, ndmin=2), np.loadtxt(rp, ndmin=2)
        grid = {}
        same_grid = np.array_equal(a[:, 0], b[:, 0])
        if not same_grid:
            common, ia, ib = np.intersect1d(a[:, 0], b[:, 0], return_indices=True)
            if len(common) < 2:
                raise ValueError(
                    f"Phase grids need at least two shared native energies: {left}, {right}"
                )
            grid = {
                "comparison_grid": "intersection of printed native energies; no interpolation",
                "left_energy_points": len(a),
                "right_energy_points": len(b),
                "compared_energy_points": len(common),
            }
            a, b = a[ia], b[ib]
        phase = np.remainder(a[:, 1] - b[:, 1] + np.pi / 2, np.pi) - np.pi / 2
        fits = [records[n]["result"]["native_resonances"]["B1"] for n in (left, right)]
        comparison = {
            **grid,
            "left": left,
            "right": right,
            "energy_range_ev": [float(a[0, 0]), float(a[-1, 0])],
            "max_phase_difference_rad_mod_pi": float(np.max(np.abs(phase))),
            "fit_counts": [len(f) for f in fits],
        }
        if all(len(f) == 1 for f in fits):
            e = [f[0]["energy_hartree"] * HARTREE_TO_EV for f in fits]
            g = [f[0]["width_hartree"] * HARTREE_TO_EV for f in fits]
            comparison |= {
                "position_difference_ev": abs(e[1] - e[0]),
                "width_difference_ev": abs(g[1] - g[0]),
                "width_fractional_difference_relative_to_left": abs(g[1] - g[0]) / g[0],
            }
        comparisons.append(comparison)
    return {
        **(provenance or {}),
        "hartree_to_ev": HARTREE_TO_EV,
        "status": "Calibration evidence; electronic model convergence not established",
        "batches": batches,
        "runs": records,
        "comparisons": comparisons,
    }


def main() -> None:
    """Write a reproducible JSON calibration snapshot from copied host artifacts."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--pairs", type=Path, help="JSON list of [left_name, right_name] pairs")
    parser.add_argument(
        "--provenance", type=Path, help="JSON object with campaign date/host metadata"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    pairs = json.loads(args.pairs.read_text()) if args.pairs else []
    provenance = json.loads(args.provenance.read_text()) if args.provenance else None
    result = collect(args.root, pairs, provenance)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Recorded {len(result['runs'])} jobs and {len(result['comparisons'])} comparisons")


if __name__ == "__main__":
    main()
