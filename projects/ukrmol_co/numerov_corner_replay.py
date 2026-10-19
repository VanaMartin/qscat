"""Independently reconstruct finite corner decisions from saved/public bytes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from projects.ukrmol_co.native_full_reconstruction import digest, load_model
from projects.ukrmol_co.native_reconstruction_replay import compare, scattering
from projects.ukrmol_co.saved_matrices import read_k


def decisions(rows, reference, steps, budget):
    """Separate expressions: never calls the producer evaluator or propagator."""
    expected = {(h, 400, 8) for h in steps}
    expected.update((steps[2], r, 8) for r in [200, 800])
    expected.add((steps[2], 400, 6))
    keys = [tuple(row[k] for k in ["step_bohr", "endpoint_bohr", "tail_order"]) for row in rows]
    if len(keys) != 6 or set(keys) != expected:
        raise ValueError("Saved configuration coverage differs")
    k = {key: np.asarray(row["k"]) for key, row in zip(keys, rows, strict=True)}
    for key, row in zip(keys, rows, strict=True):
        if k[key].shape != reference.shape or not np.isfinite(k[key]).all():
            raise ValueError("Invalid saved K")
        s = scattering(k[key])
        stored = np.asarray(row["s_real"]) + 1j * np.asarray(row["s_imag"])
        if np.max(abs(s - stored)) > 1e-12:
            raise ValueError("Stored S differs from independent Cayley transform")
    step = [compare(k[steps[i], 400, 8], k[steps[i + 1], 400, 8]) for i in [0, 1]]
    endpoint = [compare(k[steps[2], a, 8], k[steps[2], b, 8]) for a, b in [(200, 400), (400, 800)]]
    order = compare(k[steps[2], 400, 6], k[steps[2], 400, 8])
    checks = {}
    metrics = dict(
        step_changes=step, endpoint_changes=endpoint, order_change=order, finest_endpoint_native={}
    )
    for name, comparison in [("step", step[1]), ("endpoint", endpoint[1]), ("order", order)]:
        for key, value in comparison.items():
            checks[name + "_" + key] = value <= budget
    for radius in [200, 400, 800]:
        comparison = compare(k[steps[2], radius, 8], reference)
        metrics["finest_endpoint_native"][str(radius)] = comparison
        for key, value in comparison.items():
            checks[f"native_{radius}_{key}"] = value <= budget
    checks["spatial_trend"] = step[1]["s"] <= max(step[0]["s"] / 8, 1e-8)
    matrices = [scattering(matrix) for matrix in k.values()]
    checks["unitarity"] = all(
        np.max(abs(s.conj().T @ s - np.eye(len(s)))) <= budget for s in matrices
    )
    checks["reciprocity"] = all(np.max(abs(s - s.T)) <= budget for s in matrices)
    return dict(verdict=all(checks.values()), checks=checks, metrics=metrics)


def verify(root):
    """Replay complete bundles, original exits, provenance and stage stopping verdicts."""
    report = json.loads((root / "collected/result/report.json").read_text())
    contract = report["contract"]
    execution = json.loads((root / "collected/execution.json").read_text())
    for folder, key in [
        ("raw", "input_sha256"),
        ("saved", "saved_sha256"),
        ("prior", "prior_sha256"),
        ("increment", "increment_sha256"),
    ]:
        for name, want in contract[key].items():
            if digest(root / folder / name) != want:
                raise ValueError(f"Input hash differs: {folder}/{name}")
    manifest = json.loads((root / "collected/collection-manifest.json").read_text())
    for name, want in manifest.items():
        if digest(root / "collected" / name) != want:
            raise ValueError(f"Collection hash differs: {name}")
    results, identities, new_rows, reused_rows = [], set(), 0, 0
    for stage in ["endpoint", "corners"]:
        for item in report[stage]:
            key = item["case"], item["energy_index"]
            if key in identities:
                raise ValueError("Duplicate saved identity")
            identities.add(key)
            spec = next(x for x in contract["identities"] if (x["case"], x["energy_index"]) == key)
            if spec["stage"] != stage:
                raise ValueError("Identity appears in wrong stage")
            folder = root / "collected/result" / key[0] / f"{key[1]:03d}"
            if json.loads((folder / "report.json").read_text()) != item:
                raise ValueError("Identity disk/report mismatch")
            if len(item["rows"]) != 6:
                raise ValueError("Incomplete bundle requires separate failure diagnosis")
            if [
                json.loads(line) for line in (folder / "rows.jsonl").read_text().splitlines()
            ] != item["rows"]:
                raise ValueError("Stream/report mismatch")
            model = load_model(root / "raw" / key[0])
            energy = float(model["energies"][key[1]])
            if energy.hex() != spec["energy_rydberg_hex"] or item["energy_rydberg"] != energy:
                raise ValueError("Exact energy differs")
            native_root = root / ("prior" if spec["reference"] == "64" else "saved")
            energies, reference, _ = read_k(
                native_root / "result/native" / key[0] / spec["reference"] / model["unit"]
            )
            if spec["reference"] == "64":
                if energies.shape != (1,) or energies[0] != energy:
                    raise ValueError("Native64 identity differs")
                reference = reference[0]
            else:
                if not np.array_equal(energies, model["energies"]):
                    raise ValueError("Native32 energy array differs")
                e16, k16, _ = read_k(root / "saved/result/native" / key[0] / "16" / model["unit"])
                if not np.array_equal(e16, energies):
                    raise ValueError("Native16 energy array differs")
                controls = {
                    k: v <= contract["absolute_budget"]
                    for k, v in compare(k16[key[1]], reference[key[1]]).items()
                }
                if controls != item["reference_checks"]:
                    raise ValueError("Native reference control decision differs")
                reference = reference[key[1]]
            replay = decisions(
                item["rows"], reference, spec["steps_bohr"], contract["absolute_budget"]
            )
            if replay["checks"] != item["checks"] or replay["verdict"] != item["verdict"]:
                raise ValueError("Independent case decisions differ")
            for row in item["rows"]:
                if row["reused"]:
                    previous = json.loads((root / "increment/result/report.json").read_text())
                    old = next(
                        x
                        for x in previous["rows"]
                        if x["case"] == key[0]
                        and all(
                            x[n] == row[n] for n in ["step_bohr", "endpoint_bohr", "tail_order"]
                        )
                    )
                    if any(row[n] != old[n] for n in old):
                        raise ValueError("Reused record changed")
                    reused_rows += 1
                else:
                    new_rows += 1
            results.append(dict(case=key[0], energy_index=key[1], **replay))
    for stage, verdict in report["stage_verdicts"].items():
        rows = report[stage]
        expected = [x for x in contract["identities"] if x["stage"] == stage]
        observed = (
            False
            if any(x["verdict"] is False for x in rows)
            else (
                True
                if len(rows) == len(expected) and all(x["verdict"] is True for x in rows)
                else None
            )
        )
        if observed is not verdict:
            raise ValueError("Stage decision differs")
    overall = (
        False
        if any(x["verdict"] is False for x in results)
        else (True if len(identities) == 14 else None)
    )
    if overall is not report["experiment_verdict"]:
        raise ValueError("Aggregate decision differs")
    if execution["original_command_exit"] != (0 if overall is True else 1) or (
        execution["original_command_exit"] != execution["original_container_exit"]
    ):
        raise ValueError("Original producer exits differ")
    return dict(
        saved_decision_replay_passed=True,
        experiment_verdict=overall,
        identities_verified=len(results),
        new_rows_verified=new_rows,
        reused_rows_verified=reused_rows,
        cases=results,
        execution=execution,
        original_reconstruction_verdict=False,
        original_reconstruction_exits=[1, 1],
        native_full_k_reconstruction_verdict=False,
        co_complex_pole_verdict=None,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(verify(args.root), indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
