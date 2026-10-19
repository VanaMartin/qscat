"""Independent saved-decision replay of the rejected analytic interior prerequisite."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from projects.ukrmol_co.native_reconstruction_replay import compare, scattering
from projects.ukrmol_co.numerov_increment_study import exact_discrete_scalar, exact_free_k
from projects.ukrmol_co.numerov_interior_analytic import boundary
from projects.ukrmol_co.numerov_interior_selection import select


def verify(root):
    report = json.loads((root / "collected/result/report.json").read_text())
    execution = json.loads((root / "collected/execution.json").read_text())
    selection = json.loads((root / "selection.json").read_text())
    pins = [(root / "saved" / n, w) for n, w in selection["saved_sha256"].items()]
    pins += [
        (root / "saved/result/report.json", selection["history_report_sha256"]),
        (root / "corner-qualified.json", selection["qualified_report_sha256"]),
    ]
    manifest = json.loads((root / "collected/collection-manifest.json").read_text())
    pins += [(root / "collected" / n, w) for n, w in manifest.items()]
    for path, want in pins:
        if hashlib.sha256(path.read_bytes()).hexdigest() != want:
            raise ValueError(f"Pinned replay input differs: {path}")
    cases = [x["case"] for x in selection["peaks"]]
    inventory = select(
        root / "saved", root / "saved/result/report.json", root / "corner-qualified.json", cases
    )
    for key in [
        "new_identities",
        "reused_qualified_peaks",
        "new_unique_identities_if_all_computed",
    ]:
        if inventory[key] != selection[key]:
            raise ValueError("Pinned saved-native selection differs")
    expected = set()
    for name, steps in [
        ("compressed", [0.00625, 0.003125, 0.0015625]),
        ("equilibrium", [0.025, 0.0125, 0.00625]),
    ]:
        for coupled in [False, True]:
            expected.update((name, coupled, h, 400.0, 8) for h in steps)
            expected.update((name, coupled, steps[-1], r, 8) for r in [200.0, 800.0])
            expected.add((name, coupled, steps[-1], 400.0, 6))
    found, checks, rows = set(), {}, []
    for row in report["rows"]:
        name, coupled = row["model"], row["coupled_closed"]
        h, r, p = row["step_bohr"], row["endpoint_bohr"], row["tail_order"]
        key = name, coupled, h, r, p
        if key in found or key not in expected:
            raise ValueError("Duplicate or undeclared analytic row")
        found.add(key)
        bits, target, fine = (
            ("0x1.4eecbfb15b574p-2", 320.0, 0.0015625)
            if name == "compressed"
            else ("0x1.9df3b645a1cb3p-3", -55.0, 0.00625)
        )
        if row["energy_rydberg_hex"] != bits or row["target_tangent"] != target:
            raise ValueError("Analytic model definition differs")
        energy = float.fromhex(bits)
        q, thresholds = boundary(energy, target, coupled)
        low = exact_free_k(q, thresholds, energy, 18.0, 80)
        high = exact_free_k(q, thresholds, energy, 18.0, 110)
        with mp.workdps(80):
            precision = float(max(abs(a - b) for a, b in zip(low, high, strict=True)))
        k = np.asarray(row["k"])
        errors = compare(k, np.array(high.tolist(), float))
        for metric, value in errors.items():
            if abs(value - row["continuous_errors"][metric]) > 1e-9:
                raise ValueError("Independent continuous error differs")
        prefix = f"{name}_{coupled}"
        checks[prefix + "_oracle_precision"] = precision < 1e-20
        label = f"{prefix}_{h}_{r}_{p}"
        s = scattering(k)
        unitarity = float(np.max(abs(s.conj().T @ s - np.eye(len(s)))))
        reciprocity = float(np.max(abs(s - s.T)))
        checks[label + "_invariants"] = max(unitarity, reciprocity) <= 1e-6
        if h == fine:
            checks[label + "_continuous"] = all(v <= 1e-6 for v in errors.values())
        item = dict(
            model=name, coupled_closed=coupled, configuration=[h, r, p], continuous_errors=errors
        )
        if not coupled:
            discrete = exact_discrete_scalar(float(q[0, 0]), energy, 18.0, r, h, 80)
            refined = exact_discrete_scalar(float(q[0, 0]), energy, 18.0, r, h, 110)
            with mp.workdps(80):
                rounding = float(abs(mp.mpf(float(k[0, 0])) - refined))
                oracle_drift = float(abs(discrete - refined))
            if abs(rounding - row["discrete_error"]) > 1e-9:
                raise ValueError("Independent discrete error differs")
            checks[label + "_roundoff"] = rounding <= 1e-7
            checks[label + "_discrete_precision"] = oracle_drift < 1e-20
            item["discrete_error"] = rounding
        rows.append(item)
    if found != expected or checks != report["checks"]:
        raise ValueError("Saved coverage or independent decisions differ")
    verdict = all(checks.values())
    if verdict is not report["analytic_interior_verdict"]:
        raise ValueError("Aggregate differs")
    if (
        execution["original_command_exit"] != (0 if verdict else 1)
        or execution["original_container_exit"] != execution["original_command_exit"]
    ):
        raise ValueError("Original exits differ")
    return dict(
        saved_decision_replay_passed=True,
        analytic_interior_verdict=verdict,
        checks_passed=sum(checks.values()),
        checks_total=len(checks),
        checks=checks,
        rows=rows,
        original_exits=[execution["original_command_exit"], execution["original_container_exit"]],
        execution=execution,
        co_successor_submitted=False,
        new_co_independent_solves=0,
        new_unique_energy_identities=0,
        candidate_qualified_identities=14,
        unique_energy_coverage=132,
        uncomputed_identities=462,
        native_full_k_reconstruction_verdict=False,
        co_complex_pole_verdict=None,
        native_invocations=0,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(verify(args.root), indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
