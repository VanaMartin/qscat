"""Independently express saved precision-attribution decisions without propagation."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from projects.ukrmol_co.native_reconstruction_replay import compare, scattering
from projects.ukrmol_co.numerov_discrete_free import discrete_free_k
from projects.ukrmol_co.numerov_increment_study import exact_discrete_scalar, exact_free_k


def verify(root):
    report = json.loads((root / "collected/result/report.json").read_text())
    execution = json.loads((root / "collected/execution.json").read_text())
    contract = report["contract"]
    prior = root / "prior-report.json"
    if hashlib.sha256(prior.read_bytes()).hexdigest() != contract["prior_report_sha256"]:
        raise ValueError("Prior report hash differs")
    for name, want in json.loads((root / "collected/collection-manifest.json").read_text()).items():
        if hashlib.sha256((root / "collected" / name).read_bytes()).hexdigest() != want:
            raise ValueError("Collected source or result hash differs")
    if (
        len(report["baseline"]) != 12
        or len(report["reexecutions"]) != 2
        or len(report["variants"]) != 8
    ):
        raise ValueError("Attribution coverage is incomplete")
    energy = float.fromhex(contract["energy_rydberg_hex"])
    old = json.loads(prior.read_text())
    archived = [x for x in old["rows"] if x["model"] == "compressed"]
    checks = {}
    errors = {}
    for category in ["baseline", "variants"]:
        for i, row in enumerate(report[category]):
            q = np.asarray(row["q"])
            thresholds = np.asarray(row["thresholds"])
            endpoint, step = row["endpoint_bohr"], row["step_bohr"]
            low = discrete_free_k(q, thresholds, energy, 18.0, endpoint, step, 80)
            high = discrete_free_k(q, thresholds, energy, 18.0, endpoint, step, 110)
            cont = exact_free_k(q, thresholds, energy, 18.0, 80)
            cont_high = exact_free_k(q, thresholds, energy, 18.0, 110)
            with mp.workdps(80):
                precision = max(float(abs(a - b)) for a, b in zip(low, high, strict=True))
                continuous_precision = max(
                    float(abs(a - b)) for a, b in zip(cont, cont_high, strict=True)
                )
                stored = mp.matrix(row["discrete_k_decimal"])
                stored_cont = mp.matrix(row["continuous_k_decimal"])
                if max(abs(a - b) for a, b in zip(stored, high, strict=True)) > mp.mpf("1e-60"):
                    raise ValueError("Stored discrete oracle differs")
                if max(abs(a - b) for a, b in zip(stored_cont, cont_high, strict=True)) > mp.mpf(
                    "1e-60"
                ):
                    raise ValueError("Stored continuous oracle differs")
                difference = mp.matrix(row["k"]) - high
                truncation = high - cont_high
                roundoff = float(max(abs(x) for x in difference))
                discretisation = float(max(abs(x) for x in truncation))
                if (
                    abs(roundoff - row["roundoff_k"]) > 1e-15
                    or abs(discretisation - row["truncation_k"]) > 1e-15
                ):
                    raise ValueError("Error decomposition differs")
            if category == "baseline":
                checks[f"baseline_{i}_precision"] = max(precision, continuous_precision) < 1e-20
                oldrow = next(
                    x
                    for x in archived
                    if (x["coupled_closed"], x["step_bohr"], x["endpoint_bohr"], x["tail_order"])
                    == (row["coupled_closed"], step, endpoint, row["tail_order"])
                )
                if row["k"] != oldrow["k"]:
                    raise ValueError("Archived binary64 K differs")
                if not row["coupled_closed"]:
                    scalar = exact_discrete_scalar(
                        float(q[0, 0]), energy, 18.0, endpoint, step, 110
                    )
                    with mp.workdps(80):
                        checks[f"baseline_{i}_scalar_consistency"] = abs(
                            scalar - high[0, 0]
                        ) < mp.mpf("1e-20")
            else:
                name = row["representation"]
                index = contract["variant_coordinates"].index([step, endpoint])
                checks[f"{name}_{index}_precision"] = max(precision, continuous_precision) < 1e-20
                # The full boundary is saved in each coupled baseline row.
                full = next(x for x in report["baseline"] if x["coupled_closed"])
                original = exact_free_k(
                    np.asarray(full["q"]), np.asarray(full["thresholds"]), energy, 18.0, 110
                )
                transport = (
                    np.asarray(report["eigenvectors"]) if name == "reduced_diagonal" else np.eye(2)
                )
                equivalent = transport @ np.array(cont_high.tolist(), float) @ transport.T
                checks[f"{name}_{index}_equivalence"] = all(
                    v <= 1e-8
                    for v in compare(equivalent, np.array(original.tolist(), float)).values()
                )
                k = np.asarray(row["k"])
                s = scattering(k)
                checks[f"{name}_{index}_invariants"] = (
                    max(
                        float(np.max(abs(s.conj().T @ s - np.eye(len(s))))),
                        float(np.max(abs(s - s.T))),
                    )
                    <= 1e-6
                )
                transported = transport @ np.array(difference.tolist(), float) @ transport.T
                errors[name, step, endpoint] = float(np.max(abs(transported)))
    for row in report["reexecutions"]:
        coupled = row["coupled_closed"]
        endpoint = row["endpoint_bohr"]
        original = next(
            x
            for x in archived
            if (x["coupled_closed"], x["step_bohr"], x["endpoint_bohr"], x["tail_order"])
            == (coupled, 0.0015625, endpoint, 8)
        )
        checks[f"reexecution_{coupled}_exact"] = row["k"] == original["k"]
        s = scattering(np.asarray(row["k"]))
        checks[f"reexecution_{coupled}_invariants"] = (
            max(float(np.max(abs(s.conj().T @ s - np.eye(len(s))))), float(np.max(abs(s - s.T))))
            <= 1e-6
        )
    if checks != report["checks"]:
        raise ValueError("Independent diagnostic checks differ")
    points = [
        x
        for x in report["baseline"]
        if x["coupled_closed"]
        and x["step_bohr"] == 0.0015625
        and x["endpoint_bohr"] in [200.0, 400.0]
        and x["tail_order"] == 8
    ]
    dominance = all(x["roundoff_k"] >= 5 * x["truncation_k"] for x in points)
    representation = all(
        errors["reduced_diagonal", 0.0015625, r]
        <= min(0.1 * errors["reduced_rotated", 0.0015625, r], 1e-7)
        for r in [200.0, 400.0, 800.0]
    )
    hypotheses = report["hypotheses"]
    if (
        dominance is not hypotheses["coupled_rejected_points_roundoff_dominates_5fold"]
        or representation
        is not hypotheses["all_finest_diagonal_roundoff_below_tenth_rotated_and_1em7"]
    ):
        raise ValueError("Independent hypothesis decisions differ")
    if report["diagnostic_verdict"] is not all(checks.values()):
        raise ValueError("Diagnostic aggregate differs")
    if execution["original_command_exit"] != execution["original_container_exit"] or execution[
        "original_command_exit"
    ] != (0 if report["diagnostic_verdict"] else 1):
        raise ValueError("Original diagnostic exits differ")
    return dict(
        saved_decision_replay_passed=True,
        diagnostic_verdict=report["diagnostic_verdict"],
        checks_passed=sum(checks.values()),
        checks_total=len(checks),
        hypotheses=hypotheses,
        baseline_rows_verified=12,
        new_solves_verified=10,
        original_analytic_verdict=False,
        original_analytic_exits=[1, 1],
        diagnostic_exits=[execution["original_command_exit"], execution["original_container_exit"]],
        co_successor_submitted=False,
        new_co_independent_solves=0,
        new_unique_energy_identities=0,
        candidate_qualified_identities=14,
        unique_energy_coverage=132,
        uncomputed_identities=462,
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
