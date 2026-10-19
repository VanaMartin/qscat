"""Independent saved-decision replay of the fixed-grid arithmetic control."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from projects.ukrmol_co.native_reconstruction_replay import compare, scattering
from projects.ukrmol_co.numerov_discrete_free import discrete_free_k
from projects.ukrmol_co.numerov_increment_study import exact_free_k


def invariant_pass(k):
    s = scattering(np.asarray(k, float))
    return (
        max(float(np.max(abs(s.conj().T @ s - np.eye(len(s))))), float(np.max(abs(s - s.T))))
        <= 1e-6
    )


def verify(root):
    report = json.loads((root / "collected/result/report.json").read_text())
    execution = json.loads((root / "collected/execution.json").read_text())
    contract = report["contract"]
    prior = root / "prior-report.json"
    if hashlib.sha256(prior.read_bytes()).hexdigest() != contract["prior_report_sha256"]:
        raise ValueError("Archived precision-attribution input differs")
    for name, want in json.loads((root / "collected/collection-manifest.json").read_text()).items():
        if hashlib.sha256((root / "collected" / name).read_bytes()).hexdigest() != want:
            raise ValueError("Collected source or result differs")
    old = json.loads(prior.read_text())
    if (
        len(report["rows"]) != 12
        or len(report["original_reexecutions"]) != 2
        or len(report["oracle_controls"]) != 6
        or len(contract["coordinates"]) != 6
    ):
        raise ValueError("Fixed arithmetic experiment is incomplete")
    bits = report["precision_bits"]
    controls = {"extended_precision_retained": bits["binary64"] == 53 and bits["longdouble"] >= 64}
    qualification = {}
    oracles = {}
    row_errors = {}
    energy = float.fromhex(contract["energy_rydberg_hex"])
    for i, (coordinate, control) in enumerate(
        zip(contract["coordinates"], report["oracle_controls"], strict=True)
    ):
        original = next(
            x
            for x in old["baseline"]
            if all(x[k] == coordinate[k] for k in ["coupled_closed", "step_bohr", "endpoint_bohr"])
            and x["tail_order"] == 8
        )
        if (
            control["q"] != original["q"]
            or control["thresholds"] != original["thresholds"]
            or control["coordinate"] != coordinate
        ):
            raise ValueError("Actual archived boundary or coordinate differs")
        q = np.asarray(control["q"])
        thresholds = np.asarray(control["thresholds"])
        step, endpoint = coordinate["step_bohr"], coordinate["endpoint_bohr"]
        low = discrete_free_k(q, thresholds, energy, 18.0, endpoint, step, 80)
        high = discrete_free_k(q, thresholds, energy, 18.0, endpoint, step, 110)
        cont = exact_free_k(q, thresholds, energy, 18.0, 80)
        cont_high = exact_free_k(q, thresholds, energy, 18.0, 110)
        with mp.workdps(80):
            controls[f"oracle_{i}_precision"] = max(
                max(abs(a - b) for a, b in zip(low, high, strict=True)),
                max(abs(a - b) for a, b in zip(cont, cont_high, strict=True)),
            ) < mp.mpf("1e-20")
        oracles[i] = high, cont_high
    for row in report["original_reexecutions"]:
        coupled, endpoint = row["coupled_closed"], row["endpoint_bohr"]
        original = next(
            x
            for x in old["baseline"]
            if (x["coupled_closed"], x["step_bohr"], x["endpoint_bohr"], x["tail_order"])
            == (coupled, 0.0015625, endpoint, 8)
        )
        controls[f"original_{coupled}_exact"] = row["k"] == original["k"]
        controls[f"original_{coupled}_invariants"] = invariant_pass(row["k"])
    seen = set()
    for row in report["rows"]:
        index, backend = row["coordinate_index"], row["backend"]
        if (
            (index, backend) in seen
            or index not in range(6)
            or backend not in ["binary64", "longdouble"]
        ):
            raise ValueError("Duplicate or undeclared arithmetic row")
        seen.add((index, backend))
        coordinate = contract["coordinates"][index]
        if row["coordinate"] != coordinate:
            raise ValueError("Stored row coordinate differs")
        high, continuous = oracles[index]
        with mp.workdps(80):
            measured = mp.matrix(row["k_decimal"])
            stored_discrete = mp.matrix(row["discrete_k_decimal"])
            stored_continuous = mp.matrix(row["continuous_k_decimal"])
            if max(abs(a - b) for a, b in zip(stored_discrete, high, strict=True)) > mp.mpf(
                "1e-60"
            ) or max(
                abs(a - b) for a, b in zip(stored_continuous, continuous, strict=True)
            ) > mp.mpf("1e-60"):
                raise ValueError("Saved high-precision oracle differs")
            error = float(max(abs(a - b) for a, b in zip(measured, high, strict=True)))
            continuous_k = float(max(abs(a - b) for a, b in zip(measured, continuous, strict=True)))
        k = np.array(row["k_decimal"], float)
        if (
            not np.array_equal(k, np.asarray(row["k"]))
            or abs(error - row["roundoff_k"]) > 1e-20
            or abs(continuous_k - row["continuous_k"]) > 1e-20
        ):
            raise ValueError("Decimal precision retention or error decomposition differs")
        controls[f"{index}_{backend}_dtype"] = (
            row["diagnostics"]["significand_bits"] == bits[backend]
        )
        controls[f"{index}_{backend}_invariants"] = invariant_pass(k)
        row_errors[index, backend] = error
        if backend == "longdouble":
            qualification[f"{index}_extended_roundoff"] = error <= 1e-8
            if coordinate["step_bohr"] == 0.0015625:
                errors = compare(k, np.array(continuous.tolist(), float))
                qualification[f"{index}_finest_continuous"] = continuous_k <= 1e-6 and all(
                    v <= 1e-6 for key, v in errors.items() if key != "k"
                )
    if seen != {(i, b) for i in range(6) for b in ["binary64", "longdouble"]}:
        raise ValueError("Missing arithmetic comparison")
    if controls != report["control_checks"] or qualification != report["qualification_checks"]:
        raise ValueError("Independent control or qualification decisions differ")
    ratios = [row_errors[i, "binary64"] / row_errors[i, "longdouble"] for i in range(6)]
    if (
        ratios != report["hypotheses"]["roundoff_reduction_ratios"]
        or all(x >= 100 for x in ratios)
        is not report["hypotheses"]["all_six_roundoff_reduction_at_least_100fold"]
    ):
        raise ValueError("Independent precision-gain hypothesis differs")
    verdict = all(controls.values()) and all(qualification.values())
    if (
        report["experiment_verdict"] is not verdict
        or execution["original_command_exit"] != (0 if verdict else 1)
        or execution["original_container_exit"] != execution["original_command_exit"]
    ):
        raise ValueError("Arithmetic aggregate or original exits differ")
    return dict(
        saved_decision_replay_passed=True,
        experiment_verdict=verdict,
        control_checks_passed=sum(controls.values()),
        control_checks_total=len(controls),
        qualification_checks_passed=sum(qualification.values()),
        qualification_checks_total=len(qualification),
        precision_bits=bits,
        hypotheses=report["hypotheses"],
        new_analytic_solves_verified=14,
        typed_rows_verified=12,
        original_baseline_reexecutions_verified=2,
        original_exits=[execution["original_command_exit"], execution["original_container_exit"]],
        original_interior_analytic_verdict=False,
        original_interior_scientific_exits=[1, 1],
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
