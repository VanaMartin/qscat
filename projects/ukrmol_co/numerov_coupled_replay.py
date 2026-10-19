"""Saved propagation decisions and explicit preflight-diagnostic audit; no propagation."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from projects.ukrmol_co.coupled_tail import exact_tail
from projects.ukrmol_co.native_reconstruction_replay import compare, scattering
from projects.ukrmol_co.numerov_coupled_oracle import continuous_k, direct_discrete_k
from projects.ukrmol_co.numerov_discrete_free import discrete_free_k
from projects.ukrmol_co.numerov_increment_study import exact_free_k


def mp_max_difference(a, b):
    with mp.workdps(100):
        return float(max(abs(x - y) for x, y in zip(a, b, strict=True)))


def invariants(k):
    s = scattering(np.asarray(k, float))
    return dict(
        unitarity=float(np.max(abs(s.conj().T @ s - np.eye(len(s))))),
        reciprocity=float(np.max(abs(s - s.T))),
    )


def retained(row):
    with mp.workdps(100):
        high = mp.matrix(row["k_decimal"])
    array = np.array(row["k_decimal"], float)
    if not np.array_equal(array, np.asarray(row["k"])):
        raise ValueError("Decimal components and binary64 presentation differ")
    return high, array


def verify(root):
    r = json.loads((root / "collected/result/report.json").read_text())
    execution = json.loads((root / "collected/execution.json").read_text())
    c = r["contract"]
    prior_path = root / "prior-report.json"
    if hashlib.sha256(prior_path.read_bytes()).hexdigest() != c["prior_report_sha256"]:
        raise ValueError("Prior arithmetic report differs")
    for n, w in json.loads((root / "collected/collection-manifest.json").read_text()).items():
        if hashlib.sha256((root / "collected" / n).read_bytes()).hexdigest() != w:
            raise ValueError("Collected source/result/resource hash differs")
    prior = json.loads(prior_path.read_text())
    checks = {}
    bits = r.get("precision_bits", {})
    completed = {stage["stage"] for stage in r["stages"]}
    if "preflight" in completed:
        if bits.get("binary64") != 53 or bits.get("longdouble", 0) < 64:
            raise ValueError("Declared precision is unsupported")
        for p in r["preflight"]:
            backend = p["backend"]
            if p["kind"] == "rank_one":
                name = f"preflight_algebra_{backend}_{p['channels']}"
                checks[name] = p["forward_error"] <= (2e-15 if backend == "binary64" else 2e-17)
                checks[name + "_residual"] = p["diagnostics"][
                    "maximum_solve_backward_error"
                ] <= 32 * 2.0 ** (1 - bits[backend])
            else:
                checks[f"preflight_tail_{backend}_mp"] = max(p["mp_errors"]) <= (
                    2e-13 if backend == "binary64" else 2e-15
                )
                checks[f"preflight_tail_{backend}_bits"] = (
                    p["diagnostics"]["complex_significand_bits"] == bits[backend]
                )
                checks[f"preflight_tail_{backend}_hankel"] = p["hankel_error"] <= 1e-11
    if "short" in completed:
        short = c["short"]
        for model in short["models"]:
            name = model["name"]
            stored_ref = next((x for x in r["references"] if x["model"] == name), None)
            if stored_ref is None:
                break
            q, ts, c2, c3 = (np.asarray(model[k]) for k in ["q", "thresholds", "c2", "c3"])
            ang = model["angular"]
            args = (
                q,
                ts,
                ang,
                c2,
                c3,
                short["energy_rydberg"],
                short["radius_bohr"],
                short["endpoint_bohr"],
                short["tail_order"],
            )
            refs = [continuous_k(*args, *t)[0] for t in short["continuous_tolerances"]]
            stability = float(np.max(abs(refs[0] - refs[1])))
            if np.max(abs(refs[1] - np.asarray(stored_ref["k"]))) > 1e-12:
                raise ValueError("Continuous reference reconstruction differs")
            checks[f"short_{name}_continuous_reference"] = (
                stability <= short["continuous_reference_stability_budget"]
            )
            if not checks[f"short_{name}_continuous_reference"]:
                break
            for step in short["steps_bohr"]:
                low, high = [
                    direct_discrete_k(
                        q,
                        ts,
                        ang,
                        c2,
                        c3,
                        short["energy_rydberg"],
                        short["radius_bohr"],
                        short["endpoint_bohr"],
                        step,
                        short["tail_order"],
                        d,
                    )
                    for d in short["oracle_digits"]
                ]
                key = f"short_{name}_{step}"
                checks[key + "_oracle_precision"] = (
                    mp_max_difference(low, high) < short["oracle_precision_budget"]
                )
                if not checks[key + "_oracle_precision"]:
                    break
                for backend in c["backends"]:
                    row = next(
                        x
                        for x in r["rows"]
                        if x["stage"] == "short"
                        and x["model"] == name
                        and x["step_bohr"] == step
                        and x["backend"] == backend
                    )
                    high_k, k = retained(row)
                    error = mp_max_difference(high_k, high)
                    if abs(error - row["discrete_arithmetic_error"]) > 1e-20:
                        raise ValueError("Short discrete error differs")
                    metrics = compare(k, refs[1])
                    checks[key + "_" + backend + "_arithmetic"] = (
                        error <= short["discrete_arithmetic_budgets"][backend]
                    )
                    if step == short["steps_bohr"][-1]:
                        checks[key + "_" + backend + "_continuous"] = all(
                            v <= short["finest_continuous_K_S_phase_budget"]
                            for v in metrics.values()
                        )
                        actual, want = invariants(k), invariants(refs[1])
                        checks[key + "_" + backend + "_invariant_reference"] = all(
                            abs(actual[k] - want[k])
                            <= short["finest_invariant_reference_difference_budget"]
                            for k in want
                        )
            for backend in c["backends"]:
                family = [
                    x
                    for x in r["rows"]
                    if x["stage"] == "short" and x["model"] == name and x["backend"] == backend
                ]
                if len(family) != 3:
                    continue
                a, b = [compare(np.asarray(x["k"]), refs[1])["k"] for x in family[-2:]]
                checks[f"short_{name}_{backend}_trend"] = (
                    b <= a / short["raw_K_error_reduction_minimum"]
                    or max(a, b) <= short["both_errors_trend_floor"]
                )
    if "physical" in completed:
        p = c["physical"]
        with mp.workdps(80):
            expected = np.array(exact_tail(1.0, p["dipole"]).tolist(), complex)
        for backend in c["backends"]:
            rows = [x for x in r["rows"] if x["stage"] == "physical" and x["backend"] == backend]
            errors = [float(np.max(abs(scattering(retained(x)[1]) - expected))) for x in rows]
            checks[backend + "_physical_S"] = errors[-1] <= p["finest_S_and_invariant_budget"]
            checks[backend + "_physical_invariants"] = all(
                v <= p["finest_S_and_invariant_budget"] for v in invariants(rows[-1]["k"]).values()
            )
            checks[backend + "_physical_trend"] = (
                errors[-1] <= errors[-2] / p["raw_S_error_reduction_minimum"]
                or max(errors[-2:]) <= p["both_errors_trend_floor"]
            )
        checks = {
            ("physical_" + k if "_physical_" in k and not k.startswith("physical_") else k): v
            for k, v in checks.items()
        }
    if "stress" in completed:
        s = c["stress"]
        energy = float.fromhex(s["energy_rydberg_hex"])
        ratios = []
        for endpoint in s["endpoints_bohr"]:
            src = next(
                x
                for x in prior["oracle_controls"]
                if x["coordinate"]["coupled_closed"]
                and x["coordinate"]["endpoint_bohr"] == endpoint
                and x["coordinate"]["step_bohr"] == s["step_bohr"]
            )
            q, ts = np.asarray(src["q"]), np.asarray(src["thresholds"])
            discrete = discrete_free_k(
                q, ts, energy, s["radius_bohr"], endpoint, s["step_bohr"], 100
            )
            continuous = exact_free_k(q, ts, energy, s["radius_bohr"], 100)
            old = next(
                x
                for x in prior["rows"]
                if x["backend"] == "longdouble"
                and x["coordinate"]["coupled_closed"]
                and x["coordinate"]["endpoint_bohr"] == endpoint
                and x["coordinate"]["step_bohr"] == s["step_bohr"]
            )
            errors = {}
            for backend in c["backends"]:
                row = next(
                    x
                    for x in r["rows"]
                    if x["stage"] == "stress"
                    and x["endpoint_bohr"] == endpoint
                    and x["backend"] == backend
                )
                high, k = retained(row)
                errors[backend] = mp_max_difference(high, discrete)
                checks[f"stress_{endpoint}_{backend}_invariants"] = all(
                    v <= s["continuous_K_S_phase_and_invariant_budget"]
                    for v in invariants(k).values()
                )
                if backend == "longdouble":
                    continuous_error = mp_max_difference(high, continuous)
                    metrics = compare(k, np.array(continuous.tolist(), float))
                    with mp.workdps(100):
                        previous = mp.matrix(old["k_decimal"])
                    checks[f"stress_{endpoint}_extended_arithmetic"] = (
                        errors[backend] <= s["extended_roundoff_K_budget"]
                    )
                    checks[f"stress_{endpoint}_extended_continuous"] = continuous_error <= s[
                        "continuous_K_S_phase_and_invariant_budget"
                    ] and all(
                        v <= s["continuous_K_S_phase_and_invariant_budget"]
                        for v in metrics.values()
                    )
                    checks[f"stress_{endpoint}_previous_extended"] = (
                        mp_max_difference(high, previous)
                        <= s["previous_extended_K_difference_budget"]
                    )
            ratios.append(errors["binary64"] / errors["longdouble"])
        hypothesis = dict(
            all_three_precision_gain_at_least_100fold=all(
                x >= s["hypothesis_all_three_precision_gain_minimum"] for x in ratios
            ),
            roundoff_reduction_ratios=ratios,
        )
        if hypothesis != r["hypotheses"]:
            raise ValueError("Precision-gain hypothesis differs")
    if "dense" in completed:
        d = c["dense"]
        energy, radius, endpoint, step = (
            d[k] for k in ["energy_rydberg", "radius_bohr", "endpoint_bohr", "step_bohr"]
        )
        with mp.workdps(100):
            expected = mp.zeros(d["open_channels"])
            block = discrete_free_k(
                d["active_q"],
                [0.0, 0.0, d["closed_threshold_rydberg"]],
                energy,
                radius,
                endpoint,
                step,
                100,
            )
            for i in range(2):
                for j in range(2):
                    expected[i, j] = block[i, j]
            for i, q in enumerate(d["other_open_q"], 2):
                expected[i, i] = discrete_free_k([[q]], [0.0], energy, radius, endpoint, step, 100)[
                    0, 0
                ]
        for backend in c["backends"]:
            row = next(x for x in r["rows"] if x["stage"] == "dense" and x["backend"] == backend)
            high, k = retained(row)
            checks[f"dense_{backend}_dense_arithmetic"] = (
                mp_max_difference(high, expected) <= d["discrete_K_budgets"][backend]
            )
            checks[f"dense_{backend}_dense_invariants"] = all(
                v <= d["invariant_budget"] for v in invariants(k).values()
            )
            checks[f"dense_{backend}_dense_backward_residual"] = row["diagnostics"][
                "maximum_solve_backward_error"
            ] <= 32 * 2.0 ** (1 - bits[backend])
    checks = {k: bool(v) for k, v in checks.items()}
    if checks != r["checks"]:
        missing = set(checks) ^ set(r["checks"])
        changed = [k for k in checks if checks[k] != r["checks"].get(k)]
        raise ValueError(f"Saved decisions differ: {missing}; vetoes {changed}")
    verdict = all(checks.values()) and len(completed) == 5
    if (
        r["experiment_verdict"] is not verdict
        or execution["original_command_exit"] != (0 if verdict else 1)
        or execution["original_container_exit"] != execution["original_command_exit"]
    ):
        raise ValueError("Aggregate verdict or original exits differ")
    return dict(
        saved_decision_replay_passed=True,
        experiment_verdict=verdict,
        checks_passed=sum(checks.values()),
        checks_total=len(checks),
        propagation_rows_verified=len(r["rows"]),
        preflight_role=(
            "Audits stored algebra/tail diagnostic values; "
            "complete public numerical replay reexecutes those controls."
        ),
        original_exits=[execution["original_command_exit"], execution["original_container_exit"]],
        hypotheses=r["hypotheses"],
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
