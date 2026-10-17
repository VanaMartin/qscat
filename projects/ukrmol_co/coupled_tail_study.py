"""Frozen matrix tail, threshold elimination and outgoing-sheet validation."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from qscat.units import HARTREE_TO_EV
from scipy.optimize import root

from projects.ukrmol_co.coupled_tail import (
    exact_jost,
    exact_tail,
    exact_threshold,
    exact_threshold_pole,
    tail_scattering,
    threshold_jost,
    threshold_scattering,
)
from projects.ukrmol_co.physical_pole import exact_pole


def as_array(value):
    return np.array(value.tolist(), complex)


def reference_data(contract):
    previous = None
    precision = 0.0
    for digits in contract["oracle_decimal_precisions"]:
        with mp.workdps(digits):
            qs = [complex(q) for q in contract["tail_real_q"]]
            qs += [complex(*q) for q in contract["tail_complex_q"]]
            tails = {(d, q): exact_tail(q, d) for d in contract["tail_dipoles"] for q in qs}
            shells = {
                (c, e): exact_threshold(e, c)
                for c in contract["threshold_couplings"]
                for e in contract["threshold_energies_epsilon"]
            }
            poles = {
                c: exact_threshold_pole(c, contract["oracle_pole_secant_seeds_q"])
                for c in contract["threshold_couplings"]
            }
            current = [
                v for values in [tails, shells] for matrix in values.values() for v in matrix
            ]
            current += list(poles.values())
            if previous is not None:
                precision = float(max(abs(v - w) for v, w in zip(current, previous, strict=True)))
            previous = current
            jost = max(float(abs(exact_jost(q, c))) for c, q in poles.items())
            incoming = min(float(abs(exact_jost(q, c, incoming=True))) for c, q in poles.items())
            lambert = float(abs(poles[0] - exact_pole(10, 1)))
    return {
        "tails": {key: as_array(v) for key, v in tails.items()},
        "shells": {key: as_array(v) for key, v in shells.items()},
        "poles": {c: complex(q) for c, q in poles.items()},
        "precision_error": precision,
        "jost_residual": jost,
        "minimum_incoming_jost": incoming,
        "decoupled_lambert_error": lambert,
    }


def matrix_metrics(result, reference, condition, real_axis):
    if not np.isfinite(result).all():
        raise ValueError("Nonfinite S matrix")
    metrics = {
        "s_real": result.real.tolist(),
        "s_imag": result.imag.tolist(),
        "s_entry_error": float(np.max(abs(result - reference))),
        "matching_condition_number": condition,
        "open_channels": len(result),
        "real_axis": bool(real_axis),
    }
    if real_axis:
        metrics.update(
            {
                "unitarity_error": float(
                    np.max(abs(result.conj().T @ result - np.eye(len(result))))
                ),
                "reciprocity_error": float(np.max(abs(result - result.T))),
            }
        )
    return metrics


def pole_solve(coupling, radius, rtol, truth, contract):
    def residual(pair):
        value = threshold_jost(complex(*pair), coupling, radius, rtol)
        return [value.real, value.imag]

    solution = root(
        residual,
        contract["pole_seed_q"],
        tol=contract["root_tolerance"],
        options={"maxfev": contract["root_maxfev"]},
    )
    q = complex(*solution.x)
    energy, exact = (
        q * q / (2 * contract["radius_bohr"] ** 2),
        truth * truth / (2 * contract["radius_bohr"] ** 2),
    )
    return {
        "solver_success": bool(solution.success),
        "solver_status": int(solution.status),
        "solver_message": str(solution.message),
        "evaluations": int(solution.nfev),
        "seed_q": contract["pole_seed_q"],
        "q": [q.real, q.imag],
        "oracle_q": [truth.real, truth.imag],
        "q_error": abs(q - truth),
        "outgoing_residual": abs(threshold_jost(q, coupling, radius, rtol)),
        "energy_hartree": [energy.real, energy.imag],
        "position_ev": float(energy.real * HARTREE_TO_EV),
        "full_width_ev": float(-2 * energy.imag * HARTREE_TO_EV),
        "position_error_ev": float(abs(energy.real - exact.real) * HARTREE_TO_EV),
        "full_width_error_ev": float(2 * abs(energy.imag - exact.imag) * HARTREE_TO_EV),
        "proper_outgoing_sheet": bool(
            q.real > 0 and q.imag < 0 and 0 < (q * q).real < 16 and energy.imag < 0
        ),
    }


def study(contract):
    started = time.monotonic()
    refs = reference_data(contract)
    tails, shells, poles = [], [], []
    for tolerance in contract["ode_relative_tolerances"]:
        for d, q in refs["tails"]:
            for radius in contract["tail_matching_radius_ratios"]:
                row = {
                    "dipole": d,
                    "q": [q.real, q.imag],
                    "radius_ratio": radius,
                    "ode_rtol": tolerance,
                    "completed": False,
                }
                try:
                    s, condition = tail_scattering(q, d, radius, tolerance)
                    row.update(matrix_metrics(s, refs["tails"][d, q], condition, q.imag == 0))
                    row["completed"] = True
                except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
                    row["failure"] = f"{type(error).__name__}: {error}"
                tails.append(row)
        for c, e in refs["shells"]:
            for radius in contract["threshold_matching_radius_ratios"]:
                row = {
                    "coupling": c,
                    "epsilon": e,
                    "radius_ratio": radius,
                    "ode_rtol": tolerance,
                    "completed": False,
                }
                try:
                    s, condition = threshold_scattering(e, c, radius, tolerance)
                    row.update(matrix_metrics(s, refs["shells"][c, e], condition, True))
                    row["completed"] = True
                except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
                    row["failure"] = f"{type(error).__name__}: {error}"
                shells.append(row)
        for c, truth in refs["poles"].items():
            for radius in contract["threshold_matching_radius_ratios"]:
                row = {
                    "coupling": c,
                    "radius_ratio": radius,
                    "ode_rtol": tolerance,
                    "completed": False,
                }
                try:
                    row.update(pole_solve(c, radius, tolerance, truth, contract))
                    row["completed"] = True
                except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
                    row["failure"] = f"{type(error).__name__}: {error}"
                poles.append(row)
    finest = contract["ode_relative_tolerances"][-1]
    cutoff = []
    for radius in contract["cutoff_radius_ratios"]:
        q, d = contract["cutoff_q"], contract["cutoff_dipole"]
        s, condition = tail_scattering(q, d, radius, finest, matching_dipole=0)
        cutoff.append(
            {"radius_ratio": radius, **matrix_metrics(s, refs["tails"][d, q], condition, True)}
        )
    with mp.workdps(contract["oracle_decimal_precisions"][-1]):
        c, e = contract["negative_coupling"], contract["negative_closed_epsilon"]
        negative_closed = float(
            np.max(abs(as_array(exact_threshold(e, c, -1)) - refs["shells"][c, e]))
        )
        e = contract["negative_flux_epsilon"]
        negative_flux = float(
            np.max(abs(as_array(exact_threshold(e, c, flux=False)) - refs["shells"][c, e]))
        )
    complete = all(r["completed"] for r in tails + shells + poles)
    checks = {
        "complete_fixed_matrix": complete
        and len(tails) == contract["expected_tail_samples"]
        and len(shells) == contract["expected_threshold_samples"]
        and len(poles) == contract["expected_pole_solves"],
        "oracle_precision": refs["precision_error"] <= contract["oracle_precision_absolute_budget"],
        "oracle_jost": refs["jost_residual"] <= contract["oracle_precision_absolute_budget"],
        "decoupled_lambert_pole": refs["decoupled_lambert_error"]
        <= contract["oracle_precision_absolute_budget"],
        "outgoing_poles_not_cancelled": refs["minimum_incoming_jost"]
        > contract["incoming_oracle_jost_minimum"],
        "wrong_closed_sheet_detected": negative_closed > contract["negative_s_difference_minimum"],
        "missing_flux_normalization_detected": negative_flux
        > contract["negative_s_difference_minimum"],
        "tail_omission_is_visible": cutoff[0]["s_entry_error"]
        > contract["cutoff_first_error_minimum"],
        "tail_cutoff_radius_improves": cutoff[-1]["s_entry_error"]
        < cutoff[0]["s_entry_error"] / contract["cutoff_required_improvement_factor"],
    }
    convergence, spread = [], None
    if complete:
        fine = [r for r in tails + shells if r["ode_rtol"] == finest]
        pf = [p for p in poles if p["ode_rtol"] == finest]
        checks["finest_matrix_accuracy"] = (
            max(r["s_entry_error"] for r in fine) <= contract["finest_s_entry_absolute_budget"]
        )
        for key in ["unitarity_error", "reciprocity_error"]:
            checks["finest_" + key] = (
                max(r[key] for r in fine if r["real_axis"])
                <= contract["finest_unitarity_reciprocity_absolute_budget"]
            )
        checks["all_root_statuses_success"] = all(p["solver_success"] for p in poles)
        checks["all_outgoing_sheets"] = all(p["proper_outgoing_sheet"] for p in poles)
        for key, budget in [
            ("q_error", "finest_q_absolute_budget"),
            ("outgoing_residual", "finest_jost_residual_absolute_budget"),
            ("position_error_ev", "finest_position_full_width_budget_ev"),
            ("full_width_error_ev", "finest_position_full_width_budget_ev"),
        ]:
            checks["finest_" + key] = max(p[key] for p in pf) <= contract[budget]
        for tolerance in contract["ode_relative_tolerances"]:
            convergence.append(
                {
                    "ode_rtol": tolerance,
                    "tail_s_entry_error": max(
                        r["s_entry_error"] for r in tails if r["ode_rtol"] == tolerance
                    ),
                    "threshold_s_entry_error": max(
                        r["s_entry_error"] for r in shells if r["ode_rtol"] == tolerance
                    ),
                    "threshold_q_error": max(
                        p["q_error"] for p in poles if p["ode_rtol"] == tolerance
                    ),
                }
            )
        checks["three_rung_convergence"] = all(
            convergence[-1][key]
            <= max(
                convergence[0][key] / contract["convergence_improvement_factor"],
                contract["convergence_roundoff_floor"],
            )
            for key in ["tail_s_entry_error", "threshold_s_entry_error", "threshold_q_error"]
        )
        spread = 0.0
        for c in contract["threshold_couplings"]:
            qs = np.array([complex(*p["q"]) for p in pf if p["coupling"] == c])
            spread = max(spread, float(np.max(abs(qs[:, None] - qs))))
        checks["finest_radius_q_invariance"] = spread <= contract["finest_radius_q_spread_budget"]
    return {
        "contract": contract,
        "tail_samples": tails,
        "threshold_samples": shells,
        "pole_solves": poles,
        "oracle_precision_error": refs["precision_error"],
        "oracle_jost_residual": refs["jost_residual"],
        "oracle_minimum_incoming_jost": refs["minimum_incoming_jost"],
        "oracle_decoupled_lambert_error": refs["decoupled_lambert_error"],
        "convergence": convergence,
        "finest_radius_q_spread": spread,
        "cutoff_controls": cutoff,
        "negative_closed_sheet_s_difference": negative_closed,
        "negative_flux_s_difference": negative_flux,
        "checks": checks,
        "separate_tail_threshold_benchmark_verdict": bool(all(checks.values())),
        "wall_seconds": time.monotonic() - started,
        "native_data_read": False,
        "native_solves": 0,
        "combined_unequal_threshold_polar_tail_verdict": None,
        "native_boundary_normalization_verdict": None,
        "co_complex_pole_verdict": None,
        "full_sweep_released": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = args.contract.read_bytes()
    result = study(json.loads(raw))
    result["contract_sha256"] = hashlib.sha256(raw).hexdigest()
    with args.output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if result["separate_tail_threshold_benchmark_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
