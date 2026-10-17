"""Run the frozen physical shell real-axis/pole matrix, preserving every outcome."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from qscat.units import HARTREE_TO_EV

from projects.ukrmol_co.physical_pole import (
    boundary_s,
    exact_pole,
    exact_s,
    shell_boundary,
    shoot_pole,
)


def oracle(q_grid: np.ndarray, contract: dict) -> dict:
    """Evaluate closed forms at two precisions, retaining a precision/error audit."""
    previous = None
    precision_error = 0.0
    residual_error = 0.0
    for digits in contract["oracle_decimal_precisions"]:
        with mp.workdps(digits):
            values = {s: [exact_s(q, s) for q in q_grid] for s in contract["shell_lambdas"]}
            poles = {
                (s, n): exact_pole(s, n)
                for s in contract["shell_lambdas"]
                if s > 0
                for n in contract["pole_indices"]
            }
            for (s, _), q in poles.items():
                residual_error = max(
                    residual_error, float(abs(1 + s * mp.exp(1j * q) * mp.sin(q) / q))
                )
            if previous is not None:
                differences = [
                    abs(v - w)
                    for s in values
                    for v, w in zip(values[s], previous[0][s], strict=True)
                ]
                differences += [abs(v - previous[1][key]) for key, v in poles.items()]
                precision_error = float(max(differences))
            previous = values, poles
    return {
        "values": {s: np.array([complex(v) for v in a]) for s, a in values.items()},
        "poles": {key: complex(v) for key, v in poles.items()},
        "precision_error": precision_error,
        "jost_residual_error": residual_error,
    }


def real_profile(strength, radius, tolerance, grid, reference, contract):
    """Full real-axis grid, including boundary representations and threshold controls."""
    samples = []
    representation = 0.0
    for q in grid:
        state = shell_boundary(q, strength, radius, tolerance)
        value = boundary_s(q, state, radius)
        for b in contract["bloch_coefficients"]:
            alternative = boundary_s(q, state * contract["boundary_pair_scale"], radius, b)
            representation = max(representation, abs(alternative - value))
        samples.append(value)
    samples = np.asarray(samples)
    threshold = np.isin(grid, contract["threshold_q"])
    return {
        "s_values": [[v.real, v.imag] for v in samples],
        "maximum_s_error": float(max(abs(samples - reference))),
        "maximum_phase_error_rad_mod_pi": float(max(abs(np.angle(samples * reference.conj()))) / 2),
        "maximum_unitarity_error": float(max(abs(abs(samples) ** 2 - 1))),
        "maximum_bloch_rescaling_error": float(representation),
        "maximum_scattering_length_ratio_error": float(
            max(
                abs(
                    -np.angle(samples[threshold]) / (2 * grid[threshold])
                    - strength / (1 + strength)
                )
            )
        ),
    }


def pole_profile(strength, index, radius, tolerance, reference, contract):
    result = shoot_pole(strength, index, radius, tolerance, contract)
    q = complex(*result["q"])
    energy = q * q / (2 * contract["radius_bohr"] ** 2)
    truth = reference * reference / (2 * contract["radius_bohr"] ** 2)
    result.update(
        {
            "q_error": abs(q - reference),
            "oracle_q": [reference.real, reference.imag],
            "energy_hartree": [energy.real, energy.imag],
            "oracle_energy_hartree": [truth.real, truth.imag],
            "position_error_ev": abs(energy.real - truth.real) * HARTREE_TO_EV,
            "full_width_error_ev": 2 * abs(energy.imag - truth.imag) * HARTREE_TO_EV,
            "position_ev": energy.real * HARTREE_TO_EV,
            "full_width_ev": -2 * energy.imag * HARTREE_TO_EV,
            "outgoing_quadrant_positive_position_width": bool(
                q.real > 0 and q.imag < 0 and energy.real > 0 and energy.imag < 0
            ),
        }
    )
    return result


def study(contract: dict) -> dict:
    start = time.monotonic()
    low, high, points = contract["real_q_linear_grid"]
    grid = np.unique(
        np.r_[
            np.linspace(low, high, int(points)),
            np.pi * np.array(contract["real_q_node_multiples_pi"]),
            contract["threshold_q"],
        ]
    )
    reference = oracle(grid, contract)
    real, poles = [], []
    for tolerance in contract["ode_relative_tolerances"]:
        for strength in contract["shell_lambdas"]:
            for radius in contract["matching_radius_ratios"]:
                row = {
                    "lambda": strength,
                    "matching_radius_ratio": radius,
                    "ode_rtol": tolerance,
                    "reading_completed": False,
                }
                try:
                    row.update(
                        real_profile(
                            strength,
                            radius,
                            tolerance,
                            grid,
                            reference["values"][strength],
                            contract,
                        )
                    )
                    row["reading_completed"] = True
                except (ValueError, RuntimeError, FloatingPointError) as error:
                    row["failure"] = f"{type(error).__name__}: {error}"
                real.append(row)
                if strength <= 0:
                    continue
                for index in contract["pole_indices"]:
                    row = {
                        "lambda": strength,
                        "pole_index": index,
                        "matching_radius_ratio": radius,
                        "ode_rtol": tolerance,
                        "reading_completed": False,
                    }
                    try:
                        row.update(
                            pole_profile(
                                strength,
                                index,
                                radius,
                                tolerance,
                                reference["poles"][strength, index],
                                contract,
                            )
                        )
                        row["reading_completed"] = True
                    except (ValueError, RuntimeError, FloatingPointError) as error:
                        row["failure"] = f"{type(error).__name__}: {error}"
                    poles.append(row)
    complete = all(r["reading_completed"] for r in real + poles)
    checks = {
        "complete_declared_matrix": complete
        and len(grid) == contract["expected_real_q_points"]
        and len(real) == contract["expected_real_profiles"]
        and len(poles) == contract["expected_pole_solves"],
        "oracle_precision": reference["precision_error"]
        <= contract["oracle_precision_absolute_budget"],
        "oracle_unmultiplied_jost_residual": reference["jost_residual_error"]
        <= contract["oracle_jost_absolute_budget"],
    }
    convergence, radius_spread = [], None
    if complete:
        finest = contract["ode_relative_tolerances"][-1]
        rf = [r for r in real if r["ode_rtol"] == finest]
        pf = [r for r in poles if r["ode_rtol"] == finest]
        for key, budget in [
            ("maximum_s_error", "finest_s_absolute_budget"),
            ("maximum_phase_error_rad_mod_pi", "finest_phase_budget_rad_mod_pi"),
            (
                "maximum_scattering_length_ratio_error",
                "threshold_scattering_length_ratio_absolute_budget",
            ),
        ]:
            checks[key] = max(r[key] for r in rf) <= contract[budget]
        for key, budget in [
            ("maximum_unitarity_error", "unitarity_absolute_budget"),
            ("maximum_bloch_rescaling_error", "representation_absolute_budget"),
        ]:
            checks[key] = max(r[key] for r in real) <= contract[budget]
        checks["all_root_solver_statuses"] = all(p["solver_success"] for p in poles)
        checks["all_outgoing_quadrants_positive_widths"] = all(
            p["outgoing_quadrant_positive_position_width"] for p in poles
        )
        for key, budget in [
            ("q_error", "finest_q_absolute_budget"),
            ("position_error_ev", "finest_position_budget_ev"),
            ("full_width_error_ev", "finest_full_width_budget_ev"),
            ("outgoing_residual_absolute", "finest_root_residual_absolute_budget"),
        ]:
            checks[key] = max(p[key] for p in pf) <= contract[budget]
        checks["no_numerator_cancellation"] = (
            min(p["incoming_residual_absolute"] for p in pf)
            >= contract["noncancelling_incoming_residual_minimum"]
        )
        for tolerance in contract["ode_relative_tolerances"]:
            convergence.append(
                {
                    "ode_rtol": tolerance,
                    "maximum_s_error": max(
                        r["maximum_s_error"] for r in real if r["ode_rtol"] == tolerance
                    ),
                    "maximum_q_error": max(
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
            for key in ["maximum_s_error", "maximum_q_error"]
        )
        radius_spread = 0.0
        for s, n in reference["poles"]:
            values = np.array(
                [complex(*p["q"]) for p in pf if p["lambda"] == s and p["pole_index"] == n]
            )
            radius_spread = max(radius_spread, float(np.max(abs(values[:, None] - values))))
        checks["finest_matching_radius_invariance"] = (
            radius_spread <= contract["finest_radius_q_spread_budget"]
        )
    q, strength = contract["negative_control_q"], contract["negative_control_lambda"]
    with mp.workdps(contract["oracle_decimal_precisions"][-1]):
        exact = complex(exact_s(q, strength))
        negative = {
            "reversed_shell_jump": abs(complex(exact_s(q, -strength)) - exact),
            "missing_factor_two_shell_jump": abs(complex(exact_s(q, strength / 2)) - exact),
        }
    checks["mandatory_negative_controls"] = (
        min(negative.values()) > contract["negative_control_minimum_s_difference"]
    )
    # A finite value at q=0 rules out the removable Lambert-W principal solution.
    state = shell_boundary(0, strength, 1, contract["ode_relative_tolerances"][-1])
    zero_error = float(abs(state[1] - (1 + strength)))
    checks["threshold_not_a_pole"] = zero_error <= contract["representation_absolute_budget"]
    return {
        "contract": contract,
        "q_grid": grid.tolist(),
        "real_profiles": real,
        "pole_solves": poles,
        "oracle_precision_difference": reference["precision_error"],
        "oracle_maximum_jost_residual": reference["jost_residual_error"],
        "convergence": convergence,
        "finest_radius_q_spread": radius_spread,
        "negative_controls_s_differences": negative,
        "threshold_jost_error": zero_error,
        "checks": checks,
        "physical_shell_benchmark_verdict": bool(all(checks.values())),
        "wall_seconds": time.monotonic() - start,
        "native_data_read": False,
        "native_solves": 0,
        "coupled_polar_tail_verdict": None,
        "native_boundary_normalization_verdict": None,
        "co_complex_pole_verdict": None,
        "full_sweep_released": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = args.contract.read_bytes()
    result = study(json.loads(raw))
    result["contract_sha256"] = hashlib.sha256(raw).hexdigest()
    with args.output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if result["physical_shell_benchmark_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
