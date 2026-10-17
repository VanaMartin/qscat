"""Fixed noncommuting threshold/polar-tail boundary-operator experiment."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
import warnings
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.linalg import LinAlgWarning

from projects.ukrmol_co.combined_exterior import (
    exact_map,
    outgoing_momenta,
    riccati_map,
    spectral_map,
)


def unpack(row, method):
    return np.array(row[method + "_real"]) + 1j * np.array(row[method + "_imag"])


def entry_error(left, right):
    return float(np.max(abs(left - right)))


def evaluate(energy, dipole, length, angle, rtol, degree, separation):
    row = {
        "epsilon": [energy.real, energy.imag],
        "dipole": dipole,
        "endpoint": length,
        "angle_degrees": angle,
        "rtol": rtol,
        "degree": degree,
        "separation": separation,
        "completed": False,
    }
    try:
        y, evaluations = riccati_map(energy, dipole, length, angle, rtol, separation)
        with warnings.catch_warnings():
            warnings.simplefilter("error", LinAlgWarning)
            spectral, condition, residual = spectral_map(
                energy, dipole, length, angle, degree, separation
            )
        if not np.isfinite(y).all() or not np.isfinite(spectral).all():
            raise ValueError("Nonfinite boundary operator")
        for name, matrix in [("ode", y), ("spectral", spectral)]:
            row[name + "_real"] = matrix.real.tolist()
            row[name + "_imag"] = matrix.imag.tolist()
            row[name + "_symmetry_error"] = entry_error(matrix, matrix.T)
            if energy.imag == 0:
                flux = (matrix - matrix.conj().T) / (2j)
                row[name + "_minimum_flux_eigenvalue"] = float(np.linalg.eigvalsh(flux).min())
        row.update(
            {
                "completed": True,
                "ode_evaluations": evaluations,
                "spectral_condition_inf": condition,
                "spectral_residual": residual,
                "method_difference": entry_error(y, spectral),
            }
        )
    except (ValueError, RuntimeError, np.linalg.LinAlgError, LinAlgWarning) as error:
        row["failure"] = f"{type(error).__name__}: {error}"
    return row


def reference_maps(contract):
    previous = None
    precision_error = None
    cases = [(complex(*e), 0.0, 16.0) for e in contract["energies_epsilon"]]
    cases += [
        (complex(*e), d, 0.0)
        for e in contract["degenerate_control_energies"]
        for d in contract["degenerate_control_dipoles"]
    ]
    for digits in contract["oracle_decimal_precisions"]:
        with mp.workdps(digits):
            current = {key: exact_map(*key) for key in cases}
            if previous is not None:
                precision_error = float(
                    max(
                        abs(v - w)
                        for key in current
                        for v, w in zip(current[key], previous[key], strict=True)
                    )
                )
            previous = current
    return {
        key: np.array(value.tolist(), complex) for key, value in current.items()
    }, precision_error


def summarize(rows, controls, refs, precision, contract):
    checks = {
        "complete_matrix": len(rows) == contract["expected_cases"]
        and len(controls) == contract["expected_degenerate_controls"]
        and all(r["completed"] for r in rows + controls),
        "analytic_oracle_precision": precision <= contract["oracle_precision_absolute_budget"],
    }
    summary = {"checks": checks, "reference_precision_error": precision}
    if not checks["complete_matrix"]:
        return summary
    finest = contract["ode_relative_tolerances"][-1]
    fine = [r for r in rows if r["rtol"] == finest]
    checks["finest_method_agreement"] = (
        max(r["method_difference"] for r in fine + controls)
        <= contract["matrix_entry_absolute_budget"]
    )
    analytic_errors = []
    for row in fine + controls:
        key = (complex(*row["epsilon"]), row["dipole"], row["separation"])
        if key in refs:
            analytic_errors.extend(
                entry_error(unpack(row, method), refs[key]) for method in ["ode", "spectral"]
            )
    checks["analytic_limits"] = max(analytic_errors) <= contract["matrix_entry_absolute_budget"]
    refinement, domain = {}, {}
    for method in ["ode", "spectral"]:
        groups = {}
        for row in rows:
            key = (tuple(row["epsilon"]), row["dipole"], row["endpoint"], row["angle_degrees"])
            groups.setdefault(key, []).append(row)
        differences = [[], []]
        for group in groups.values():
            for i in range(2):
                differences[i].append(
                    entry_error(unpack(group[i], method), unpack(group[i + 1], method))
                )
        coarse, final = map(max, differences)
        refinement[method] = {"coarse_to_middle": coarse, "middle_to_finest": final}
        checks[method + "_refinement"] = final <= contract[
            "refinement_absolute_budget"
        ] and final <= max(
            coarse / contract["convergence_improvement_factor"],
            contract["convergence_roundoff_floor"],
        )
        groups = {}
        for row in fine:
            groups.setdefault((tuple(row["epsilon"]), row["dipole"]), []).append(row)
        spread, endpoint_differences = 0.0, [0.0, 0.0]
        for group in groups.values():
            for left in group:
                for right in group:
                    spread = max(spread, entry_error(unpack(left, method), unpack(right, method)))
            for angle in contract["angles_degrees"]:
                ordered = [r for r in group if r["angle_degrees"] == angle]
                for i in range(2):
                    endpoint_differences[i] = max(
                        endpoint_differences[i],
                        entry_error(unpack(ordered[i], method), unpack(ordered[i + 1], method)),
                    )
        domain[method] = {
            "all_endpoint_angle_spread": spread,
            "successive_endpoint_differences": endpoint_differences,
        }
        checks[method + "_domain_contour"] = spread <= contract["domain_contour_absolute_budget"]
        real = [r for r in fine if r["epsilon"][1] == 0]
        checks[method + "_real_transpose_symmetry"] = (
            max(r[method + "_symmetry_error"] for r in real)
            <= contract["real_symmetry_flux_absolute_budget"]
        )
        checks[method + "_real_flux_positivity"] = (
            min(r[method + "_minimum_flux_eigenvalue"] for r in real)
            >= -contract["real_symmetry_flux_absolute_budget"]
        )
    target = next(
        r
        for r in fine
        if r["epsilon"] == [8.0, 0.0]
        and r["dipole"] == 1.2
        and r["endpoint"] == 80
        and r["angle_degrees"] == 35
    )
    correct = unpack(target, "ode")
    wrong_incoming = refs[8 + 0j, 0.0, 16.0].copy()
    wrong_incoming[0, 0] = wrong_incoming[0, 0].conjugate()
    negatives = {
        "omitted_dipole": entry_error(correct, refs[8 + 0j, 0.0, 16.0]),
        "omitted_derivative_rotation": entry_error(correct, np.exp(1j * np.deg2rad(35)) * correct),
        "wrong_incoming_open_sheet": entry_error(wrong_incoming, refs[8 + 0j, 0.0, 16.0]),
    }
    for name, error in negatives.items():
        checks[name + "_detected"] = error > contract["negative_difference_minimum"]
    checks["all_rotated_momenta_decay"] = all(
        np.min(
            (
                outgoing_momenta(complex(*r["epsilon"]), r["separation"])
                * np.exp(1j * np.deg2rad(r["angle_degrees"]))
            ).imag
        )
        > 0
        for r in rows + controls
    )
    summary.update(
        {
            "maximum_finest_method_difference": max(
                r["method_difference"] for r in fine + controls
            ),
            "maximum_finest_analytic_error": max(analytic_errors),
            "maximum_spectral_condition_inf": max(
                r["spectral_condition_inf"] for r in rows + controls
            ),
            "maximum_spectral_residual": max(r["spectral_residual"] for r in rows + controls),
            "refinement": refinement,
            "domain_contour": domain,
            "negative_controls": negatives,
        }
    )
    return summary


def study(contract):
    start = time.monotonic()
    refs, precision = reference_maps(contract)
    rows = []
    for rtol, degree in zip(
        contract["ode_relative_tolerances"], contract["chebyshev_degrees"], strict=True
    ):
        for dipole in contract["dipoles"]:
            for energy in contract["energies_epsilon"]:
                for length in contract["endpoint_lengths"]:
                    for angle in contract["angles_degrees"]:
                        rows.append(
                            evaluate(complex(*energy), dipole, length, angle, rtol, degree, 16.0)
                        )
    controls = [
        evaluate(complex(*energy), dipole, 80, 35, 1e-11, 256, 0.0)
        for energy in contract["degenerate_control_energies"]
        for dipole in contract["degenerate_control_dipoles"]
    ]
    summary = summarize(rows, controls, refs, precision, contract)
    return {
        "contract": contract,
        "cases": rows,
        "degenerate_controls": controls,
        **summary,
        "combined_exterior_dtn_verdict": bool(all(summary["checks"].values())),
        "wall_seconds": time.monotonic() - start,
        "native_data_read": False,
        "native_solves": 0,
        "native_full_k_reconstruction_verdict": None,
        "co_complex_pole_verdict": None,
        "full_sweep_released": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = args.contract.read_bytes()
    report = study(json.loads(raw))
    report["contract_sha256"] = hashlib.sha256(raw).hexdigest()
    with args.output.open("x") as stream:
        stream.write(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if report["combined_exterior_dtn_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
