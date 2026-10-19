"""Prospective high-momentum/large-tangent precision control of the fixed grids."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np

from projects.ukrmol_co.native_exterior_attribution import algebra, differences
from projects.ukrmol_co.numerov_corner_control import configurations
from projects.ukrmol_co.numerov_increment import increment_numerov_k
from projects.ukrmol_co.numerov_increment_study import exact_discrete_scalar, exact_free_k


def boundary(energy, target, coupled):
    k = np.sqrt(energy)
    q = np.tan(k * 18 + np.arctan(target)) / k
    if not coupled:
        return np.array([[q]]), np.array([0.0])
    frame = np.array([[0.8, -0.6], [0.6, 0.8]])
    effective = frame @ np.diag([q, np.tan(k * 18 + np.arctan(0.3)) / k]) @ frame.T
    coupling, kap = np.array([0.12, -0.08]), np.sqrt(1.2 - energy)
    matrix = np.zeros((3, 3))
    matrix[:2, :2] = effective + np.outer(coupling, coupling) * kap / (1 + 0.3 * kap)
    matrix[:2, 2] = matrix[2, :2] = coupling
    matrix[2, 2] = 0.3
    return matrix, np.array([0.0, 0.0, 1.2])


def study():
    rows, checks = [], {}
    begin = time.monotonic()
    models = [
        ("compressed", "0x1.4eecbfb15b574p-2", 320.0, [0.00625, 0.003125, 0.0015625]),
        ("equilibrium", "0x1.9df3b645a1cb3p-3", -55.0, [0.025, 0.0125, 0.00625]),
    ]
    for name, bits, target, steps in models:
        energy = float.fromhex(bits)
        for coupled in [False, True]:
            q, thresholds = boundary(energy, target, coupled)
            continuous = exact_free_k(q, thresholds, energy, 18.0, 80)
            refined = exact_free_k(q, thresholds, energy, 18.0, 110)
            with mp.workdps(80):
                drift = float(max(abs(a - b) for a, b in zip(continuous, refined, strict=True)))
            prefix = f"{name}_{coupled}"
            checks[prefix + "_oracle_precision"] = drift < 1e-20
            reference = np.array(refined.tolist(), float)
            for step, endpoint, order in configurations(steps):
                zero = np.zeros_like(q)
                k, diagnostics = increment_numerov_k(
                    q,
                    thresholds,
                    np.zeros(len(q), int),
                    zero,
                    zero,
                    energy,
                    18.0,
                    endpoint,
                    step,
                    order,
                )
                errors, invariants = differences(k, reference), algebra(k)
                key = f"{prefix}_{step}_{endpoint}_{order}"
                checks[key + "_invariants"] = all(v <= 1e-6 for v in invariants.values())
                if step == steps[-1]:
                    checks[key + "_continuous"] = all(v <= 1e-6 for v in errors.values())
                row = dict(
                    model=name,
                    coupled_closed=coupled,
                    energy_rydberg_hex=bits,
                    target_tangent=target,
                    step_bohr=step,
                    endpoint_bohr=endpoint,
                    tail_order=order,
                    k=k.tolist(),
                    continuous_errors=errors,
                    invariants=invariants,
                    continuous_precision_drift=drift,
                    diagnostics=diagnostics,
                )
                if not coupled:
                    discrete = exact_discrete_scalar(
                        float(q[0, 0]), energy, 18.0, endpoint, step, 80
                    )
                    high = exact_discrete_scalar(float(q[0, 0]), energy, 18.0, endpoint, step, 110)
                    with mp.workdps(80):
                        error = float(abs(mp.mpf(float(k[0, 0])) - high))
                        precision = float(abs(discrete - high))
                        truncation = float(abs(high - refined[0, 0]))
                    row.update(
                        discrete_error=error,
                        discrete_precision_drift=precision,
                        discrete_truncation_error=truncation,
                    )
                    checks[key + "_roundoff"] = error <= 1e-7
                    checks[key + "_discrete_precision"] = precision < 1e-20
                rows.append(row)
    return dict(
        analytic_interior_verdict=all(checks.values()),
        checks=checks,
        rows=rows,
        wall_seconds=time.monotonic() - begin,
        native_surface_or_potential_used=False,
        native_full_k_reconstruction_verdict=False,
        co_complex_pole_verdict=None,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir()
    report = study()
    (args.output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if report["analytic_interior_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
