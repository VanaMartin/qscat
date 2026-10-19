"""Fixed analytic large-K qualification before any saved CO matrix is used."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np

from projects.ukrmol_co.exterior_numerov import numerov_k
from projects.ukrmol_co.native_exterior_attribution import algebra, differences
from projects.ukrmol_co.numerov_increment import increment_numerov_k


def free_boundary(case):
    """Binary64 boundary inputs constructed independently of either propagator."""
    energy, radius = 0.0735, 18.0
    k = np.sqrt(energy)
    if case == "scalar":
        q = np.array([[np.tan(k * radius + np.arctan(-185.0)) / k]])
        return q, np.array([0.0]), energy, radius
    if case != "coupled_closed":
        raise ValueError("Only two frozen free-boundary cases are defined")
    frame = np.array([[0.8, -0.6], [0.6, 0.8]])
    effective = frame @ np.diag(np.tan(k * radius + np.arctan([-185.0, 0.3])) / k) @ frame.T
    coupling, closed_q, kap = np.array([0.12, -0.08]), 0.3, np.sqrt(1.2 - energy)
    q = np.zeros((3, 3))
    q[:2, :2] = effective + np.outer(coupling, coupling) * kap / (1 + closed_q * kap)
    q[:2, 2] = q[2, :2] = coupling
    q[2, 2] = closed_q
    return q, np.array([0.0, 0.0, 1.2]), energy, radius


def exact_free_k(q, thresholds, energy, radius, digits):
    """High-precision direct boundary matching; no numerical propagation or tails."""
    with mp.workdps(digits):
        q = mp.matrix(q.tolist())
        n = len(thresholds)
        opened = [i for i, t in enumerate(thresholds) if t < energy]
        closed = [i for i, t in enumerate(thresholds) if t > energy]
        f, g, fp, gp = [mp.zeros(n, len(opened)) for _ in range(4)]
        decaying, dp = mp.zeros(n, len(closed)), mp.zeros(n, len(closed))
        for j, i in enumerate(opened):
            k = mp.sqrt(mp.mpf(energy) - mp.mpf(float(thresholds[i])))
            f[i, j], g[i, j] = mp.sin(k * radius) / mp.sqrt(k), mp.cos(k * radius) / mp.sqrt(k)
            fp[i, j], gp[i, j] = mp.sqrt(k) * mp.cos(k * radius), -mp.sqrt(k) * mp.sin(k * radius)
        for j, i in enumerate(closed):
            decaying[i, j], dp[i, j] = 1, -mp.sqrt(mp.mpf(float(thresholds[i])) - mp.mpf(energy))
        lhs, rhs = mp.zeros(n), -(f - q * fp)
        irregular, decaying = g - q * gp, decaying - q * dp
        for i in range(n):
            for j in range(len(opened)):
                lhs[i, j] = irregular[i, j]
            for j in range(len(closed)):
                lhs[i, len(opened) + j] = decaying[i, j]
        answer = mp.zeros(len(opened))
        for j in range(len(opened)):
            column = mp.lu_solve(lhs, rhs[:, j])
            for i in range(len(opened)):
                answer[i, j] = column[i]
        return answer


def exact_discrete_scalar(q, energy, radius, endpoint, step, digits=80):
    """Closed-form discrete transfer with the original sixth-order Taylor maps."""
    with mp.workdps(digits):
        q, e, b, h = map(mp.mpf, [q, energy, endpoint, step])
        n = round((endpoint - radius) / step)
        w, k = -e, mp.sqrt(e)
        theta = mp.acos((12 / (1 - h**2 * w / 12) - 10) / 2)

        def taylor(value, derivative, s):
            derivatives = [value, derivative]
            for j in range(5):
                derivatives.append(w * derivatives[j])
            return mp.fsum(s**j * value / mp.factorial(j) for j, value in enumerate(derivatives))

        u0, u1 = q, taylor(q, mp.mpf(1), h)
        un = (mp.sin(n * theta) * u1 - mp.sin((n - 1) * theta) * u0) / mp.sin(theta)
        before = (mp.sin((n - 1) * theta) * u1 - mp.sin((n - 2) * theta) * u0) / mp.sin(theta)
        coefficient = taylor(mp.mpf(1), mp.mpf(0), -h) - before / un
        derivative = taylor(mp.mpf(0), mp.mpf(1), -h)
        f, g, fp, gp = mp.sin(k * b), mp.cos(k * b), k * mp.cos(k * b), -k * mp.sin(k * b)
        return -(coefficient * f + derivative * fp) / (coefficient * g + derivative * gp)


def study():
    rows, checks = [], {}
    begin = time.monotonic()
    for case in ["scalar", "coupled_closed"]:
        q, thresholds, energy, radius = free_boundary(case)
        oracle = exact_free_k(q, thresholds, energy, radius, 80)
        refined = exact_free_k(q, thresholds, energy, radius, 110)
        with mp.workdps(80):
            drift = float(
                max(
                    abs(oracle[i, j] - refined[i, j])
                    for i in range(oracle.rows)
                    for j in range(oracle.cols)
                )
            )
        checks[case + "_oracle_precision"] = drift < 1e-20
        reference = np.array(oracle.tolist(), float)
        zero = np.zeros_like(q)
        for step in [0.025, 0.0125, 0.00625]:
            for name, function in [("original", numerov_k), ("increment", increment_numerov_k)]:
                start = time.monotonic()
                measured, diagnostics = function(
                    q, thresholds, np.zeros(len(q), int), zero, zero, energy, radius, 400.0, step, 8
                )
                row = {
                    "case": case,
                    "formulation": name,
                    "step_bohr": step,
                    "endpoint_bohr": 400,
                    "oracle_precision_drift": drift,
                    "k": measured.tolist(),
                    "continuous_errors": differences(measured, reference),
                    "invariants": algebra(measured),
                    "wall_seconds": time.monotonic() - start,
                    "diagnostics": diagnostics,
                }
                if case == "scalar":
                    discrete = exact_discrete_scalar(float(q[0, 0]), energy, radius, 400.0, step)
                    discrete110 = exact_discrete_scalar(
                        float(q[0, 0]), energy, radius, 400.0, step, 110
                    )
                    with mp.workdps(80):
                        discrete_error = float(abs(mp.mpf(float(measured[0, 0])) - discrete))
                        truncation_error = float(abs(discrete - oracle[0, 0]))
                        discrete_drift = float(abs(discrete - discrete110))
                    row.update(
                        discrete_error=discrete_error,
                        discrete_truncation_error=truncation_error,
                        discrete_oracle_precision_drift=discrete_drift,
                    )
                    if name == "increment":
                        checks[f"scalar_discrete_precision_{step}"] = discrete_drift < 1e-20
                        checks[f"scalar_roundoff_{step}"] = discrete_error <= 1e-7
                if name == "increment":
                    checks[f"{case}_invariants_{step}"] = all(
                        x <= 1e-6 for x in row["invariants"].values()
                    )
                    if step == 0.00625:
                        checks[case + "_finest_continuous"] = all(
                            x <= 1e-6 for x in row["continuous_errors"].values()
                        )
                rows.append(row)
    return {
        "analytic_qualification_verdict": all(checks.values()),
        "checks": checks,
        "rows": rows,
        "wall_seconds": time.monotonic() - begin,
        "native_data_used": False,
        "original_reconstruction_verdict": False,
        "original_reconstruction_exits": [1, 1],
        "new_inner_solves": 0,
        "full_sweep_released": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir()
    result = study()
    (args.output / "report.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if result["analytic_qualification_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
