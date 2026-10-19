"""Frozen endpoint/spacing Cartesian control against high-precision free oracles."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import mpmath as mp
import numpy as np

from projects.ukrmol_co.native_exterior_attribution import algebra, differences
from projects.ukrmol_co.numerov_increment import increment_numerov_k
from projects.ukrmol_co.numerov_increment_study import (
    exact_discrete_scalar,
    exact_free_k,
    free_boundary,
)


def study():
    rows, checks = [], {}
    begin = time.monotonic()
    for case in ["scalar", "coupled_closed"]:
        q, thresholds, energy, radius = free_boundary(case)
        continuous = exact_free_k(q, thresholds, energy, radius, 80)
        refined = exact_free_k(q, thresholds, energy, radius, 110)
        with mp.workdps(80):
            drift = float(max(abs(a - b) for a, b in zip(continuous, refined, strict=True)))
        checks[case + "_continuous_precision"] = drift < 1e-20
        zero = np.zeros_like(q)
        reference = np.array(refined.tolist(), float)
        for step in [0.0125, 0.00625]:
            for endpoint in [200.0, 400.0, 800.0]:
                k, diagnostics = increment_numerov_k(
                    q,
                    thresholds,
                    np.zeros(len(q), int),
                    zero,
                    zero,
                    energy,
                    radius,
                    endpoint,
                    step,
                    8,
                )
                errors, invariants = differences(k, reference), algebra(k)
                key = f"{case}_{step}_{endpoint}"
                checks[key + "_invariants"] = all(v <= 1e-6 for v in invariants.values())
                if step == 0.00625:
                    checks[key + "_continuous"] = all(v <= 5e-7 for v in errors.values())
                row = dict(
                    case=case,
                    step_bohr=step,
                    endpoint_bohr=endpoint,
                    k=k.tolist(),
                    continuous_errors=errors,
                    invariants=invariants,
                    diagnostics=diagnostics,
                    continuous_oracle_precision_drift=drift,
                )
                if case == "scalar":
                    discrete = exact_discrete_scalar(float(q[0, 0]), energy, radius, endpoint, step)
                    refined_discrete = exact_discrete_scalar(
                        float(q[0, 0]), energy, radius, endpoint, step, 110
                    )
                    with mp.workdps(80):
                        rounding = float(abs(mp.mpf(float(k[0, 0])) - refined_discrete))
                        truncation = float(abs(refined_discrete - refined[0, 0]))
                        discrete_drift = float(abs(discrete - refined_discrete))
                    row.update(
                        discrete_error=rounding,
                        discrete_truncation_error=truncation,
                        discrete_precision_drift=discrete_drift,
                    )
                    checks[key + "_roundoff"] = rounding <= 1e-7
                    checks[key + "_discrete_precision"] = discrete_drift < 1e-20
                rows.append(row)
    return dict(
        analytic_endpoint_verdict=all(checks.values()),
        checks=checks,
        rows=rows,
        wall_seconds=time.monotonic() - begin,
        native_data_used=False,
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
    if report["analytic_endpoint_verdict"] is not True:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
