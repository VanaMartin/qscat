"""Finite exact dual-row transport control; never imports or reads a CO model."""

from __future__ import annotations

import argparse
import itertools
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


def exact_row(k, kappa, angle, length, epsilon):
    z = length * np.exp(1j * np.deg2rad(angle))
    return np.asarray([1, epsilon * np.exp(-(kappa + 1j * k) * z)], complex)


def integrate(k, kappa, angle, length, epsilon, rtol, contract):
    rotation = np.exp(1j * np.deg2rad(angle))
    z = length * rotation
    initial = np.exp(-1j * k * z) * np.asarray([1, epsilon], complex)
    started = time.monotonic()
    calls = 0

    def rhs(t, a):
        nonlocal calls
        calls += 1
        if (
            calls > contract["rhs_cap"]
            or time.monotonic() - started > contract["per_integration_seconds"]
        ):
            raise TimeoutError("Frozen toy integration cap exhausted")
        return rotation * np.asarray([-1j * k * a[0], kappa * a[1]])

    result = solve_ivp(
        rhs,
        (length, 0),
        initial,
        method="DOP853",
        rtol=rtol,
        atol=rtol * contract["atol_factor"],
        t_eval=[0],
    )
    if not result.success or not np.isfinite(result.y).all():
        raise ValueError("Toy integration failed: " + result.message)
    return result.y[:, -1], calls


def observables(a, k, kappa):
    a = np.asarray(a)[None, :]
    wall = np.linalg.lstsq(a.conj().T, a.T, rcond=None)[0]
    flux = (np.diag([1j * k, -kappa]) - np.diag([1j * k, -kappa]).conj().T) / (2j) - a.conj().T @ a
    return dict(s_identity_error=float(abs(wall[0, 0] - 1)), flux_defect=float(abs(flux).max()))


def run(contract, output):
    output.mkdir(exist_ok=False)
    started = time.monotonic()
    report = dict(
        scientific_exit=2,
        toy_verdict=None,
        rows=[],
        co_maps=0,
        native_calls=0,
        pole_searches=0,
        automatic_successor=False,
        full_sweep_released=False,
    )

    def save():
        report["wall_seconds"] = time.monotonic() - started
        for name in ["memory.peak", "cpu.stat"]:
            path = Path("/sys/fs/cgroup") / name
            if path.exists():
                report["cgroup_" + name] = path.read_text()
        partial = output / "result.partial"
        partial.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        partial.replace(output / "result.json")

    try:
        k = np.sqrt(contract["open_energy_rydberg"])
        kappa = np.sqrt(contract["closed_threshold_rydberg"] - contract["open_energy_rydberg"])
        finest = []
        for angle, length, epsilon, rtol in itertools.product(
            contract["angles_degrees"],
            contract["lengths_bohr"],
            contract["endpoint_contamination"],
            contract["rtols"],
        ):
            if time.monotonic() - started > contract["total_seconds"]:
                raise TimeoutError("Frozen combined toy cap exhausted")
            report["attempted_case"] = dict(
                angle_degrees=angle, length_bohr=length, epsilon=epsilon, rtol=rtol
            )
            save()
            actual, calls = integrate(k, kappa, angle, length, epsilon, rtol, contract)
            expected = exact_row(k, kappa, angle, length, epsilon)
            scaled = float(
                np.max(abs(actual - expected) / (contract["scaled_error_floor"] + abs(expected)))
            )
            row = dict(
                angle_degrees=angle,
                length_bohr=length,
                epsilon=epsilon,
                rtol=rtol,
                actual=[[float(v.real), float(v.imag)] for v in actual],
                exact=[[float(v.real), float(v.imag)] for v in expected],
                contamination_ratio=float(abs(actual[1]) / epsilon) if epsilon else None,
                scaled_error=scaled,
                rhs_evaluations=calls,
                **observables(actual, k, kappa),
            )
            report["rows"].append(row)
            if rtol == contract["rtols"][-1]:
                finest.append(row)
                if scaled > contract["finest_scaled_error_budget"] or (
                    epsilon == 0 and abs(actual[1]) > contract["zero_injection_closed_budget"]
                ):
                    raise ValueError("Exact toy oracle gate failed")
            save()
        negative = [
            r
            for r in finest
            if r["epsilon"]
            and r["s_identity_error"] <= contract["s_negative_control_budget"]
            and r["flux_defect"] > contract["flux_negative_control_budget"]
        ]
        if time.monotonic() - started > contract["total_seconds"]:
            raise TimeoutError("Frozen combined toy cap exhausted")
        if len(report["rows"]) != contract["expected_integrations"] or not negative:
            raise ValueError("Missing exact count or intended S-pass/flux-fail control")
        report.update(
            scientific_exit=0,
            toy_verdict=True,
            negative_control_cases=len(negative),
            finest_maximum_scaled_error=max(r["scaled_error"] for r in finest),
        )
    except Exception as error:
        report["error"] = repr(error)
    finally:
        save()
    return report["scientific_exit"]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("contract", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    raise SystemExit(run(json.loads(args.contract.read_text()), args.output))
