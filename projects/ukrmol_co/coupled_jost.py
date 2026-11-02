"""Exact coupled endpoint control; no CO data or production-method import."""

from __future__ import annotations

import argparse
import itertools
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


def superpotential(z, g):
    k = np.asarray([1, np.sqrt(3)])
    b = np.exp(-k * z)
    denominator = 1 / g + np.sum((1 - np.exp(-2 * k * z)) / (2 * k))
    correction = np.outer(b, b) / denominator
    w = np.diag(k) + correction
    derivative = -(k[:, None] + k[None, :]) * correction - correction @ correction
    return w, derivative


def exact_maps(z, g):
    w, derivative = superpotential(z, g)
    s = np.asarray([1j, -1])
    matrix = np.diag(s) - w
    y = np.linalg.solve(matrix.T, (-derivative + matrix * s).T).T
    a = (1j - 1) * np.exp(-1j * z) * np.linalg.solve(matrix.T, np.asarray([1, 0]))
    return y, a[None, :]


def potential(z, g):
    return np.diag([0, 2]) - 2 * superpotential(z, g)[1]


def metrics(y, a, reference):
    wall = np.linalg.lstsq(a.conj().T, a.T, rcond=None)[0]
    return dict(
        map_error=float(max(abs(y - reference[0]).max(), abs(a - reference[1]).max())),
        flux_defect=float(abs((y - y.conj().T) / (2j) - a.conj().T @ a).max()),
        wronskian_residual=float(abs(a.conj().T @ wall - a.T).max()),
        s_identity_error=float(abs(wall[0, 0] - 1)),
    )


def integrate(contract, angle, length, rtol, initialization):
    rotation = np.exp(1j * np.deg2rad(angle))
    radius, g = contract["radius_bohr"], contract["g"]
    endpoint = radius + length * rotation
    if initialization == "exact":
        y, a = exact_maps(endpoint, g)
    elif initialization == "free":
        y = np.diag([1j, -1])
        a = np.asarray([[np.exp(-1j * endpoint), 0]])
    else:
        raise ValueError("Unknown endpoint prescription")
    calls, started = 0, time.monotonic()

    def rhs(t, state):
        nonlocal calls
        calls += 1
        if (
            calls > contract["rhs_cap"]
            or time.monotonic() - started > contract["per_integration_seconds"]
        ):
            raise TimeoutError("Frozen coupled integration cap")
        y, a = state[:4].reshape(2, 2), state[4:].reshape(1, 2)
        dy = potential(radius + t * rotation, g) - np.eye(2) - y @ y
        return rotation * np.concatenate([dy.ravel(), (-a @ y).ravel()])

    result = solve_ivp(
        rhs,
        (length, 0),
        np.concatenate([y.ravel(), a.ravel()]),
        method="DOP853",
        rtol=rtol,
        atol=rtol * contract["atol_factor"],
        t_eval=[0],
    )
    if not result.success or not np.isfinite(result.y).all():
        raise ValueError("Coupled integration failed: " + result.message)
    state = result.y[:, -1]
    return state[:4].reshape(2, 2), state[4:].reshape(1, 2), calls


def encoded(array):
    return np.stack([array.real, array.imag], axis=-1).tolist()


def run(contract, output):
    output.mkdir(exist_ok=False)
    started = time.monotonic()
    result = dict(
        scientific_exit=2,
        coupled_control_verdict=None,
        rows=[],
        co_maps=0,
        native_calls=0,
        automatic_successor=False,
        full_sweep_released=False,
    )

    def save():
        result["science_wall_seconds"] = time.monotonic() - started
        for name in ["memory.peak", "cpu.stat"]:
            path = Path("/sys/fs/cgroup") / name
            if path.exists():
                result["cgroup_" + name] = path.read_text()
        partial = output / "result.partial"
        partial.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        partial.replace(output / "result.json")

    try:
        reference = exact_maps(contract["radius_bohr"], contract["g"])
        for angle, length, rtol, initialization in itertools.product(
            contract["angles_degrees"],
            contract["lengths_bohr"],
            contract["rtols"],
            contract["initializations"],
        ):
            if time.monotonic() - started > contract["science_seconds"]:
                raise TimeoutError("Frozen combined science cap")
            case = dict(
                angle_degrees=angle, length_bohr=length, rtol=rtol, initialization=initialization
            )
            result["attempted_case"] = case
            save()
            y, a, calls = integrate(contract, angle, length, rtol, initialization)
            row = dict(
                **case,
                y=encoded(y),
                a=encoded(a),
                rhs_evaluations=calls,
                **metrics(y, a, reference),
            )
            result["rows"].append(row)
            if rtol == contract["rtols"][-1] and (
                initialization == "exact" or length == contract["lengths_bohr"][-1]
            ):
                if any(
                    row[name] > contract[budget]
                    for name, budget in [
                        ("map_error", "map_error_budget"),
                        ("flux_defect", "flux_budget"),
                        ("wronskian_residual", "wronskian_budget"),
                    ]
                ):
                    raise ValueError("Frozen coupled map/flux/Wronskian gate failed")
            save()
        reductions = []
        for angle in contract["angles_degrees"]:
            rows = [
                row
                for row in result["rows"]
                if row["angle_degrees"] == angle
                and row["rtol"] == contract["rtols"][-1]
                and row["initialization"] == "free"
            ]
            short, long = rows[0]["map_error"], rows[-1]["map_error"]
            if (
                short <= contract["minimum_resolved_tail_error"]
                or short < contract["minimum_tail_error_reduction"] * long
            ):
                raise ValueError("No resolved convergent omitted-tail error")
            reductions.append(
                dict(
                    angle_degrees=angle,
                    short_error=short,
                    long_error=long,
                    reduction=short / long if long else None,
                )
            )
        if len(result["rows"]) != contract["expected_integrations"]:
            raise ValueError("Wrong integration count")
        result.update(
            scientific_exit=0, coupled_control_verdict=True, tail_error_reductions=reductions
        )
    except Exception as error:
        result["error"] = repr(error)
    finally:
        save()
    return result["scientific_exit"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("contract", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    raise SystemExit(run(json.loads(args.contract.read_text()), args.output))
