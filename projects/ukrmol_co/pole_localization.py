"""Finite sufficient-decrease safeguards for the existing local boundary zero."""

from __future__ import annotations

import warnings

import numpy as np
from scipy.linalg import LinAlgWarning, solve, svd

from projects.ukrmol_co.pole_reference import bordered_search, matching


def configured_search(contract):
    """Only an explicit frozen selection changes historical local-search behavior."""
    name = contract.get("root_search", "frozen-secant")
    if name == "frozen-secant":
        if "safeguard_controls" in contract:
            raise ValueError("Safeguard controls require explicit safeguarded search")
        return bordered_search, {}
    if name == "safeguarded-newton-v1":
        controls = contract.get("safeguard_controls")
        if not isinstance(controls, dict) or set(controls) != {
            "max_step",
            "max_backtracks",
            "armijo",
        }:
            raise ValueError("Complete explicit safeguard controls are required")
        return safeguarded_search, controls.copy()
    raise ValueError("Unsupported frozen root-search selection")


def safeguarded_search(
    evaluate,
    seeds,
    box,
    limit,
    *,
    max_step=0.025,
    max_backtracks=8,
    armijo=1e-4,
    observer=None,
):
    """Refresh a local derivative and require decrease; never claim global localization.

    Both seeds define one starting center and a positive real probe increment.
    Every unique model evaluation, including rejected trials, consumes the cap.
    Model/ODE failures propagate unchanged; safeguards do not repair failed maps.
    """
    seeds = np.asarray(seeds, dtype=complex)
    box = np.asarray(box, dtype=float)
    if (
        seeds.shape != (2,)
        or box.shape != (4,)
        or not np.isfinite(seeds).all()
        or not np.isfinite(box).all()
        or box[0] >= box[1]
        or box[2] >= box[3]
        or type(limit) is not int
        or limit < 1
        or type(max_backtracks) is not int
        or max_backtracks < 0
        or not np.isfinite([max_step, armijo]).all()
        or max_step <= 0
        or not 0 < armijo < 1
    ):
        raise ValueError("Finite seeds/box and positive bounded search controls required")
    delta = seeds[1] - seeds[0]
    if delta.real <= 0 or abs(delta.imag) > 1e-15:
        raise ValueError("The second seed must define a positive real derivative probe")
    h = float(delta.real)
    cache = {}
    attempted = set()
    report = dict(evaluations=0, trace=[], centers=[], probes=[], trials=[])
    b = c = None

    def notify():
        if observer is not None:
            observer(report)

    def inside(z):
        return box[0] <= z.real <= box[1] and box[2] <= z.imag <= box[3]

    def coordinates(z):
        return [float(z.real), float(z.imag)]

    def point(z):
        nonlocal b, c
        if not inside(z):
            raise ValueError("Candidate left the frozen search box")
        if z in cache:
            return cache[z]
        if len(attempted) >= limit:
            raise ValueError("Frozen root evaluation cap exhausted")
        attempted.add(z)
        report["evaluations"] = len(attempted)
        row = dict(z=coordinates(z), completed=False)
        report["trace"].append(row)
        notify()
        data = evaluate(z)
        q, y = (np.asarray(data[key]) for key in ("q", "y"))
        if (
            q.ndim != 2
            or q.shape[0] == 0
            or q.shape[0] != q.shape[1]
            or y.shape != q.shape
            or not np.isfinite(q).all()
            or not np.isfinite(y).all()
        ):
            raise ValueError("Finite compatible Q/Y matrices required")
        m = matching(q, y)
        u, singular, vh = svd(m)
        if b is None:
            b, c = u[:, -1], vh[-1]
        n = len(m)
        border = np.zeros((n + 1, n + 1), complex)
        border[:n, :n], border[:n, n], border[n, :n] = m, b, c
        rhs = np.zeros(n + 1, complex)
        rhs[-1] = 1
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", LinAlgWarning)
                alpha = solve(border, rhs)[-1]
        except (np.linalg.LinAlgError, LinAlgWarning) as error:
            raise ValueError("Singular or ill-conditioned frozen border") from error
        if not np.isfinite(alpha):
            raise ValueError("Nonfinite frozen-border scalar")
        accepted = abs(alpha) <= 1e-10 and singular[-1] <= 1e-8 and (n == 1 or singular[-2] >= 1e-4)
        cache[z] = alpha, accepted
        row.update(
            dict(
                completed=True,
                alpha=coordinates(alpha),
                alpha_abs=float(abs(alpha)),
                sigma_min=float(singular[-1]),
                sigma_second=float(singular[-2]) if n > 1 else None,
                candidate=bool(accepted),
            )
        )
        notify()
        return alpha, accepted

    def finish(z):
        report["candidate"] = coordinates(z)
        notify()
        return z, report

    center = complex(seeds[0])
    value, accepted = point(center)
    report["centers"].append(dict(z=coordinates(center), alpha_abs=float(abs(value))))
    if accepted:
        return finish(center)
    while True:
        increment = h if inside(center + h) else -h
        probe = center + increment
        if not inside(probe) or probe == center:
            raise ValueError("No resolved in-box derivative probe")
        report["probes"].append(
            dict(
                center=coordinates(center),
                z=coordinates(probe),
                direction="forward" if increment > 0 else "backward",
            )
        )
        probe_value, accepted = point(probe)
        if accepted:
            return finish(probe)
        derivative = (probe_value - value) / increment
        if not np.isfinite(derivative) or abs(derivative) == 0:
            raise ValueError("Zero or nonfinite local derivative")
        full_step = -value / derivative
        if not np.isfinite(full_step):
            raise ValueError("Nonfinite local Newton step")
        clipped = min(1.0, max_step / abs(full_step)) if abs(full_step) else 1.0
        moved = False
        for backtrack in range(max_backtracks + 1):
            fraction = clipped * 2.0**-backtrack
            trial = center + fraction * full_step
            row = dict(
                z=coordinates(trial), backtrack=backtrack, fraction=float(fraction), accepted=False
            )
            report["trials"].append(row)
            if not inside(trial) or trial == center:
                row["reason"] = "outside box" if not inside(trial) else "stagnation"
                notify()
                continue
            trial_value, accepted = point(trial)
            row["alpha_abs"] = float(abs(trial_value))
            if accepted or (
                abs(trial_value) < abs(value)
                and abs(trial_value) <= (1 - armijo * fraction) * abs(value)
            ):
                row.update(accepted=True, reason="candidate" if accepted else "sufficient decrease")
                center, value = trial, trial_value
                report["centers"].append(dict(z=coordinates(center), alpha_abs=float(abs(value))))
                notify()
                if accepted:
                    return finish(center)
                moved = True
                break
            row["reason"] = "insufficient decrease"
            notify()
        if not moved:
            raise ValueError("No decreasing in-box step within the frozen backtracking cap")
