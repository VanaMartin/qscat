"""Bounded boundary zeros and independent resolvent/channel contour witnesses."""

from __future__ import annotations

import numpy as np
from qscat.units import hartree_to_ev
from scipy.linalg import solve, svd


def matching(q, y):
    return np.eye(len(q)) - q @ y


def pole_difference(left, right):
    """Physical-eV position and FULL-width differences, never requested-grid eV."""
    delta = left - right
    return dict(
        position_ev=float(hartree_to_ev(abs(delta.real))),
        full_width_ev=float(hartree_to_ev(2 * abs(delta.imag))),
    )


def bordered_search(evaluate, seeds, box, limit, observer=None):
    """Complex secant on a frozen border; callback returns Q,Y, not a phase fit."""
    cache = {}
    trace = []

    def point(z):
        if not (box[0] <= z.real <= box[1] and box[2] <= z.imag <= box[3]):
            raise ValueError("Candidate left the frozen search box")
        if z not in cache:
            if len(cache) >= limit:
                raise ValueError("Frozen root evaluation cap exhausted")
            cache[z] = evaluate(z)
        return cache[z]

    z0, z1 = map(complex, seeds)
    first = point(z0)
    m = matching(first["q"], first["y"])
    u, _, vh = svd(m)
    b = u[:, -1]
    c = vh[-1]

    def value(z):
        data = point(z)
        m = matching(data["q"], data["y"])
        n = len(m)
        border = np.zeros((n + 1, n + 1), complex)
        border[:n, :n] = m
        border[:n, n] = b
        border[n, :n] = c
        rhs = np.zeros(n + 1, complex)
        rhs[-1] = 1
        alpha = solve(border, rhs)[-1]
        singular = svd(m, compute_uv=False)
        trace.append(
            dict(
                z=[float(z.real), float(z.imag)],
                alpha_abs=float(abs(alpha)),
                sigma_min=float(singular[-1]),
                sigma_second=float(singular[-2]) if n > 1 else None,
            )
        )
        if observer is not None:
            observer(dict(evaluations=len(cache), trace=list(trace)))
        accepted = abs(alpha) <= 1e-10 and singular[-1] <= 1e-8 and (n == 1 or singular[-2] >= 1e-4)
        return alpha, accepted

    f0, ok = value(z0)
    if ok:
        return z0, dict(evaluations=len(cache), trace=trace)
    f1, ok = value(z1)
    while not ok:
        if f1 == f0:
            raise ValueError("Singular frozen secant step")
        step = -f1 * (z1 - z0) / (f1 - f0)
        if abs(step) > 0.025:
            step *= 0.025 / abs(step)
        z0, f0 = z1, f1
        z1 = z1 + step
        f1, ok = value(z1)
    return z1, dict(evaluations=len(cache), trace=trace)


def projected_resolvent(poles, w, total, y, a, x, q):
    """Projected interior (E-H_eff)^-1 and analytic outgoing T via a Schur solve."""
    p = 1 / (poles - total)
    b = w.T @ (p[:, None] * x)
    rhs = np.column_stack([b, q @ a.T])
    solution = solve(matching(q, y), rhs)
    count = x.shape[1]
    interior = -(x.T @ (p[:, None] * x) + b.T @ y @ solution[:, :count])
    outgoing = a @ solution[:, count:]
    return interior, outgoing


def _ratio(value, reference):
    norm = float(np.linalg.norm(reference))
    return (
        float(np.linalg.norm(value - reference) / norm)
        if norm
        else (0.0 if not np.any(value) else 1e300)
    )


def contour_witness(center, points, interior, outgoing, phases, budget_ev=1e-4):
    """Two nested contour rules; zero interior/open signals must not pass."""
    points = np.asarray(points)
    interior = np.asarray(interior)
    outgoing = np.asarray(outgoing)
    phases = np.asarray(phases)
    if len(points) != 32:
        raise ValueError("Exactly32 ordered contour points are required")
    moments = []
    for indices in [np.arange(0, 32, 2), np.arange(32)]:
        delta = points[indices] - center
        c0 = np.mean(delta[:, None, None] * interior[indices], axis=0)
        c1 = np.mean((delta * points[indices])[:, None, None] * interior[indices], axis=0)
        residue = 2j * np.mean(delta[:, None, None] * outgoing[indices], axis=0)
        norm = float(np.linalg.norm(c0))
        pole = np.vdot(c0, c1) / np.vdot(c0, c0) if norm > 1e-20 else None
        phase = phases[indices]
        jumps = np.angle(np.exp(1j * (np.roll(phase, -1) - phase)))
        moments.append(
            dict(
                c0=c0,
                c1=c1,
                residue=residue,
                pole=pole,
                winding=float(np.sum(jumps) / (2 * np.pi)),
                max_phase_jump=float(np.max(abs(jumps))),
            )
        )
    low, high = moments
    sr = svd(high["c0"], compute_uv=False)
    ss = svd(high["residue"], compute_uv=False)
    rnorm = float(np.linalg.norm(high["c0"]))
    snorm = float(np.linalg.norm(high["residue"]))
    width = -2 * center.imag
    ratios = dict(
        interior_rank_ratio=float(sr[1] / sr[0]) if len(sr) > 1 and sr[0] > 0 else 0.0,
        open_rank_ratio=float(ss[1] / ss[0]) if len(ss) > 1 and ss[0] > 0 else 0.0,
        interior_rule_change=_ratio(low["c0"], high["c0"]),
        open_rule_change=_ratio(low["residue"], high["residue"]),
    )
    errors = pole_difference(high["pole"], center) if high["pole"] is not None else None
    rule_errors = (
        pole_difference(low["pole"], high["pole"])
        if low["pole"] is not None and high["pole"] is not None
        else None
    )
    checks = dict(
        interior_signal=rnorm > 1e-8,
        interior_rank=ratios["interior_rank_ratio"] <= 1e-3,
        open_signal=width > 0 and snorm / width > 1e-6,
        open_rank=ratios["open_rank_ratio"] <= 1e-3,
        interior_rule=ratios["interior_rule_change"] <= 1e-5,
        open_rule=ratios["open_rule_change"] <= 1e-5,
        root_agreement=errors is not None and max(errors.values()) <= budget_ev,
        pole_rule_agreement=rule_errors is not None and max(rule_errors.values()) <= budget_ev,
        winding=all(abs(value["winding"] - 1) <= 1e-6 for value in moments),
        winding_resolution=all(value["max_phase_jump"] <= np.pi / 2 for value in moments),
    )
    report = dict(
        interior_residue_norm=rnorm,
        open_residue_norm_hartree=snorm,
        pole=[float(high["pole"].real), float(high["pole"].imag)]
        if high["pole"] is not None
        else None,
        root_difference=errors,
        rule_pole_difference=rule_errors,
        winding=[value["winding"] for value in moments],
        max_phase_jump=[value["max_phase_jump"] for value in moments],
        **ratios,
        checks=checks,
    )
    return report, moments


def derivative_residue(q, y, a, derivative):
    """Bilinear simple-pole residue; derivative is with respect to incident Hartree."""
    m = matching(q, y)
    u, s, vh = svd(m)
    v = vh[-1].conj()
    ell = u[:, -1].conj()
    denominator = ell @ derivative @ v
    if not np.isfinite(denominator) or abs(denominator) == 0:
        raise ValueError("Vanishing or nonfinite simple-pole derivative")
    boundary = np.outer(v, ell @ q) / denominator
    scattering = 2j * a @ boundary @ a.T
    return scattering, dict(
        derivative_denominator_abs=float(abs(denominator)),
        sigma_min=float(s[-1]),
        boundary_transpose_relative_error=_ratio(boundary.T, boundary),
    )
