"""Small finite-model bound matching; no molecular readers or execution paths.

Energies are incident/total Hartree; exterior thresholds and C2/C3 use the
native Rydberg radial equation. A returned root is not a tail-error enclosure,
electronic-accuracy certificate or qualified molecular bound state.
"""

from dataclasses import dataclass
from time import monotonic

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh, solve
from scipy.optimize import brentq

from projects.ukrmol_co.outgoing_map import _tail_series


@dataclass(frozen=True)
class BoundCaps:
    """Cooperative numerical caps, deliberately sized for synthetic controls."""

    timeout_seconds: float = 30.0
    max_rhs: int = 20000
    max_evaluations: int = 96
    max_channels: int = 32
    max_poles: int = 256

    def __post_init__(self):
        if not np.isfinite(self.timeout_seconds) or self.timeout_seconds <= 0:
            raise ValueError("A finite positive timeout cap is required")
        for value in [self.max_rhs, self.max_evaluations, self.max_channels, self.max_poles]:
            if type(value) is not int or value < 1:
                raise ValueError("Positive integer numerical caps are required")


_DEFAULT_CAPS = BoundCaps()


def _real(value, name):
    if np.iscomplexobj(value):
        raise ValueError(f"{name} must be real and finite")
    result = np.asarray(value, dtype=float)
    if not np.isfinite(result).all():
        raise ValueError(f"{name} must be real and finite")
    return result


def _negative_energy(value):
    energy = _real(value, "Incident energy")
    if energy.ndim != 0 or energy >= 0:
        raise ValueError("A strictly negative incident Hartree energy is required")
    return float(energy)


def _symmetric(matrix, name):
    if not np.allclose(matrix, matrix.T, atol=1e-8, rtol=0):
        raise ValueError(f"{name} must be real symmetric within numerical tolerance")


def decaying_dtn(
    incident_hartree,
    thresholds_rydberg,
    c2,
    c3,
    radius_bohr,
    *,
    endpoint_bohr,
    tail_order=8,
    rtol=1e-11,
    caps=_DEFAULT_CAPS,
):
    """Integrate an all-closed real-axis DtN inward from finite tail data.

    The tail recurrence is shared with the outgoing implementation; independent
    analytic/BVP controls, not agreement of those two routes, validate it.
    Diagnostics expose the raw transpose defect and last tail term, without
    projecting away defects or interpreting them as a rigorous error bound.
    Riccati/endpoint chart failures stop; no representation-switch fallback.
    """
    started = monotonic()
    energy = _negative_energy(incident_hartree)
    thresholds = _real(thresholds_rydberg, "Thresholds").copy()
    c2, c3 = _real(c2, "C2").copy(), _real(c3, "C3").copy()
    n = len(thresholds) if thresholds.ndim == 1 else 0
    if not 1 <= n <= caps.max_channels:
        raise ValueError("Channel dimension exceeds the declared cap or is empty")
    if np.any(thresholds < 0) or any(m.shape != (n, n) for m in [c2, c3]):
        raise ValueError("Nonnegative thresholds and compatible square multipoles required")
    _symmetric(c2, "C2")
    _symmetric(c3, "C3")
    scalars = _real([radius_bohr, endpoint_bohr, rtol], "Exterior settings")
    radius, endpoint, rtol = scalars
    if not 0 < radius < endpoint or not 0 < rtol < 1:
        raise ValueError("A positive real exterior interval and 0 < rtol < 1 are required")
    if type(tail_order) is not int or not 0 <= tail_order <= 16:
        raise ValueError("Tail order must be an integer from zero through sixteen")
    kappa = np.sqrt(thresholds - 2 * energy)
    coefficients = _tail_series(thresholds, c2, c3, -kappa, tail_order).real
    powers = endpoint ** -np.arange(tail_order + 1, dtype=float)
    series = np.einsum("p,pij->ij", powers, coefficients)
    derivative = -np.einsum(
        "p,pij->ij", np.arange(tail_order + 1) * powers / endpoint, coefficients
    )
    if not np.isfinite(series).all() or np.linalg.cond(series) > 1e10:
        raise ValueError("Singular or unresolved decaying tail chart")
    initial = solve(series.T, (derivative - series * kappa).T).T
    diagonal = np.diag(thresholds - 2 * energy)
    calls = 0

    def check_cap():
        if monotonic() - started >= caps.timeout_seconds:
            raise TimeoutError("Decaying-map time cap reached")

    def rhs(radius, flat):
        nonlocal calls
        calls += 1
        check_cap()
        if calls > caps.max_rhs:
            raise TimeoutError("Decaying-map RHS cap reached")
        y = flat.reshape(n, n)
        if not np.isfinite(y).all() or np.max(abs(y)) > 1e8:
            raise ValueError("Unresolved DtN chart or divergent integration")
        return (diagonal + c2 / radius**2 + c3 / radius**3 - y @ y).ravel()

    integrated = solve_ivp(
        rhs,
        (endpoint, radius),
        initial.ravel(),
        method="DOP853",
        rtol=rtol,
        atol=rtol * 0.01,
        t_eval=[radius],
    )
    check_cap()
    if not integrated.success or not np.isfinite(integrated.y).all():
        raise ValueError("Decaying-map integration failed; bound result remains unset")
    y = integrated.y[:, -1].reshape(n, n)
    _symmetric(y, "Decaying DtN")
    return y, dict(
        decaying=True,
        rhs_evaluations=calls,
        wall_seconds=monotonic() - started,
        transpose_defect=float(np.max(abs(y - y.T))),
        endpoint_bohr=float(endpoint),
        tail_order=tail_order,
        tail_last_term_norm=float(np.linalg.norm(coefficients[-1] * powers[-1])),
        total_error_bound=None,
    )


def find_bound_state(
    poles_hartree,
    amplitudes,
    ground_hartree,
    negative_interval_hartree,
    decaying_map,
    *,
    mode_index=0,
    energy_tolerance_hartree=1e-12,
    residual_tolerance_hartree=1e-8,
    isolation_gap_hartree=1e-6,
    finite_box_boundary=None,
    caps=_DEFAULT_CAPS,
):
    """Bracket one preselected ordered eigenvalue of D-T-W Y W.T.

    The pure callback returns (Y, diagnostics) and must enforce its own finite
    work limits. Time checks occur before/after it, not inside arbitrary caller
    code; an external owner is still needed to preempt a blocked callback.
    Only a local finite-model root is returned. Search tolerance is not a total
    energy uncertainty or a proof of global root count/state assignment.
    An explicit finite_box_boundary opts into a tagged Dirichlet/Neumann map;
    that result is finite-box-root-only, not an infinite-domain bound claim.
    """
    started = monotonic()
    if finite_box_boundary not in [None, "dirichlet", "neumann"]:
        raise ValueError("Finite-box boundary must be dirichlet or neumann")
    poles = _real(poles_hartree, "Inner poles").copy()
    w = _real(amplitudes, "Boundary amplitudes").copy()
    ground = _real(ground_hartree, "Ground total")
    if ground.ndim != 0 or poles.ndim != 1 or not 1 <= len(poles) <= caps.max_poles:
        raise ValueError("Ground scalar and nonempty inner spectrum within its cap required")
    if w.ndim != 2 or w.shape[0] != len(poles) or not 1 <= w.shape[1] <= caps.max_channels:
        raise ValueError("Boundary amplitude/channel dimensions exceed declared caps")
    if type(mode_index) is not int or not 0 <= mode_index < len(poles):
        raise ValueError("A valid fixed inner eigenvalue index is required")
    interval = _real(negative_interval_hartree, "Search interval")
    if interval.shape != (2,) or not interval[0] < interval[1] < 0:
        raise ValueError("An ordered strictly negative energy interval is required")
    tolerances = _real(
        [energy_tolerance_hartree, residual_tolerance_hartree, isolation_gap_hartree],
        "Root tolerances",
    )
    if np.any(tolerances <= 0) or interval[1] >= -energy_tolerance_hartree:
        raise ValueError("Positive tolerances and a resolved negative search interval required")
    evaluations = 0
    cache = {}

    def cap_and_inputs():
        if monotonic() - started >= caps.timeout_seconds:
            raise TimeoutError("Bound-root time cap reached")
        if not np.array_equal(np.asarray(poles_hartree), poles) or not np.array_equal(
            np.asarray(amplitudes), w
        ):
            raise ValueError("Inner inputs changed during bound search")

    def evaluate(energy, *, repeat=False):
        nonlocal evaluations
        cap_and_inputs()
        if not interval[0] <= energy <= interval[1]:
            raise ValueError("Bound candidate left the frozen negative interval")
        if energy in cache and not repeat:
            return cache[energy]
        if evaluations >= caps.max_evaluations:
            raise TimeoutError("Bound-root evaluation cap reached")
        evaluations += 1
        y, report = decaying_map(energy)
        cap_and_inputs()
        y = _real(y, "Decaying DtN")
        if y.shape != (w.shape[1], w.shape[1]):
            raise ValueError("Decaying DtN channel dimension differs")
        if finite_box_boundary is None:
            if (
                not report.get("decaying", True)
                or report.get("finite_box", False)
                or report.get("boundary") in ["dirichlet", "neumann"]
            ):
                raise ValueError("A physical decaying map is required")
        elif (
            report.get("boundary") != finite_box_boundary
            or report.get("finite_box") is not True
            or report.get("decaying") is not False
        ):
            raise ValueError("Finite-box callback boundary tags differ from explicit opt-in")
        _symmetric(y, "Decaying DtN")
        if repeat and not np.allclose(y, cache[energy][3], atol=1e-12, rtol=0):
            raise ValueError("Decaying map changed on identical-input replay")
        matrix = np.diag(poles - float(ground) - energy) - w @ y @ w.T
        _symmetric(matrix, "Inner matching matrix")
        values, vectors = eigh(matrix)
        cap_and_inputs()
        data = (values, vectors, matrix, y.copy(), dict(report))
        cache[energy] = data
        return data

    def scalar(energy):
        return float(evaluate(energy)[0][mode_index])

    left, right = map(float, interval)
    fleft, fright = scalar(left), scalar(right)
    if fleft != 0 and fright != 0 and np.signbit(fleft) == np.signbit(fright):
        raise ValueError("The selected finite-model eigenvalue has no sign bracket")
    energy = brentq(
        scalar,
        left,
        right,
        xtol=energy_tolerance_hartree,
        rtol=4 * np.finfo(float).eps,
        maxiter=caps.max_evaluations,
    )
    values, vectors, matrix, _, report = evaluate(energy, repeat=True)
    vector = vectors[:, mode_index]
    residual = float(np.linalg.norm(matrix @ vector))
    if (
        residual > residual_tolerance_hartree
        or abs(values[mode_index]) > residual_tolerance_hartree
    ):
        raise ValueError("Bound root residual fails; possible DtN chart pole or false zero")
    peers = np.delete(values, mode_index)
    gap = float(np.min(abs(peers - values[mode_index]))) if len(peers) else float("inf")
    if gap < isolation_gap_hartree:
        raise ValueError("Bound candidate is not a simple isolated finite-model root")
    return dict(
        energy_hartree=float(energy),
        total_energy_hartree=float(ground) + energy,
        inner_coefficients=vector,
        residual_hartree=residual,
        matching_transpose_defect_hartree=float(np.max(abs(matrix - matrix.T))),
        isolation_gap_hartree=gap,
        search_tolerance_hartree=energy_tolerance_hartree,
        total_uncertainty_hartree=None,
        qualification=(
            "finite-model-root-only" if finite_box_boundary is None else "finite-box-root-only"
        ),
        map_diagnostics=report,
        map_evaluations=evaluations,
        wall_seconds=monotonic() - started,
    )
