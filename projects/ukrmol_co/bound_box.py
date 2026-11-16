"""Finite-radius closed-channel maps, not infinite-domain bound certificates."""

from time import monotonic

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve

from projects.ukrmol_co.bound_matching import (
    _DEFAULT_CAPS,
    _negative_energy,
    _real,
    _symmetric,
)


def _exterior_inputs(thresholds, c2, c3, radius, caps):
    thresholds = _real(thresholds, "Thresholds").copy()
    c2, c3 = _real(c2, "C2").copy(), _real(c3, "C3").copy()
    n = len(thresholds) if thresholds.ndim == 1 else 0
    if not 1 <= n <= caps.max_channels:
        raise ValueError("Channel dimension exceeds the declared cap or is empty")
    if np.any(thresholds < 0) or any(m.shape != (n, n) for m in [c2, c3]):
        raise ValueError("Nonnegative thresholds and compatible square multipoles required")
    for matrix in [c2, c3]:
        if not np.array_equal(matrix, matrix.T):
            raise ValueError("Finite-box multipoles must be exactly real symmetric")
    radius = _real(radius, "Radius")
    if radius.ndim != 0 or radius <= 0:
        raise ValueError("A positive scalar radius is required")
    return thresholds, c2, c3, float(radius)


def nonnegative_exterior_witness(thresholds_rydberg, c2, c3, radius_bohr, *, caps=_DEFAULT_CAPS):
    """Sufficient endpoint-matrix test for V(r)>=0 at every r>=a.

    In exact arithmetic, r**2 V(r) >= a**2 diag(thresholds)+C2+C3/r.
    The latter is a convex combination of its endpoints at 1/r=0 and 1/a.
    Nonnegative Gershgorin bounds at both endpoints therefore suffice. A
    floating-point safety allowance is applied, but is not interval certification.
    Rejection is inconclusive, not evidence for a negative exterior or no binding.
    """
    started = monotonic()
    thresholds, c2, c3, radius = _exterior_inputs(thresholds_rydberg, c2, c3, radius_bohr, caps)
    n = len(thresholds)
    base = radius**2 * np.diag(thresholds) + c2
    lower_bounds, allowances = [], []
    for matrix in [base, base + c3 / radius]:
        if monotonic() - started >= caps.timeout_seconds:
            raise TimeoutError("Exterior-witness time cap reached")
        if not np.isfinite(matrix).all():
            raise ValueError("Exterior-witness matrix is nonfinite")
        off_diagonal = matrix.copy()
        np.fill_diagonal(off_diagonal, 0)
        lower = np.diag(matrix) - np.sum(abs(off_diagonal), axis=1)
        scale = np.sum(radius**2 * np.diag(thresholds) + abs(c2) + abs(c3) / radius, axis=1)
        allowance = 64 * (n + 1) * np.finfo(float).eps * scale
        if np.any(lower < allowance):
            raise ValueError("Nonnegative exterior not established by sufficient endpoint test")
        lower_bounds.append(float(np.min(lower)))
        allowances.append(float(np.max(allowance)))
    if monotonic() - started >= caps.timeout_seconds:
        raise TimeoutError("Exterior-witness time cap reached")
    return dict(
        nonnegative_established=True,
        radius_bohr=radius,
        endpoint_gershgorin_lower_bounds=lower_bounds,
        floating_safety_allowances=allowances,
        roundoff_certified=False,
        scope="Sufficient finite-input algebraic witness, not an electronic or ODE error bound",
        wall_seconds=monotonic() - started,
    )


def finite_box_dtn(
    incident_hartree,
    thresholds_rydberg,
    c2,
    c3,
    radius_bohr,
    *,
    endpoint_bohr,
    boundary,
    rtol=1e-11,
    caps=_DEFAULT_CAPS,
):
    """Return Y(a) with u(R)=0 or u'(R)=0; never switch a failed chart.

    Dirichlet starts the inverse chart Z=U(U')^-1 at zero and integrates
    Z'=I-Z V Z. Neumann starts Y=U'U^-1 at zero and integrates Y'=V-Y Y.
    The final Dirichlet Z is inverted only after condition/symmetry checks.
    Both use V=diag(thresholds-2*epsilon)+C2/r**2+C3/r**3 in native units.
    """
    started = monotonic()
    energy = _negative_energy(incident_hartree)
    thresholds, second, third, radius = _exterior_inputs(
        thresholds_rydberg, c2, c3, radius_bohr, caps
    )
    if boundary not in ["dirichlet", "neumann"]:
        raise ValueError("An explicit dirichlet or neumann boundary is required")
    settings = _real([endpoint_bohr, rtol], "Finite-box settings")
    endpoint, rtol = settings
    if not radius < endpoint or not 0 < rtol < 1:
        raise ValueError("An exterior endpoint beyond a and 0 < rtol < 1 are required")
    n = len(thresholds)
    diagonal = np.diag(thresholds - 2 * energy)
    identity = np.eye(n)
    calls = 0

    def check_cap():
        if monotonic() - started >= caps.timeout_seconds:
            raise TimeoutError("Finite-box map time cap reached")

    def rhs(r, flat):
        nonlocal calls
        calls += 1
        check_cap()
        if calls > caps.max_rhs:
            raise TimeoutError("Finite-box map RHS cap reached")
        chart = flat.reshape(n, n)
        if not np.isfinite(chart).all() or np.max(abs(chart)) > 1e8:
            raise ValueError("Unresolved finite-box chart or divergent integration")
        potential = diagonal + second / r**2 + third / r**3
        derivative = (
            identity - chart @ potential @ chart
            if boundary == "dirichlet"
            else potential - chart @ chart
        )
        return derivative.ravel()

    result = solve_ivp(
        rhs,
        (endpoint, radius),
        np.zeros(n * n),
        method="DOP853",
        rtol=rtol,
        atol=rtol * 0.01,
        t_eval=[radius],
    )
    check_cap()
    if not result.success or not np.isfinite(result.y).all():
        raise ValueError("Finite-box integration failed; no alternate chart attempted")
    if any(
        not np.array_equal(original, copied)
        for original, copied in [
            (thresholds_rydberg, thresholds),
            (c2, second),
            (c3, third),
        ]
    ):
        raise ValueError("Exterior inputs changed during finite-box propagation")
    chart = result.y[:, -1].reshape(n, n)
    _symmetric(chart, "Finite-box chart")
    condition = float(np.linalg.cond(chart)) if boundary == "dirichlet" else None
    if condition is not None and (not np.isfinite(condition) or condition > 1e10):
        raise ValueError("Singular or unresolved finite-box inverse chart")
    y = solve(chart, identity) if boundary == "dirichlet" else chart
    if not np.isfinite(y).all() or np.max(abs(y)) > 1e8:
        raise ValueError("Unresolved finite-box DtN chart")
    _symmetric(y, "Finite-box DtN")
    check_cap()
    return y, dict(
        finite_box=True,
        boundary=boundary,
        decaying=False,
        endpoint_bohr=float(endpoint),
        radius_bohr=radius,
        rhs_evaluations=calls,
        transpose_defect=float(np.max(abs(y - y.T))),
        inverse_chart_condition=condition,
        total_error_bound=None,
        wall_seconds=monotonic() - started,
    )
