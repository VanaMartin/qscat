"""Experimental cancellation-reduced ratio algebra; original Numerov oracle retained."""

from __future__ import annotations

import math

import numpy as np
from scipy.linalg import solve

from projects.ukrmol_co.exterior_numerov import (
    _potential_derivatives,
    _taylor_value,
    match_ratio,
    tail_coefficients,
    tail_columns,
)


def taylor_increment(value, derivative, potential, step):
    """Sixth-order value change, omitting the order-zero value before summation."""
    derivatives = [value, derivative]
    for n in range(5):
        derivatives.append(
            sum(math.comb(n, j) * potential[j] @ derivatives[n - j] for j in range(n + 1))
        )
    return sum(step**n / math.factorial(n) * derivatives[n] for n in range(1, 7))


def increment_numerov_k(
    q, thresholds, angular, c2, c3, energy, radius, endpoint, step, order, degeneracy=1e-10
):
    """Same real-axis Numerov problem using D=I-R without identity subtraction.

    See docs/physics/co-numerov-increment-control.md for ordered matrix algebra.
    This experimental function does not symmetrize K or modify the saved method.
    """
    q, c2, c3 = (np.asarray(x, float) for x in (q, c2, c3))
    thresholds, angular = np.asarray(thresholds, float), np.asarray(angular)
    count = len(thresholds)
    if (
        thresholds.ndim != 1
        or angular.shape != thresholds.shape
        or any(x.shape != (count, count) for x in [q, c2, c3])
        or not all(np.isfinite(x).all() for x in [q, thresholds, angular, c2, c3])
        or not all(np.isfinite(x) for x in [energy, radius, endpoint, step, degeneracy])
        or radius <= 0
        or endpoint <= radius
        or step <= 0
        or degeneracy < 0
        or not isinstance(order, int)
        or order < 0
        or not np.any(thresholds < energy)
    ):
        raise ValueError("Finite compatible real arrays and a positive radial grid required")
    if any(not np.allclose(x, x.T, atol=1e-10, rtol=0) for x in [q, c2, c3]):
        raise ValueError("Real symmetric boundary and potential coefficients required")
    intervals = round((endpoint - radius) / step)
    if intervals < 2 or abs(radius + intervals * step - endpoint) > 1e-9:
        raise ValueError("The endpoint must lie on the declared uniform grid")
    identity = np.eye(count)
    diagonal = np.diag(thresholds - energy)

    def w(r):
        return diagonal + c2 / r**2 + c3 / r**3

    def b(r):
        return identity - step**2 * w(r) / 12

    potential = _potential_derivatives(diagonal, c2, c3, radius, 5)
    change = taylor_increment(q, identity, potential, step)
    b1 = b(radius + step)
    b_difference = -(step**2) * (w(radius + step) - w(radius)) / 12
    # Right division: (z1-z0) z1^-1. Inverting Q is neither needed nor valid.
    delta = solve((b1 @ (q + change)).T, (b_difference @ q + b1 @ change).T).T
    for i in range(1, intervals):
        current_w = w(radius + i * step)
        current_b = identity - step**2 * current_w / 12
        x = delta + step**2 * solve(current_b, current_w, assume_a="sym")
        delta = solve(identity + x, x)
    coefficients, s, spread = tail_coefficients(thresholds, energy, c2, c3, order, degeneracy)
    final = tail_columns(thresholds, angular, energy, coefficients, s, endpoint, endpoint)
    derivative = tail_columns(
        thresholds, angular, energy, coefficients, s, endpoint, endpoint, derivative=True
    )
    potential = _potential_derivatives(diagonal, c2, c3, endpoint, 5)
    zero = np.zeros_like(identity)
    b_previous = b(endpoint - step)
    b_difference = -(step**2) * (w(endpoint - step) - w(endpoint)) / 12
    value_coefficient = b_previous @ taylor_increment(identity, zero, potential, -step)
    value_coefficient += b_difference + delta @ b(endpoint)
    derivative_coefficient = b_previous @ _taylor_value(zero, identity, potential, -step)
    k, diagnostics = match_ratio(value_coefficient, derivative_coefficient, final, derivative)
    if not np.isfinite(k).all():
        raise ValueError("Nonfinite reconstructed K")
    diagnostics.update(
        radial_intervals=intervals,
        threshold_group_maximum_spread_rydberg=spread,
        unsymmetrized_k_transpose_error=float(np.max(abs(k - k.T))),
        formulation="identity-defect increment",
    )
    return k, diagnostics
