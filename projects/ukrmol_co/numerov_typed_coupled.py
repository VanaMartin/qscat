"""Experimental real-axis coupled Numerov with retained extended arithmetic."""

from __future__ import annotations

import math
import warnings

import numpy as np
from scipy.linalg import LinAlgWarning, lu_factor, lu_solve

from projects.ukrmol_co.numerov_typed_free import small_solve


def refined_solve(a, b, stats):
    """Binary64 LU preconditioning with bounded same-dtype residual refinement."""
    if len(a) <= 3:
        return small_solve(a, b)
    dtype = a.dtype.type
    if b.dtype != a.dtype or a.dtype.kind != "f":
        raise ValueError("Compatible real arithmetic required")

    def norm(x):
        return np.max(np.sum(abs(x), axis=1))

    with warnings.catch_warnings():
        warnings.simplefilter("error", LinAlgWarning)
        factors = lu_factor(np.asarray(a, float), check_finite=True)
        x = np.asarray(lu_solve(factors, np.asarray(b, float)), dtype=dtype)
        for iteration in range(4):
            residual = b - a @ x
            denominator = norm(a) * norm(x) + norm(b)
            error = norm(residual) / denominator if denominator else dtype(0)
            if not np.isfinite(error):
                raise ValueError("Nonfinite mixed-precision residual")
            if error <= 32 * np.finfo(dtype).eps:
                stats["maximum_solve_backward_error"] = max(
                    stats["maximum_solve_backward_error"], float(error)
                )
                stats["maximum_refinement_corrections"] = max(
                    stats["maximum_refinement_corrections"], iteration
                )
                stats["large_system_solves"] += 1
                return x
            if iteration < 3:
                x += np.asarray(lu_solve(factors, np.asarray(residual, float)), dtype=dtype)
    raise np.linalg.LinAlgError("Mixed-precision refinement exhausted three corrections")


def typed_tail(thresholds, angular, energy, c2, c3, order, radius, degeneracy):
    """Outgoing series and analytic derivative, without complex128 coercion."""
    real = thresholds.dtype.type
    complex_dtype = np.result_type(thresholds.dtype, np.complex64).type
    if np.finfo(complex_dtype).nmant != np.finfo(real).nmant:
        raise ValueError("Complex arithmetic lost real significand bits")
    n = len(thresholds)
    s = np.sqrt(np.asarray(thresholds - energy, dtype=complex_dtype))
    # Principal sqrt gives positive real closed exponents; choose decaying sign.
    s = np.where(thresholds > energy, -s, s)
    gaps = thresholds[:, None] - thresholds
    grouped = abs(gaps) <= degeneracy
    coefficients = [np.eye(n, dtype=complex_dtype)]
    zero = np.zeros_like(coefficients[0])
    for p in range(1, order + 1):
        previous = coefficients[p - 1]
        before = coefficients[p - 2] if p >= 2 else zero
        earlier = coefficients[p - 3] if p >= 3 else zero
        rhs = -2 * s * (p - 1) * previous + (p - 2) * (p - 1) * before - c2 @ before - c3 @ earlier
        current = np.zeros_like(zero)
        np.divide(rhs, gaps, out=current, where=~grouped)
        degenerate = (p * (p - 1) * previous - c2 @ previous - c3 @ before) / (2 * s * p)
        current[grouped] = degenerate[grouped]
        coefficients.append(current)
    series = sum(value / radius**p for p, value in enumerate(coefficients))
    derivative = series * s - sum(
        p * value / radius ** (p + 1) for p, value in enumerate(coefficients)
    )
    opened = thresholds < energy
    momenta = s[opened].imag
    weights = np.ones(n, dtype=complex_dtype)
    pi = np.arccos(real(-1))
    weights[opened] = np.exp(
        complex_dtype(1j) * (momenta * radius - angular[opened] * pi / 2)
    ) / np.sqrt(momenta)
    values, derivatives = series * weights, derivative * weights
    if (
        max(
            np.max(abs(values[:, ~opened].imag), initial=0),
            np.max(abs(derivatives[:, ~opened].imag), initial=0),
        )
        > 128 * np.finfo(real).eps
    ):
        raise ValueError("Closed real-axis series is not real")

    def columns(x):
        return x[:, opened].imag, x[:, opened].real, x[:, ~opened].real

    return (
        columns(values),
        columns(derivatives),
        dict(
            complex_significand_bits=np.finfo(complex_dtype).nmant + 1,
            threshold_group_maximum_spread_rydberg=float(np.max(abs(gaps[grouped]))),
        ),
    )


def potential_derivatives(diagonal, c2, c3, radius):
    values = []
    for j in range(5):
        value = (-1) ** j * (
            math.factorial(j + 1) * c2 / radius ** (j + 2)
            + math.factorial(j + 2) * c3 / (2 * radius ** (j + 3))
        )
        values.append(value + diagonal if j == 0 else value)
    return values


def taylor_change(value, derivative, potential, step):
    derivatives = [value, derivative]
    for j in range(5):
        derivatives.append(
            sum(math.comb(j, k) * potential[k] @ derivatives[j - k] for k in range(j + 1))
        )
    return sum(step**j / math.factorial(j) * derivatives[j] for j in range(1, 7))


def typed_coupled_k(
    q,
    thresholds,
    angular,
    c2,
    c3,
    energy,
    radius,
    endpoint,
    step,
    order,
    degeneracy=1e-10,
    dtype=np.float64,
):
    """Unsymmetrized open K for the declared finite-order real-axis exterior."""
    dtype = np.dtype(dtype).type
    q, c2, c3 = (np.asarray(x, dtype=dtype) for x in [q, c2, c3])
    thresholds = np.asarray(thresholds, dtype=dtype)
    angular = np.asarray(angular)
    energy, radius, endpoint, h, degeneracy = map(
        dtype, [energy, radius, endpoint, step, degeneracy]
    )
    n = len(thresholds)
    if (
        thresholds.ndim != 1
        or angular.shape != thresholds.shape
        or any(x.shape != (n, n) for x in [q, c2, c3])
        or not all(np.isfinite(x).all() for x in [q, thresholds, angular, c2, c3])
        or not np.isfinite([energy, radius, endpoint, h, degeneracy]).all()
        or radius <= 0
        or endpoint <= radius
        or h <= 0
        or degeneracy < 0
        or not isinstance(order, int)
        or order < 0
        or np.any(thresholds == energy)
        or not np.any(thresholds < energy)
    ):
        raise ValueError("Finite real compatible arrays and positive uniform grid required")
    if any(np.max(abs(x - x.T)) > dtype(1e-10) for x in [q, c2, c3]):
        raise ValueError("Symmetric real boundary/potential required")
    intervals = round(float((endpoint - radius) / h))
    if intervals < 2 or abs(radius + intervals * h - endpoint) > dtype(1e-9):
        raise ValueError("Endpoint is not on the declared grid")
    stats = dict(
        maximum_solve_backward_error=0.0, maximum_refinement_corrections=0, large_system_solves=0
    )
    identity = np.eye(n, dtype=dtype)
    diagonal = np.diag(thresholds - energy)

    def w(r):
        return diagonal + c2 / r**2 + c3 / r**3

    def b(r):
        return identity - h * h * w(r) / 12

    change = taylor_change(q, identity, potential_derivatives(diagonal, c2, c3, radius), h)
    b1 = b(radius + h)
    difference = -h * h * (w(radius + h) - w(radius)) / 12
    delta = refined_solve((b1 @ (q + change)).T, (difference @ q + b1 @ change).T, stats).T
    constant = not np.any(c2) and not np.any(c3)
    cached = h * h * refined_solve(b(radius), w(radius), stats) if constant else None
    for i in range(1, intervals):
        if constant:
            increment = cached
        else:
            current_w = w(radius + i * h)
            increment = h * h * refined_solve(identity - h * h * current_w / 12, current_w, stats)
        x = delta + increment
        delta = refined_solve(identity + x, x, stats)
    potential = potential_derivatives(diagonal, c2, c3, endpoint)
    zero = np.zeros_like(identity)
    previous_b = b(endpoint - h)
    difference = -h * h * (w(endpoint - h) - w(endpoint)) / 12
    a = previous_b @ taylor_change(identity, zero, potential, -h) + difference + delta @ b(endpoint)
    c = previous_b @ taylor_change(zero, identity, potential, -h)
    columns, derivatives, tail_stats = typed_tail(
        thresholds, angular, energy, c2, c3, order, endpoint, degeneracy
    )
    regular, irregular, closed = [
        a @ value + c @ derivative for value, derivative in zip(columns, derivatives, strict=True)
    ]
    matching = np.column_stack([irregular, closed])
    scales = np.max(abs(matching), axis=0)
    if np.any(scales == 0):
        raise ValueError("Vanishing matching column")
    solution = refined_solve(matching / scales, -regular, stats) / scales[:, None]
    k = solution[: regular.shape[1]]
    if not np.isfinite(k).all():
        raise ValueError("Nonfinite coupled K")
    stats.update(
        tail_stats,
        radial_intervals=intervals,
        significand_bits=np.finfo(dtype).nmant + 1,
        dtype=np.dtype(dtype).name,
        unsymmetrized_k_transpose_error=float(np.max(abs(k - k.T))),
        matching_residual=float(np.max(abs(matching @ solution + regular))),
    )
    return k, stats
