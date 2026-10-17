"""Real-axis matrix Numerov matching with a coupled inverse-power exterior tail."""

from __future__ import annotations

import math

import numpy as np
from scipy.linalg import solve


def tail_coefficients(thresholds, energy, c2, c3, order, degeneracy=1e-10):
    """Outgoing/decaying series coefficients, with threshold grouping explicit."""
    thresholds = np.asarray(thresholds, float)
    momenta = np.sqrt((energy - thresholds).astype(complex))
    if np.any(momenta == 0):
        raise ValueError("An exact channel threshold is outside this method")
    s = 1j * momenta
    gaps = thresholds[:, None] - thresholds
    grouped = abs(gaps) <= degeneracy
    coefficients = [np.eye(len(thresholds), dtype=complex)]
    zero = np.zeros_like(coefficients[0])
    for p in range(1, order + 1):
        previous = coefficients[p - 1]
        before = coefficients[p - 2] if p >= 2 else zero
        earlier = coefficients[p - 3] if p >= 3 else zero
        rhs = -2 * s * (p - 1) * previous
        rhs += (p - 2) * (p - 1) * before - c2 @ before - c3 @ earlier
        current = np.zeros_like(zero)
        np.divide(rhs, gaps, out=current, where=~grouped)
        degenerate = (p * (p - 1) * previous - c2 @ previous - c3 @ before) / (2 * s * p)
        current[grouped] = degenerate[grouped]
        coefficients.append(current)
    return np.asarray(coefficients), s, float(np.max(abs(gaps[grouped])))


def tail_columns(thresholds, angular, energy, coefficients, s, radius, anchor, *, derivative=False):
    """Return native-sign f, g and closed-decaying columns at a finite radius.

    Closed columns omit their common exponential at anchor; all columns still
    carry their coupled algebraically forced components in the other channels.
    """
    thresholds, angular = np.asarray(thresholds), np.asarray(angular)
    opened = thresholds < energy
    powers = radius ** -np.arange(len(coefficients), dtype=float)
    series = np.einsum("p,pij->ij", powers, coefficients)
    if derivative:
        series = series * s - np.einsum(
            "p,pij->ij", np.arange(len(coefficients)) * powers / radius, coefficients
        )
    weights = np.empty(len(thresholds), complex)
    k = s[opened].imag
    weights[opened] = np.exp(1j * (k * radius - angular[opened] * np.pi / 2)) / np.sqrt(k)
    weights[~opened] = np.exp(s[~opened] * (radius - anchor))
    columns = series * weights
    closed = columns[:, ~opened]
    if np.max(abs(closed.imag), initial=0) > 1e-12:
        raise ValueError("Closed columns are not real on this real-energy sheet")
    return columns[:, opened].imag, columns[:, opened].real, closed.real


def match_ratio(value_coefficient, derivative_coefficient, columns, derivatives):
    """Match unsymmetrized open K and arbitrary closed column amplitudes."""
    residuals = [
        value_coefficient @ value + derivative_coefficient @ derivative
        for value, derivative in zip(columns, derivatives, strict=True)
    ]
    regular, irregular, closed = residuals
    matching = np.column_stack([irregular, closed])
    scales = np.max(abs(matching), axis=0)
    if np.any(scales == 0):
        raise ValueError("A matching column vanishes")
    equilibrated = matching / scales
    solution = solve(equilibrated, -regular) / scales[:, None]
    residual = float(np.max(abs(matching @ solution + regular)))
    return solution[: regular.shape[1]], {
        "matching_condition_inf": float(np.linalg.cond(matching, p=np.inf)),
        "equilibrated_matching_condition_inf": float(np.linalg.cond(equilibrated, p=np.inf)),
        "matching_residual": residual,
    }


def _potential_derivatives(diagonal, c2, c3, radius, count):
    derivatives = []
    for n in range(count):
        value = (-1) ** n * (
            math.factorial(n + 1) * c2 / radius ** (n + 2)
            + math.factorial(n + 2) * c3 / (2 * radius ** (n + 3))
        )
        if n == 0:
            value = value + diagonal
        derivatives.append(value)
    return derivatives


def _taylor_value(value, derivative, potential, step):
    derivatives = [value, derivative]
    for n in range(5):
        derivatives.append(
            sum(math.comb(n, j) * potential[j] @ derivatives[n - j] for j in range(n + 1))
        )
    return sum(step**n / math.factorial(n) * value for n, value in enumerate(derivatives))


def numerov_k(
    q, thresholds, angular, c2, c3, energy, radius, endpoint, step, order, degeneracy=1e-10
):
    """Propagate all channels and return K for the declared finite-order tail.

    Energies and the radial equation use native Rydberg coefficients; q is the
    length-valued zero-Bloch surface map. Thresholds need not be sorted, but
    at least one channel must be open. No K symmetrization or native fitting.
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

    def b(r):
        return identity - step**2 * (diagonal + c2 / r**2 + c3 / r**3) / 12

    potential = _potential_derivatives(diagonal, c2, c3, radius, 5)
    next_value = _taylor_value(q, identity, potential, step)
    # Right division constructs z_0 z_1^-1 without inverting Q.
    ratio = solve((b(radius + step) @ next_value).T, (b(radius) @ q).T).T
    for i in range(1, intervals):
        current = b(radius + i * step)
        recurrence = 12 * solve(current, identity, assume_a="sym") - 10 * identity - ratio
        ratio = solve(recurrence, identity)
    coefficients, s, spread = tail_coefficients(thresholds, energy, c2, c3, order, degeneracy)
    final = tail_columns(thresholds, angular, energy, coefficients, s, endpoint, endpoint)
    final_derivative = tail_columns(
        thresholds, angular, energy, coefficients, s, endpoint, endpoint, derivative=True
    )
    potential = _potential_derivatives(diagonal, c2, c3, endpoint, 5)
    zero = np.zeros_like(identity)
    value_coefficient = b(endpoint - step) @ _taylor_value(identity, zero, potential, -step)
    value_coefficient -= ratio @ b(endpoint)
    derivative_coefficient = b(endpoint - step) @ _taylor_value(zero, identity, potential, -step)
    k, diagnostics = match_ratio(value_coefficient, derivative_coefficient, final, final_derivative)
    if not np.isfinite(k).all():
        raise ValueError("Nonfinite reconstructed K")
    diagnostics.update(
        {
            "radial_intervals": intervals,
            "threshold_group_maximum_spread_rydberg": spread,
            "unsymmetrized_k_transpose_error": float(np.max(abs(k - k.T))),
        }
    )
    return k, diagnostics
