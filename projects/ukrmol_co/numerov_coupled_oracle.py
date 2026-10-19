"""Independent small-model fundamental-value and tail-series reference calculations."""

from __future__ import annotations

import math

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp


def mp_tail(thresholds, angular, energy, c2, c3, order, endpoint):
    """Entrywise inverse-power coefficient equations in the current mp context."""
    n = len(thresholds)
    exponents = [
        mp.j * mp.sqrt(energy - t) if t < energy else -mp.sqrt(t - energy) for t in thresholds
    ]
    terms = [mp.eye(n)]
    zero = mp.zeros(n)
    for power in range(1, order + 1):
        out = mp.zeros(n)
        previous = terms[power - 1]
        before = terms[power - 2] if power >= 2 else zero
        earlier = terms[power - 3] if power >= 3 else zero
        for i in range(n):
            for j in range(n):
                if thresholds[i] == thresholds[j]:
                    out[i, j] = (
                        power * (power - 1) * previous[i, j]
                        - sum(c2[i, k] * previous[k, j] + c3[i, k] * before[k, j] for k in range(n))
                    ) / (2 * exponents[j] * power)
                else:
                    out[i, j] = (
                        -2 * exponents[j] * (power - 1) * previous[i, j]
                        + (power - 2) * (power - 1) * before[i, j]
                        - sum(c2[i, k] * before[k, j] + c3[i, k] * earlier[k, j] for k in range(n))
                    ) / (thresholds[i] - thresholds[j])
        terms.append(out)
    values, derivatives = mp.zeros(n), mp.zeros(n)
    opened, closed = [], []
    for j, t in enumerate(thresholds):
        if t < energy:
            opened.append(j)
            k = mp.sqrt(energy - t)
            phase = mp.exp(mp.j * (k * endpoint - angular[j] * mp.pi / 2)) / mp.sqrt(k)
        else:
            closed.append(j)
            phase = 1
        for i in range(n):
            values[i, j] = phase * sum(terms[p][i, j] / endpoint**p for p in range(order + 1))
            derivatives[i, j] = phase * sum(
                (exponents[j] - mp.mpf(p) / endpoint) * terms[p][i, j] / endpoint**p
                for p in range(order + 1)
            )

    def select(matrix, indices, part):
        return (
            mp.matrix([[part(matrix[i, j]) for j in indices] for i in range(n)])
            if indices
            else mp.zeros(n, 0)
        )

    def columns(matrix):
        return (
            select(matrix, opened, mp.im),
            select(matrix, opened, mp.re),
            select(matrix, closed, mp.re),
        )

    return columns(values), columns(derivatives)


def _inputs(q, thresholds, angular, c2, c3, energy, radius, endpoint, step):
    def matrix(x):
        return mp.matrix(np.asarray(x).tolist())

    return (
        matrix(q),
        [mp.mpf(float(t)) for t in thresholds],
        list(angular),
        matrix(c2),
        matrix(c3),
        *map(lambda x: mp.mpf(float(x)), [energy, radius, endpoint, step]),
    )


def _match(values, derivative, tail, tail_prime):
    logarithmic = derivative * values**-1
    residuals = [d - logarithmic * v for v, d in zip(tail, tail_prime, strict=True)]
    regular, irregular, closed = residuals
    combined = mp.matrix(
        [
            [
                irregular[i, j] if j < irregular.cols else closed[i, j - irregular.cols]
                for j in range(values.rows)
            ]
            for i in range(values.rows)
        ]
    )
    answer = combined**-1 * (-regular)
    return answer[: regular.cols, :]


def direct_discrete_k(
    q, thresholds, angular, c2, c3, energy, radius, endpoint, step, order, digits=70
):
    """Propagate fundamental values, without ratio or identity-defect recurrences."""
    intervals = round((endpoint - radius) / step)
    with mp.workdps(digits):
        q, thresholds, angular, c2, c3, energy, radius, endpoint, h = _inputs(
            q, thresholds, angular, c2, c3, energy, radius, endpoint, step
        )
        n = len(thresholds)
        identity = mp.eye(n)

        def w(r):
            return mp.diag([t - energy for t in thresholds]) + c2 / r**2 + c3 / r**3

        def taylor(value, derivative, at, offset):
            potential = [
                mp.matrix(
                    [
                        [mp.diff(lambda r, i=i, j=j: w(r)[i, j], at, k) for j in range(n)]
                        for i in range(n)
                    ]
                )
                for k in range(5)
            ]
            derivatives = [value, derivative]
            for k in range(5):
                derivatives.append(
                    sum(
                        (math.comb(k, j) * potential[j] * derivatives[k - j] for j in range(k + 1)),
                        mp.zeros(n),
                    )
                )
            return sum(
                (offset**k / mp.factorial(k) * v for k, v in enumerate(derivatives)), mp.zeros(n)
            )

        previous, current = q, taylor(q, identity, radius, h)

        def b(r):
            return identity - h * h * w(r) / 12

        for i in range(1, intervals):
            at = radius + i * h
            after = b(at + h) ** -1 * (
                (2 * identity + 5 * h * h * w(at) / 6) * current - b(at - h) * previous
            )
            previous, current = current, after
        p = taylor(identity, mp.zeros(n), endpoint, -h)
        d = taylor(mp.zeros(n), identity, endpoint, -h)
        derivative = d**-1 * (previous - p * current)
        columns, derivatives = mp_tail(thresholds, angular, energy, c2, c3, order, endpoint)
        return +_match(current, derivative, columns, derivatives)


def continuous_k(q, thresholds, angular, c2, c3, energy, radius, endpoint, order, rtol, atol):
    """DOP853 fundamental columns, matched with independently evaluated mp tails."""
    n = len(thresholds)

    def rhs(r, state):
        value, derivative = state.reshape(2, n, n)
        w = np.diag(thresholds - energy) + c2 / r**2 + c3 / r**3
        return np.stack([derivative, w @ value]).ravel()

    solution = solve_ivp(
        rhs,
        (radius, endpoint),
        np.stack([q, np.eye(n)]).ravel(),
        method="DOP853",
        rtol=rtol,
        atol=atol,
    )
    if not solution.success:
        raise ValueError("Continuous fundamental integration failed")
    value, derivative = solution.y[:, -1].reshape(2, n, n)
    with mp.workdps(70):
        _, ts, ang, c2m, c3m, e, _, end, _ = _inputs(
            q, thresholds, angular, c2, c3, energy, radius, endpoint, 1.0
        )
        columns, derivatives = mp_tail(ts, ang, e, c2m, c3m, order, end)
        k = _match(mp.matrix(value.tolist()), mp.matrix(derivative.tolist()), columns, derivatives)
        return np.array(k.tolist(), float), dict(function_evaluations=solution.nfev)
