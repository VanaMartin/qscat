"""High-precision backward discrete Jost matching for a constant diagonal exterior."""

from __future__ import annotations

import mpmath as mp
import numpy as np


def discrete_free_k(q, thresholds, energy, radius, endpoint, step, digits=80):
    """Solve the same sixth-order-endpoint Numerov difference problem in closed form.

    Only a free diagonal exterior and zero angular momenta are represented.
    Closed Jost columns are independently scaled; no ratio recurrence is used.
    """
    intervals = round((endpoint - radius) / step)
    if intervals < 2 or abs(radius + intervals * step - endpoint) > 1e-9:
        raise ValueError("The endpoint must lie on a positive uniform grid")
    with mp.workdps(digits):
        q = mp.matrix(np.asarray(q).tolist())
        thresholds = [mp.mpf(float(t)) for t in thresholds]
        energy, endpoint, step = map(mp.mpf, [energy, endpoint, step])
        n = len(thresholds)
        opened = [i for i, t in enumerate(thresholds) if t < energy]
        closed = [i for i, t in enumerate(thresholds) if t > energy]
        if len(opened) + len(closed) != n or not opened:
            raise ValueError("Exact thresholds and closed-only problems are unsupported")
        columns = [mp.zeros(n, len(opened)), mp.zeros(n, len(opened)), mp.zeros(n, len(closed))]
        for i, t in enumerate(thresholds):
            w = t - energy
            p = mp.fsum(step ** (2 * j) * w**j / mp.factorial(2 * j) for j in range(4))
            d = mp.fsum(step ** (2 * j + 1) * w**j / mp.factorial(2 * j + 1) for j in range(3))
            a = 2 + step**2 * w / (1 - step**2 * w / 12)
            if t < energy:
                theta, k = mp.acos(a / 2), mp.sqrt(energy - t)
                values = {
                    j: mp.sin(j * theta) for j in [1, intervals, intervals - 1, intervals - 2]
                }
                j = opened.index(i)
                endpoint_pairs = [
                    (0, j, mp.sin(k * endpoint) / mp.sqrt(k), mp.sqrt(k) * mp.cos(k * endpoint)),
                    (1, j, mp.cos(k * endpoint) / mp.sqrt(k), -mp.sqrt(k) * mp.sin(k * endpoint)),
                ]
            else:
                gamma, kap = mp.acosh(a / 2), mp.sqrt(t - energy)
                values = {
                    j: mp.sinh(j * gamma) for j in [1, intervals, intervals - 1, intervals - 2]
                }
                scale = mp.exp(-intervals * gamma)
                endpoint_pairs = [(2, closed.index(i), scale, -kap * scale)]
            for kind, j, value, derivative in endpoint_pairs:
                before = p * value - d * derivative
                v0 = (values[intervals] * before - values[intervals - 1] * value) / values[1]
                v1 = (values[intervals - 1] * before - values[intervals - 2] * value) / values[1]
                # Assemble boundary residuals by columns, avoiding inversion of Q.
                initial_derivative = (v1 - p * v0) / d
                for row in range(n):
                    columns[kind][row, j] += (v0 if row == i else 0) - q[
                        row, i
                    ] * initial_derivative
        regular, irregular, decaying = columns
        lhs = mp.matrix(n)
        for i in range(n):
            for j in range(len(opened)):
                lhs[i, j] = irregular[i, j]
            for j in range(len(closed)):
                lhs[i, len(opened) + j] = decaying[i, j]
        result = mp.matrix(len(opened))
        for j in range(len(opened)):
            solution = mp.lu_solve(lhs, -regular[:, j])
            for i in range(len(opened)):
                result[i, j] = solution[i]
        return result
