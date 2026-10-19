"""Independent scalar and direct mixed-channel recurrence checks of the Jost oracle."""

import mpmath as mp
import numpy as np
import pytest

from projects.ukrmol_co.numerov_discrete_free import discrete_free_k
from projects.ukrmol_co.numerov_increment_study import exact_discrete_scalar


@pytest.mark.parametrize(
    "endpoint,step", [(400.0, 0.003125), (200.0, 0.0015625), (800.0, 0.0015625)]
)
def test_matrix_oracle_agrees_with_independent_scalar_discrete_transfer(endpoint, step):
    actual = discrete_free_k([[0.4]], [0.0], 0.327075, 18.0, endpoint, step, 90)
    expected = exact_discrete_scalar(0.4, 0.327075, 18.0, endpoint, step, 90)
    with mp.workdps(90):
        assert abs(actual[0, 0] - expected) < mp.mpf("1e-65")


@pytest.mark.parametrize("thresholds", [[0.0, 0.2], [0.0, 1.2], [0.0, 0.0, 1.2]])
def test_backward_oracle_matches_direct_forward_matrix_recurrence(thresholds):
    n = len(thresholds)
    q = np.full((n, n), 0.08)
    np.fill_diagonal(q, [0.4, 0.3, 0.2][:n])
    with mp.workdps(90):
        energy, h, b = mp.mpf(0.7), mp.mpf(0.125), mp.mpf(3.0)
        w = mp.diag([mp.mpf(t) - energy for t in thresholds])
        p = mp.eye(n)
        d = mp.zeros(n)
        for j in range(1, 4):
            p += h ** (2 * j) * w**j / mp.factorial(2 * j)
        for j in range(3):
            d += h ** (2 * j + 1) * w**j / mp.factorial(2 * j + 1)
        a = 12 * (mp.eye(n) - h * h * w / 12) ** -1 - 10 * mp.eye(n)
        previous = mp.matrix(q.tolist())
        current = p * previous + d
        for _ in range(1, 8):
            previous, current = current, a * current - previous
        coef = p - previous * current**-1
        opens = [i for i, t in enumerate(thresholds) if t < energy]
        closed = [i for i, t in enumerate(thresholds) if t > energy]
        f, g = mp.zeros(n, len(opens)), mp.zeros(n, len(opens))
        fp, gp = mp.zeros(n, len(opens)), mp.zeros(n, len(opens))
        decay, dp = mp.zeros(n, len(closed)), mp.zeros(n, len(closed))
        for j, i in enumerate(opens):
            k = mp.sqrt(energy - thresholds[i])
            f[i, j], g[i, j] = mp.sin(k * b) / mp.sqrt(k), mp.cos(k * b) / mp.sqrt(k)
            fp[i, j], gp[i, j] = mp.sqrt(k) * mp.cos(k * b), -mp.sqrt(k) * mp.sin(k * b)
        for j, i in enumerate(closed):
            decay[i, j], dp[i, j] = 1, -mp.sqrt(thresholds[i] - energy)
        rhs = -(coef * f - d * fp)
        open_lhs = coef * g - d * gp
        closed_lhs = coef * decay - d * dp
        lhs = mp.matrix(n)
        for i in range(n):
            for j in range(len(opens)):
                lhs[i, j] = open_lhs[i, j]
            for j in range(len(closed)):
                lhs[i, len(opens) + j] = closed_lhs[i, j]
        expected = mp.matrix(len(opens))
        for j in range(len(opens)):
            solution = mp.lu_solve(lhs, rhs[:, j])
            for i in range(len(opens)):
                expected[i, j] = solution[i]
        actual = discrete_free_k(q, thresholds, 0.7, 2.0, 3.0, 0.125, 90)
        assert max(abs(x - y) for x, y in zip(actual, expected, strict=True)) < mp.mpf("1e-65")
