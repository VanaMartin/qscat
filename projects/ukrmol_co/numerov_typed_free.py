"""Small dtype-controlled free Numerov diagnostic; no LAPACK precision downcast."""

from __future__ import annotations

import numpy as np


def small_solve(a, b, stats=None):
    """Partial-pivot elimination in the input dtype, for one to three channels."""
    a, b = np.array(a, copy=True), np.array(b, copy=True)
    n = len(a)
    if (
        a.shape != (n, n)
        or b.ndim != 2
        or b.shape[0] != n
        or not 1 <= n <= 3
        or a.dtype != b.dtype
        or a.dtype.kind != "f"
    ):
        raise ValueError("Compatible one-to-three-channel matrices of one dtype required")
    if n == 1:
        if a[0, 0] == 0:
            raise np.linalg.LinAlgError("Singular one-channel solve")
        return b / a[0, 0]
    for column in range(n):
        pivot = column + int(np.argmax(abs(a[column:, column])))
        scale = np.max(abs(a[column:, column:]))
        if a[pivot, column] == 0:
            raise np.linalg.LinAlgError("Singular small matrix solve")
        if stats is not None:
            stats["minimum_relative_pivot"] = min(
                stats["minimum_relative_pivot"], float(abs(a[pivot, column]) / scale)
            )
        if pivot != column:
            a[[column, pivot]] = a[[pivot, column]]
            b[[column, pivot]] = b[[pivot, column]]
        factors = a[column + 1 :, column] / a[column, column]
        a[column + 1 :, column + 1 :] -= factors[:, None] * a[column, column + 1 :]
        b[column + 1 :] -= factors[:, None] * b[column]
    answer = np.zeros_like(b)
    for row in range(n - 1, -1, -1):
        answer[row] = (b[row] - a[row, row + 1 :] @ answer[row + 1 :]) / a[row, row]
    return answer


def typed_free_k(q, thresholds, energy, radius, endpoint, step, dtype=np.float64):
    """Constant diagonal exterior, l=0, up to three channels; all arithmetic typed."""
    dtype = np.dtype(dtype).type
    q = np.asarray(q, dtype=dtype)
    thresholds = np.asarray(thresholds, dtype=dtype)
    n = len(thresholds)
    energy, radius, endpoint, h = map(dtype, [energy, radius, endpoint, step])
    intervals = round(float((endpoint - radius) / h))
    if (
        q.shape != (n, n)
        or not 1 <= n <= 3
        or h <= 0
        or intervals < 2
        or abs(radius + intervals * h - endpoint) > dtype(1e-9)
        or not np.isfinite(q).all()
        or not np.isfinite(thresholds).all()
        or not np.isfinite([energy, radius, endpoint, h]).all()
        or np.any(thresholds == energy)
        or not np.any(thresholds < energy)
    ):
        raise ValueError("Finite compatible free boundary and uniform grid required")
    identity = np.eye(n, dtype=dtype)
    w = thresholds - energy
    b = 1 - h * h * w / 12
    p = h * h * w / 2 + h**4 * w * w / 24 + h**6 * w**3 / 720
    d = h + h**3 * w / 6 + h**5 * w * w / 120
    change = p[:, None] * q + np.diag(d)
    stats = {"minimum_relative_pivot": 1.0}
    delta = small_solve((b[:, None] * (q + change)).T, (b[:, None] * change).T, stats).T
    recurrence = np.diag(h * h * w / b)
    for _ in range(1, intervals):
        x = delta + recurrence
        delta = small_solve(identity + x, x, stats)
    value_coefficient = np.diag(b * p) + delta * b
    derivative_coefficient = np.diag(-b * d)
    opened = np.flatnonzero(thresholds < energy)
    closed = np.flatnonzero(thresholds > energy)
    regular = np.zeros((n, len(opened)), dtype=dtype)
    irregular = np.zeros_like(regular)
    regular_prime = np.zeros_like(regular)
    irregular_prime = np.zeros_like(regular)
    decay = np.zeros((n, len(closed)), dtype=dtype)
    decay_prime = np.zeros_like(decay)
    for j, i in enumerate(opened):
        k = np.sqrt(energy - thresholds[i])
        angle = k * endpoint
        regular[i, j], irregular[i, j] = np.sin(angle) / np.sqrt(k), np.cos(angle) / np.sqrt(k)
        regular_prime[i, j], irregular_prime[i, j] = (
            np.sqrt(k) * np.cos(angle),
            -np.sqrt(k) * np.sin(angle),
        )
    for j, i in enumerate(closed):
        decay[i, j], decay_prime[i, j] = 1, -np.sqrt(thresholds[i] - energy)
    lhs = np.column_stack(
        [
            value_coefficient @ irregular + derivative_coefficient @ irregular_prime,
            value_coefficient @ decay + derivative_coefficient @ decay_prime,
        ]
    )
    rhs = -(value_coefficient @ regular + derivative_coefficient @ regular_prime)
    scales = np.max(abs(lhs), axis=0)
    if np.any(scales == 0):
        raise np.linalg.LinAlgError("Vanishing free matching column")
    solution = small_solve(lhs / scales, rhs, stats) / scales[:, None]
    k = solution[: len(opened)]
    if not np.isfinite(k).all():
        raise ValueError("Nonfinite typed K matrix")
    stats.update(
        radial_intervals=intervals,
        significand_bits=np.finfo(dtype).nmant + 1,
        dtype=np.dtype(dtype).name,
        unsymmetrized_k_transpose_error=float(np.max(abs(k - k.T))),
    )
    return k, stats
