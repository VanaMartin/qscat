"""Rectangular loader-aware Gaussian cross metrics on one common physical ball.

No orbital renormalization, native export, state identity or missing-vector
recovery is performed. See docs/physics/co-loader-aware-self-metric.md.
"""

from __future__ import annotations

import numpy as np

from projects.ukrmol_co.native_orbital_metric import (
    _basis,
    _scaled_tail_roots,
    radial_overlap,
    shell_overlap,
)


def all_space_cross_metric(left, right, *, method="moments", extra=0, max_work=100_000_000):
    """Independent rectangular shell-pair assembly, without mirrored reversal."""
    left, a = _basis(left)
    right, b = _basis(right)
    if (
        method not in ("moments", "hermite")
        or type(extra) is not int
        or not 0 <= extra <= 2
        or type(max_work) is not int
        or not 1 <= max_work <= 100_000_000
    ):
        raise ValueError("Known method, bounded Hermite increment and work cap required")
    work = sum(
        ((x.angular + y.angular + 2) // 2 + extra) ** 3
        * (2 * x.angular + 2 * y.angular + 2)
        * len(x.exponents)
        * len(y.exponents)
        for x in left
        for y in right
    )
    if work > max_work:
        raise ValueError("Rectangular Gaussian-product point-function work cap exceeded")
    result = np.empty((a[-1], b[-1]), complex)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[a[i] : a[i + 1], b[j] : b[j + 1]] = shell_overlap(
                x, y, method=method, extra=extra
            )
    return result


def cross_ball_witness(left, right, all_space, radius, left_scales, right_scales):
    """Restored hybrid cross Gram plus an entrywise uniform-ball tail bound.

    Only origin-centered continuum pairs use explicit ball integrals. Every
    target-containing exterior term has a contraction-aware Cauchy bound. The
    two generator-scale vectors are fixed independently of the physical radius.
    """
    left, a = _basis(left)
    right, b = _basis(right)
    matrix = np.asarray(all_space, complex)
    da, db = np.asarray(left_scales), np.asarray(right_scales)
    if (
        matrix.shape != (a[-1], b[-1])
        or da.shape != (a[-1],)
        or db.shape != (b[-1],)
        or not np.isfinite(matrix).all()
        or not np.isfinite(radius)
        or radius <= 0
        or any(not np.isrealobj(d) or not np.isfinite(d).all() or np.any(d <= 0) for d in (da, db))
    ):
        raise ValueError(
            "Matching finite cross Gram, positive AO scales and physical radius required"
        )
    hybrid = matrix.copy()
    ca = np.concatenate([np.full(2 * s.angular + 1, s.continuum, bool) for s in left])
    cb = np.concatenate([np.full(2 * s.angular + 1, s.continuum, bool) for s in right])
    for i, x in enumerate(left):
        if not x.continuum:
            continue
        for j, y in enumerate(right):
            if not y.continuum:
                continue
            block = np.zeros((2 * x.angular + 1, 2 * y.angular + 1), complex)
            if x.angular == y.angular:
                block = np.eye(2 * x.angular + 1) * radial_overlap(x, y, radius)
            hybrid[a[i] : a[i + 1], b[j] : b[j + 1]] = block
    error = np.outer(_scaled_tail_roots(left, radius, da), _scaled_tail_roots(right, radius, db))
    error[np.logical_and.outer(ca, cb)] = 0
    return da[:, None] * hybrid * db[None, :], error
