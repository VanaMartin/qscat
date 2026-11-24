"""Common-real-region Gram matrix for two nuclear-basis projected DVR fields."""

from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss
from qscat.dvr import dvr_interpolation_matrix

_FIELD_POINT_BUDGET = 2_000_000


def _quadrature(grids, lo, hi, order):
    boundaries = [lo, hi]
    for grid in grids:
        edges = grid.spec.x_min + np.r_[0, np.cumsum([e.length for e in grid.spec.elements])]
        boundaries.extend(edges[(edges > lo) & (edges < hi)])
    edges = np.unique(boundaries)
    nodes, weights = leggauss(order)
    half = np.diff(edges) / 2
    middle = (edges[1:] + edges[:-1]) / 2
    return (middle[:, None] + half[:, None] * nodes).ravel(), (half[:, None] * weights).ravel()


def projected_real_region_gram(left, right, *, radial_stop: float, quadrature: int = 16):
    """Return a 2x2 Hermitian Gram matrix on the common UNSCALED rectangle.

    Each input is `(radial_grid, real_nuclear_grid, basis, coefficients)` with
    real nuclear DVR coefficient columns B and radial-major coefficients V.
    Fields are `(P_r @ V) @ (P_R @ B).T`, not interpolated coefficient arrays.
    This is a local shape witness, not the global bilinear ECS c-product.
    """
    states = [left, right]
    for radial, nuclear, basis, coefficients in states:
        B, V = np.asarray(basis), np.asarray(coefficients)
        if (
            B.ndim != 2
            or B.shape[0] != nuclear.n
            or B.shape[1] < 1
            or V.shape != (radial.n, B.shape[1])
            or not np.isrealobj(B)
            or not np.isfinite(B).all()
            or not np.isfinite(V).all()
            or np.max(abs(nuclear.points.imag)) > 1e-13
        ):
            raise ValueError(
                "Real nuclear basis and matching finite projected coefficients required"
            )
    radial_grids, nuclear_grids = [s[0] for s in states], [s[1] for s in states]
    rlo = max(g.spec.x_min for g in radial_grids)
    Rlo = max(g.spec.x_min for g in nuclear_grids)
    Rhi = min(g.R0 for g in nuclear_grids)
    if (
        not np.isfinite(radial_stop)
        or not rlo < radial_stop < min(g.R0 for g in radial_grids)
        or Rlo >= Rhi
        or not isinstance(quadrature, int)
        or quadrature < max(g.nq for g in radial_grids + nuclear_grids)
    ):
        raise ValueError("Common real rectangle and sufficient integer quadrature required")
    r, wr = _quadrature(radial_grids, rlo, radial_stop, quadrature)
    R, wR = _quadrature(nuclear_grids, Rlo, Rhi, quadrature)
    if len(R) > _FIELD_POINT_BUDGET:
        raise ValueError("A common nuclear row exceeds the field memory budget")
    radial_values, nuclear_values = [], []
    for radial, nuclear, B, V in states:
        radial_values.append(dvr_interpolation_matrix(radial, r) @ V)
        nuclear_values.append(dvr_interpolation_matrix(nuclear, R) @ B)
    chunk = max(1, _FIELD_POINT_BUDGET // len(R))
    gram = np.zeros((2, 2), dtype=complex)
    for start in range(0, len(r), chunk):
        end = min(start + chunk, len(r))
        fields = []
        for rv, nv in zip(radial_values, nuclear_values, strict=True):
            field = rv[start:end] @ nv.T
            field *= np.sqrt(wr[start:end])[:, None]
            field *= np.sqrt(wR)[None, :]
            fields.append(field.ravel())
        gram[0, 0] += np.vdot(fields[0], fields[0])
        gram[1, 1] += np.vdot(fields[1], fields[1])
        gram[0, 1] += np.vdot(fields[0], fields[1])
    gram[1, 0] = gram[0, 1].conjugate()
    return gram
