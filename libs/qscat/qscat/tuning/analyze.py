"""Pure potential analysis for the FEM-DVR-ECS discretisation tuner.

`analyze_potential` samples a plain potential callable `V(x)` on a dense
real grid and produces a `PotentialProfile`: the local-wavenumber /
forbidden-region-decay profile, the classical turning points (where
`e_max - Re(V(x))` changes sign), and any boundary singularities (e.g. the
`-1/r` origin of a Coulomb tail). It takes no models -- just a callable and
floats -- so the mesh/ECS generators built on top of it can consume the
profile without depending on `qscat.model` or `qscat.core`.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]

# Growth test for boundary singularities: an edge sample is flagged as
# singular when |V| there exceeds both a floor (guards against an all-zero
# or tiny-scale potential) and a multiple of the "bulk" (interior) scale of
# |V| over the rest of the domain.
_SINGULARITY_GROWTH_FACTOR = 50.0
_SINGULARITY_ABS_FLOOR = 10.0


@dataclass(frozen=True)
class PotentialProfile:
    """The potential-analysis output consumed by the mesh/ECS generators."""

    x: FloatArray
    V: FloatArray
    k: FloatArray
    kappa: FloatArray
    turning_points: FloatArray
    singularities: FloatArray


def analyze_potential(
    V: Callable[[NDArray[np.float64]], NDArray[np.complexfloating | np.floating]],
    x_min: float,
    x_max: float,
    m: float,
    e_max: float,
    *,
    n_sample: int = 4000,
) -> PotentialProfile:
    """Analyze `V` on `[x_min, x_max]`, returning its `PotentialProfile`.

    `V` may return complex values (potentials do, on ECS tails); the
    analysis works on `Re(V)` and is only ever sampled at real `x`.
    """
    singularities: list[float] = []

    x_lo = float(x_min)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        v_lo = float(np.real(np.asarray(V(np.array([x_lo]))))[0])
    if not np.isfinite(v_lo):
        # Non-finite right at the boundary (e.g. -1/r at r=0): record the
        # singularity there and start sampling a hair above it.
        singularities.append(x_lo)
        span = float(x_max) - x_lo
        x_lo = x_lo + max(span * 1e-6, 1e-12)

    x = np.linspace(x_lo, float(x_max), n_sample)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        v_sample = np.asarray(V(x))
    v_x = np.real(v_sample).astype(np.float64)

    # Growth test: does |V| blow up near the low edge relative to the bulk
    # scale elsewhere in the domain?
    bulk = np.abs(v_x[v_x.size // 4 :])
    bulk_finite = bulk[np.isfinite(bulk)]
    bulk_scale = float(np.median(bulk_finite)) if bulk_finite.size else 0.0
    threshold = max(_SINGULARITY_GROWTH_FACTOR * bulk_scale, _SINGULARITY_ABS_FLOOR)
    edge_n = max(n_sample // 20, 1)
    edge_abs = np.abs(v_x[:edge_n])
    if edge_abs.size and np.isfinite(edge_abs).any() and np.nanmax(edge_abs) > threshold:
        singularities.append(float(x[0]))

    k = np.sqrt(2.0 * m * np.maximum(e_max - v_x, 0.0))
    kappa = np.sqrt(2.0 * m * np.maximum(v_x - e_max, 0.0))

    turning_points = _turning_points(x, e_max - v_x)

    return PotentialProfile(
        x=x,
        V=v_x,
        k=k,
        kappa=kappa,
        turning_points=turning_points,
        singularities=np.array(sorted(set(singularities)), dtype=np.float64),
    )


def _turning_points(x: FloatArray, f: FloatArray) -> FloatArray:
    """Locate sign changes of `f = e_max - V` and interpolate the crossing."""
    sign = np.sign(f)
    crossings = np.nonzero(np.diff(sign))[0]
    points = np.empty(crossings.size, dtype=np.float64)
    for out_i, i in enumerate(crossings):
        x0, x1 = x[i], x[i + 1]
        f0, f1 = f[i], f[i + 1]
        if f1 == f0:
            points[out_i] = x0
        else:
            t = -f0 / (f1 - f0)
            points[out_i] = x0 + t * (x1 - x0)
    return points
