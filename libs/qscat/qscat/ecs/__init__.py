"""Exterior Complex Scaling (ECS) utilities.

The ECS contour maps the real radial coordinate onto a path that runs
straight out to a pivot `R0` and then bends into the complex plane at a
fixed angle `theta`. Rotating the outgoing (continuum) coordinate this way
turns divergent scattering states into decaying ones and exposes resonance
poles, while leaving bound-state energies unchanged (see
`docs/physics/femdvr-ecs.md` and `.superpowers/sdd/femdvr-ecs-extraction.md`).

This module is the single source of that transform: `qscat.dvr.grid` uses
`ecs_map` to place its complex-tail quadrature points.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = ["ecs_map"]


def ecs_map(
    x: npt.ArrayLike, R0: float, theta_deg: float
) -> npt.NDArray[np.complexfloating]:
    """Exterior-complex-scaling coordinate map.

    ``z(x) = x`` for ``x <= R0``, and ``z(x) = R0 + (x - R0) * exp(i*theta)``
    for ``x > R0``. `x` may be a scalar or array; `theta_deg` is in degrees.
    """
    xa = np.asarray(x, dtype=np.float64)
    eit = np.exp(1j * np.deg2rad(theta_deg))
    return np.where(xa <= R0, xa.astype(np.complex128), R0 + (xa - R0) * eit)
