"""ECS-tail tuner: max stable rotation angle + exp-growth absorbing tail.

`max_stable_angle` scans the ECS rotation angle `theta` from small values up
to `angle_cap` (the double-ECS bound, ~35 deg) and evaluates the potential
`V` on the rotated contour `z(x) = R0 + (x - R0) * exp(i*theta)` (`qscat.ecs.
ecs_map`) out to `R0 + tail_extent`. Some interactions (e.g. a Gaussian
`exp(-alpha r**2)`) analytically continue to a GROWING function past a
critical angle (45 deg for a pure Gaussian, since `Re(z**2)` itself turns
over and heads to `-infinity` once `cos(2*theta) < 0`) -- rotating past that
angle would blow up the ECS tail rather than absorb it. `max_stable_angle`
rejects any `theta` where `|V|` grows anywhere along the tail relative to its
running minimum, and returns the largest angle that does not, capped at
`angle_cap` regardless (a bare `-1/r`, or anything else non-diverging, is
limited only by that cap).

`tune_ecs_tail` sizes an exponentially-growing ECS tail (the same
`base * exp(alpha * (i - skip + 1))` growth as `qscat.core.grids._ecs_tail`)
long enough to absorb the fastest outgoing wave `exp(-K * (x - R0) *
sin(theta))` down to `decay_target`. The required tail length is `L =
-ln(decay_target) / (K * sin(theta))`; `tail_skip` flat elements at `base`
followed by exponentially growing ones (rate `tail_alpha`) are appended
until their cumulative length reaches `L`.
"""

from __future__ import annotations

import math
from collections.abc import Callable

import numpy as np
import numpy.typing as npt

from qscat.ecs import ecs_map

__all__ = ["max_stable_angle", "tune_ecs_tail"]

# Number of points sampled along the tail `[R0, R0 + tail_extent]` when
# probing a single angle for growth -- independent of `n_probe` (which
# controls the angle scan resolution instead).
_N_TAIL_SAMPLES = 400

# Relative tolerance for the "|V| grows along the tail" divergence check --
# guards against float round-off flagging a flat/monotonically-decaying
# profile as diverging.
_GROWTH_RELTOL = 1e-6
# Absolute floor (Hartree) under which |V| is "zero" for the growth check: a
# potential that has already decayed to round-off on the tail (O2's EMO
# neutral, beta(inf) = 3.9, is at 1e-17 by the pivot) otherwise shows
# noise-level "growth" against a running minimum of 1e-18 and reads as
# diverging at EVERY angle, returning 0 degrees (measured).
_GROWTH_ABS_FLOOR = 1e-12


def max_stable_angle(
    V: Callable[[npt.NDArray[np.complexfloating]], npt.ArrayLike],
    R0: float,
    tail_extent: float,
    *,
    angle_cap: float = 35.0,
    n_probe: int = 40,
) -> float:
    """Largest ECS rotation angle (deg) for which `V` stays bounded on the tail.

    Scans `n_probe` angles from `angle_cap / n_probe` up to `angle_cap`
    (ascending); at each, evaluates `V` on the rotated contour over
    `_N_TAIL_SAMPLES` points spanning `[R0, R0 + tail_extent]`. The first
    angle at which `|V|` grows anywhere along the tail (relative to its
    running minimum, beyond `_GROWTH_RELTOL`; a minimum below
    `_GROWTH_ABS_FLOOR` counts as that floor, so a potential already at
    round-off cannot "grow") stops the scan; the last angle that did NOT
    diverge is returned. Never exceeds `angle_cap`.
    """
    if n_probe < 1:
        raise ValueError(f"n_probe must be >= 1, got {n_probe}")

    x = R0 + np.linspace(0.0, tail_extent, _N_TAIL_SAMPLES)
    angles = np.linspace(angle_cap / n_probe, angle_cap, n_probe)

    best = 0.0
    for theta in angles:
        z = ecs_map(x, R0, float(theta))
        mag = np.abs(np.asarray(V(z), dtype=np.complex128))
        running_min = np.maximum(np.minimum.accumulate(mag), _GROWTH_ABS_FLOOR)
        diverges = bool(np.any(mag > running_min * (1.0 + _GROWTH_RELTOL)))
        if diverges:
            break
        best = float(theta)
    return min(best, angle_cap)


def tune_ecs_tail(
    K: float,
    R0: float,
    *,
    angle: float,
    order: int,
    tail_alpha: float = 0.2,
    tail_skip: int = 2,
    decay_target: float = 1e-12,
    base: float | None = None,
) -> list[float]:
    """Exp-growth ECS-tail element lengths absorbing wavenumber `K` at `angle`.

    Returns element LENGTHS (not a grid) spanning `[R0, R0 + L]` where `L =
    -ln(decay_target) / (K * sin(angle))` is the distance over which the
    fastest outgoing wave `exp(-K * (x - R0) * sin(angle))` decays to
    `decay_target`. The first `tail_skip` elements are flat at `base`; the
    rest grow as `base * exp(tail_alpha * (i - tail_skip + 1))` (the same
    growth law as `qscat.core.grids._ecs_tail`), appended until the running
    sum reaches `L`.

    `order` is the DVR quadrature order the caller will attach to each tail
    element when building the actual grid; it does not affect the LENGTHS
    returned here, but is validated (`order >= 2`, matching `GridSpec`) since
    `tune_ecs_tail` is meant to be a pure function of `(K, R0, angle, order,
    growth)` per the ECS-tail tuner design.

    `base`, if not given, defaults to `L / 8` -- a small fraction of the
    target tail span, chosen so `tail_skip` flat elements plus a handful of
    growing ones reach `L` without either an excessively coarse first
    element or an excessively long list.
    """
    if order < 2:
        raise ValueError(f"order must be >= 2, got {order}")
    if tail_skip < 0:
        raise ValueError(f"tail_skip must be >= 0, got {tail_skip}")

    sin_theta = math.sin(math.radians(angle))
    if sin_theta <= 0.0:
        raise ValueError(f"angle must be in (0, 180) so sin(angle) > 0, got {angle}")
    if K <= 0.0:
        raise ValueError(f"K must be > 0, got {K}")

    tail_length = -math.log(decay_target) / (K * sin_theta)
    if base is None:
        base = tail_length / 8.0
    if base <= 0.0:
        raise ValueError(f"base must be > 0, got {base}")

    lengths: list[float] = []
    total = 0.0
    i = 0
    while total < tail_length:
        h = base if i < tail_skip else base * math.exp(tail_alpha * (i - tail_skip + 1))
        lengths.append(h)
        total += h
        i += 1
    return lengths
