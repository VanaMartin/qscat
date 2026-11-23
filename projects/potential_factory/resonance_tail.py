"""Project-level complex-resonance adapter to the real-channel ECS tail tuner."""

from __future__ import annotations

import cmath
import math

from qscat.tuning import tune_ecs_tail

__all__ = ["tune_resonance_tail"]


def tune_resonance_tail(
    momentum: complex,
    R0: float,
    *,
    angle: float,
    order: int,
    decay_target: float,
    base: float | None = None,
    tail_alpha: float = 0.2,
    tail_skip: int = 2,
) -> list[float]:
    """Tail lengths giving global plane-wave damping at an exposed resonance.

    Potential-continuation stability and actual pole convergence are separate
    checks. This does not change the existing real-channel tuner or defaults.
    """
    if (
        not all(math.isfinite(x) for x in (momentum.real, momentum.imag, R0, angle, decay_target))
        or momentum.real <= 0
        or momentum.imag > 0
        or R0 < 0
        or not 0 < angle < 90
        or not 0 < decay_target < 1
        or not math.isfinite(tail_alpha)
        or tail_alpha < 0
        or (base is not None and (not math.isfinite(base) or base <= 0))
    ):
        raise ValueError(
            "Finite outgoing momentum, nonnegative pivot and valid angle/budget required"
        )
    theta = math.radians(angle)
    rate = (momentum * cmath.exp(1j * theta)).imag
    if rate <= 0:
        raise ValueError("ECS angle does not expose this resonance momentum")
    pivot_adjusted_target = decay_target * math.exp(momentum.imag * R0)
    if pivot_adjusted_target <= 0:
        raise ValueError("Global attenuation budget is not representable at this pivot")
    return tune_ecs_tail(
        rate / math.sin(theta),
        R0,
        angle=angle,
        order=order,
        decay_target=pivot_adjusted_target,
        base=base,
        tail_alpha=tail_alpha,
        tail_skip=tail_skip,
    )
