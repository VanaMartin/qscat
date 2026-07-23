"""Parametrized electronic FEM-DVR-ECS grid for the exact 2-D N2 model.

Layout follows eMoScat's `input/experimental/N2-model.json` `grids.electronic`:
a finely-resolved region near the origin where the interaction
`-lambda(R) exp(-alpha_c r^2)` lives, coarsening outward, then an ECS tail of
exponentially growing elements at a single angle.

EVERY parameter is exposed because sub-project #6's convergence study (Task 4)
varies all of them -- eMoScat asserted 35 degrees and a 98-bohr box without
ever documenting the study that justified them, so we redo it.
"""

from __future__ import annotations

import numpy as np
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec

__all__ = ["n2_electronic_grid"]

# (segment_end_bohr, element_length_bohr); the final segment runs to r_max.
_INNER_SEGMENTS: tuple[tuple[float, float], ...] = ((1.0, 0.2), (5.0, 1.0), (7.0, 2.0), (10.0, 3.0))
_OUTER_LENGTH = 4.0


def _ecs_tail(base: float, n: int, *, skip: int, alpha: float) -> list[float]:
    """eMoScat `uniform_increment`/`exp`: `skip` elements at `base`, then growing."""
    return [base if i < skip else base * float(np.exp(alpha * (i - skip + 1))) for i in range(n)]


def n2_electronic_grid(
    *,
    r_max: float = 30.0,
    angle_deg: float = 35.0,
    order: int = 8,
    n_complex: int = 8,
    tail_alpha: float = 0.2,
    tail_skip: int = 2,
) -> FemDvrEcsGrid:
    """Electronic radial grid: real region [0, r_max] + an ECS tail at `angle_deg`.

    The ECS pivot is `R0 == r_max` by construction.
    """
    if r_max <= _INNER_SEGMENTS[-1][0]:
        raise ValueError(f"r_max must exceed {_INNER_SEGMENTS[-1][0]} bohr, got {r_max}")

    elements: list[ElementSpec] = []
    start = 0.0
    for end, length in _INNER_SEGMENTS:
        k = round((end - start) / length)
        elements += [ElementSpec((end - start) / k) for _ in range(k)]
        start = end

    k_out = max(1, round((r_max - start) / _OUTER_LENGTH))
    elements += [ElementSpec((r_max - start) / k_out) for _ in range(k_out)]

    base = (r_max - start) / k_out
    elements += [
        ElementSpec(h, angle_deg)
        for h in _ecs_tail(base, n_complex, skip=tail_skip, alpha=tail_alpha)
    ]
    return FemDvrEcsGrid(GridSpec(quadrature=order, elements=elements, x_min=0.0))
