"""Thin shim: `n2_electronic_grid` now delegates to `qscat.core.grids.electronic_grid`.

The FEM-DVR-ECS electronic-grid element layout was promoted verbatim into
`qscat.core.grids` (sub-project #A, Task 3) -- see
`docs/superpowers/specs/2026-07-27-diatomic-ve-scattering-library-design.md`.
This module is kept only so existing callers (`n2_2d_cross_section` and
sibling projects/tests) don't need to change their imports; it carries no
logic of its own.
"""

from __future__ import annotations

import qscat.core.grids
from qscat.dvr import FemDvrEcsGrid

__all__ = ["n2_electronic_grid"]


def n2_electronic_grid(
    *,
    r_max: float = 30.0,
    angle_deg: float = 35.0,
    order: int = 8,
    n_complex: int = 8,
    tail_alpha: float = 0.2,
    tail_skip: int = 2,
) -> FemDvrEcsGrid:
    """See `qscat.core.grids.electronic_grid` for the implementation."""
    return qscat.core.grids.electronic_grid(
        r_max=r_max,
        angle_deg=angle_deg,
        order=order,
        n_complex=n_complex,
        tail_alpha=tail_alpha,
        tail_skip=tail_skip,
    )
