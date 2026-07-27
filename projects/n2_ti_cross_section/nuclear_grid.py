"""Thin shim: `n2_nuclear_grid` now delegates to `qscat.core.grids.nuclear_grid`.

The FEM-DVR-ECS nuclear-grid element layout was promoted verbatim into
`qscat.core.grids` (sub-project #A, Task 3) -- see
`docs/superpowers/specs/2026-07-27-diatomic-ve-scattering-library-design.md`.
This module is kept only so existing callers (`n2_ti_cross_section` and
sibling projects/tests) don't need to change their imports; it carries no
logic of its own.
"""

from __future__ import annotations

import qscat.core.grids
from qscat.dvr import FemDvrEcsGrid

__all__ = ["n2_nuclear_grid"]


def n2_nuclear_grid(
    *,
    angle_deg: float = 35.0,
    r_max: float = 40.0,
    n_complex: int = 10,
    quadrature: int = 14,
) -> FemDvrEcsGrid:
    """See `qscat.core.grids.nuclear_grid` for the implementation."""
    return qscat.core.grids.nuclear_grid(
        angle_deg=angle_deg,
        r_max=r_max,
        n_complex=n_complex,
        quadrature=quadrature,
    )
