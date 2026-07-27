"""The model-independent VE-scattering engine.

`qscat.core` holds everything the electron-diatomic vibrational-excitation
(VE) solver stack needs that does NOT depend on which molecule is being
solved: the FEM-DVR-ECS grid builders (`grids`) and the neutral-molecule
vibrational-states solver (`vibrational`), promoted from the N2 projects
(sub-project #3) as the first slice of the generalized library -- see
`docs/superpowers/specs/2026-07-27-diatomic-ve-scattering-library-design.md`.

**Hard boundary: `qscat.core` must never import `qscat.model` (nor any
`projects.*`).** Anything molecule-specific (a potential-energy surface, a
parameter set) is passed in by the caller -- e.g. `vibrational_states` takes
`v0` as a callable rather than importing a hardcoded potential. This keeps
`core` reusable for NO, F2, and any future model the `qscat.model`
`ResonanceModel` protocol admits.

Public API:
  - `electronic_grid`, `nuclear_grid` -- FEM-DVR-ECS radial grid builders,
    parameterized (extents/orders are config, not baked in).
  - `vibrational_states` -- the `n` lowest bound eigenpairs of
    `T_nuc(mu) + diag(v0(R))` on a nuclear grid.
"""

from __future__ import annotations

from .grids import electronic_grid, nuclear_grid
from .vibrational import vibrational_states

__all__ = ["electronic_grid", "nuclear_grid", "vibrational_states"]
