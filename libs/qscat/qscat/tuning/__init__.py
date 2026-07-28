"""Automatic FEM-DVR-ECS discretisation tuner.

Public API:
  - `analyze_potential` -- pure potential analysis (no models): samples a
    plain callable `V(x)` on a dense real grid and returns the local-
    wavenumber / forbidden-region-decay profile, classical turning points,
    and boundary singularities (e.g. the `-1/r` origin). This is the sole
    input the mesh/ECS generators consume.
  - `PotentialProfile` -- the frozen dataclass `analyze_potential` returns.
  - `equidistribution_elements` -- adaptive real-region element lengths
    equidistributing de Broglie phase per element (capped by kappa-decay
    length in classically forbidden runs, refined near turning points /
    singularities).
  - `optimal_real_mesh` -- h/p sweep over DVR orders, returning the
    `(mesh, order)` combination with the fewest DVR points.
"""

from __future__ import annotations

from .analyze import PotentialProfile, analyze_potential
from .mesh import equidistribution_elements, optimal_real_mesh

__all__ = [
    "PotentialProfile",
    "analyze_potential",
    "equidistribution_elements",
    "optimal_real_mesh",
]
