"""Automatic FEM-DVR-ECS discretisation tuner.

Public API:
  - `analyze_potential` -- pure potential analysis (no models): samples a
    plain callable `V(x)` on a dense real grid and returns the local-
    wavenumber / forbidden-region-decay profile, classical turning points,
    and boundary singularities (e.g. the `-1/r` origin). This is the sole
    input the mesh/ECS generators (later tasks) consume.
  - `PotentialProfile` -- the frozen dataclass `analyze_potential` returns.
"""

from __future__ import annotations

from .analyze import PotentialProfile, analyze_potential

__all__ = ["PotentialProfile", "analyze_potential"]
