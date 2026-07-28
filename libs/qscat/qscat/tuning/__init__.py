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
  - `max_stable_angle` -- the largest ECS rotation angle (capped at the
    double-ECS bound, ~35 deg) for which a potential `V` stays bounded on
    the rotated tail contour.
  - `tune_ecs_tail` -- exp-growth ECS-tail element lengths sized to absorb
    a given outgoing wavenumber `K` down to a target decay.
  - `ProbeResult`, `refine`, `probe_nuclear`, `probe_electronic`,
    `probe_channel_representation` -- the decoupled 1-D convergence probes:
    empirical validators that tell the tuner whether a candidate grid
    resolves the physics (nuclear vibrational levels, the electronic
    bound-state energy, and the cheap/diagnostic channel-representation
    check that catches an unresolved fast outgoing wave) -- see
    `qscat.tuning.probes`.
"""

from __future__ import annotations

from .analyze import PotentialProfile, analyze_potential
from .ecs import max_stable_angle, tune_ecs_tail
from .mesh import equidistribution_elements, optimal_real_mesh
from .probes import (
    ProbeResult,
    probe_channel_representation,
    probe_electronic,
    probe_nuclear,
    refine,
)

__all__ = [
    "PotentialProfile",
    "ProbeResult",
    "analyze_potential",
    "equidistribution_elements",
    "max_stable_angle",
    "optimal_real_mesh",
    "probe_channel_representation",
    "probe_electronic",
    "probe_nuclear",
    "refine",
    "tune_ecs_tail",
]
