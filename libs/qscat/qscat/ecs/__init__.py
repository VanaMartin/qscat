"""Exterior Complex Scaling (ECS) utilities.

The ECS contour maps the real radial coordinate onto a path that runs
straight out to a pivot `R0` and then bends into the complex plane at a
fixed angle `theta`. Rotating the outgoing (continuum) coordinate this way
turns divergent scattering states into decaying ones and exposes resonance
poles, while leaving bound-state energies unchanged (see
`docs/physics/femdvr-ecs.md` and `.superpowers/sdd/femdvr-ecs-extraction.md`).

This module is the single source of that transform: `qscat.dvr.grid` uses
`ecs_map` to place its complex-tail quadrature points.

It also carries `find_resonance_pole`, the general two-spectrum resonance-pole
matcher promoted from the N2 resonance project (`projects/n2_resonance/pole.py`,
sub-project #2) -- see `docs/physics/n2-resonance.md`.
"""

from __future__ import annotations

from .map import ecs_map
from .pole import find_resonance_pole

__all__ = ["ecs_map", "find_resonance_pole"]
