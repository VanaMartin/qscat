"""Wavefunction visualisation: project a 2-D FEM-DVR-ECS state to a uniform grid
and domain-colour the complex field.

Two composable pieces, ported from eMoScat's ``EquidistantProjector2d`` +
``display_wf.py``:

  - `EquidistantProjector` -- cached sparse projection of a 2-D `TensorGrid`
    state onto an equidistant sampling grid (build once, apply per frame).
  - `complex_to_rgb` / `complex_to_hsv` -- domain colouring (phase -> hue,
    magnitude -> brightness), pure numpy.
  - `plot_wavefunction_2d` -- project + colour + draw (matplotlib, the optional
    ``qscat[plot]`` extra; imported lazily).

This is the building block for time-dependent wavefunction animations: propagate
`psi(t)`, project + colour each frame, stitch. matplotlib is NOT imported at
module load, so ``import qscat.viz`` works without the plotting extra.
"""

from __future__ import annotations

from .coloring import complex_to_hsv, complex_to_rgb, hsv_to_rgb
from .levels import energy_contour_levels
from .plot import plot_wavefunction_2d
from .projector import EquidistantProjector

__all__ = [
    "EquidistantProjector",
    "complex_to_hsv",
    "complex_to_rgb",
    "hsv_to_rgb",
    "energy_contour_levels",
    "plot_wavefunction_2d",
]
