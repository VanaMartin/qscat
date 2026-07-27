"""The model layer: the `ResonanceModel` protocol, the shared diatomic
resonance-model form, and the per-molecule parameter registry.

Public API:
  - `ResonanceModel` -- the structural protocol `qscat.core`'s solvers depend
    on (never the concrete class below).
  - `DiatomicResonanceModel` -- the shared Morse + sigmoid + Gaussian-in-r
    form (N2/NO/F2 differ only in parameters).
  - `N2`, `NO`, `F2` -- the registry instances.

See `docs/superpowers/specs/2026-07-27-diatomic-ve-scattering-library-design.md`.
"""

from __future__ import annotations

from .diatomic import DiatomicResonanceModel, ResonanceModel
from .library import F2, N2, NO

__all__ = ["ResonanceModel", "DiatomicResonanceModel", "N2", "NO", "F2"]
