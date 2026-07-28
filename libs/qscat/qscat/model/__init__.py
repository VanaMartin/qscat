"""The model layer: the `ResonanceModel` protocol, the shared diatomic
resonance-model form, the H2+ ionic resonance-model form, and the
per-molecule parameter registry.

Public API:
  - `ResonanceModel` -- the structural protocol `qscat.core`'s solvers depend
    on (never a concrete class below).
  - `DiatomicResonanceModel` -- the shared Morse + sigmoid + Gaussian-in-r
    form (N2/NO/F2 differ only in parameters).
  - `IonicResonanceModel` -- the H2+ Morse + sigma-capture + Coulomb-tail form.
  - `N2`, `NO`, `F2`, `H2P` -- the registry instances.

See `docs/superpowers/specs/2026-07-27-diatomic-ve-scattering-library-design.md`.
"""

from __future__ import annotations

from .diatomic import DiatomicResonanceModel, ResonanceModel
from .ionic import IonicResonanceModel
from .library import F2, H2P, N2, NO

__all__ = [
    "ResonanceModel",
    "DiatomicResonanceModel",
    "IonicResonanceModel",
    "N2",
    "NO",
    "F2",
    "H2P",
]
