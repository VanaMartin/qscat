"""The model layer: the `ResonanceModel` protocol, the shared diatomic
resonance-model form, the H2+ ionic resonance-model form, and the
per-molecule parameter registry.

Public API:
  - `ResonanceModel` -- the structural protocol `qscat.core`'s solvers depend
    on (never a concrete class below).
  - `DiatomicResonanceModel` -- the shared Morse + sigmoid + Gaussian-in-r
    form (N2/NO/F2 differ only in parameters).
  - `IonicResonanceModel` -- the H2+ Morse + sigma-capture + Coulomb-tail form.
  - `FlexibleDiatomicModel` (+ `SmoothR`/`TailR`) -- the potential factory's
    fitted form: EMO neutral + Gaussian well with `lam(R)`, `alpha(R)` and an
    optional shell; embeds the `DiatomicResonanceModel`s exactly
    (`from_diatomic`) -- see docs/physics/potential-factory.md.
  - `N2`, `NO`, `F2`, `H2P`, `O2` -- the registry instances (`O2` is the
    first FITTED model, not a published parameter set), plus `O2_SO12` /
    `O2_SO32`, its two spin-orbit components (statistical factor 1/3 each).

Adding a molecule means adding a registry entry (parameters) plus its
validation -- never solver code, which lives model-free in `qscat.core`. See
`docs/physics/qscat-core-scattering.md`.
"""

from __future__ import annotations

from .diatomic import DiatomicResonanceModel, ResonanceModel
from .flexible import FlexibleDiatomicModel, SmoothR, TailR, from_diatomic
from .ionic import IonicResonanceModel
from .library import F2, H2P, N2, NO, O2, O2_SO12, O2_SO32

__all__ = [
    "F2",
    "H2P",
    "N2",
    "NO",
    "O2",
    "O2_SO12",
    "O2_SO32",
    "DiatomicResonanceModel",
    "FlexibleDiatomicModel",
    "IonicResonanceModel",
    "ResonanceModel",
    "SmoothR",
    "TailR",
    "from_diatomic",
]
