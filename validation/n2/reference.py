"""Reference values: literature resonance + Houfek golden-data anchor coordinates.

Anchor *values* are looked up from CSVE.V00.J00, never hardcoded. RTOL applies when a
future time-independent solver's output is compared at the anchor coordinates.
"""

from __future__ import annotations

import numpy as np

import loader

# Electron–N₂ ²Π_g shape-resonance acceptance windows (eV) for the FUTURE B1 check.
# Deliberately generous *plausibility bands*, not tight literature centres: the literature
# centre is ~2.3–2.4 eV / Γ~0.4–0.5 eV (Schulz; Berman/Domcke), while the port-scout ECS
# prototype gave ~2.44 eV / ~0.46 eV. The windows below bracket BOTH so the eigensolver
# check (once ported) is not tripped by the method/grid-dependent offset. Tighten when a
# converged ECS solver exists.
LITERATURE = {"E_res_eV": (2.3, 2.5), "Gamma_eV": (0.35, 0.55)}

# Anchor coordinates: (energy_Ha, channel_index). channel 0 = elastic, j = v=0->j.
# Chosen near E=0.2 Ha (resonance region), one mid-range, one near-threshold.
ANCHOR_COORDS = [
    (0.2, 0), (0.2, 1), (0.2, 2), (0.2, 3),
    (0.1, 1), (0.02, 1),
]
RTOL = 0.05  # 5% — tune when the TI solver lands


def anchors() -> list[tuple[float, int, float]]:
    """Resolve (energy, channel) -> reference sigma (bohr²) from the golden data."""
    d = loader.load()
    out: list[tuple[float, int, float]] = []
    for e, ch in ANCHOR_COORDS:
        i = int(np.argmin(np.abs(d.energy - e)))
        out.append((float(d.energy[i]), ch, float(d.sigma[i, ch])))
    return out
