"""Reference values: literature resonance + Houfek golden-data anchor coordinates.

Anchor *values* are looked up from CSVE.V00.J00, never hardcoded. RTOL applies when a
future time-independent solver's output is compared at the anchor coordinates.
"""

from __future__ import annotations

import numpy as np

import loader

# Literature electron–N₂ ²Π_g shape resonance (eV) — Schulz; Berman/Domcke.
LITERATURE = {"E_res_eV": (2.3, 2.4), "Gamma_eV": (0.35, 0.55)}

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
