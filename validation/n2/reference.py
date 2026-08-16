"""Reference values: literature resonance + Houfek golden-data anchor coordinates.

Anchor *values* are looked up from CSVE.V00.J00, never hardcoded. `ANCHOR_FACTOR` /
`ANCHOR_MARGIN_HA` are the C5 gating constants used by `cross_section.py` to compare
the TI solver's output at the anchor coordinates against those looked-up values (see
that module's docstring for the GATED-vs-DOCUMENTED-LIMITED classification).
"""

from __future__ import annotations

import numpy as np

from validation.n2 import loader

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

# C5 gating constants (see `cross_section.py` for how they classify/gate each anchor).
#
# ANCHOR_FACTOR: loose, documented cross-model bound between our 1D LCP-derived TI
# formula and Houfek's independent, explicit 2D time-independent calculation -- an
# anchor is a real PASS iff 1/ANCHOR_FACTOR <= sigma_computed/sigma_houfek <=
# ANCHOR_FACTOR. Matches `projects/n2_ti_cross_section/test_cross_section.py`'s
# `ANCHOR_FACTOR` (kept in lockstep).
ANCHOR_FACTOR = 3.0
#
# ANCHOR_MARGIN_HA: an anchor's VE channel (v'>=1) is only GATED (subject to
# ANCHOR_FACTOR) if it sits clear of its own threshold by more than this margin,
# i.e. `E_tot - eps[channel] > ANCHOR_MARGIN_HA` (E_tot = E + eps[0], v_init=0
# throughout). Below this margin the LCP's energy-INDEPENDENT local width Gamma(R)
# gives the model the wrong (non-Wigner) threshold law, so it diverges as ~1/E
# purely from the model's structure, not from a solver defect (see Task 3 / the
# physics docs). Chosen to be approximately one full vibrational quantum of the
# model's neutral N2 ladder (eps1-eps0 ~= 0.0124 Ha, see
# `projects/n2_ti_cross_section/test_vibrational.py`): comfortably excludes the
# (E=0.02 Ha, v'=1) anchor (only ~0.0076 Ha above its own threshold) while
# comfortably including every anchor in the E=0.1-0.2 Ha resonance region (excess
# >= 0.088 Ha there).
ANCHOR_MARGIN_HA = 0.0124


def anchors() -> list[tuple[float, int, float]]:
    """Resolve (energy, channel) -> reference sigma (bohr²) from the golden data."""
    d = loader.load()
    out: list[tuple[float, int, float]] = []
    for e, ch in ANCHOR_COORDS:
        i = int(np.argmin(np.abs(d.energy - e)))
        out.append((float(d.energy[i]), ch, float(d.sigma[i, ch])))
    return out
