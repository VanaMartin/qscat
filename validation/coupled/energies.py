"""The energy mesh for the coupled VE sweep.

Two scales have to be resolved and they differ by an order of magnitude.
Threshold cusps are non-analytic AT each channel opening, so they need points
bracketing the threshold tightly rather than a fine mesh everywhere. The
overlapping-resonance interference is spread across the whole range instead,
on the scale of a resonance width -- at the campaign's own (s, kappa) =
(0.3, 0.5) point, over R where each curve is genuinely resonant (Re E_res >
0), the N_l = 4 width runs 33.7-134.7 mHa and the fixed-l (N_l = 1) width
58.1-328 mHa, while the 24 vibrational levels this module computes are
spaced 5.1-9.1 mHa apart (widest near threshold, narrowing with v) -- so the
resonances genuinely overlap.

Hence a background grid plus clusters. Measured, this gives 1008 energies
against 2961 for a uniform mesh at the same 0.05 mHa resolution, which would
cost three times as much to resolve twenty places.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.core.vibrational import vibrational_states
from qscat.model import NO

from validation.diatomic.config import CONFIGS

__all__ = [
    "BACKGROUND",
    "CLUSTER_HALF",
    "CLUSTER_STEP",
    "DEDUP_TOL",
    "E_HI",
    "E_LO",
    "sweep_energies",
    "vibrational_thresholds",
]

# Hartree. The window covers the near-threshold region, the 20 vibrational
# thresholds that fall inside it, and every measured cross-section peak
# (E_max 0.051-0.084 Ha across v'=0..4, both models) -- NOT either model's
# full E_res(R) curve, which runs far wider (0.008-0.338 Ha for N_l=4,
# 0.002-0.342 Ha for N_l=1, at (s,kappa)=(0.3,0.5)) because E_res(R) keeps
# rising at the small-R points nobody scatters at. NO's DA threshold is
# +0.172 Ha, above all of it.
E_LO = 0.002
E_HI = 0.150
# 0.25 mHa background: 130+ points across the narrowest width used in the
# campaign (33.7 mHa, N_l = 4 at (s, kappa) = (0.3, 0.5)).
BACKGROUND = 2.5e-4
# 21 points at 0.05 mHa spanning +-0.5 mHa around each threshold.
CLUSTER_HALF = 5.0e-4
CLUSTER_STEP = 5.0e-5
# Two energies closer than this are one energy. Dedup by TOLERANCE, not exact
# equality: rounding alone leaves pairs a ten-thousandth of a mHa apart, and
# each one is a wasted ~15 s solve on the production deck.
DEDUP_TOL = 1e-6


def vibrational_thresholds(n_levels: int = 24) -> npt.NDArray[np.float64]:
    """`eps_v - eps_0` for each neutral level: where channel `v'` opens for VE
    out of `v = 0`."""
    nuclear = CONFIGS["NO"].da_grid().grids[1]
    eps, _chi = vibrational_states(nuclear, NO.mu, n_levels, NO.v0)
    return np.asarray(eps - eps[0], dtype=np.float64)


def sweep_energies() -> npt.NDArray[np.float64]:
    """The mesh: background grid plus a dense cluster at every threshold."""
    parts = [np.arange(E_LO, E_HI + 0.5 * BACKGROUND, BACKGROUND)]
    for t in vibrational_thresholds():
        if E_LO < t < E_HI:
            parts.append(
                np.arange(t - CLUSTER_HALF, t + CLUSTER_HALF + 0.5 * CLUSTER_STEP, CLUSTER_STEP)
            )
    E = np.sort(np.concatenate(parts))
    E = E[(E >= E_LO) & (E <= E_HI)]
    keep = np.concatenate(([True], np.diff(E) >= DEDUP_TOL))
    return np.asarray(E[keep], dtype=np.float64)
