"""Energy-unit conversions.

**Every quantity in qscat is in atomic units** (hbar = m_e = e = 1): energies
in Hartree, lengths in bohr, times in hbar/Hartree. Nothing in the library
converts units internally -- these helpers exist only for presenting results.
"""

from __future__ import annotations

import numpy as np

__all__ = ["HARTREE_TO_EV", "EV_TO_HARTREE", "hartree_to_ev", "ev_to_hartree"]

# CODATA 2018
HARTREE_TO_EV: float = 27.211386245988
EV_TO_HARTREE: float = 1.0 / HARTREE_TO_EV


def hartree_to_ev(x: float | np.ndarray) -> float | np.ndarray:
    """Convert energy in Hartree to electron-volts."""
    return x * HARTREE_TO_EV


def ev_to_hartree(x: float | np.ndarray) -> float | np.ndarray:
    """Convert energy in electron-volts to Hartree."""
    return x * EV_TO_HARTREE
