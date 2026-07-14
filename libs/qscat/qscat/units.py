"""Atomic-unit conversions. All physics in qModeling uses atomic units."""

from __future__ import annotations

import numpy as np

# CODATA 2018
HARTREE_TO_EV: float = 27.211386245988
EV_TO_HARTREE: float = 1.0 / HARTREE_TO_EV


def hartree_to_ev(x: float | np.ndarray) -> float | np.ndarray:
    """Convert energy in Hartree to electron-volts."""
    return x * HARTREE_TO_EV


def ev_to_hartree(x: float | np.ndarray) -> float | np.ndarray:
    """Convert energy in electron-volts to Hartree."""
    return x * EV_TO_HARTREE
