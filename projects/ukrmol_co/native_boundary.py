"""Zero-Bloch native surface sum: stored amplitudes already include 1/sqrt(2)."""

from __future__ import annotations

import numpy as np


def surface_matrix(poles: np.ndarray, amplitudes: np.ndarray, total_hartree: complex) -> np.ndarray:
    """Return length-valued Q mapping u' to u from an uncorrected pole subset.

    Amplitudes are pole-major, channel-minor as decoded by boundary_spectrum.
    This has no missing-spectrum correction and uses no complex conjugation.
    """
    poles = np.asarray(poles, dtype=float)
    amplitudes = np.asarray(amplitudes, dtype=float)
    if (
        poles.ndim != 1
        or amplitudes.ndim != 2
        or len(poles) != len(amplitudes)
        or len(poles) == 0
        or amplitudes.shape[1] == 0
        or not np.isfinite(poles).all()
        or not np.isfinite(amplitudes).all()
        or not np.isfinite(total_hartree)
    ):
        raise ValueError(
            "Finite pole energies, pole-major amplitudes and total Hartree energy required"
        )
    denominators = poles - total_hartree
    if np.any(denominators == 0):
        raise ValueError("Energy lies on an exported inner-region pole")
    return amplitudes.T @ (amplitudes / denominators[:, None])
