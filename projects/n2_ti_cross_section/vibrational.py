"""Thin shim: `vibrational_states` now delegates to
`qscat.core.vibrational.vibrational_states`, binding N2's neutral potential.

The model-independent solver (`T_nuc(mu) + diag(v0(R))` eigendecomposition)
was promoted verbatim into `qscat.core.vibrational` (sub-project #A, Task 3)
-- see
`docs/superpowers/specs/2026-07-27-diatomic-ve-scattering-library-design.md`.
This module keeps the original 3-arg signature (`grid, mu, n`) that existing
callers use, binding `qscat.model.N2.v0` as the potential -- verified
value-identical to the old `projects.n2_resonance.potential.v0` import
(`libs/qscat/tests/test_model.py`).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import qscat.core.vibrational
from qscat.dvr import FemDvrEcsGrid
from qscat.model import N2

__all__ = ["vibrational_states"]


def vibrational_states(
    grid: FemDvrEcsGrid, mu: float, n: int
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.complex128]]:
    """See `qscat.core.vibrational.vibrational_states` for the implementation."""
    return qscat.core.vibrational.vibrational_states(grid, mu, n, N2.v0)
