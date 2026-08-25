"""Numerical config for the H2+ dissociative-recombination (DR) validation.

The full deck transcribes eMoScat's `input/experimental/H2p.json` grid
layout verbatim (electronic real region to 1300 bohr + a 5-degree
exp-growth ECS tail; nuclear real region to 14 bohr + a 22-degree exp-growth
ECS tail; order-8 quadrature throughout), built via
`qscat.core.grids.fem_grid_exp_tail`.
`full_grid()` is Docker/MUMPS-sized (the electronic grid alone runs to 1300
bohr); `proxy_grid()` shrinks the electronic real region to ~60 bohr with a
smaller ECS tail so a laptop SuperLU solve is feasible, keeping the nuclear
grid unchanged.

`N_CHANNELS = 3` is the Rydberg exit-channel cutoff -- eMoScat's
`cross_sections.testfunctions.dissociative_recombination.channels` in
H2p.json (its `rydbergs.txt` deck listing individual Rydberg state data is
not present in this repo's `reference/` snapshot, but the JSON's `channels: 3`
pins the count independent of that file). The energy grid brackets eMoScat's
`cross_sections.range` ([0.0, 0.05] Ha) at a coarser step than its TD-FFT
`dE=1e-5` (that resolution is for the correlation-function transform, not a
practical TI energy-sweep step).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.core.grids import fem_grid_exp_tail
from qscat.dvr import TensorGrid

__all__ = ["E_HI", "E_LO", "E_STEP", "N_CHANNELS", "energy_grid", "full_grid", "proxy_grid"]

N_CHANNELS = 3

E_LO = 0.001
E_HI = 0.050
E_STEP = 0.001


def energy_grid() -> npt.NDArray[np.float64]:
    return np.round(np.arange(E_LO, E_HI + 0.5 * E_STEP, E_STEP), 6)


def full_grid() -> TensorGrid:
    """The real H2+ deck: electronic real->1300 + 5-degree ECS tail; nuclear
    real->14 + 22-degree ECS tail; order 8. Docker/MUMPS-sized.
    """
    electronic = fem_grid_exp_tail(
        [(10, 1.0), (10, 4.0), (16, 20.0), (20, 100.0), (120, 1300.0)],
        angle_deg=5.0,
        quadrature=8,
        tail_n=25,
    )
    nuclear = fem_grid_exp_tail(
        [(5, 1.0), (20, 4.0), (67, 14.0)],
        angle_deg=22.0,
        quadrature=8,
        tail_n=25,
    )
    return TensorGrid([electronic, nuclear])


def proxy_grid() -> TensorGrid:
    """The reduced laptop grid: electronic real region truncated to ~60 bohr
    (dropping the 100/1300 segments) with a smaller ECS tail; nuclear
    unchanged in extent but with a proportionally smaller tail.
    """
    electronic = fem_grid_exp_tail(
        [(10, 1.0), (10, 4.0), (16, 20.0), (10, 60.0)],
        angle_deg=5.0,
        quadrature=8,
        tail_n=8,
    )
    nuclear = fem_grid_exp_tail(
        [(5, 1.0), (20, 4.0), (40, 14.0)],
        angle_deg=22.0,
        quadrature=8,
        tail_n=8,
    )
    return TensorGrid([electronic, nuclear])
