"""Thin N2-binding shim over the model-agnostic, promoted implementation in
`qscat.core.time_dependent` (sub-project #A, Task 5) -- `propagate` here
just supplies N2's `build_h2d(tgrid)` as the default Hamiltonian (the
model-agnostic core has no default; the engine gets `H` from the caller) and
delegates. See `qscat.core.time_dependent`'s module docstring for the full
propagation/sampling writeup (the two cadences, the Pade evolution operator,
norm/c-product conventions) -- unchanged by this promotion.

Kept as a module (not deleted) so existing callers/imports in this project
(and its tests) are unaffected by the move; no new physics or numerics live
here.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from qscat.core.time_dependent import PropagationResult, Snapshot
from qscat.core.time_dependent import propagate as _propagate
from qscat.dvr import TensorGrid

from projects.n2_2d_cross_section.hamiltonian2d import build_h2d

__all__ = ["PropagationResult", "Snapshot", "propagate"]


def propagate(
    tgrid: TensorGrid,
    psi0: npt.NDArray[np.complex128],
    out_channels: list[npt.NDArray[np.complex128]],
    *,
    dt: float,
    n_steps: int,
    sample_period: int = 0,
    snapshot_times: list[float] | None = None,
    keep_psi_at: list[float] | None = None,
    hamiltonian: sp.spmatrix | None = None,
    order: int = 3,
) -> PropagationResult:
    """N2-binding shim over `qscat.core.time_dependent.propagate`; see there
    for the full docstring.

    `hamiltonian` overrides the propagation Hamiltonian; when `None` (the
    default) `build_h2d(tgrid)` is used -- the N2-specific default the
    model-agnostic core dropped.
    """
    h = build_h2d(tgrid) if hamiltonian is None else hamiltonian
    return _propagate(
        tgrid,
        psi0,
        out_channels,
        dt=dt,
        n_steps=n_steps,
        sample_period=sample_period,
        snapshot_times=snapshot_times,
        keep_psi_at=keep_psi_at,
        hamiltonian=h,
        order=order,
    )
