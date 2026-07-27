"""The exact 2-D electron-N2 potential surface and Hamiltonian.

    H_2D = -(1/2) d^2/dr^2 - (1/2 mu) d^2/dR^2
           + v0(R) + l(l+1)/(2 r^2) - lambda(R) exp(-alpha_c r^2)

verbatim eMoScat `Neutral2dPotential` (`source/Model2d/Potentials2d.cpp:18`).

Thin N2-binding shim over `qscat.model.N2`, the single source of truth for
the N2 potential-surface/Hamiltonian model (sub-project #A, Task 6) -- every
function here just delegates to the corresponding `N2` method/attribute.
Kept as a module (not deleted) so existing callers/imports in this project
(and its tests) are unaffected; no potential physics lives here anymore.

TWO potentials live here and confusing them is a physics error, not a
convention choice:

- `potential_2d` is the FULL surface, which goes into `H_2D`.
- `interaction_2d` is `V_int = -lambda(R) exp(-alpha_c r^2)` ALONE -- the only
  perturbation relative to the entrance channel. `v0(R)` (the neutral
  molecule's own potential) and `l(l+1)/2r^2` (the centrifugal barrier) are
  CHANNEL potentials: they survive as `r -> infinity`, and the asymptotic
  channel function is an eigenfunction of the Hamiltonian containing them.
  Sweeping them into the driving term would produce a plausible-looking but
  wrong T-matrix.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from qscat.dvr import TensorGrid
from qscat.model import N2

__all__ = [
    "MU",
    "ELL",
    "potential_2d",
    "interaction_2d",
    "build_h2d",
    "interaction_diag",
]

MU: float = N2.mu  # N2 nuclear reduced mass (a.u.)
ELL: int = N2.ell  # fixed partial wave, l = 2


def interaction_2d(r: npt.ArrayLike, R: npt.ArrayLike) -> npt.ArrayLike:
    """`V_int(r,R) = -lambda(R) exp(-alpha_c r^2)` -- the perturbation ALONE."""
    return N2.v_int(r, R)


def potential_2d(r: npt.ArrayLike, R: npt.ArrayLike) -> npt.ArrayLike:
    """The full surface `v0(R) + l(l+1)/(2 r^2) + V_int(r,R)`.

    Must not coerce to a real dtype: `r`/`R` are complex on the ECS tails.
    """
    return N2.surface(r, R)


def interaction_diag(tgrid: TensorGrid) -> npt.NDArray[np.complex128]:
    """`V_int` evaluated on the tensor grid, flattened (C order)."""
    return N2.interaction_diag(tgrid)


def build_h2d(tgrid: TensorGrid) -> sp.csr_matrix:
    """`H_2D` on `tgrid` (axis 0 = electronic r, axis 1 = nuclear R)."""
    return N2.hamiltonian(tgrid)
