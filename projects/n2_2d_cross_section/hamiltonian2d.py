"""The exact 2-D electron-N2 potential surface and Hamiltonian.

    H_2D = -(1/2) d^2/dr^2 - (1/2 mu) d^2/dR^2
           + v0(R) + l(l+1)/(2 r^2) - lambda(R) exp(-alpha_c r^2)

verbatim eMoScat `Neutral2dPotential` (`source/Model2d/Potentials2d.cpp:18`),
built from the same model functions already ported and verified in
sub-project #2. No new potential physics.

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
from qscat.dvr import TensorGrid, hamiltonian_nd, potential_nd

from projects.n2_resonance.potential import PARAMS, v0, v_int

__all__ = [
    "MU",
    "ELL",
    "potential_2d",
    "interaction_2d",
    "build_h2d",
    "interaction_diag",
]

MU: float = 12766.36                       # N2 nuclear reduced mass (a.u.)
ELL: int = int(PARAMS["impulsemomentum"])  # fixed partial wave, l = 2


def interaction_2d(r: npt.ArrayLike, R: npt.ArrayLike) -> npt.ArrayLike:
    """`V_int(r,R) = -lambda(R) exp(-alpha_c r^2)` -- the perturbation ALONE."""
    return v_int(r, R)


def potential_2d(r: npt.ArrayLike, R: npt.ArrayLike) -> npt.ArrayLike:
    """The full surface `v0(R) + l(l+1)/(2 r^2) + V_int(r,R)`.

    Must not coerce to a real dtype: `r`/`R` are complex on the ECS tails.
    """
    rr = np.asarray(r)
    return v0(R) + ELL * (ELL + 1) / (2.0 * rr**2) + v_int(rr, R)


def interaction_diag(tgrid: TensorGrid) -> npt.NDArray[np.complex128]:
    """`V_int` evaluated on the tensor grid, flattened (C order)."""
    return potential_nd(tgrid, interaction_2d)


def build_h2d(tgrid: TensorGrid) -> sp.csr_matrix:
    """`H_2D` on `tgrid` (axis 0 = electronic r, axis 1 = nuclear R)."""
    return hamiltonian_nd(tgrid, [1.0, MU], potential_2d)
