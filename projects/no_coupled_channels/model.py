"""The coupled-channel model built on a two-centre well.

`CoupledModel` assembles the Lambda-block Hamiltonian for l = Lambda ...
Lambda + n_channels - 1, on a fixed-R electronic grid or on the full 2-D
tensor grid. `n_channels = 1` is not a degenerate case: it IS the fixed-l
model, the approximation under test, and it runs through exactly the same
code as the coupled one so that the comparison is differential.

`DiagonalChannelModel` presents a single channel as a plain `ResonanceModel`,
which lets the SHIPPED `qscat.core.lcp.local_complex_potential` be run on it
as an independent cross-check of the pole walk.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from qscat.dvr import TensorGrid, hamiltonian_nd, kinetic_sparse, potential_nd
from qscat.dvr.grid import FemDvrEcsGrid

from projects.no_coupled_channels.anisotropy import TwoCentreWell
from projects.no_coupled_channels.blocks import assemble_coupled

__all__ = ["CoupledModel", "DiagonalChannelModel"]


@dataclass(frozen=True)
class DiagonalChannelModel:
    """One channel of the coupled problem, as a `ResonanceModel`."""

    well: TwoCentreWell
    l: int

    @property
    def mu(self) -> float:
        """Nuclear reduced mass (a.u.), from the wrapped model."""
        return self.well.base.mu

    @property
    def ell(self) -> int:
        """This channel's partial wave."""
        return self.l

    @property
    def charge(self) -> int:
        """Coulomb charge of the residual channel, from the wrapped model."""
        return self.well.base.charge

    def v0(self, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """Neutral curve (Hartree), from the wrapped model."""
        return self.well.base.v0(R)

    def v_int(self, r: npt.ArrayLike, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """This channel's diagonal interaction `V_ll(r, R)`."""
        return self.well.v_block(self.l, self.l, r, R)

    def surface(self, r: npt.ArrayLike, R: npt.ArrayLike) -> npt.NDArray[np.complex128]:
        """`v0(R) + l(l+1)/(2 r^2) + V_ll(r, R)`."""
        rr = np.asarray(r, dtype=np.complex128)
        out = self.v0(R) + self.l * (self.l + 1) / (2.0 * rr**2) + self.v_int(rr, R)
        return np.asarray(out, dtype=np.complex128)

    def hamiltonian(self, tgrid: TensorGrid) -> sp.csr_matrix:
        """`H_2D` for this single channel."""
        return hamiltonian_nd(tgrid, [1.0, self.mu], self.surface)

    def interaction_diag(self, tgrid: TensorGrid) -> npt.NDArray[np.complex128]:
        """`V_ll` on the tensor grid, flattened."""
        return potential_nd(tgrid, self.v_int)


@dataclass(frozen=True)
class CoupledModel:
    """The Lambda-block coupled-channel model. `n_channels = 1` is the
    fixed-l model."""

    well: TwoCentreWell
    n_channels: int = 1

    @property
    def mu(self) -> float:
        """Nuclear reduced mass (a.u.), from the wrapped model."""
        return self.well.base.mu

    @property
    def ell(self) -> int:
        """The lowest channel's partial wave, `Lambda`."""
        return self.well.Lambda

    @property
    def charge(self) -> int:
        """Coulomb charge of the residual channel, from the wrapped model."""
        return self.well.base.charge

    def channel_ells(self) -> tuple[int, ...]:
        """The partial waves in the block: `Lambda, Lambda+1, ...`."""
        return tuple(self.well.Lambda + i for i in range(self.n_channels))

    def _coupling_table(
        self, r: npt.ArrayLike, R: npt.ArrayLike
    ) -> list[list[npt.NDArray[np.complex128] | None]]:
        """Off-diagonal `V_{ll'}` for every channel pair, flattened."""
        ells = self.channel_ells()
        n_ch = len(ells)
        table: list[list[npt.NDArray[np.complex128] | None]] = [[None] * n_ch for _ in range(n_ch)]
        for i, l in enumerate(ells):
            for j, lp in enumerate(ells):
                if i == j:
                    continue
                vals = np.broadcast_arrays(self.well.v_block(l, lp, r, R))[0]
                table[i][j] = np.asarray(vals, dtype=np.complex128).ravel()
        return table

    def electronic_hamiltonian(self, grid: FemDvrEcsGrid, R: complex) -> sp.csr_matrix:
        """The fixed-`R` coupled electronic Hamiltonian (mass-1 electron)."""
        r = np.asarray(grid.points, dtype=np.complex128)
        T = kinetic_sparse(grid, 1.0)
        diagonal = [
            sp.csr_matrix(
                T + sp.diags(DiagonalChannelModel(self.well, l).surface(r, R), format="csr")
            )
            for l in self.channel_ells()
        ]
        return assemble_coupled(diagonal, self._coupling_table(r, R))

    def interaction_matrix(self, tgrid: TensorGrid) -> sp.csr_matrix:
        """The coupled interaction `V_{ll'}` on `tgrid`, sparse, channel-outermost.

        The PERTURBATION alone: no kinetic energy, no `v0(R)` and no
        centrifugal term. Those belong to the free channel Hamiltonian, whose
        solutions `channel_vector` already supplies -- putting them here would
        drive the Lippmann-Schwinger equation with the wrong operator and
        silently produce a plausible wrong T-matrix.

        The single-channel sibling is `DiagonalChannelModel.interaction_diag`,
        which returns a flat array because one channel's interaction really is
        diagonal. Coupled it is not: the off-diagonal blocks ARE the coupling.
        """
        r, R = tgrid.points()
        diagonal = [
            sp.diags(
                np.asarray(self.well.v_block(l, l, r, R), dtype=np.complex128).ravel(),
                format="csr",
            )
            for l in self.channel_ells()
        ]
        return assemble_coupled(diagonal, self._coupling_table(r, R))

    def hamiltonian(self, tgrid: TensorGrid) -> sp.csr_matrix:
        """The coupled `H_2D` on `tgrid` (axis 0 electronic `r`, axis 1 `R`)."""
        diagonal = [
            hamiltonian_nd(tgrid, [1.0, self.mu], DiagonalChannelModel(self.well, l).surface)
            for l in self.channel_ells()
        ]
        r, R = tgrid.points()
        return assemble_coupled(diagonal, self._coupling_table(r, R))
