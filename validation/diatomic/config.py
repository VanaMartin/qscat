"""Per-molecule numerical config for the diatomic VE-scattering validation.

The *physics* (potential form + parameters) lives in `qscat.model` (`N2`/`NO`/
`F2`); this module holds the *numerical* choices that are a convergence/study
decision rather than a model property: the FEM-DVR-ECS grid extents and the
energy grid over which the exact-2D σ(E) oracle is computed. Adding a molecule
= one `MoleculeConfig` entry here (plus a `qscat.model` registry entry).

Grid: the N₂-style FEM-DVR-ECS layout (electronic r_max=16 / nuclear r_max=22)
is converged for NO and F₂ as well -- verified directly: the exact-2D σ(E) is
unchanged (<1%) at electronic r_max = 16/24/32 for NO. The resonances sit LOW
(NO ~0.02-0.05 Ha, a P-wave shape resonance; F₂ ~0.01-0.04 Ha, near threshold,
weakly bound), so the energy grids are correspondingly low and dense.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.core.grids import electronic_grid, segmented_grid
from qscat.dvr import FemDvrEcsGrid, TensorGrid
from qscat.model import F2, NO, DiatomicResonanceModel

__all__ = ["MoleculeConfig", "CONFIGS"]


@dataclass(frozen=True)
class MoleculeConfig:
    """The model + the numerical grid/energy config for one molecule's oracle."""

    name: str
    model: DiatomicResonanceModel
    # electronic FEM-DVR-ECS grid
    e_r_max: float
    e_order: int
    e_n_complex: int
    # nuclear FEM-DVR-ECS grid
    n_r_max: float
    n_quadrature: int
    n_n_complex: int
    n_vib: int  # number of neutral vibrational states to resolve
    e_lo: float  # energy grid (Hartree)
    e_hi: float
    e_step: float
    # eMoScat per-molecule nuclear deck for DA (transcribed verbatim from
    # reference/eMoScat/input/{NO,F2}/grids.txt, 2nd/nuclear declaration):
    # (n_elements, endpoint) segment pairs, per qscat.core.grids.segmented_grid.
    nuc_real: tuple[tuple[int, float], ...]
    nuc_complex: tuple[tuple[int, float], ...]
    nuc_angle: float
    nuc_quad: int
    # LCP fixed-R electronic ECS angles (two-angle pole-matching, sub-project B)
    lcp_angle_a: float = 35.0
    lcp_angle_b: float = 44.0

    def energy_grid(self) -> npt.NDArray[np.float64]:
        return np.round(np.arange(self.e_lo, self.e_hi + 0.5 * self.e_step, self.e_step), 6)

    def da_grid(self) -> TensorGrid:
        """Electronic (VE-validated, r_max=e_r_max) x eMoScat nuclear deck.

        A NEW grid path used only by DA: the shared N2-style nuclear grid
        (`n_r_max`/`n_quadrature`/`n_n_complex`, used by `energy_grid`'s VE
        callers via `curves.build_grid`) under-resolves the K_R~58
        dissociation wave. eMoScat's per-molecule nuclear deck
        (`nuc_real`/`nuc_complex`/`nuc_angle`/`nuc_quad`) resolves it.
        """
        return TensorGrid(
            [
                electronic_grid(r_max=self.e_r_max, order=self.e_order, n_complex=self.e_n_complex),
                segmented_grid(
                    self.nuc_real,
                    self.nuc_complex,
                    angle_deg=self.nuc_angle,
                    quadrature=self.nuc_quad,
                ),
            ]
        )

    def lcp_elec_grids(self) -> tuple[FemDvrEcsGrid, FemDvrEcsGrid]:
        """The two ECS-angle-matched fixed-R electronic grids for LCP pole finding."""
        return (
            electronic_grid(
                r_max=self.e_r_max, order=self.e_order, n_complex=self.e_n_complex,
                angle_deg=self.lcp_angle_a,
            ),
            electronic_grid(
                r_max=self.e_r_max, order=self.e_order, n_complex=self.e_n_complex,
                angle_deg=self.lcp_angle_b,
            ),
        )

    def lcp_nuclear_grid(self) -> FemDvrEcsGrid:
        """The eMoScat per-molecule nuclear deck (same fine grid as `da_grid()`)."""
        return segmented_grid(
            self.nuc_real,
            self.nuc_complex,
            angle_deg=self.nuc_angle,
            quadrature=self.nuc_quad,
        )


# The N2-style grid (converged for NO/F2, see module docstring). Energy ranges
# bracket each molecule's low-lying resonance.
CONFIGS: dict[str, MoleculeConfig] = {
    "NO": MoleculeConfig(
        name="NO",
        model=NO,
        e_r_max=16.0,
        e_order=8,
        e_n_complex=6,
        n_r_max=22.0,
        n_quadrature=10,
        n_n_complex=5,
        n_vib=4,
        e_lo=0.004,
        e_hi=0.120,
        e_step=0.004,
        nuc_real=((1, 1.0), (1, 1.6), (37, 9.0)),
        nuc_complex=((1, 9.25), (1, 10.0), (1, 12.0), (4, 42.0)),
        nuc_angle=45.0,
        nuc_quad=14,
    ),
    "F2": MoleculeConfig(
        name="F2",
        model=F2,
        e_r_max=16.0,
        e_order=8,
        e_n_complex=6,
        n_r_max=22.0,
        n_quadrature=10,
        n_n_complex=5,
        n_vib=4,
        e_lo=0.004,
        e_hi=0.100,
        e_step=0.004,
        nuc_real=((9, 1.8), (1, 2.0), (5, 2.5), (4, 2.596908), (4, 2.7), (40, 10.7)),
        nuc_complex=(
            (1, 10.8),
            (1, 11.0),
            (1, 11.5),
            (1, 12.5),
            (1, 14.0),
            (1, 18.0),
            (4, 30.0),
            (2, 101.0),
        ),
        nuc_angle=35.0,
        nuc_quad=14,
    ),
}
