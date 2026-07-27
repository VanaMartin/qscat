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

    def energy_grid(self) -> npt.NDArray[np.float64]:
        return np.round(np.arange(self.e_lo, self.e_hi + 0.5 * self.e_step, self.e_step), 6)


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
    ),
}
