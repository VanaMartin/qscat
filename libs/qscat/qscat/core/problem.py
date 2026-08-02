"""`ScatteringProblem` — the high-level, object API over the core solvers.

The functional solvers (`ve_cross_section`, `da_cross_section`, …) all share the
same `(tgrid, model, eps, chi, v_init, …)` leading argument group. `ScatteringProblem`
bundles that group once — you give it a grid, a model, and how many vibrational
states to solve — and exposes each observable as a method taking only the
per-call arguments (`vprimes`, `E`, and the keyword options):

    from qscat.core import ScatteringProblem
    from qscat.dvr import TensorGrid
    from qscat.model import N2
    from qscat.core.grids import electronic_grid, nuclear_grid

    prob = ScatteringProblem(
        grid=TensorGrid([electronic_grid(r_max=16.0, order=7, n_complex=5),
                         nuclear_grid(r_max=22.0, quadrature=10, n_complex=5)]),
        model=N2,
        n_vib=4,
    )
    sigma = prob.ve_cross_section(vprimes=[0, 1, 2], E=[0.10, 0.15, 0.20])

This is the recommended entry point. The functional solvers remain public (they
are the low-level layer this delegates to, and are marked provisional pending any
future signature change — see ADR 0004); `ScatteringProblem` is stable API.

The vibrational basis (`eps`, `chi`) is solved once, on construction, from
`model.mu`/`model.v0` on the nuclear grid, and reused across every observable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np
import numpy.typing as npt

from .dissociation import da_cross_section, dr_cross_section
from .driven import ve_cross_section
from .time_dependent import (
    td_da_cross_section,
    td_da_cross_sections_all,
    td_ve_cross_section,
    td_ve_cross_sections_all,
)
from .vibrational import VibrationalBasis, vibrational_states

if TYPE_CHECKING:
    from qscat.dvr import TensorGrid
    from qscat.model import ResonanceModel

__all__ = ["ScatteringProblem"]


@dataclass(frozen=True)
class ScatteringProblem:
    """A fully-specified electron-diatomic scattering problem.

    Bundles the grid, model, and vibrational basis shared by every observable.
    `grid` is the electronic x nuclear `TensorGrid`; `model` is any
    `qscat.model.ResonanceModel`; `n_vib` is the number of vibrational states to
    solve; `v_init` is the initial vibrational level (default 0). The basis is
    solved once on construction and exposed as `.eps` / `.chi` / `.basis`.
    """

    grid: TensorGrid
    model: ResonanceModel
    n_vib: int
    v_init: int = 0
    basis: VibrationalBasis = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        basis = vibrational_states(self.grid.grids[1], self.model.mu, self.n_vib, self.model.v0)
        object.__setattr__(self, "basis", basis)  # frozen dataclass: set via object

    @property
    def eps(self) -> npt.NDArray[np.float64]:
        """Vibrational energies (Ha), ascending."""
        return self.basis.eps

    @property
    def chi(self) -> npt.NDArray[np.complex128]:
        """Vibrational eigenvectors, one row per level."""
        return self.basis.chi

    # --- time-independent observables ---------------------------------------

    def ve_cross_section(
        self, vprimes: list[int], E: float | npt.ArrayLike, **kwargs: Any
    ) -> Any:
        """Vibrational-excitation cross section; see `qscat.core.ve_cross_section`."""
        return ve_cross_section(
            self.grid, self.model, self.eps, self.chi, self.v_init, vprimes, E, **kwargs
        )

    def da_cross_section(self, E: float | npt.ArrayLike, **kwargs: Any) -> Any:
        """Dissociative-attachment cross section; see `qscat.core.da_cross_section`."""
        return da_cross_section(
            self.grid, self.model, self.eps, self.chi, self.v_init, E, **kwargs
        )

    def dr_cross_section(self, E: float | npt.ArrayLike, **kwargs: Any) -> Any:
        """Dissociative-recombination cross section; see `qscat.core.dr_cross_section`."""
        return dr_cross_section(
            self.grid, self.model, self.eps, self.chi, self.v_init, E, **kwargs
        )

    # --- time-dependent observables -----------------------------------------

    def td_ve_cross_section(
        self, vprimes: list[int], E: float | npt.ArrayLike, **kwargs: Any
    ) -> Any:
        """Time-dependent VE cross section; see `qscat.core.td_ve_cross_section`."""
        return td_ve_cross_section(
            self.grid, self.model, self.eps, self.chi, self.v_init, vprimes, E, **kwargs
        )

    def td_ve_cross_sections_all(
        self, vprimes: list[int], E: float | npt.ArrayLike, **kwargs: Any
    ) -> Any:
        """All three TD-VE extractors from one propagation; see the functional twin."""
        return td_ve_cross_sections_all(
            self.grid, self.model, self.eps, self.chi, self.v_init, vprimes, E, **kwargs
        )

    def td_da_cross_section(self, E: float | npt.ArrayLike, **kwargs: Any) -> Any:
        """Time-dependent DA cross section; see `qscat.core.td_da_cross_section`."""
        return td_da_cross_section(
            self.grid, self.model, self.eps, self.chi, self.v_init, E, **kwargs
        )

    def td_da_cross_sections_all(self, E: float | npt.ArrayLike, **kwargs: Any) -> Any:
        """All three TD-DA extractors from one propagation; see the functional twin."""
        return td_da_cross_sections_all(
            self.grid, self.model, self.eps, self.chi, self.v_init, E, **kwargs
        )
