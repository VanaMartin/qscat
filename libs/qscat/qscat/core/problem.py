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
from typing import TYPE_CHECKING, Literal, overload

import numpy as np
import numpy.typing as npt

from .dissociation import da_cross_section, dr_cross_section
from .driven import ve_cross_section
from .time_dependent import (
    Method,
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

# Mirrors the private copies in driven.py / dissociation.py / lcp.py (scipy
# splu's permc_spec). The library-structure pass (lib-M12) consolidates all of
# them into a public `qscat.linalg.Ordering`; if that has already landed,
# import that name here instead of redefining.
_Ordering = Literal["NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD"]

# Return/parameter conventions, identical to the functional solvers mirrored
# below (see driven.py / dissociation.py / time_dependent.py).
_Sigma = npt.NDArray[np.float64]
_Psi = npt.NDArray[np.complex128] | None
_PsiOut = _Psi | list[_Psi]
_Amp = npt.NDArray[np.complex128]
_WpIn = dict[str, float]
_WpOut = dict[str, float]


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

    @overload
    def ve_cross_section(
        self,
        vprimes: list[int],
        E: float | npt.ArrayLike,
        *,
        ordering: _Ordering = ...,
        lam_scale: float = ...,
        return_wavefunction: Literal[False] = ...,
    ) -> _Sigma: ...

    @overload
    def ve_cross_section(
        self,
        vprimes: list[int],
        E: float | npt.ArrayLike,
        *,
        ordering: _Ordering = ...,
        lam_scale: float = ...,
        return_wavefunction: Literal[True],
    ) -> tuple[_Sigma, _PsiOut]: ...

    def ve_cross_section(
        self,
        vprimes: list[int],
        E: float | npt.ArrayLike,
        *,
        ordering: _Ordering = "COLAMD",
        lam_scale: float = 1.0,
        return_wavefunction: bool = False,
    ) -> _Sigma | tuple[_Sigma, _PsiOut]:
        """Vibrational-excitation cross section; same parameters, defaults and
        return convention as `qscat.core.ve_cross_section` (which this
        delegates to with the bundled grid/model/basis)."""
        if return_wavefunction:
            return ve_cross_section(
                self.grid,
                self.model,
                self.eps,
                self.chi,
                self.v_init,
                vprimes,
                E,
                ordering=ordering,
                lam_scale=lam_scale,
                return_wavefunction=True,
            )
        return ve_cross_section(
            self.grid,
            self.model,
            self.eps,
            self.chi,
            self.v_init,
            vprimes,
            E,
            ordering=ordering,
            lam_scale=lam_scale,
        )

    @overload
    def da_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = ...,
        ordering: _Ordering = ...,
        return_wavefunction: Literal[False] = ...,
    ) -> _Sigma: ...

    @overload
    def da_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = ...,
        ordering: _Ordering = ...,
        return_wavefunction: Literal[True],
    ) -> tuple[_Sigma, _PsiOut]: ...

    def da_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = 1,
        ordering: _Ordering = "COLAMD",
        return_wavefunction: bool = False,
    ) -> _Sigma | tuple[_Sigma, _PsiOut]:
        """Dissociative-attachment cross section; see `qscat.core.da_cross_section`."""
        if return_wavefunction:
            return da_cross_section(
                self.grid,
                self.model,
                self.eps,
                self.chi,
                self.v_init,
                E,
                n_channels=n_channels,
                ordering=ordering,
                return_wavefunction=True,
            )
        return da_cross_section(
            self.grid,
            self.model,
            self.eps,
            self.chi,
            self.v_init,
            E,
            n_channels=n_channels,
            ordering=ordering,
        )

    @overload
    def dr_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = ...,
        ordering: _Ordering = ...,
        return_wavefunction: Literal[False] = ...,
        return_amplitude: Literal[False] = ...,
    ) -> _Sigma: ...

    @overload
    def dr_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = ...,
        ordering: _Ordering = ...,
        return_wavefunction: Literal[True],
        return_amplitude: Literal[False] = ...,
    ) -> tuple[_Sigma, _PsiOut]: ...

    @overload
    def dr_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = ...,
        ordering: _Ordering = ...,
        return_wavefunction: Literal[False] = ...,
        return_amplitude: Literal[True],
    ) -> tuple[_Sigma, _Amp]: ...

    @overload
    def dr_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = ...,
        ordering: _Ordering = ...,
        return_wavefunction: Literal[True],
        return_amplitude: Literal[True],
    ) -> tuple[_Sigma, _PsiOut, _Amp]: ...

    def dr_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = 3,
        ordering: _Ordering = "COLAMD",
        return_wavefunction: bool = False,
        return_amplitude: bool = False,
    ) -> _Sigma | tuple[_Sigma, _PsiOut] | tuple[_Sigma, _Amp] | tuple[_Sigma, _PsiOut, _Amp]:
        """Dissociative-recombination cross section (ionic target); see
        `qscat.core.dr_cross_section`."""
        if return_wavefunction and return_amplitude:
            return dr_cross_section(
                self.grid,
                self.model,
                self.eps,
                self.chi,
                self.v_init,
                E,
                n_channels=n_channels,
                ordering=ordering,
                return_wavefunction=True,
                return_amplitude=True,
            )
        if return_wavefunction:
            return dr_cross_section(
                self.grid,
                self.model,
                self.eps,
                self.chi,
                self.v_init,
                E,
                n_channels=n_channels,
                ordering=ordering,
                return_wavefunction=True,
            )
        if return_amplitude:
            return dr_cross_section(
                self.grid,
                self.model,
                self.eps,
                self.chi,
                self.v_init,
                E,
                n_channels=n_channels,
                ordering=ordering,
                return_amplitude=True,
            )
        return dr_cross_section(
            self.grid,
            self.model,
            self.eps,
            self.chi,
            self.v_init,
            E,
            n_channels=n_channels,
            ordering=ordering,
        )

    # --- time-dependent observables -----------------------------------------

    def td_ve_cross_section(
        self,
        vprimes: list[int],
        E: float | npt.ArrayLike,
        *,
        dt: float,
        n_steps: int,
        wp_in: _WpIn,
        wp_out: _WpOut | None = None,
        order: int = 3,
        subtract_free_reference: bool = True,
        method: Method = "tw",
        position: int | None = None,
        surface: int | None = None,
    ) -> npt.NDArray[np.float64]:
        """Time-dependent VE cross section; see `qscat.core.td_ve_cross_section`
        (same method/`wp_out`/`position`/`surface` contract)."""
        return td_ve_cross_section(
            self.grid,
            self.model,
            self.eps,
            self.chi,
            self.v_init,
            vprimes,
            E,
            dt=dt,
            n_steps=n_steps,
            wp_in=wp_in,
            wp_out=wp_out,
            order=order,
            subtract_free_reference=subtract_free_reference,
            method=method,
            position=position,
            surface=surface,
        )

    def td_ve_cross_sections_all(
        self,
        vprimes: list[int],
        E: float | npt.ArrayLike,
        *,
        dt: float,
        n_steps: int,
        wp_in: _WpIn,
        wp_out: _WpOut,
        position: int,
        surface: int,
        order: int = 3,
        subtract_free_reference: bool = True,
    ) -> dict[str, npt.NDArray[np.float64]]:
        """All three TD-VE extractors from ONE propagation; see
        `qscat.core.td_ve_cross_sections_all`."""
        return td_ve_cross_sections_all(
            self.grid,
            self.model,
            self.eps,
            self.chi,
            self.v_init,
            vprimes,
            E,
            dt=dt,
            n_steps=n_steps,
            wp_in=wp_in,
            wp_out=wp_out,
            position=position,
            surface=surface,
            order=order,
            subtract_free_reference=subtract_free_reference,
        )

    def td_da_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        dt: float,
        n_steps: int,
        wp_in: _WpIn,
        method: Method = "flow",
        surface: int | None = None,
        position: int | None = None,
        wp_out: _WpOut | None = None,
        n_channels: int = 1,
        order: int = 3,
    ) -> npt.NDArray[np.float64]:
        """Time-dependent DA cross section; see `qscat.core.td_da_cross_section`."""
        return td_da_cross_section(
            self.grid,
            self.model,
            self.eps,
            self.chi,
            self.v_init,
            E,
            dt=dt,
            n_steps=n_steps,
            wp_in=wp_in,
            method=method,
            surface=surface,
            position=position,
            wp_out=wp_out,
            n_channels=n_channels,
            order=order,
        )

    def td_da_cross_sections_all(
        self,
        E: float | npt.ArrayLike,
        *,
        dt: float,
        n_steps: int,
        wp_in: _WpIn,
        surface: int,
        position: int,
        wp_out: _WpOut,
        n_channels: int = 1,
        order: int = 3,
    ) -> dict[str, npt.NDArray[np.float64]]:
        """All three TD-DA extractors from ONE propagation; see
        `qscat.core.td_da_cross_sections_all`."""
        return td_da_cross_sections_all(
            self.grid,
            self.model,
            self.eps,
            self.chi,
            self.v_init,
            E,
            dt=dt,
            n_steps=n_steps,
            wp_in=wp_in,
            surface=surface,
            position=position,
            wp_out=wp_out,
            n_channels=n_channels,
            order=order,
        )
