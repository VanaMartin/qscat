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

This is the recommended entry point for every observable it carries: the TI
cross sections (VE/DA/DR) and TD cross sections (VE/DA), the LCP DA cross
section, the BO resonance levels, the exact 2-D resonance states, and the
NRM VE/DA cross sections. The
functional solvers remain public (they are the low-level layer this delegates
to, and each carries ADR 0004's *provisional* marker pending the pre-1.0
signature freeze); `ScatteringProblem` is the stable API. Deliberately NOT on
the facade: `lcp_resonance_levels` (its inputs are hand-built curves on
angle-paired grids -- `resonance_levels` here is the model-first route that
computes them internally), the `td_nrm_*` solvers (knob surface still
settling), and the curve/state builders (`local_complex_potential`,
`resonance_pole_walk`, `qscat.core.bo.*` -- ingredients, not observables).

The vibrational basis (`eps`, `chi`) is solved once, on construction, from
`model.mu`/`model.v0` on the nuclear grid, and reused across every observable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import numpy as np
import numpy.typing as npt

from .dissociation import da_cross_section, dr_cross_section
from .driven import ve_cross_section
from .lcp import lcp_da_cross_section, resonance_levels
from .resonance import exact_resonance_states
from .time_dependent import (
    Method,
    td_da_cross_section,
    td_da_cross_sections_all,
    td_ve_cross_section,
    td_ve_cross_sections_all,
)
from .vibrational import VibrationalBasis, vibrational_states

if TYPE_CHECKING:
    from qscat.dvr import FemDvrEcsGrid, TensorGrid
    from qscat.model import ResonanceModel

    from .lcp import ResonanceLevels
    from .nrm import DiscreteState, NrmIngredients
    from .resonance import ExactResonanceStates

__all__ = ["ScatteringProblem"]

# Mirrors the private copies in driven.py / dissociation.py / lcp.py (scipy
# splu's permc_spec). The library-structure pass (lib-M12) consolidates all of
# them into a public `qscat.linalg.Ordering`; if that has already landed,
# import that name here instead of redefining.
_Ordering = Literal["NATURAL", "MMD_ATA", "MMD_AT_PLUS_A", "COLAMD"]
_Window = tuple[float, float, float, float]

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

    @property
    def _bundle(
        self,
    ) -> tuple[
        TensorGrid,
        ResonanceModel,
        npt.NDArray[np.float64],
        npt.NDArray[np.complex128],
        int,
    ]:
        """The `(grid, model, eps, chi, v_init)` group every solver takes first."""
        return (self.grid, self.model, self.eps, self.chi, self.v_init)

    # --- time-independent observables ---------------------------------------

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
        delegates to with the bundled grid/model/basis). Returns the plain
        sigma array unless the `return_*` flag is set (see the union in the
        signature); literal-flag callers who want the narrowed type can call
        the functional `qscat.core.ve_cross_section` directly."""
        return ve_cross_section(
            *self._bundle,
            vprimes,
            E,
            ordering=ordering,
            lam_scale=lam_scale,
            return_wavefunction=return_wavefunction,
        )

    def da_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        n_channels: int = 1,
        ordering: _Ordering = "COLAMD",
        return_wavefunction: bool = False,
    ) -> _Sigma | tuple[_Sigma, _PsiOut]:
        """Dissociative-attachment cross section; see `qscat.core.da_cross_section`.
        Returns the plain sigma array unless the `return_*` flag is set (see
        the union in the signature); literal-flag callers who want the
        narrowed type can call the functional `qscat.core.da_cross_section`
        directly."""
        return da_cross_section(
            *self._bundle,
            E,
            n_channels=n_channels,
            ordering=ordering,
            return_wavefunction=return_wavefunction,
        )

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
        `qscat.core.dr_cross_section`. Returns the plain sigma array unless
        the `return_*` flag is set (see the union in the signature);
        literal-flag callers who want the narrowed type can call the
        functional `qscat.core.dr_cross_section` directly."""
        return dr_cross_section(
            *self._bundle,
            E,
            n_channels=n_channels,
            ordering=ordering,
            return_wavefunction=return_wavefunction,
            return_amplitude=return_amplitude,
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
            *self._bundle,
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
            *self._bundle,
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
            *self._bundle,
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
            *self._bundle,
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

    # --- LCP / resonance observables ----------------------------------------

    def lcp_da_cross_section(
        self,
        E: float | npt.ArrayLike,
        *,
        Vd: npt.NDArray[np.complex128],
        Gamma: npt.NDArray[np.float64],
        ordering: _Ordering = "COLAMD",
        return_wavefunction: bool = False,
    ) -> _Sigma | tuple[_Sigma, _PsiOut]:
        """LCP dissociative-attachment cross section on this problem's NUCLEAR
        grid; see `qscat.core.lcp_da_cross_section`. The curve `(Vd, Gamma)`
        is per-call (compute it with `resonance_levels(..., return_curve=True)`
        -- see that docstring for why not `local_complex_potential` directly);
        `mu`/`eps`/`chi`/`v_init` come from the bundle, which is what pays down
        the functional signature's documented argument-order exception
        (docs/adr/0007). The LCP magnitude needs the FINE per-molecule nuclear
        deck -- construct the problem on it for physical numbers. Returns the
        plain sigma array unless `return_wavefunction` is set (see the union
        in the signature); literal-flag callers who want the narrowed type
        can call the functional `qscat.core.lcp_da_cross_section` directly."""
        return lcp_da_cross_section(
            self.grid.grids[1],
            self.model.mu,
            Vd,
            Gamma,
            self.eps,
            self.chi,
            self.v_init,
            E,
            ordering=ordering,
            return_wavefunction=return_wavefunction,
        )

    # EXTENSION POINT: when `qscat.core.lcp.lcp_ve_cross_section` (the LCP
    # VIBRATIONAL-EXCITATION route, currently a validation-layer driver)
    # graduates into qscat.core.lcp, add the matching typed method here,
    # mirroring `lcp_da_cross_section` above: the bundle supplies
    # grid/mu/eps/chi/v_init, the curve arrives as `Vd`/`Gamma` keywords.

    def resonance_levels(
        self,
        nuclear_grid_b: FemDvrEcsGrid,
        elec_grid_b: FemDvrEcsGrid,
        *,
        re_half_width: float = 0.05,
        im_half_width: float = 0.05,
        resid_tol: float = 1e-3,
        window: _Window | None = None,
        n_levels: int | None = None,
        rel_tol: float = 1e-4,
        atol: float = 1e-8,
        golden_rule: bool = True,
        return_curve: bool = False,
    ) -> (
        ResonanceLevels
        | tuple[ResonanceLevels, npt.NDArray[np.complex128], npt.NDArray[np.float64]]
    ):
        """Born-Oppenheimer quasi-bound levels of this problem's anion; see
        `qscat.core.resonance_levels`. This problem's own electronic/nuclear
        grids are the `_a` partners; `nuclear_grid_b`/`elec_grid_b` are the
        angle-moved partners (share every real node, differ only in the ECS
        tail angle -- `qscat.core.grids.ecs_angle_family` builds a valid
        family). `return_curve=True` also returns the `(Vd, Gamma)` curve the
        levels were computed in -- the input `lcp_da_cross_section` needs.
        Returns the plain levels unless `return_curve` is set (see the union
        in the signature); literal-flag callers who want the narrowed type
        can call the functional `qscat.core.resonance_levels` directly."""
        return resonance_levels(
            self.model,
            self.grid.grids[1],
            nuclear_grid_b,
            self.grid.grids[0],
            elec_grid_b,
            re_half_width=re_half_width,
            im_half_width=im_half_width,
            resid_tol=resid_tol,
            window=window,
            n_levels=n_levels,
            rel_tol=rel_tol,
            atol=atol,
            golden_rule=golden_rule,
            return_curve=return_curve,
        )

    def exact_resonance_states(
        self,
        grid_electronic: TensorGrid,
        grid_nuclear: TensorGrid,
        *,
        shifts: npt.ArrayLike,
        window: _Window,
        k: int = 8,
        rel_tol: float = 1e-4,
        atol: float = 1e-8,
    ) -> ExactResonanceStates:
        """Exact 2-D resonance states by two-angle ECS stability; see
        `qscat.core.exact_resonance_states`. This problem's grid is the base;
        `grid_electronic`/`grid_nuclear` are the one-angle-moved partner
        TensorGrids (`ecs_angle_family` builds all three consistently).
        Seeds (`shifts`) are passed in -- typically `resonance_levels`'s
        output -- so the exact solver never depends on the approximation it
        measures."""
        return exact_resonance_states(
            self.model,
            self.grid,
            grid_electronic,
            grid_nuclear,
            shifts=shifts,
            window=window,
            k=k,
            rel_tol=rel_tol,
            atol=atol,
        )

    # --- nonlocal resonance model (NRM) observables --------------------------

    def nrm_ve_cross_section(
        self,
        phi_d: DiscreteState,
        vprimes: list[int],
        E: float | npt.ArrayLike,
        *,
        ingredients: NrmIngredients | None = None,
        n_states: int | None = None,
        include_background: bool = True,
    ) -> npt.NDArray[np.float64]:
        """NRM vibrational-excitation cross section; see
        `qscat.core.nrm.nrm_ve_cross_section` (this problem's nuclear and
        electronic grids fill its leading NUCLEAR-grid-first pair)."""
        # Deferred import: `import qscat.core` must never pull `nrm` in
        # (the hard boundary documented in qscat.core.__init__).
        from .nrm import nrm_ve_cross_section

        return nrm_ve_cross_section(
            self.grid.grids[1],
            self.grid.grids[0],
            self.model,
            phi_d,
            self.eps,
            self.chi,
            self.v_init,
            vprimes,
            E,
            ingredients=ingredients,
            n_states=n_states,
            include_background=include_background,
        )

    def nrm_da_cross_section(
        self,
        phi_d: DiscreteState,
        E: float | npt.ArrayLike,
        *,
        ingredients: NrmIngredients | None = None,
        n_states: int | None = None,
    ) -> npt.NDArray[np.float64]:
        """NRM dissociative-attachment cross section; see
        `qscat.core.nrm.nrm_da_cross_section`."""
        from .nrm import nrm_da_cross_section  # deferred: see nrm_ve_cross_section

        return nrm_da_cross_section(
            self.grid.grids[1],
            self.grid.grids[0],
            self.model,
            phi_d,
            self.eps,
            self.chi,
            self.v_init,
            E,
            ingredients=ingredients,
            n_states=n_states,
        )
