"""NRM vs LCP vs exact-2D VIBRATIONAL EXCITATION for N2 and F2.

The sibling of `validation.diatomic.nrm` (dissociative attachment), and the
comparison that has published curves behind it. PRA 77 plots a VE cross
section for **every** molecule in the study -- N2 in Fig. 4 (choice A) and
Fig. 8 (choice B), F2 0->1 in Fig. 6 (A) and Fig. 8 (B) -- whereas its only
published DA panel is F2's. VE is also O(1) rather than exponentially
sensitive near threshold, so a ratio against the exact oracle means what it
appears to mean.

WHICH PAIRS HAVE AN EXTERNAL ANCHOR. From the paper's own figure panels
(inventory in `reference/literature/houfek-2008-pra77-012710.md`): N2 0->0 is
in Figs. 4 and 8; F2 0->1 is in Figs. 6 and 8; **N2 0->1 appears only in
Fig. 4** -- Fig. 8 omits it, "because results of all calculations are
practically the same in this particular case" (p. 012710-10, Fig. 8 caption,
quoted in that note), which is itself a testable claim; and **F2 0->0 is not
plotted at all**. All four pairs are gated here against OUR exact 2-D solver,
which exists for every transition. What F2 0->0 lacks is external
corroboration, not an oracle.

THAT CAPTION IS A STATEMENT ABOUT A LINEAR AXIS, NOT ABOUT A RATIO. Figs. 4-6
and 8 are plotted linearly in units of a0^2 (p. 012710-8, Fig. 4), so two
curves are "practically the same" when their ABSOLUTE difference is small on
the panel -- which it can be while their ratio is far from 1, wherever sigma
itself is small. `test_ve_nrm.py` measures both and keeps them apart.

FOUR ROUTES per (molecule, v'):
  exact  -- `qscat.core.driven.ve_cross_section`, the 2-D oracle
  lcp    -- `qscat.core.lcp.local_complex_potential` +
            `projects.n2_ti_cross_section.ve_cross_section`
  A      -- NRM, `PhysicalDiscreteState` (the scattering function at
            Re E_res(R)), with and without the Eq. (37) background
  B      -- NRM, `AsymptoticDiscreteState` (the R-independent bound state),
            with and without the background

ENERGY WINDOWS COME FROM THE PAPER: N2 VE is plotted over 0.05-0.17 Ha
(Fig. 4 top), F2 VE 0->1 over 0-0.10 Ha (Fig. 6 top). `test_ve_nrm.py` stays
inside those.

GRIDS. N2's deck is `qscat_run.presets`' `N2:emoscat` TI grid (107 electronic
x 251 nuclear); F2's is `validation.diatomic.config`'s `da_grid()` (132 x 974),
which `test_da_grid.py` locks byte-identical to the same preset's TI grid. The
NRM's second electronic ECS angle is 40 deg against the decks' own 35, the
pairing `validation/diatomic/nrm.py` and `qscat_run.presets.resolve_nrm_grids`
both use.

`validation/` may import `qscat` and `projects`; the reverse is forbidden.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.core.driven import ve_cross_section as exact_ve_cross_section
from qscat.core.grids import electronic_grid
from qscat.core.lcp import local_complex_potential
from qscat.core.nrm import (
    AsymptoticDiscreteState,
    DiscreteState,
    NrmIngredients,
    PhysicalDiscreteState,
    nrm_ingredients,
    nrm_ve_cross_section,
)
from qscat.core.vibrational import vibrational_states
from qscat.dvr import FemDvrEcsGrid, TensorGrid
from qscat.model import ResonanceModel
from qscat_run.presets import MODELS, PRESETS

from projects.n2_ti_cross_section.cross_section import (
    ve_cross_section as lcp_ve_cross_section,
)

from .config import CONFIGS

__all__ = ["VeComparison", "VeSetup", "compare", "discrete_states", "nrm_sigma", "setup"]

# The second electronic grid (a different ECS angle) both the LCP pole walk
# and `PhysicalDiscreteState` need for two-angle pole stability. 40 deg against
# the decks' own 35 deg -- the pairing `validation/diatomic/nrm.py` measured
# its numbers on and `qscat_run.presets.resolve_nrm_grids` ships.
_ANGLE_B_DEG = 40.0

# Initial vibrational level for every route.
_V_INIT = 0


@dataclass(frozen=True)
class VeSetup:
    """The per-molecule inputs every route shares.

    Built once by `setup` so an energy sweep, a state-sum ladder and the two
    discrete-state choices all run against exactly the same grids and
    vibrational basis -- the comparison is only meaningful if they do.

    Attributes
    ----------
    molecule : str
        `"N2"` or `"F2"`.
    model : ResonanceModel
        The molecule's model.
    tgrid : TensorGrid
        Electronic x nuclear product grid (the exact-2D solver's grid).
    elec, nuc : FemDvrEcsGrid
        The two factors of `tgrid`.
    elec_b : FemDvrEcsGrid
        The second-angle electronic grid used for pole matching.
    eps, chi : ndarray
        Neutral vibrational energies and states on `nuc`.
    R_desc : ndarray
        `nuc`'s real DVR nodes, strictly descending -- what both the discrete
        states and `nrm_ingredients` are built on (`nonlocal_operator`
        requires the ingredient nodes to be exactly these).
    """

    molecule: str
    model: ResonanceModel
    tgrid: TensorGrid
    elec: FemDvrEcsGrid
    nuc: FemDvrEcsGrid
    elec_b: FemDvrEcsGrid
    eps: npt.NDArray[np.float64]
    chi: npt.NDArray[np.complex128]
    R_desc: npt.NDArray[np.float64]


@dataclass(frozen=True)
class VeComparison:
    """One molecule's `sigma_{0->v'}(E)` by all six curves.

    Every `sigma_*` array has shape `(len(energies), len(vprimes))`, in
    `qscat.core.driven.ve_cross_section`'s own convention (bohr^2).
    """

    molecule: str
    vprimes: list[int]
    energies: npt.NDArray[np.float64]
    sigma_exact: npt.NDArray[np.float64]
    sigma_lcp: npt.NDArray[np.float64]
    sigma_nrm_a: npt.NDArray[np.float64]
    sigma_nrm_a_nobg: npt.NDArray[np.float64]
    sigma_nrm_b: npt.NDArray[np.float64]
    sigma_nrm_b_nobg: npt.NDArray[np.float64]


# The ELECTRONIC deck parameters per molecule, transcribed from the decks
# below so the second-ECS-angle grid can be rebuilt with the same shape.
# `setup` asserts the rebuild at the deck's own angle reproduces the deck's
# electronic factor node-for-node, so a drift here cannot pass silently.
_ELEC_PARAMS: dict[str, dict[str, float | int]] = {
    # `qscat_run.presets._n2_ti_grid`
    "N2": {"r_max": 16.0, "order": 7, "n_complex": 5, "angle_deg": 35.0},
    # `validation.diatomic.config.CONFIGS["F2"]` (== `presets._f2_ti_grid`)
    "F2": {"r_max": 16.0, "order": 8, "n_complex": 6, "angle_deg": 35.0},
}


def _deck(molecule: str) -> tuple[TensorGrid, ResonanceModel, int]:
    """`(tensor grid, model, n_vib)` for `molecule` -- the deck each route runs on.

    N2 comes from `qscat_run.presets` (it has no `validation/` deck of its
    own); F2 from `validation.diatomic.config`, whose `da_grid()` is locked
    byte-identical to the same preset's TI grid by `test_da_grid.py`.
    """
    if molecule == "N2":
        preset = PRESETS["N2:emoscat"]
        return preset.ti_grid(), MODELS["N2"], preset.n_vib
    cfg = CONFIGS[molecule]
    return cfg.da_grid(), cfg.model, cfg.n_vib


def setup(molecule: str) -> VeSetup:
    """Build the shared grids and vibrational basis for `molecule`.

    Parameters
    ----------
    molecule : str
        `"N2"` or `"F2"`.

    Returns
    -------
    VeSetup

    Raises
    ------
    KeyError
        If `molecule` is neither `"N2"` nor `"F2"`.
    ValueError
        If `_ELEC_PARAMS[molecule]` no longer reproduces the deck's own
        electronic grid -- the second-angle grid would then be built on a
        different discretisation from the one every other route uses.
    """
    params = dict(_ELEC_PARAMS[molecule])
    tgrid, model, n_vib = _deck(molecule)
    elec, nuc = tgrid.grids

    check = electronic_grid(**params)  # type: ignore[arg-type]
    if check.n != elec.n or not np.allclose(check.points, elec.points):
        raise ValueError(f"_ELEC_PARAMS[{molecule!r}] no longer matches the deck's electronic grid")
    elec_b = electronic_grid(**{**params, "angle_deg": _ANGLE_B_DEG})  # type: ignore[arg-type]

    eps, chi = vibrational_states(nuc, model.mu, n_vib, model.v0)
    real = nuc.points.imag == 0.0
    R_desc = np.sort(nuc.points[real].real)[::-1]
    return VeSetup(
        molecule=molecule,
        model=model,
        tgrid=tgrid,
        elec=elec,
        nuc=nuc,
        elec_b=elec_b,
        eps=eps,
        chi=chi,
        R_desc=R_desc,
    )


def nrm_sigma(
    s: VeSetup,
    phi_d: DiscreteState,
    energies: npt.NDArray[np.float64],
    vprimes: list[int],
    *,
    n_states: int | None,
    include_background: bool,
    ingredients: NrmIngredients | None = None,
) -> npt.NDArray[np.float64]:
    """`sigma_{0->v\'}(E)` in the NRM for one discrete-state choice.

    Parameters
    ----------
    s : VeSetup
        The shared per-molecule setup.
    phi_d : DiscreteState
        Choice A (`PhysicalDiscreteState`) or B (`AsymptoticDiscreteState`).
    energies : ndarray
        Incident electron kinetic energies (hartree).
    vprimes : list of int
        Final vibrational levels.
    n_states : int, optional
        Eq. (60) state-sum truncation; `None` uses every state.
    include_background : bool
        Add the Eq. (37) background T-matrix. `True` is the paper's
        "nonlocal + bg" curve, `False` its bare "nonlocal" one.
    ingredients : NrmIngredients, optional
        Prebuilt ingredients for `phi_d`; built here if omitted. They are the
        expensive, energy-INDEPENDENT half, so pass them in when sweeping
        `n_states` or toggling the background.

    Returns
    -------
    ndarray
        `sigma`, shape `(len(energies), len(vprimes))`, bohr^2.
    """
    ing = (
        ingredients
        if ingredients is not None
        else nrm_ingredients(s.elec, s.model, phi_d, s.R_desc)
    )
    return np.asarray(
        nrm_ve_cross_section(
            s.nuc,
            s.elec,
            s.model,
            phi_d,
            s.eps,
            s.chi,
            _V_INIT,
            vprimes,
            energies,
            ingredients=ing,
            n_states=n_states,
            include_background=include_background,
        ),
        dtype=np.float64,
    )


def discrete_states(s: VeSetup) -> tuple[DiscreteState, DiscreteState]:
    """`(choice A, choice B)` for `s`.

    `R_inf` is a NUCLEAR coordinate, so it is `nuc.R0` (N2 12.0, F2 10.7) and
    never `elec.R0` (16.0, the electronic grid's ECS radius). It must also be
    the outermost node the ingredients are built on, or `phi_d` is not an
    eigenvector of `H_el` there and `nonlocal_operator`'s tail-coupling guard
    rejects the set.
    """
    ds_a = PhysicalDiscreteState(s.elec, s.model, s.R_desc, s.elec_b)
    ds_b = AsymptoticDiscreteState(s.elec, s.model, R_inf=s.nuc.R0)
    return ds_a, ds_b


def compare(
    molecule: str,
    energies: npt.ArrayLike,
    vprimes: list[int],
    *,
    n_states: int | None,
) -> VeComparison:
    """Run all six `sigma_{0->v\'}(E)` curves for `molecule` on its own deck.

    Parameters
    ----------
    molecule : str
        `"N2"` or `"F2"`.
    energies : array_like
        Incident electron kinetic energies (hartree), inside the molecule's
        published window (see the module docstring).
    vprimes : list of int
        Final vibrational levels.
    n_states : int, optional
        The Eq. (60) state-sum truncation, measured per molecule by the
        ladders recorded in `test_ve_nrm.py`. `None` uses every state.

    Returns
    -------
    VeComparison
    """
    s = setup(molecule)
    e = np.asarray(energies, dtype=np.float64)

    sigma_exact = np.asarray(
        exact_ve_cross_section(s.tgrid, s.model, s.eps, s.chi, _V_INIT, vprimes, e),
        dtype=np.float64,
    )

    vd, gamma = local_complex_potential(s.model, s.nuc, s.elec, s.elec_b)
    sigma_lcp = np.asarray(
        lcp_ve_cross_section(s.nuc, s.model.mu, vd, gamma, s.eps, s.chi, _V_INIT, vprimes, e),
        dtype=np.float64,
    )

    ds_a, ds_b = discrete_states(s)
    curves: dict[str, npt.NDArray[np.float64]] = {}
    for label, phi_d in (("a", ds_a), ("b", ds_b)):
        ing = nrm_ingredients(s.elec, s.model, phi_d, s.R_desc)
        for suffix, bg in (("", True), ("_nobg", False)):
            curves[label + suffix] = nrm_sigma(
                s, phi_d, e, vprimes, n_states=n_states, include_background=bg, ingredients=ing
            )

    return VeComparison(
        molecule=molecule,
        vprimes=list(vprimes),
        energies=e,
        sigma_exact=sigma_exact,
        sigma_lcp=sigma_lcp,
        sigma_nrm_a=curves["a"],
        sigma_nrm_a_nobg=curves["a_nobg"],
        sigma_nrm_b=curves["b"],
        sigma_nrm_b_nobg=curves["b_nobg"],
    )
