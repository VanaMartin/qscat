"""NRM vs LCP vs exact-2D dissociative attachment for F2 and NO.

The research question spec 1 exists to answer: how much of the LCP's
documented error does nonlocality alone buy back, and where does the nonlocal
model itself break down? The exact 2-D solver is the oracle; both the LCP and
the NRM are approximations under test.

Both discrete-state choices of PRA 77 Sec. VI are run:
  A -- `PhysicalDiscreteState`, the scattering function at Re E_res(R)
  B -- `AsymptoticDiscreteState`, the R-independent bound state
The paper predicts B near-exact for DA and A degraded by a Born-Oppenheimer
breakdown. Measured here, that holds on F2: B lands within 1.9% of the exact
2-D oracle at the anchor nearest threshold and within 0.06-0.33% at the other
four, while A falls to 0.29 of exact.

PRA 77 PUBLISHES A DA CROSS SECTION FOR F2 ONLY. NO's DA channel opens at
+0.1719 Ha, above every NO window the paper plots (0.01-0.08 Ha), and N2's at
+0.5016 Ha -- both are energetically shut throughout the published data, so
Sec. VI B's "gives exact results" for DA rests on the single F2 panel. NO is
therefore an UNTESTED regime, not a counter-example: what the four routes do
there is recorded by `test_nrm.py` as an observation, and choice B's absolute
scale is unexplained.

ENERGIES ARE PER MOLECULE. The DA channel opens at `E = eps_e - eps[v_init]`,
measured here at **-0.069 Ha for F2** (open at every positive `E`) but at
**+0.1719 Ha for NO**. Below its threshold every route returns exactly `0.0`
and every ratio is `0/0`, so a shared energy grid cannot compare both
molecules; `validation/diatomic/test_nrm.py` carries the per-molecule anchors.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.core.dissociation import da_cross_section
from qscat.core.grids import electronic_grid
from qscat.core.lcp import lcp_da_cross_section, local_complex_potential
from qscat.core.nrm import (
    AsymptoticDiscreteState,
    DiscreteState,
    NrmIngredients,
    PhysicalDiscreteState,
    nrm_da_cross_section,
    nrm_ingredients,
)
from qscat.core.vibrational import vibrational_states
from qscat.dvr import FemDvrEcsGrid, TensorGrid
from qscat.model import ResonanceModel

from .config import CONFIGS

__all__ = ["NrmComparison", "NrmSetup", "compare", "nrm_sigma", "setup"]

# The second electronic grid (a different ECS angle) both the LCP pole walk and
# `PhysicalDiscreteState` need for two-angle pole stability. 40 deg against the
# decks' own 30 deg, matching `qscat.core.lcp`'s own usage.
_ANGLE_B_DEG = 40.0

# Initial vibrational level for every route.
_V_INIT = 0


@dataclass(frozen=True)
class NrmSetup:
    """The per-molecule inputs every route shares.

    Built once by `setup` so an energy sweep, a state-sum ladder and the two
    discrete-state choices all run against exactly the same grids and
    vibrational basis -- the comparison is only meaningful if they do.

    Attributes
    ----------
    molecule : str
        The `CONFIGS` key.
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
class NrmComparison:
    """One molecule's sigma_DA(E) by all four routes."""

    molecule: str
    energies: npt.NDArray[np.float64]
    sigma_exact: npt.NDArray[np.float64]
    sigma_lcp: npt.NDArray[np.float64]
    sigma_nrm_a: npt.NDArray[np.float64]
    sigma_nrm_b: npt.NDArray[np.float64]


def setup(molecule: str, *, e_r_max: float | None = None) -> NrmSetup:
    """Build the shared grids and vibrational basis for `molecule`.

    Parameters
    ----------
    molecule : str
        A `validation.diatomic.config.CONFIGS` key (`"F2"` or `"NO"`).
    e_r_max : float, optional
        Electronic real-region extent (bohr), overriding the deck's own
        `cfg.e_r_max`. BOTH electronic grids move together -- the pole walk
        needs its two ECS angles on grids that differ only in angle, so
        overriding one and not the other would silently compare different
        discretisations. Used by `da_figure.py` to measure how far the LCP's
        pole walk moves with the box; leave it `None` for the shipped deck.

    Returns
    -------
    NrmSetup
    """
    cfg = CONFIGS[molecule]
    r_max = cfg.e_r_max if e_r_max is None else e_r_max
    tgrid = (
        cfg.da_grid()
        if e_r_max is None
        else TensorGrid(
            [
                electronic_grid(r_max=r_max, order=cfg.e_order, n_complex=cfg.e_n_complex),
                cfg.da_grid().grids[1],
            ]
        )
    )
    elec, nuc = tgrid.grids
    elec_b = electronic_grid(
        r_max=r_max,
        order=cfg.e_order,
        n_complex=cfg.e_n_complex,
        angle_deg=_ANGLE_B_DEG,
    )
    eps, chi = vibrational_states(nuc, cfg.model.mu, cfg.n_vib, cfg.model.v0)
    real = nuc.points.imag == 0.0
    R_desc = np.sort(nuc.points[real].real)[::-1]
    return NrmSetup(
        molecule=molecule,
        model=cfg.model,
        tgrid=tgrid,
        elec=elec,
        nuc=nuc,
        elec_b=elec_b,
        eps=eps,
        chi=chi,
        R_desc=R_desc,
    )


def nrm_sigma(
    s: NrmSetup,
    phi_d: DiscreteState,
    energies: npt.NDArray[np.float64],
    *,
    n_states: int,
    ingredients: NrmIngredients | None = None,
) -> npt.NDArray[np.float64]:
    """`sigma_DA(E)` in the NRM for one discrete-state choice.

    Parameters
    ----------
    s : NrmSetup
        The shared per-molecule setup.
    phi_d : DiscreteState
        Choice A (`PhysicalDiscreteState`) or B (`AsymptoticDiscreteState`).
    energies : ndarray
        Incident electron kinetic energies (hartree).
    n_states : int
        Eq. (60) state-sum truncation.
    ingredients : NrmIngredients, optional
        Prebuilt ingredients for `phi_d`; built here if omitted. They are the
        expensive, energy-INDEPENDENT half, so pass them in when sweeping
        `n_states`.

    Returns
    -------
    ndarray
        `sigma_DA` per energy, bohr^2.
    """
    ing = (
        ingredients
        if ingredients is not None
        else nrm_ingredients(s.elec, s.model, phi_d, s.R_desc)
    )
    return np.asarray(
        nrm_da_cross_section(
            s.nuc,
            s.elec,
            s.model,
            phi_d,
            s.eps,
            s.chi,
            _V_INIT,
            energies,
            ingredients=ing,
            n_states=n_states,
        ),
        dtype=np.float64,
    )


def compare(
    molecule: str,
    energies: npt.ArrayLike,
    *,
    n_states: int,
) -> NrmComparison:
    """Run all four sigma_DA routes for `molecule` on the per-molecule deck.

    Parameters
    ----------
    molecule : str
        A `CONFIGS` key (`"F2"` or `"NO"`).
    energies : array_like
        Incident electron kinetic energies (hartree). Must lie above the
        molecule's own DA threshold -- see the module docstring.
    n_states : int
        The Eq. (60) state-sum truncation, measured per molecule and per
        discrete-state choice by the ladders in `test_nrm.py`'s docstrings.

    Returns
    -------
    NrmComparison
    """
    s = setup(molecule)
    e = np.asarray(energies, dtype=np.float64)

    sigma_exact = np.asarray(
        da_cross_section(s.tgrid, s.model, s.eps, s.chi, _V_INIT, e), dtype=np.float64
    )[:, 0]

    vd, gamma = local_complex_potential(s.model, s.nuc, s.elec, s.elec_b)
    sigma_lcp = np.asarray(
        lcp_da_cross_section(s.nuc, s.model.mu, vd, gamma, s.eps, s.chi, _V_INIT, e),
        dtype=np.float64,
    )

    # R_inf is a NUCLEAR coordinate, so it is `nuc.R0` (F2 10.7, NO 9.0) and
    # never `elec.R0` (16.0, the electronic grid's ECS radius). It must also be
    # the outermost node the ingredients are built on, or `phi_d` is not an
    # eigenvector of `H_el` there and `nonlocal_operator`'s tail-coupling guard
    # rejects the set.
    ds_b = AsymptoticDiscreteState(s.elec, s.model, R_inf=s.nuc.R0)
    sigma_b = nrm_sigma(s, ds_b, e, n_states=n_states)

    ds_a = PhysicalDiscreteState(s.elec, s.model, s.R_desc, s.elec_b)
    sigma_a = nrm_sigma(s, ds_a, e, n_states=n_states)

    return NrmComparison(
        molecule=molecule,
        energies=e,
        sigma_exact=sigma_exact,
        sigma_lcp=sigma_lcp,
        sigma_nrm_a=sigma_a,
        sigma_nrm_b=sigma_b,
    )
