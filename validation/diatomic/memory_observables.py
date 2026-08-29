"""The three-molecule memory-observable campaign of the time-dependent NRM.

    uv run --no-sync python -m validation.diatomic.memory_observables N2
    uv run --no-sync python -m validation.diatomic.memory_observables F2 --outdir <dir>
    uv run --no-sync python -m validation.diatomic.memory_observables NO --converge
    uv run --no-sync python -m validation.diatomic.memory_observables figure

`qscat.core.nrm.memory` records what the time-dependent nonlocal resonance
model is DOING while it runs: where the amplitude goes (`arm_norm`), the rate
at which the coupling feeds the discrete state (`exchange`, against its
Markovian limit `exchange_local`), and the residual that says how much of the
first of those may be read as a transfer (`imbalance`). This module runs those
observables on N2, F2 and NO and writes one `.npz` per molecule plus the
three-column figure assembled from them.

THE QUESTION IT EXISTS TO ANSWER. In the ENERGY domain the three molecules'
local-complex-potential failures are ordered: N2's LCP is mild, F2's ratio to
the exact oracle sweeps 0.263 -> 1.736 across 0.010-0.050 Ha, and NO's is
UNDETERMINED because its pole walk does not converge in the electronic box
(3.98e4 spread, `da_figure.py`). Does the return flux -- the thing the LCP
cannot represent at all -- reproduce that ordering?

THE ANSWER IS YES, ON THE INTEGRAL AND NOT ON THE FLUX, AND IT IS CONDITIONAL.
The returning flux itself is resolved on NO alone (`COARSE_GRAINED_RETURN`), so
it cannot carry a three-way comparison. `nonlocality` can: it converges under
refinement on all three and orders them N2 < NO < F2 without overlap -- but
only on the rungs where it measures the kernel rather than the launch
transient. Near a threshold the Markovian REFERENCE in its denominator
collapses and the ratio inflates for reasons unrelated to the kernel.
`ENERGY_LADDER` carries the criterion, the thirteen rungs it was measured on,
and the record of this module having retracted the ordering once on a ladder
that included the invalid ones.

NORMALIZATION IS PART OF EVERY NUMBER. `exchange` ships UNNORMALIZED, and the
two normalizations in use differ by orders: on the N2 gate deck a raw positive
maximum of +8.776e-7 is +2.420e-4 divided by `S_d(0)` and +1.64e-3 divided by
`S_d(t)`. Every summary field below carries its normalization in its NAME, and
every axis label in `write_figure` carries it in the label. The two exchange
curves are compared as a DIFFERENCE or with BOTH divided by `S_d`; never by
each other, and never anything divided by `Gamma_loc`, which decays to ~1e-10
outside the open-channel window rather than vanishing (`local_width`'s
docstring).

`arm_norm` IS NOT A POPULATION. Under ECS `H_ext` is complex symmetric, the
conjugating norm is conserved by nothing, and Task 1 measured the coupling's
two one-sided rates disagreeing by median 0.822 of the larger. `arm_norm` and
`arm_norm_by_channel` are a RELATIVE CHANNEL DECOMPOSITION -- read across
channels and against themselves over time -- and the partition panel's two
curves do not sum to anything conserved. `imbalance_median_rel` reports that
residual per molecule rather than as a footnote.

EVERY PROPAGATION USES THE COMPLETE ARM SET (`n_states=None`). A truncated arm
set leaves `H_ext` with growing eigenmodes; the arms ARE the dissipation
(`docs/superpowers/specs/2026-08-26-nrm-memory-observables-design.md` Sec. 4).

O2 IS DELIBERATELY ABSENT even though `qscat.model.library` registers `O2`,
`O2_SO12` and `O2_SO32`: those parameters are the potential factory's fit, not
a published set, so there is no external LCP comparison for them and the
comparative question above is undefined.

`validation/` may import `qscat` and `projects`; the reverse is forbidden.
"""

from __future__ import annotations

import argparse
import resource
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

import numpy as np
import numpy.typing as npt
from qscat.core.grids import electronic_grid, segmented_grid
from qscat.core.nrm.coupling import gamma_from_coupling, v_dk_plus
from qscat.core.nrm.discrete_state import AsymptoticDiscreteState
from qscat.core.nrm.extended import extended_hamiltonian, initial_packet
from qscat.core.nrm.ingredients import NrmIngredients, nrm_ingredients
from qscat.core.nrm.memory import MemorySpec, local_width
from qscat.core.nrm.nonlocal_potential import nonlocal_operator
from qscat.core.nrm.propagation import propagate_nrm
from qscat.core.vibrational import vibrational_states
from qscat.dvr import FemDvrEcsGrid
from qscat.model import F2, N2, NO, ResonanceModel

__all__ = ["DECKS", "CampaignDeck", "Summary", "run", "summarize", "write_figure"]

_FIGURE_DIR: Final[Path] = Path(__file__).resolve().parents[2] / "docs" / "physics" / "figures"

#: Initial vibrational level, every molecule. The launch packet is `chi_v0`
#: times the incident electron's amplitude; `nrm_ingredients`' choice-B
#: discrete state is the R-INDEPENDENT bound state (`AsymptoticDiscreteState`),
#: the one PRA 77 Sec. VI predicts near-exact and `validation/diatomic/nrm.py`
#: measures at 0.06-0.33% of the exact oracle on F2.
_V_INIT: Final[int] = 0

#: Neutral vibrational levels resolved on the nuclear grid. Four is what every
#: other diatomic driver here uses; only `eps[_V_INIT]` is read.
_N_VIB: Final[int] = 4


@dataclass(frozen=True)
class CampaignDeck:
    """One molecule's grids, energy and propagation window.

    The nuclear decks are the shipped ones -- N2's from the propagation gate,
    F2's and NO's the eMoScat per-molecule decks of
    `validation/diatomic/config.py`. What is NOT inherited is the ELECTRONIC
    box: see `e_r_max` below.

    Attributes
    ----------
    name : str
        Molecule label, and the `.npz` stem.
    model : ResonanceModel
        The molecule's model.
    nuc_real, nuc_complex : tuple
        `segmented_grid` `(n_elements, endpoint)` pairs for the nuclear grid.
    nuc_angle : float
        Nuclear ECS angle, degrees.
    nuc_quad : int
        Nuclear DVR quadrature order.
    e_r_max : float
        Electronic real-region extent, bohr. **Chosen by the `--converge`
        ladder, not inherited**: `local_width` is the Markovian reference every
        observable here is read against, and it is not converged on the boxes
        the existing gate decks use. See `E_BOX_LADDER` for the measurement.
    e_order, e_n_complex : int
        Electronic DVR order and number of ECS-tail elements.
    e_kin : float
        Incident electron kinetic energy, hartree. Per molecule and stated with
        every number: `F(E)` carries one total energy, so `Gamma_loc` -- and
        therefore `exchange_local` -- is the width at THIS energy's local
        electron energy `eps_loc(R) = E_tot - v0(R)` and no other. One
        propagation per energy for that reason; a multi-energy `LaunchBasis`
        would share a single `MemorySpec.gamma_local` across columns that do
        not share an energy.
    dt, n_steps : float, int
        Propagation window. `dt = 1` is `td_nrm_figures.py`'s setting; the
        order-3 Pade error falls as `dt^6`.
    order : int
        Diagonal-Pade order of the propagator. Three everywhere; `--order 4`
        exists so the RESOLUTION of the exchange rate can be measured without
        changing the number of solves, which is a different experiment from
        halving `dt` and discriminates truncation error from accumulated
        round-off.
    """

    name: str
    model: ResonanceModel
    nuc_real: tuple[tuple[int, float], ...]
    nuc_complex: tuple[tuple[int, float], ...]
    nuc_angle: float
    nuc_quad: int
    e_r_max: float
    e_order: int
    e_n_complex: int
    e_kin: float
    dt: float
    n_steps: int
    order: int = 3

    def nuclear(self) -> FemDvrEcsGrid:
        """The nuclear FEM-DVR-ECS grid."""
        return segmented_grid(
            self.nuc_real, self.nuc_complex, angle_deg=self.nuc_angle, quadrature=self.nuc_quad
        )

    def electronic(self, r_max: float | None = None) -> FemDvrEcsGrid:
        """The electronic FEM-DVR-ECS grid, or a ladder rung of it."""
        return electronic_grid(
            r_max=self.e_r_max if r_max is None else r_max,
            order=self.e_order,
            n_complex=self.e_n_complex,
        )


#: The campaign decks.
#:
#: ELECTRONIC BOXES ARE MEASURED, NOT INHERITED, and none of the three is the
#: box its molecule's existing gate deck uses. `local_width` is 22% high on the
#: N2 propagation gate deck (`r_max = 11`, `n_complex = 3`), which was sized for
#: an algebraic transform identity that cannot see an under-converged electronic
#: box at all; `td_nrm_figures.py`'s F2 decks use `r_max = 13, order 5,
#: n_complex 2`, which reads 1.54x the Eq. (68) width. See `E_BOX_LADDER`.
#:
#: ENERGIES. N2 at 0.10 Ha is the propagation gate's own and sits in the
#: vibrational-excitation window (`docs/physics/nonlocal-resonance-model.md`
#: gates N2 VE over 0.06-0.16 Ha); N2's DA channel opens only at +0.5016 Ha.
#: F2 at 0.030 Ha sits at the LCP/exact DA crossing (the ratio passes through
#: unity near 0.032), so the two flanking energies 0.010 and 0.050 -- where the
#: LCP reads 0.263 and 1.736 -- are the natural extra rungs and `--energy`
#: takes them. NO at 0.200 Ha is inside `da_figure.py`'s 0.150-0.300 grid and
#: above NO's +0.1719 Ha DA threshold, so the dissociation channel is open.
DECKS: Final[dict[str, CampaignDeck]] = {
    "N2": CampaignDeck(
        name="N2",
        model=N2,
        # The nuclear half of `test_nrm_propagation.py`'s gate deck, verbatim.
        nuc_real=((3, 1.5), (8, 3.0), (2, 4.0), (4, 8.0)),
        nuc_complex=((3, 20.0),),
        nuc_angle=35.0,
        nuc_quad=10,
        e_r_max=16.0,
        e_order=6,
        e_n_complex=4,
        e_kin=0.10,
        dt=1.0,
        n_steps=4000,
    ),
    "F2": CampaignDeck(
        name="F2",
        model=F2,
        # `validation/diatomic/config.py`'s eMoScat F2 nuclear deck, verbatim.
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
        e_r_max=20.0,
        e_order=6,
        e_n_complex=4,
        e_kin=0.030,
        dt=1.0,
        n_steps=4000,
    ),
    "NO": CampaignDeck(
        name="NO",
        model=NO,
        # `validation/diatomic/config.py`'s eMoScat NO nuclear deck, verbatim.
        nuc_real=((1, 1.0), (1, 1.6), (37, 9.0)),
        nuc_complex=((1, 9.25), (1, 10.0), (1, 12.0), (4, 42.0)),
        nuc_angle=45.0,
        nuc_quad=14,
        e_r_max=16.0,
        e_order=6,
        e_n_complex=4,
        e_kin=0.200,
        dt=1.0,
        n_steps=4000,
    ),
}

#: MEASURED 2026-08-27 by `--converge` (refresh with
#: `python -m validation.diatomic.memory_observables <mol> --converge`).
#: Per rung: `(r_max, order, n_complex) -> (peak Gamma_loc, median
#: Gamma_loc/Gamma_Eq68 over the open window, min `_sign_align` overlap)`.
#:
#: READ THE FIRST COLUMN, NOT THE SECOND. The median ratio against Eq. (68)
#: compares two quantities built on the SAME electronic grid, so both move
#: together and it can read 1.000 on a box that is 15% from converged: N2 at
#: `r_max = 11, n_complex = 4` scores 1.0004 while its peak `Gamma_loc` is
#: 2.890e-2 against a converged 2.504e-2. The criterion used to pick each
#: deck's box is therefore the box-to-box change in `Gamma_loc` ITSELF, with
#: the ratio as a sanity check and `min_overlap` as the veto -- below 0.5,
#: `nrm_ingredients` warns that its adiabatic tracking has paired the wrong
#: P-space state between adjacent nuclear nodes, which would corrupt every
#: ingredient downstream.
#:
#: N2: peak `Gamma_loc` settles at 2.504e-2 from `r_max = 14` on, but 14 warns
#: (`min_overlap` 0.374) and 20/24 track worse still (0.216 / 0.015). 16 is the
#: one converged rung with clean tracking, and it is what the deck uses.
#: F2: 13:5:2 (the shipped TD figure deck) is 1.10x the converged peak and
#: reads 1.54x Eq. (68); 16:6:4 is 4% off pointwise; 20:6:4 agrees with the
#: 128568-square 16:8:6 production box to 1.9e-3 at the same 81816 size.
#: NO: converged from `r_max = 14` on -- peak `Gamma_loc` stable to four
#: figures over 14-20 and both orders, `min_overlap` >= 0.998 throughout.
E_BOX_LADDER: Final[dict[str, dict[tuple[float, int, int], tuple[float, float, float]]]] = {
    "N2": {
        (11.0, 6, 3): (3.828e-2, 1.2181, 0.973),
        (11.0, 6, 4): (2.890e-2, 1.0004, 0.996),
        (14.0, 6, 4): (2.503e-2, 0.9975, 0.374),
        (16.0, 6, 4): (2.509e-2, 0.9958, 0.820),
        (20.0, 6, 4): (2.504e-2, 0.9966, 0.216),
        (24.0, 6, 4): (2.504e-2, 0.9963, 0.015),
        (16.0, 8, 6): (2.504e-2, 0.9968, 0.022),
    },
    "F2": {
        (13.0, 5, 2): (8.442e-3, 1.5447, 1.000),
        (11.0, 6, 4): (4.313e-3, 0.8713, 1.000),
        (14.0, 6, 4): (7.691e-3, 1.0276, 0.993),
        (16.0, 6, 4): (7.711e-3, 0.9799, 1.000),
        (20.0, 6, 4): (7.700e-3, 1.0102, 1.000),
        (16.0, 8, 6): (7.698e-3, 1.0155, 0.996),
    },
    "NO": {
        (14.0, 6, 4): (6.9375e-2, 0.9787, 0.998),
        (16.0, 6, 4): (6.9374e-2, 0.9774, 1.000),
        (20.0, 6, 4): (6.9377e-2, 0.9852, 1.000),
        (16.0, 8, 6): (6.9376e-2, 0.9856, 1.000),
        (20.0, 8, 6): (6.9372e-2, 0.9833, 1.000),
    },
}

#: WHAT EACH CONVERGED DECK COSTS, so "run these on a big machine" is a number
#: someone can plan against rather than an instruction.
#: `(N_R, n_arm, H_ext order, H_ext nnz, peak RSS in GB, wall-clock minutes for
#: 4000 order-3 steps)`.
#:
#: MEASURED 2026-08-27 on a 12-core / 68.7 GB laptop with the **SuperLU**
#: fallback (no system MUMPS outside Docker on a Mac), `OMP_NUM_THREADS=4`.
#: Peak RSS is `resource.getrusage(RUSAGE_SELF).ru_maxrss`, printed by `run`.
#: The wall-clock figures are from a run sharing the machine with up to four
#: others, so they are an upper bound; N2 alone took 3.7 min.
#:
#: THE HEADLINE IS THAT THIS FITS ON A LAPTOP. `td_nrm_figures.py` used to
#: state that an `H_ext` of 53570 "DOES NOT FIT ON A LAPTOP" on the strength of
#: three OS kills; F2 here is 81816 and peaks at 5.82 GB. That docstring is
#: corrected -- the kills were concurrency, not size. MUMPS would cut these
#: further (~9x at 143k unknowns) and is worth having, but is not a
#: prerequisite at this size.
DECK_COST: Final[dict[str, tuple[int, int, int, int, float, float]]] = {
    "N2": (179, 83, 15036, 192926, 0.76, 3.7),
    "NO": (597, 83, 50148, 848130, 3.4, 0.0),
    "F2": (974, 83, 81816, 1385732, 5.82, 0.0),
}

#: HOW FAR THE EXCHANGE RATE IS RESOLVED ON N2, AND THE RETRACTION THAT FOLLOWS.
#: MEASURED 2026-08-27, four propagations of the campaign N2 deck (E = 0.10 Ha,
#: T = 4000): order 3 at dt = 1, order 4 at dt = 1 (truncation only, same number
#: of solves), and order 3 at dt = 0.5 and dt = 0.25. Refresh with
#: `--order 4` / `--dt 0.5 --steps 8000` / `--dt 0.25 --steps 16000` into a
#: scratch `--outdir`, then `resolution --against N2` over it.
#:
#: `(label) -> (max positive / max|exchange|, onset t, nonlocality,
#: arm peak / S_d(0))`:
N2_RESOLUTION: Final[dict[str, tuple[float, float, float, float]]] = {
    "order 3, dt=1": (6.402e-6, 1453.00, 0.5070, 0.1364),
    "order 4, dt=1": (3.890e-6, 1472.00, 0.5097, 0.1365),
    "order 3, dt=0.5": (1.192e-4, 236.00, 0.5092, 0.1365),
    "order 3, dt=0.25": (1.649e-4, 246.75, 0.5068, 0.1365),
}

#: The pointwise floor between those runs, `max|d(exchange)| / max|exchange|`
#: over the shared times, with the fraction of steps on which the pair agrees
#: about the SIGN. The finest pair is the binding one.
N2_RESOLUTION_FLOOR: Final[dict[tuple[str, str], tuple[float, float]]] = {
    ("order 3, dt=1", "order 4, dt=1"): (2.372e-2, 0.706),
    ("order 3, dt=1", "order 3, dt=0.5"): (2.210e-2, 0.699),
    ("order 3, dt=1", "order 3, dt=0.25"): (2.425e-2, 0.715),
    ("order 4, dt=1", "order 3, dt=0.5"): (7.382e-3, 0.708),
    ("order 4, dt=1", "order 3, dt=0.25"): (8.096e-3, 0.725),
    ("order 3, dt=0.5", "order 3, dt=0.25"): (4.257e-3, 0.720),
}

#: THE DISCRIMINATOR, measured the same way on every molecule. `RESOLVED_RETURN`
#: maps `(molecule, pair)` to `(positive max / peak |exchange|, floor over the
#: returning window, fraction of the first run's returning steps that the second
#: also calls returning, median |d(exchange)|/|exchange| at those steps)`.
#:
#: THE LAST TWO ARE THE ONES TO READ. "Of the steps this run calls a RETURN, how
#: many does a better-resolved run call a return too?" is normalization-free and
#: tests exactly the claim being made. A resolved rate scores near 1; one at the
#: numerical floor scores near chance.
#:
#: The window matters and the global floor is the wrong one for this purpose:
#: `max|d(exchange)|` over a whole run is normally set in the launch transient,
#: where `|exchange|` peaks and moves fastest. On NO it lands at t = 15 and
#: reports 7.8e-2 for a run whose late structure the two propagators agree on to
#: 8.9e-3, hundreds of a.u. before the first positive step.
#:
#: THE COMPARISON WINDOW IS PART OF THE MEASUREMENT, and every row below is on
#: the FULL campaign window, T = 4000. `compare_resolution` compares whatever
#: two runs share and its docstring licenses a shorter check -- which is right
#: for bounding the floor, and wrong for the concordance, because the
#: concordance DEGRADES with window length as the two propagators drift apart.
#: Measured on NO: 0.952 / 0.945 / 0.943 / 0.908 for the order-4 pair at
#: T = 1000 / 2000 / 3000 / 4000, and 0.901 / 0.872 / 0.852 / 0.833 for the
#: dt pair. This table's NO rows previously carried the T = 1000 values while
#: its N2 rows carried T = 4000 ones, so the two molecules -- the whole point
#: of the comparison -- were being read off different windows. Corrected
#: 2026-08-28 by re-running every pair and reporting T = 4000 throughout.
#:
#: PROVENANCE. Every row re-measured 2026-08-28 on sadaharu (x86, MUMPS) except
#: where the base run is the committed campaign `.npz` (arm, SuperLU). That mix
#: is safe and was checked rather than assumed -- see `CROSS_PLATFORM`.
#: THE CONTROL THAT LICENSES READING `RESOLVED_RETURN`'S FLOORS AS TIME
#: DISCRETISATION. A floor of 1e-2 of the exchange peak is only evidence about
#: the propagator if everything else about the run is smaller than it. Each
#: campaign deck was re-propagated on sadaharu (x86, system MUMPS) and compared
#: with the committed `.npz` from the laptop (arm64, SuperLU) -- a different
#: CPU, a different BLAS and a different sparse factorisation, at identical
#: `dt` and order.
#:
#: `molecule -> (max|d(exchange)| / peak, sign agreement, concordance)`:
CROSS_PLATFORM: Final[dict[str, tuple[float, float, float]]] = {
    "N2": (1.946e-13, 1.000, 1.000),
    "NO": (1.323e-12, 1.000, 1.000),
    "F2": (8.422e-13, 1.000, 1.000),
}

#: EVERY CONCORDANCE HERE CARRIES ITS NULL, for the same reason
#: `coarse_grained_return` does: the fraction is CONDITIONAL, so the null is the
#: comparison run's own positive-step rate and not one half. This table shipped
#: without it one revision longer than the binned one did, and the correction
#: matters in both directions -- N2's 0.427 is a lift of +0.171 over a null of
#: 0.256 (a real signal, not "below chance"), and F2's larger-looking 0.622 is a
#: lift of +0.048 over 0.575. **On this metric N2 OUTRANKS F2 by 3.7x**, the
#: reverse of what the raw column suggests.
#:
#: The verdict is unchanged, because it never rested on N2-vs-F2: NO's lift is
#: +0.513 to +0.568, three to twelve times either of them, and NO is the only
#: molecule whose positive maximum also clears its own floor.
#:
#: `(molecule, pair) -> (positive max / peak |exchange|, floor over the
#: returning window, concordance, ITS NULL, median |d(exchange)|/|exchange| at
#: the returning steps)`.
RESOLVED_RETURN: Final[dict[tuple[str, str], tuple[float, float, float, float, float]]] = {
    ("N2", "order 3 vs order 4, dt=1"): (6.402e-6, 3.853e-5, 0.427, 0.256, 1.592),
    ("N2", "dt=0.5 vs dt=0.25"): (1.192e-4, 1.672e-3, 0.386, 0.219, 1.255),
    ("NO", "order 3 vs order 4, dt=1"): (1.982e-2, 8.898e-3, 0.908, 0.340, 0.314),
    ("NO", "dt=1 vs dt=0.5"): (1.982e-2, 7.694e-3, 0.833, 0.320, 0.482),
    ("F2", "order 3 vs order 4, dt=1"): (1.590e-1, 2.447e-1, 0.622, 0.575, 1.052),
    ("F2", "dt=0.5 vs dt=0.25"): (1.959e-1, 1.913e-1, 0.594, 0.549, 1.198),
}

#: THE DISCRIMINATOR THAT SURVIVES, on time-averaged bins rather than steps.
#: `(molecule, bin width) -> (min lift, max lift, min magnitude, max magnitude,
#: n returning bins)` over EVERY refinement available for that molecule -- N2
#: three (order 4, dt = 0.5, dt = 0.25), NO two, F2 three. Measured 2026-08-28;
#: `coarse_grained_return` defines all four and says why the raw concordance is
#: not among them.
#:
#: LIFT IS CONCORDANCE MINUS ITS NULL, and the null is the comparison run's own
#: positive-bin rate rather than one half, because the concordance is
#: CONDITIONAL. That correction changes the reading of this table completely and
#: it is why the raw numbers are not recorded here. F2's positive-bin rate is
#: ~0.59, so its raw 0.638-0.750 -- once quoted as "reproduces at about two
#: thirds" -- is a lift of +0.05 to +0.16, essentially chance. NO's rate is
#: ~0.35 against a raw 0.96-1.00, a lift of +0.61 to +0.65.
#:
#: MAGNITUDE IS THE OTHER HALF, because agreeing about a bin's sign says nothing
#: about its size. It separates the same way and more sharply: NO's returning
#: bins agree to 2-19%, N2's differ by 125-148%, F2's by 84-520%.
#:
#: SO THERE ARE TWO BANDS, NOT THREE, and MAGNITUDE is what settles it, since
#: that column involves no choice of metric. NO is resolved on both columns at
#: every bin width. N2 and F2 are NOT ranked against each other: both lifts are
#: small and positive, and which looks larger REVERSES between the pointwise
#: table (N2 +0.17 against F2 +0.05) and this binned one (F2 +0.05..+0.16
#: against N2 0.00..+0.12). An earlier reading put F2 in a middle band, and
#: that was the missing null rather than a finding.
COARSE_GRAINED_RETURN: Final[dict[tuple[str, float], tuple[float, float, float, float, int]]] = {
    ("N2", 5.0): (0.037, 0.122, 1.278, 1.475, 197),
    ("N2", 10.0): (0.000, 0.100, 1.323, 1.419, 96),
    ("N2", 20.0): (0.000, 0.111, 1.272, 1.333, 46),
    ("N2", 50.0): (0.000, 0.104, 1.249, 1.282, 18),
    ("NO", 5.0): (0.613, 0.618, 0.180, 0.191, 272),
    ("NO", 10.0): (0.631, 0.635, 0.145, 0.150, 137),
    ("NO", 20.0): (0.650, 0.650, 0.022, 0.022, 70),
    ("NO", 50.0): (0.650, 0.650, 0.020, 0.033, 28),
    ("F2", 5.0): (0.056, 0.096, 0.903, 1.060, 462),
    ("F2", 10.0): (0.048, 0.063, 0.929, 1.424, 243),
    ("F2", 20.0): (0.065, 0.091, 0.838, 2.120, 113),
    ("F2", 50.0): (0.058, 0.162, 0.911, 5.201, 48),
}

#: THE VERDICT, stated here because it retracts a claim this sub-project made.
#:
#: **On N2 the returning flux is NOT resolved, at any resolution tried, and the
#: design document's Sec. 2.2 headline does not hold for this molecule.** The
#: two finest runs agree to 4.257e-3 of the exchange peak; N2's largest
#: positive excursion is 1.649e-4 of that peak, THIRTY TIMES UNDER THE FLOOR.
#: The sign agreement does not improve with refinement -- 0.699 to 0.725 for
#: every pair including the finest -- so better than a quarter of the steps flip
#: sign between two runs that otherwise agree to four parts in a thousand. The
#: onset moves 1453 -> 236 -> 247 and the positive maximum by a factor of 42.
#:
#: The prototype's `+8.776e-7` is retracted twice over: it is below this floor
#: AND it was measured on the `r_max = 11` box, which `E_BOX_LADDER` shows is
#: 22% wrong in `Gamma_loc` before resolution is considered.
#:
#: WHAT IS UNAFFECTED, and it is most of the campaign. `nonlocality` spans
#: 0.5068-0.5097 (0.6%), the arm-norm peak 0.1364-0.1365, `Gamma_eff` is
#: identical to every digit printed, and the decay-law crossings agree
#: (t = 27.5-28 / 508 / 702, ratios 0.943-0.947 / 1.30e4 / 1.16e5). Those are
#: integrals and they converge. A CONVERGED INTEGRAL IS NOT EVIDENCE FOR AN
#: UNCONVERGED POINTWISE SIGN and is not offered as one. The Markovian side is
#: structural rather than numerical: `-<Psi_d|Gamma_loc|Psi_d>` is non-positive
#: by construction and measured strictly negative on every deck here.
#:
#: THIS IS A STATEMENT ABOUT N2, NOT ABOUT THE OBSERVABLE, and the tables above
#: are the proof rather than the hope. Measured the same way, at four bin widths
#: and against every refinement available, NO's returning bins reproduce at a
#: lift of +0.61 to +0.65 over their own null and agree in magnitude to 2-19%.
#: N2 scores a lift of 0.000-0.122 with magnitudes 125-148% apart. The
#: observable is sound; N2 is a molecule on which it cannot be read.
#:
#: AND SO IS F2, which this table used to place in a middle band. Its lift is
#: +0.05 to +0.16 and its magnitudes are 84-520% apart -- the same verdict as
#: N2, not an intermediate one. The middle band was the raw concordance read
#: without its null, and F2 has the highest positive-bin rate of the three, so
#: it flattered itself the most. ONE molecule of three carries a readable
#: returning flux; the other two are unreadable and are not ordered against
#: each other.
#:
#: WHICH COLUMN THE THREE-MOLECULE ORDERING IS ALLOWED TO REST ON, given all of
#: that. NOT the returning flux, which is readable on NO alone. It rests on
#: `nonlocality`, which converges under REFINEMENT on all three (N2 0.6% over
#: four runs, NO 0.04%, F2 1.9%) and, integrated from the arm-norm peak over
#: each molecule's declared energy window, orders them
#: N2 0.224-0.773 < NO 0.870-0.872 < F2 1.055-1.341 with no overlap.
#: `ENERGY_LADDER` carries the seventeen rungs, the two exclusions, and the
#: record of this module having once over-retracted the ordering.

_N2_RETRACTION = (
    "the returning-flux claim does not hold on N2: 1.649e-4 of peak against a 4.257e-3 floor"
)

#: THE CONFOUND THE ORDERING CLAIM HAD TO SURVIVE, AND HOW IT DOES.
#:
#: The three molecules propagate at three different energies -- each deck's
#: incident energy is set by where its channel is open, not by a shared choice
#: -- so a comparison across them could be a comparison in energy instead. All
#: three are laddered here, seventeen rungs. Refresh with `<MOL> --energy <E>`.
#: MEASURED 2026-08-28.
#:
#: `(molecule, E_kin) -> (nonlocality over the FULL run, nonlocality from the
#: arm-norm peak onwards, int|X_loc| -- the Markovian reference, max positive /
#: peak)`. The SECOND column is the one the ordering is read from; see
#: `nonlocality_post_peak` for why the first is the wrong window.
ENERGY_LADDER: Final[dict[tuple[str, float], tuple[float, float, float, float]]] = {
    ("N2", 0.050): (1.6334, 1.0380, 3.8472e-4, 6.3966e-6),
    ("N2", 0.060): (0.8479, 0.4156, 7.8094e-4, 6.3980e-6),
    ("N2", 0.080): (0.4020, 0.2241, 2.2780e-3, 6.4005e-6),
    ("N2", 0.100): (0.5070, 0.4872, 4.9992e-3, 6.4025e-6),
    ("N2", 0.120): (0.6036, 0.6320, 9.0746e-3, 6.4041e-6),
    ("N2", 0.150): (0.6960, 0.7484, 1.7860e-2, 6.4062e-6),
    ("N2", 0.160): (0.7188, 0.7731, 2.1461e-2, 6.4067e-6),
    ("NO", 0.175): (0.8120, 0.8709, 1.9416e-2, 1.9826e-2),
    ("NO", 0.185): (0.8128, 0.8713, 1.9614e-2, 1.9825e-2),
    ("NO", 0.200): (0.8134, 0.8717, 1.9833e-2, 1.9822e-2),
    ("NO", 0.300): (0.8113, 0.8700, 1.9734e-2, 1.9735e-2),
    ("NO", 0.400): (0.8052, 0.8655, 1.8692e-2, 1.9645e-2),
    ("F2", 0.010): (1.5092, 1.3412, 1.8555e-5, 1.5859e-1),
    ("F2", 0.020): (1.0263, 1.1400, 1.0459e-4, 1.5888e-1),
    ("F2", 0.030): (0.9463, 1.0905, 2.4989e-4, 1.5901e-1),
    ("F2", 0.040): (0.9259, 1.0675, 4.4177e-4, 1.5907e-1),
    ("F2", 0.050): (0.9246, 1.0552, 6.5549e-4, 1.5909e-1),
}

#: THE ONLY RUNGS EXCLUDED, and by a criterion that predates this ladder.
#: `test_the_deck_energies_are_where_the_lcp_comparison_lives` has declared each
#: molecule's window since before any of this work. Two rungs fall outside --
#: N2 at 0.050 and NO at 0.400, both added BY this ladder -- and both are
#: dropped. Nothing else is.
#:
#: THIS MODULE PREVIOUSLY EXCLUDED FOUR RUNGS by a share cut invented after the
#: ladder was run, to keep the launch transient out of the integral. Right
#: diagnosis, wrong remedy: the contamination is a WINDOW, so
#: `nonlocality_post_peak` removes it from every rung instead of removing
#: rungs. That cut is gone, and with it the objection that a criterion was
#: fitted to the answer. What remains is one pre-existing energy window per
#: molecule, applied to both offending rungs rather than to one of them.
ENERGY_WINDOWS: Final[dict[str, tuple[float, float]]] = {
    "N2": (0.06, 0.16),
    "F2": (0.010, 0.050),
    "NO": (0.1719, 0.300),
}

#: WHAT THE LADDER SAYS.
#:
#: **Everything the campaign reads as a RETURN is frozen in energy**, on all
#: three molecules and at every rung: over a 4-6x change in `Gamma_eff` the
#: onset does not move at all and `max positive / peak` moves in the fifth
#: figure. Those columns describe the molecule, not the collision energy.
#:
#: **`nonlocality` orders the three, on the fifteen in-window rungs:**
#:
#:   N2 0.224-0.773  <  NO 0.870-0.872  <  F2 1.055-1.341
#:                  +12.5%            +21.0%
#:
#: `NO < F2` is the strongest and needs no window argument at all: it holds on
#: the RAW full-run column too, over all seventeen rungs, NO's maximum 0.8134
#: against F2's minimum 0.9246. `N2 < F2` holds with a wide margin. `N2 < NO`
#: is the narrowest at 12.5% and is the one to re-examine first -- it is also
#: the one the discarded share-cut construction could not support, since N2's
#: in-window 0.06 rung read 0.848 there, ABOVE NO, where post-peak it reads
#: 0.416.
#:
#: NEAR-THRESHOLD BEHAVIOUR IS STILL PRESENT, just no longer in the way: N2's
#: post-peak values are non-monotone (0.416 at 0.06, dipping to 0.224 at 0.08,
#: rising to 0.773 at 0.16) and F2's fall with energy. The ordering is a
#: statement about RANGES over each declared window, not about trends.
#:
#: WHAT COLLAPSES NEAR THRESHOLD, corrected. An earlier version of this comment
#: said the open window shrinks. It barely does -- the nodes carrying
#: `Gamma_loc` go 89 -> 95 of 153 across N2's ladder and 291 -> 360 of 819
#: across F2's. What collapses is `Gamma_loc`'s MAGNITUDE over the doorway:
#: `max Gamma_loc` moves 5.5x on N2 (8.2e-3 -> 4.5e-2) and 11.8x on F2
#: (1.4e-3 -> 1.7e-2). The reference dies in size, not in extent.

#: `S_d(t)/S_d(0)` levels the decay law is read at (design Sec. 2.3).
_DECAY_LEVELS: Final[tuple[float, ...]] = (0.5, 0.1, 0.01)


@dataclass(frozen=True)
class Summary:
    """The campaign numbers for one molecule, all normalizations in the names.

    Attributes
    ----------
    molecule : str
        `DECKS` key.
    e_kin, dt, t_max : float
        Incident kinetic energy (hartree) and propagation window (a.u.).
    n_r, n_arm, n_ext : int
        Nuclear grid size, number of arms, and `H_ext`'s order.
    min_overlap : float
        `nrm_ingredients`' adiabatic-tracking diagnostic. Below 0.5 the
        ingredients are suspect and every number here with them.
    s_d0 : float
        `S_d(0)` -- the normalization the `_over_s0` fields divide by.
    gamma_eff : float
        Golden-rule rate `<Psi_d(0)|Gamma_loc|Psi_d(0)> / ||Psi_d(0)||^2`, a.u.
        This is `-exchange_local[0] / survival[0]`, i.e. the SAME `Gamma_loc`
        the Markovian curve uses, not a pole walk's.
    arm_peak_over_s0 : float
        Maximum of `sum_n ||phi_n||^2` over the run, divided by `S_d(0)`. A
        RELATIVE CHANNEL DECOMPOSITION, not a population.
    arm_peak_time : float
        When that maximum occurs.
    arm_top_blocks : tuple of int
        The four arm blocks with the largest running maximum, best first.
    arm_first_four_share : float
        Fraction of the total arm norm AT ITS PEAK carried by blocks 0-3 --
        the first four in `h_ext`'s own block order, which is what
        `MemorySpec.n_channels=4` would keep. On the N2 gate deck the first
        four were also the four largest (93.2%); this field is what says
        whether that survives on this deck.
    exchange_max_raw, exchange_min_raw : float
        Extremes of the UNNORMALIZED `exchange`, a.u.
    exchange_max_over_s0, exchange_min_over_s0 : float
        The same two divided by `S_d(0)`.
    exchange_max_over_s : float
        Maximum of `exchange(t) / S_d(t)` -- the FRACTIONAL return rate, which
        is the one that is comparable across molecules whose `S_d` differ by
        orders at the same `t`.
    n_positive, n_steps : int
        Steps at which `exchange > 0` (amplitude returning from the continuum),
        out of `n_steps + 1`. Step 0 is exactly 0 by construction -- the arms
        start empty -- and is not counted.
    t_first_positive : float
        Time of the first such step; `nan` if there is none.
    return_share : float
        `int max(exchange, 0) dt / int |exchange| dt`, dimensionless and
        normalization-free. Dominated by the launch transient's large negative
        lobe, so read it beside `n_positive`.
    nonlocality : float
        `int |exchange - exchange_local| dt / int |exchange_local| dt` -- how
        far the nonlocal rate runs from its own Markovian limit, as a
        DIFFERENCE normalized by the Markovian integral. Never a pointwise
        ratio: `Gamma_loc` is ~1e-10 outside the open-channel window.
    decay_times : tuple of float
        First `t` at which `S_d/S_d(0)` reaches 0.5, 0.1, 0.01; `nan` if not
        reached inside the window.
    decay_ratio_at_levels : tuple of float
        `S_d(t)/S_d(0)` divided by `exp(-gamma_eff t)` at those same times --
        how many times slower (or faster) than one rate constant.
    imbalance_median_rel, imbalance_max_rel : float
        Task 1's residual as a fraction of the LARGER of the two one-sided
        coupling rates, over the steps where that larger rate is nonzero. O(1)
        is the expected answer and is why `arm_norm` is not a population.
    survival_rebound : float
        `S_d(t_max) / min_t S_d(t)`. Above 1 the packet grew after its minimum;
        `propagate_nrm` warns on it. At the 1e-10 level it is the arm floor
        feeding back rather than a growing eigenmode, but it bounds how far
        into the tail the run may be read.
    """

    molecule: str
    e_kin: float
    dt: float
    t_max: float
    n_r: int
    n_arm: int
    n_ext: int
    min_overlap: float
    s_d0: float
    gamma_eff: float
    arm_peak_over_s0: float
    arm_peak_time: float
    arm_top_blocks: tuple[int, ...]
    arm_first_four_share: float
    exchange_max_raw: float
    exchange_min_raw: float
    exchange_max_over_s0: float
    exchange_min_over_s0: float
    exchange_max_over_s: float
    n_positive: int
    n_steps: int
    t_first_positive: float
    return_share: float
    nonlocality: float
    decay_times: tuple[float, ...]
    decay_ratio_at_levels: tuple[float, ...]
    imbalance_median_rel: float
    imbalance_max_rel: float
    survival_rebound: float


def _build(
    deck: CampaignDeck, *, r_max: float | None = None
) -> tuple[
    FemDvrEcsGrid,
    FemDvrEcsGrid,
    AsymptoticDiscreteState,
    NrmIngredients,
    npt.NDArray[np.float64],
    npt.NDArray[np.complex128],
]:
    """Grids, choice-B discrete state, ingredients, vibrational basis.

    Returns `(nuc, elec, phi_d, ing, eps, chi)`. `R_inf` is a NUCLEAR
    coordinate (`nuc.R0`); `nrm_ingredients` requires strictly DESCENDING `R`.
    """
    nuc = deck.nuclear()
    elec = deck.electronic(r_max)
    phi_d = AsymptoticDiscreteState(elec, deck.model, R_inf=nuc.R0)
    r_desc = np.sort(nuc.points[nuc.points.imag == 0.0].real)[::-1]
    ing = nrm_ingredients(elec, deck.model, phi_d, r_desc)
    eps, chi = vibrational_states(nuc, deck.model.mu, _N_VIB, deck.model.v0)
    return nuc, elec, phi_d, ing, eps, chi


def converge(deck: CampaignDeck, ladder: tuple[tuple[float, int, int], ...]) -> None:
    """Print the electronic-box ladder that picks `deck.e_r_max`.

    For each rung it reports the peak `Gamma_loc`, the median against Eq.
    (68)'s independently constructed width over the open-channel window, the
    box-to-box change in `Gamma_loc` itself, and `nrm_ingredients`'
    `min_overlap`. The BOX-TO-BOX CHANGE is the criterion -- the Eq. (68)
    comparison is built on the same electronic grid and moves with it, so it
    can read 1.000 on an unconverged box (see `E_BOX_LADDER`).

    Expensive (one ingredient build plus one `v_dk_plus` solve per open node
    per rung); the results are recorded in `E_BOX_LADDER` so the campaign does
    not re-run it.
    """
    nuc = deck.nuclear()
    eps, _chi = vibrational_states(nuc, deck.model.mu, _N_VIB, deck.model.v0)
    e_total = deck.e_kin + float(eps[_V_INIT])
    real = np.flatnonzero(nuc.points.imag == 0.0)
    r_real = nuc.points[real].real
    # `eps_loc > 0.02` is `test_nrm_memory.py`'s own near-threshold exclusion:
    # inside 0.02 Ha of threshold the width varies by orders across the ~10
    # nodes the kernel spans and no local limit exists to compare to.
    open_ = np.flatnonzero(e_total - np.real(deck.model.v0(r_real)) > 0.02)
    eps_loc = e_total - np.real(deck.model.v0(r_real))
    print(
        f"{deck.name}: N_R={nuc.n} real={real.size} E_kin={deck.e_kin} E_tot={e_total:.5f} "
        f"open window {open_.size} nodes, R in "
        f"[{r_real[open_].min():.4f}, {r_real[open_].max():.4f}]",
        flush=True,
    )
    prev: npt.NDArray[np.float64] | None = None
    for r_max, order, n_cx in ladder:
        elec = electronic_grid(r_max=r_max, order=order, n_complex=n_cx)
        phi_d = AsymptoticDiscreteState(elec, deck.model, R_inf=nuc.R0)
        with warnings.catch_warnings():
            # The tracking warning is the very thing being measured here;
            # `min_overlap` is reported per rung instead of raised per rung.
            warnings.simplefilter("ignore")
            ing = nrm_ingredients(elec, deck.model, phi_d, r_real[::-1])
        gamma = local_width(nonlocal_operator(ing, nuc, deck.model, e_total, n_states=None), nuc)
        eq68 = np.array(
            [
                gamma_from_coupling(
                    v_dk_plus(elec, deck.model, phi_d, np.array([r_real[j]]), float(eps_loc[j]))
                )[0]
                for j in open_
            ]
        )
        g_real = gamma[real]
        line = (
            f"  r_max={r_max:5.1f} order={order} n_complex={n_cx} n_elec={elec.n:3d} "
            f"arms={ing.E_n.shape[1]:3d} H_ext={(1 + ing.E_n.shape[1]) * nuc.n:6d} "
            f"peak Gamma_loc={g_real.max():.6e} median/Eq68={np.median(g_real[open_] / eq68):.4f} "
            f"min_overlap={ing.min_overlap:.3f}"
        )
        if prev is not None:
            change = np.abs(g_real[open_] - prev[open_]).max() / np.abs(prev).max()
            line += f" max|dGamma_loc|/peak vs previous rung={change:.3e}"
        prev = g_real
        print(line, flush=True)


def _stem(deck: CampaignDeck, *, smoke: bool = False) -> str:
    """`.npz` stem for `deck`.

    The campaign deck gets the bare `<mol>-nrm-memory-observables` name the
    figure looks for; anything overridden on the command line gets its
    override in the stem, so an energy ladder rung or a `dt` convergence check
    cannot silently overwrite the campaign run it is supposed to be checked
    against.
    """
    base = DECKS[deck.name]
    stem = f"{deck.name.lower()}-nrm-memory-observables"
    if deck.e_kin != base.e_kin:
        stem += f"-e{deck.e_kin:g}"
    if deck.dt != base.dt:
        stem += f"-dt{deck.dt:g}"
    if deck.order != base.order:
        stem += f"-order{deck.order}"
    if smoke or (deck.n_steps != base.n_steps and deck.dt == base.dt):
        stem += f"-n{20 if smoke else deck.n_steps}"
    return stem


def run(deck: CampaignDeck, *, outdir: Path, smoke: bool = False) -> Path:
    """Propagate `deck` with the memory observables on and write its `.npz`.

    `n_states=None` (the complete arm set) always: a truncated `H_ext` is not
    dissipative and `psi_d` comes back exponentially wrong.
    `MemorySpec.n_channels=None` records every channel -- measured free
    (design Sec. 5) and 2.7 MB at the largest deck here.

    Parameters
    ----------
    deck : CampaignDeck
        The molecule, grids, energy and window.
    outdir : Path
        Directory for `<mol>-nrm-memory-observables.npz`.
    smoke : bool, optional
        Twenty steps instead of `deck.n_steps` -- exercises the path, produces
        no physics.

    Returns
    -------
    Path
        The `.npz` written.
    """
    t0 = time.time()
    nuc, elec, phi_d, ing, eps, chi = _build(deck)
    e_total = deck.e_kin + float(eps[_V_INIT])
    gamma_local = local_width(nonlocal_operator(ing, nuc, deck.model, e_total, n_states=None), nuc)
    h_ext = extended_hamiltonian(ing, nuc, deck.model)
    launch = initial_packet(
        nuc, elec, deck.model, phi_d, ing, eps, chi, _V_INIT, np.array([deck.e_kin]), rank_tol=1e-10
    )
    n_arm = int(ing.E_n.shape[1])
    print(
        f"{deck.name}: N_R={nuc.n} n_elec={elec.n} arms={n_arm} H_ext={h_ext.shape[0]} "
        f"nnz={h_ext.nnz} rank={launch.rank} min_overlap={ing.min_overlap:.3f} "
        f"E_kin={deck.e_kin} E_tot={e_total:.5f} setup {time.time() - t0:.1f}s",
        flush=True,
    )

    n_steps = 20 if smoke else deck.n_steps
    t0 = time.time()
    res = propagate_nrm(
        h_ext,
        launch,
        nuc,
        dt=deck.dt,
        n_steps=n_steps,
        order=deck.order,
        memory=MemorySpec(gamma_local=gamma_local, n_channels=None),
    )
    # `ru_maxrss` is BYTES on macOS and KILOBYTES on Linux -- the one platform
    # difference in this module, and getting it wrong would misreport the
    # number `DECK_COST` exists to carry by a factor of 1024.
    max_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_gb = max_rss / 1e9 if sys.platform == "darwin" else max_rss / 1e6
    print(
        f"{deck.name}: {n_steps} steps in {time.time() - t0:.1f}s, peak RSS {peak_gb:.2f} GB",
        flush=True,
    )

    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{_stem(deck, smoke=smoke)}.npz"
    assert res.arm_norm is not None  # memory= was given, so all six are filled
    assert res.arm_norm_by_channel is not None
    assert res.arm_peak is not None
    assert res.exchange is not None
    assert res.exchange_local is not None
    assert res.imbalance is not None
    # WHICH per-channel TIME SERIES are kept, and why it is not all of them.
    # `arm_peak` (every channel's running maximum) is kept in full -- it is
    # `(n_arm,)` and it is what answers "which channels received the flux".
    # The `(n_steps+1, n_arm)` series is not: at 4001 steps and 83 arms it is
    # 2.7 MB per molecule and five sixths of it is never read. (`*.npz` is
    # gitignored, so this is load time and disk, not repository weight -- the
    # committed artifact is the `.png`.) Kept are the eight largest by
    # `arm_peak` UNION
    # blocks 0-3 -- the eight the figure and the ranking need, plus the first
    # four in `h_ext`'s block order, which are what `MemorySpec.n_channels=4`
    # would have recorded and which `arm_first_four_share` measures against.
    # `arm_channel_index` travels with the array so a consumer never has to
    # assume the mapping.
    kept = np.array(sorted(set(np.argsort(res.arm_peak[:, 0])[::-1][:8].tolist()) | {0, 1, 2, 3}))
    kept = kept[kept < n_arm]

    np.savez(
        path,
        molecule=deck.name,
        e_kin=deck.e_kin,
        dt=deck.dt,
        order=deck.order,
        e_total=e_total,
        min_overlap=ing.min_overlap,
        n_arm=n_arm,
        n_ext=h_ext.shape[0],
        e_r_max=deck.e_r_max,
        e_order=deck.e_order,
        e_n_complex=deck.e_n_complex,
        R_real=nuc.points[nuc.points.imag == 0.0].real,
        gamma_local=gamma_local,
        time=res.time,
        survival=res.survival[:, 0],
        arm_norm=res.arm_norm[:, 0],
        arm_channel_index=kept,
        arm_norm_by_channel=res.arm_norm_by_channel[:, kept, 0],
        arm_peak=res.arm_peak[:, 0],
        exchange=res.exchange[:, 0],
        exchange_local=res.exchange_local[:, 0],
        imbalance=res.imbalance[:, 0],
    )
    print(f"wrote {path}", flush=True)
    return path


def _first_crossing(t: npt.NDArray[np.float64], y: npt.NDArray[np.float64], level: float) -> float:
    """First `t` at which `y` has fallen to `level`; `nan` if it never does."""
    below = np.flatnonzero(y <= level)
    return float(t[below[0]]) if below.size else float("nan")


def summarize(data: dict[str, npt.NDArray[Any]]) -> Summary:
    """Reduce one molecule's `.npz` to the campaign numbers.

    Every ratio here is either dimensionless by construction or divides by a
    named normalization that appears in the field name. NOTHING is divided by
    `Gamma_loc` or by the other exchange curve.
    """
    t = data["time"]
    s = data["survival"]
    ex = data["exchange"]
    exl = data["exchange_local"]
    imb = data["imbalance"]
    arm = data["arm_norm"]
    arm_c = data["arm_norm_by_channel"]
    peak = data["arm_peak"]
    s0 = float(s[0])

    # Step 0 has empty arms, so `exchange[0]` is exactly 0 and is neither a
    # return nor a loss; counting it either way would be an artefact of the
    # initial condition.
    positive = np.flatnonzero(ex > 0.0)
    ex_over_s = ex / s

    gamma_eff = -float(exl[0]) / s0
    decay_t = tuple(_first_crossing(t, s / s0, lev) for lev in _DECAY_LEVELS)
    decay_ratio = tuple(
        float("nan")
        if not np.isfinite(tt)
        else float((s[int(np.flatnonzero(t == tt)[0])] / s0) / np.exp(-gamma_eff * tt))
        for tt in decay_t
    )

    # The two one-sided coupling rates: `exchange` is the discrete state's,
    # `imbalance - exchange` is the arms'. Their sum is the residual, and it is
    # reported against the LARGER of the two -- against `exchange` alone it
    # diverges wherever that rate passes through zero.
    arm_rate = imb - ex
    larger = np.maximum(np.abs(ex), np.abs(arm_rate))
    ok = larger > 0.0
    rel = np.abs(imb[ok]) / larger[ok]

    arm_peak_idx = int(arm.argmax())
    ranked = np.argsort(peak)[::-1]
    total_at_peak = float(arm[arm_peak_idx])
    # Blocks 0-3 BY BLOCK INDEX, located through the stored index map -- they
    # are not the first four columns of `arm_norm_by_channel`, which holds a
    # union of the largest and the first four (see `run`).
    idx = data["arm_channel_index"]
    first_four_cols = [int(np.flatnonzero(idx == b)[0]) for b in range(4) if (idx == b).any()]
    first_four = (
        float(arm_c[arm_peak_idx, first_four_cols].sum()) / total_at_peak if total_at_peak else 0.0
    )

    return Summary(
        molecule=str(data["molecule"]),
        e_kin=float(data["e_kin"]),
        dt=float(data["dt"]),
        t_max=float(t[-1]),
        n_r=int(data["R_real"].size),
        n_arm=int(data["n_arm"]),
        n_ext=int(data["n_ext"]),
        min_overlap=float(data["min_overlap"]),
        s_d0=s0,
        gamma_eff=gamma_eff,
        arm_peak_over_s0=float(arm.max() / s0),
        arm_peak_time=float(t[arm_peak_idx]),
        arm_top_blocks=tuple(int(b) for b in ranked[:4]),
        arm_first_four_share=first_four,
        exchange_max_raw=float(ex.max()),
        exchange_min_raw=float(ex.min()),
        exchange_max_over_s0=float(ex.max() / s0),
        exchange_min_over_s0=float(ex.min() / s0),
        exchange_max_over_s=float(ex_over_s.max()),
        n_positive=int(positive.size),
        n_steps=int(t.size - 1),
        t_first_positive=float(t[positive[0]]) if positive.size else float("nan"),
        return_share=float(np.trapezoid(np.maximum(ex, 0.0), t) / np.trapezoid(np.abs(ex), t)),
        nonlocality=float(np.trapezoid(np.abs(ex - exl), t) / np.trapezoid(np.abs(exl), t)),
        decay_times=decay_t,
        decay_ratio_at_levels=decay_ratio,
        imbalance_median_rel=float(np.median(rel)),
        imbalance_max_rel=float(rel.max()),
        survival_rebound=float(s[-1] / s.min()),
    )


def print_summary(sm: Summary) -> None:
    """Print one molecule's numbers, normalizations included."""
    print(f"\n=== {sm.molecule}  E_kin = {sm.e_kin} Ha  dt = {sm.dt}  T = {sm.t_max:g} ===")
    print(
        f"  deck: N_R(real) {sm.n_r}, arms {sm.n_arm}, H_ext {sm.n_ext}, "
        f"min_overlap {sm.min_overlap:.3f}"
    )
    print(f"  S_d(0) = {sm.s_d0:.6e}   Gamma_eff = {sm.gamma_eff:.6e} a.u.")
    print(
        f"  arm norm peak / S_d(0) = {sm.arm_peak_over_s0:.4f} at t = {sm.arm_peak_time:g} "
        "(RELATIVE channel decomposition, not a population)"
    )
    print(
        f"  arm_peak ranking: top four blocks {sm.arm_top_blocks}; blocks 0-3 carry "
        f"{100 * sm.arm_first_four_share:.1f}% of the arm norm at its peak"
    )
    print(
        f"  exchange raw: max {sm.exchange_max_raw:+.4e}, min {sm.exchange_min_raw:+.4e} a.u.\n"
        f"  exchange / S_d(0): max {sm.exchange_max_over_s0:+.4e}, "
        f"min {sm.exchange_min_over_s0:+.4e}\n"
        f"  max exchange / S_d(t) = {sm.exchange_max_over_s:+.4e} (fractional return rate)"
    )
    print(
        f"  net return: {sm.n_positive} of {sm.n_steps + 1} steps positive, first at "
        f"t = {sm.t_first_positive:g}; return share {sm.return_share:.3e}"
    )
    print(
        "  nonlocality = int|exchange - exchange_local| / int|exchange_local| = "
        f"{sm.nonlocality:.3f}"
    )
    for lev, tt, ratio in zip(_DECAY_LEVELS, sm.decay_times, sm.decay_ratio_at_levels, strict=True):
        print(f"  S/S_0 = {lev:<5} at t = {tt:>8.1f}   S / exp(-Gamma_eff t) = {ratio:.3e}")
    print(
        f"  imbalance / larger one-sided rate: median {sm.imbalance_median_rel:.3f}, "
        f"max {sm.imbalance_max_rel:.3f}"
    )
    print(f"  survival rebound S(T)/min S = {sm.survival_rebound:.3f}")


def _ordering_table(summaries: list[Summary], data: list[dict[str, npt.NDArray[Any]]]) -> None:
    """Print the campaign's comparators, one row per molecule.

    NAMED FOR A QUESTION THIS TABLE NO LONGER ANSWERS. It was written to read a
    three-molecule ordering off these columns; `ENERGY_LADDER` and
    `COARSE_GRAINED_RETURN` between them retract that, and the footer says so.
    The rows are still the campaign's summary and are still worth printing.

    In the ENERGY domain the LCP's failure is ordered N2 (mild) then F2 (its
    ratio to the exact oracle sweeps 0.263 -> 1.736 across 0.010-0.050 Ha),
    with NO UNDETERMINED -- its pole walk does not converge in the electronic
    box, so it is unranked rather than ranked last. These are the time-domain
    candidates for the same ordering.

    `max positive / peak` is the one to compare across molecules: it is each
    run's largest return as a fraction of its OWN largest exchange, so it is
    the quantity `RESOLVED_RETURN`'s floor is also a fraction of. `S_d/S_d(0)
    at first return` is here because the first thing a sceptical reader
    assumes about a positive excursion is that it is a decayed-tail artefact,
    and that column answers it without further argument.
    """
    print("\n=== the three molecules, each at one energy ===")
    header = (
        f"{'molecule':<9}{'E_kin':>7}{'max pos/peak':>14}{'first return':>14}"
        f"{'S/S_0 there':>13}{'% returning':>13}{'nonlocality':>13}"
    )
    print(header)
    rows = sorted(zip(summaries, data, strict=True), key=lambda pair: pair[0].nonlocality)
    for sm, d in rows:
        ex = d["exchange"]
        peak = float(np.abs(ex).max())
        pos = ex > 0.0
        s_at = float(d["survival"][pos][0] / d["survival"][0]) if pos.any() else float("nan")
        print(
            f"{sm.molecule:<9}{sm.e_kin:>7.3f}{sm.exchange_max_raw / peak:>14.3e}"
            f"{sm.t_first_positive:>14.1f}{s_at:>13.3f}"
            f"{100 * sm.n_positive / (sm.n_steps + 1):>12.1f}%{sm.nonlocality:>13.3f}"
        )
    print(
        "  `nonlocality` orders these three N2 < NO < F2, and ENERGY_LADDER is\n"
        "  what says whether that means anything: it holds across every rung\n"
        "  where the observable measures the kernel rather than the launch\n"
        "  transient, and fails on the near-threshold rungs where the Markovian\n"
        "  REFERENCE collapses and inflates the ratio. The return columns are\n"
        "  frozen in energy but are only RESOLVED on NO -- read them against\n"
        "  COARSE_GRAINED_RETURN before quoting them."
    )


#: Bin widths the coarse-grained discriminator is reported at. The verdict must
#: not depend on the choice, so it is measured at four and all four are printed.
_BIN_WIDTHS: Final[tuple[float, ...]] = (5.0, 10.0, 20.0, 50.0)


def nonlocality_post_peak(
    exchange: npt.NDArray[np.float64],
    exchange_local: npt.NDArray[np.float64],
    arm_norm: npt.NDArray[np.float64],
    time: npt.NDArray[np.float64],
) -> float:
    """`nonlocality` integrated from the arm-norm peak onwards.

    WHY THE FULL-RUN INTEGRAL IS THE WRONG WINDOW. While the arms are filling
    `X` is still near zero, so `|X - X_loc|` is near `|X_loc|` and the ratio is
    pinned near 1 by arithmetic rather than by the kernel. Every propagation
    passes through that, and near a threshold it dominates: the Markovian
    reference collapses (`int|X_loc|` falls 46x across N2's ladder) while that
    floor does not, and the full-run ratio inflates for reasons that are not
    about nonlocality.

    The contamination is a WINDOW, not a rung. Integrating from `t_arm_peak`
    onwards removes it directly, instead of discarding whole propagations for
    containing it -- which is what this module did first, and which cost four
    of seventeen rungs including N2's second-lowest nonlocality.

    `t_arm_peak` is not a tuned knob: it is identical at every energy within a
    molecule (18 / 55 / 40 a.u. for N2 / NO / F2), and starting at 2x or 3x it
    gives the same verdict.

    Parameters
    ----------
    exchange, exchange_local, arm_norm, time : ndarray
        One propagation's recorded series, all on `time`.

    Returns
    -------
    float
        `int|X - X_loc| / int|X_loc|` over `t >= t_arm_peak`.
    """
    start = int(np.argmax(arm_norm))
    num = np.abs(exchange - exchange_local)[start:]
    den = np.abs(exchange_local)[start:]
    t = time[start:]
    return float(np.trapezoid(num, t) / np.trapezoid(den, t))


def coarse_grained_return(
    ex_a: npt.NDArray[np.float64],
    ex_b: npt.NDArray[np.float64],
    dt_a: float,
    dt_b: float,
    width: float,
) -> tuple[float, float, int, float]:
    """How much of the RETURN survives refinement, after averaging in time.

    WHY A COARSE-GRAINED VERSION EXISTS AT ALL, since adding a second metric
    after the first one fails is exactly how a claim gets rescued instead of
    tested. The pointwise sign of `exchange` is not a converged quantity on ANY
    of the three molecules: refining `dt` does not lengthen the sign-flip period
    in atomic units, it shrinks it to stay at roughly two steps (F2
    2.09 -> 0.84 -> 0.52 a.u. at dt = 1 / 0.5 / 0.25; N2 7.33 -> 1.34 -> 0.72;
    NO 34.5 -> 3.19). A quantity whose structure sits at the step scale at every
    step size is being measured at the wrong resolution.

    The returns the campaign is about are not step-scale: on NO they are bursts
    hundreds of a.u. long, plainly visible in the figure's middle row. So the
    comparison is made on time-averaged bins, at four widths (`_BIN_WIDTHS`) so
    the answer cannot be an artefact of one. Binning two runs at DIFFERENT `dt`
    to a common width is sound -- the bin edges align, and each side is a time
    average of the same integral.

    IT REPORTS ITS OWN NULL, because the concordance is CONDITIONAL ("of the
    bins A calls returning, how many does B?") and the null for that is B's own
    positive-bin rate, NOT one half. On F2 that rate is ~0.59, so a raw
    concordance of 0.65 is barely above chance; on NO it is ~0.35 against a
    concordance of 0.96-1.00. Quoting the raw number without its null overstates
    every molecule with a high base rate -- it did exactly that for F2 here --
    which is why `lift` is what the verdict is read off.

    AND IT REPORTS MAGNITUDE, because agreeing about the SIGN of a bin says
    nothing about its SIZE, and "resolved" has to mean both. That column is what
    separates NO (2.0-19.1%) from N2 (124.9-147.5%) and F2 (83.8-520.1%).

    Parameters
    ----------
    ex_a, ex_b : ndarray
        The two `exchange` series, each on its own uniform time grid.
    dt_a, dt_b : float
        Their step sizes. The two need not match; each is binned to `width`.
    width : float
        Bin width in atomic units.

    Returns
    -------
    lift : float
        Concordance minus the null -- the fraction of `ex_a`'s returning bins
        that `ex_b` also calls returning, in excess of `ex_b`'s own positive-bin
        rate. Zero is chance. `nan` if `ex_a` has no returning bin.
    concordance : float
        The raw conditional fraction, before the null is subtracted.
    n_returning : int
        How many bins the fraction is taken over. Reported because a
        concordance over two bins and one over five hundred are not the same
        evidence, and because a `lift` of zero can mean either "B agrees at
        chance" or "B has no positive bins at all at this width".
    median_relative_magnitude : float
        Median `|ex_b - ex_a| / |ex_a|` over those same bins -- how far apart
        the two runs are in SIZE where they both see a return.
    """

    def binned(ex: npt.NDArray[np.float64], dt: float) -> npt.NDArray[np.float64]:
        per = round(width / dt)
        if per < 1:
            raise ValueError(f"bin width {width} is below the step size {dt}")
        keep = (len(ex) // per) * per
        return ex[:keep].reshape(-1, per).mean(axis=1)

    a, b = binned(ex_a, dt_a), binned(ex_b, dt_b)
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    pos = a > 0.0
    if not pos.any():
        return float("nan"), float("nan"), 0, float("nan")
    concordance = float((b[pos] > 0.0).mean())
    null = float((b > 0.0).mean())
    median_rel = float(np.median(np.abs(b - a)[pos] / np.abs(a)[pos]))
    return concordance - null, concordance, int(pos.sum()), median_rel


def compare_resolution(base: Path, other: Path) -> None:
    """How much of the exchange curve is resolved, by comparing two runs.

    WHY THIS EXISTS, AND WHAT IT FOUND. `exchange` is a rate read off a
    propagated state, and the two are not equally converged: on the N2
    campaign deck, halving `dt` moves `survival` by 3.6e-4 relative but moves
    `exchange` by O(1) RELATIVE TO ITS OWN LOCAL VALUE wherever that value is
    small. In absolute terms the two runs agree to 2.2% of the exchange PEAK,
    which is the honest floor: any feature smaller than that -- and N2's
    positive excursions are 1e-5 to 1e-4 of the peak -- is below the
    propagator's own discretization error and its SIGN is not resolved.

    So this prints both: the integral comparators, which converge, and the
    pointwise floor, which is what the returning-flux count has to be read
    against. Run it with a halved `dt` (`--dt 0.5 --steps 8000`, which changes
    the truncation AND the number of solves) and with a raised order
    (`--order 4`, which changes only the truncation; `pade_roots` stops at 4) -- the two together
    separate Pade truncation from accumulated round-off.

    Parameters
    ----------
    base, other : Path
        Two campaign `.npz` on the same deck and energy. Only the times they
        share are compared.
    """
    a = dict(np.load(base, allow_pickle=False))
    b = dict(np.load(other, allow_pickle=False))
    # The INTERSECTION of the two time grids, not the whole of either. A
    # resolution check is allowed to be shorter than the campaign run -- the
    # exchange peak is in the launch transient, so a 1000 a.u. window already
    # bounds the floor as a fraction of it -- and forcing it to cover the full
    # 4000 would double a 70-minute propagation to bound a number the first
    # quarter of it already bounds.
    in_a = np.isin(a["time"], b["time"])
    in_b = np.isin(b["time"], a["time"])
    if not in_a.any():
        raise SystemExit(f"{other.name} and {base.name} share no time points")
    ex_a, ex_b = a["exchange"][in_a], b["exchange"][in_b]
    s_a, s_b = a["survival"][in_a], b["survival"][in_b]
    peak = float(np.abs(ex_a).max())
    sm_a, sm_b = summarize(a), summarize(b)
    print(f"\n=== resolution: {base.name}  vs  {other.name} ===")
    print(
        f"  compared over the shared window t = {a['time'][in_a][0]:g} to "
        f"{a['time'][in_a][-1]:g} ({int(in_a.sum())} points); the summary rows below "
        "are each run's OWN full window"
    )
    print(
        f"  survival:   max|dS|/S_d(0) = {float(np.abs(s_b - s_a).max() / a['survival'][0]):.3e}"
        f"   S(T) ratio = {float(s_b[-1] / s_a[-1]):.4f}"
    )
    print(
        f"  exchange:   max|d(exchange)| / max|exchange| = "
        f"{float(np.abs(ex_b - ex_a).max()) / peak:.3e}   <-- the pointwise floor"
    )
    print(
        f"  sign agreement over shared steps: {float((np.sign(ex_a) == np.sign(ex_b)).mean()):.3f}"
    )
    print(
        f"  positive maximum / max|exchange|: {sm_a.exchange_max_raw / peak:.3e} "
        f"vs {sm_b.exchange_max_raw / peak:.3e}"
    )
    # THE FLOOR THAT ACTUALLY BOUNDS THE RETURNING FLUX, and why the global one
    # above does not. `max|d(exchange)|` over the WHOLE run is normally set in
    # the launch transient, where `|exchange|` is at its peak and changing
    # fastest -- on NO it lands at t = 15, hundreds of a.u. before the first
    # positive step, and reports 7.8e-2 for a run whose late structure the two
    # propagators agree on to 8.9e-3. Restricting to the window the returning
    # flux lives in compares like with like. Both are printed: the global one
    # bounds the curve, this one bounds the claim.
    pos_a = ex_a > 0.0
    if pos_a.any():
        t_shared = a["time"][in_a]
        window = t_shared >= 0.9 * float(t_shared[pos_a][0])
        local_floor = float(np.abs(ex_b - ex_a)[window].max()) / peak
        concordant = float((ex_b[pos_a] > 0.0).mean())
        median_rel = float(np.median(np.abs(ex_b - ex_a)[pos_a] / np.abs(ex_a)[pos_a]))
        print(
            f"  floor over the returning window (t >= {0.9 * float(t_shared[pos_a][0]):g}): "
            f"{local_floor:.3e}   <-- what the positive maximum must clear"
        )
        # The primary discriminator, and normalization-free: of the steps the
        # first run calls a RETURN, how many does the second call a return too?
        # A rate whose sign is resolved scores near 1; one at the numerical
        # floor scores near the fraction of steps that are positive by chance.
        print(
            f"  of {int(pos_a.sum())} returning steps in {base.name}, "
            f"{100 * concordant:.1f}% are returning in {other.name} too; "
            f"median |d(exchange)|/|exchange| there {median_rel:.3e}"
        )
    print(
        f"  returning steps: {sm_a.n_positive} (first t = {sm_a.t_first_positive:g}) "
        f"vs {sm_b.n_positive} (first t = {sm_b.t_first_positive:g})"
    )
    # The coarse-grained discriminator, at every bin width, because the
    # pointwise sign above is not a converged quantity on any molecule.
    binned = [
        coarse_grained_return(a["exchange"], b["exchange"], float(a["dt"]), float(b["dt"]), w)
        for w in _BIN_WIDTHS
    ]
    print(
        "  RETURNING BINS, lift over the null (raw concordance, n, median |d|/|x|):\n"
        + "\n".join(
            f"    bin {w:>4g}:  lift {lift:+.3f}   (raw {raw:.3f}, n = {n}, magnitude {rel:.3f})"
            for w, (lift, raw, n, rel) in zip(_BIN_WIDTHS, binned, strict=True)
        )
    )
    print(f"  nonlocality:      {sm_a.nonlocality:.4f} vs {sm_b.nonlocality:.4f}")
    print(f"  arm peak / S_d(0): {sm_a.arm_peak_over_s0:.4f} vs {sm_b.arm_peak_over_s0:.4f}")
    print(f"  Gamma_eff:        {sm_a.gamma_eff:.6e} vs {sm_b.gamma_eff:.6e}")


def write_figure(paths: list[Path], outdir: Path) -> Path:
    """The three-row, one-column-per-molecule campaign figure.

    Rows, top to bottom: the partition (`S_d` and the arm norm, BOTH divided by
    `S_d(0)`, and the three channels with the largest running maximum); the
    exchange rate against its Markovian limit (BOTH divided by `S_d(0)` -- the
    curves are never divided by each other, and nothing is divided by
    `Gamma_loc`); and the decay law against `exp(-Gamma_eff t)`.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    loaded = [dict(np.load(p, allow_pickle=False)) for p in paths]
    n = len(loaded)
    fig, axes = plt.subplots(3, n, figsize=(5.2 * n, 11.0), squeeze=False)

    for col, data in enumerate(loaded):
        sm = summarize(data)
        t = data["time"]
        s0 = sm.s_d0
        ax = axes[0][col]
        ax.plot(t, data["survival"] / s0, color="C0", lw=1.4, label=r"$S_d(t)/S_d(0)$")
        ax.plot(
            t,
            data["arm_norm"] / s0,
            color="C3",
            lw=1.4,
            label=r"$\sum_n\|\varphi_n\|^2 / S_d(0)$",
        )
        idx = list(data["arm_channel_index"])
        ranked = np.argsort(data["arm_peak"])[::-1][:3]
        for k, block in enumerate(ranked):
            ax.plot(
                t,
                data["arm_norm_by_channel"][:, idx.index(block)] / s0,
                color="C3",
                lw=0.8,
                alpha=0.55,
                ls=[":", "--", "-."][k],
                label=rf"arm block {block}",
            )
        ax.set_yscale("log")
        ax.set_ylim(1e-8, 3.0)
        ax.set_title(
            f"{sm.molecule}   $E$ = {sm.e_kin} Ha   {sm.n_arm} arms   $H_{{ext}}$ = {sm.n_ext}",
            fontsize=11,
        )
        ax.set_ylabel("relative channel decomposition,\nnormalized by $S_d(0)$ (NOT a population)")
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(alpha=0.25)

        ax = axes[1][col]
        scale = max(abs(sm.exchange_min_over_s0), abs(sm.exchange_max_over_s0))
        ax.plot(
            t,
            data["exchange"] / s0,
            color="C0",
            lw=1.2,
            label=(
                r"nonlocal  $2\,\mathrm{Im}\langle\Psi_d|"
                r"\sum_n V_{dn}\varphi_n\rangle / S_d(0)$"
            ),
        )
        ax.plot(
            t,
            data["exchange_local"] / s0,
            color="C1",
            lw=1.2,
            ls="--",
            label=r"Markovian  $-\langle\Psi_d|\Gamma_{loc}|\Psi_d\rangle / S_d(0)$",
        )
        pos = data["exchange"] > 0.0
        if pos.any():
            # The per-step markers are what the coarse-grained discriminator
            # demoted, so the legend carries that verdict rather than leaving
            # the reader to infer a resolved structure from a dense scatter.
            band = [v for (mol, _), v in COARSE_GRAINED_RETURN.items() if mol == sm.molecule]
            verdict = (
                f"; binned lift {min(lo for lo, _, _, _, _ in band):+.2f}"
                f"..{max(hi for _, hi, _, _, _ in band):+.2f}"
                if band
                else ""
            )
            ax.plot(
                t[pos],
                data["exchange"][pos] / s0,
                ls="none",
                marker=".",
                ms=2.5,
                color="C2",
                label=(
                    f"returning ({sm.n_positive} steps, first $t$ = "
                    f"{sm.t_first_positive:g}{verdict})"
                ),
            )
        ax.set_yscale("symlog", linthresh=1e-12 if scale == 0 else scale * 1e-7)
        ax.axhline(0.0, color="k", lw=0.6)
        ax.set_ylabel("exchange rate / $S_d(0)$  (a.u.$^{-1}$)")
        ax.legend(fontsize=7, loc="lower right")
        ax.grid(alpha=0.25)

        ax = axes[2][col]
        ax.plot(t, data["survival"] / s0, color="C0", lw=1.4, label=r"$S_d(t)/S_d(0)$")
        ax.plot(
            t,
            np.exp(-sm.gamma_eff * t),
            color="C1",
            lw=1.2,
            ls="--",
            label=rf"$\exp(-\Gamma_{{eff}}t)$, $\Gamma_{{eff}}$ = {sm.gamma_eff:.3e}",
        )
        for lev, tt in zip(_DECAY_LEVELS, sm.decay_times, strict=True):
            if np.isfinite(tt):
                ax.axvline(tt, color="0.6", lw=0.7, ls=":")
                ax.annotate(
                    f"{lev:g}",
                    (tt, lev),
                    fontsize=7,
                    textcoords="offset points",
                    xytext=(3, 3),
                )
        ax.set_yscale("log")
        ax.set_ylim(1e-11, 3.0)
        ax.set_xlabel("$t$ (a.u.)")
        ax.set_ylabel("$S_d(t)/S_d(0)$  vs one rate constant")
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(alpha=0.25)

    fig.suptitle(
        "TD nonlocal resonance model: memory observables\n"
        "arm norms are a RELATIVE CHANNEL DECOMPOSITION, not populations;\n"
        "both exchange curves are divided by $S_d(0)$; the per-step returning\n"
        "markers are RESOLVED ONLY ON NO (binned lift over chance, per panel)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "nrm-td-memory-observables.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


#: The energy rungs each explanatory figure reads, by molecule, as
#: `(E_kin, .npz suffix)`. `None` is the campaign run's bare stem.
_LADDER_FILES: Final[dict[str, tuple[tuple[float, str | None], ...]]] = {
    "N2": (
        (0.050, "e0.05"),
        (0.060, "e0.06"),
        (0.080, "e0.08"),
        (0.100, None),
        (0.120, "e0.12"),
        (0.150, "e0.15"),
        (0.160, "e0.16"),
    ),
    "NO": ((0.175, "e0.175"), (0.185, "e0.185"), (0.200, None), (0.300, "e0.3"), (0.400, "e0.4")),
    "F2": (
        (0.010, "e0.01"),
        (0.020, "e0.02"),
        (0.030, None),
        (0.040, "e0.04"),
        (0.050, "e0.05"),
    ),
}

#: The refinement runs the resolvability figure compares each campaign run
#: against, by molecule.
_REFINEMENT_FILES: Final[dict[str, tuple[str, ...]]] = {
    "N2": ("order4", "dt0.5", "dt0.25"),
    "NO": ("order4", "dt0.5"),
    "F2": ("order4", "dt0.5", "dt0.25"),
}

_COLOUR: Final[dict[str, str]] = {"N2": "C0", "NO": "C1", "F2": "C2"}


def _load(outdir: Path, molecule: str, suffix: str | None) -> dict[str, Any]:
    """One recorded propagation, by molecule and `.npz` suffix."""
    stem = f"{molecule.lower()}-nrm-memory-observables"
    if suffix is not None:
        stem += f"-{suffix}"
    return dict(np.load(outdir / f"{stem}.npz", allow_pickle=False))


def _draw_arrow_hamiltonian(ax: Any) -> None:
    """Schematic of `H_ext`'s arrow block structure and what reads off it."""
    from matplotlib.patches import FancyArrowPatch, Rectangle

    ax.set_xlim(-1.5, 7.6)
    ax.set_ylim(-2.6, 6.2)
    ax.set_aspect("equal")
    ax.axis("off")

    def block(i: int, j: int, colour: str, label: str) -> None:
        ax.add_patch(
            Rectangle((j, 4 - i), 0.92, 0.92, facecolor=colour, edgecolor="k", lw=1.0, alpha=0.85)
        )
        ax.text(j + 0.46, 4.46 - i, label, ha="center", va="center", fontsize=9)

    block(0, 0, "#f4a582", r"$H_{dd}$")
    for k in range(1, 5):
        block(0, k, "#92c5de", r"$V_{dn}$")
        block(k, 0, "#92c5de", r"$V_{dn}$")
        block(k, k, "#d9d9d9", r"$H_n$")
    for k in range(1, 5):
        for m in range(1, 5):
            if k != m:
                ax.text(m + 0.46, 4.46 - k, "0", ha="center", va="center", fontsize=8, color="0.6")

    ax.text(0.46, 5.25, r"$\Psi_d$", ha="center", fontsize=11)
    ax.text(3.0, 5.85, r"arms  $\varphi_n$", ha="center", fontsize=11)
    ax.annotate("", xy=(1.0, 5.45), xytext=(4.9, 5.45), arrowprops={"arrowstyle": "<->"})
    ax.text(-0.25, 4.46, r"$\Psi_d$", ha="right", va="center", fontsize=11)
    ax.text(-0.25, 2.5, r"$\varphi_n$", ha="right", va="center", fontsize=11)

    # What each part of the state is read for, to the right of the matrix.
    ax.add_patch(
        FancyArrowPatch(
            (0.92, 4.46), (5.45, 4.46), arrowstyle="-|>", mutation_scale=12, color="#b2182b"
        )
    )
    ax.text(5.6, 4.46, r"$S_d = \|\Psi_d\|^2$", va="center", fontsize=9.5, color="#b2182b")
    ax.text(
        5.6,
        3.75,
        r"$X_{loc} = -\langle\Psi_d|\Gamma_{loc}|\Psi_d\rangle$",
        va="center",
        fontsize=9.5,
        color="#b2182b",
    )
    ax.add_patch(
        FancyArrowPatch(
            (4.92, 2.5), (5.45, 2.5), arrowstyle="-|>", mutation_scale=12, color="#b2182b"
        )
    )
    ax.text(5.6, 2.5, r"$A = \sum_n \|\varphi_n\|^2$", va="center", fontsize=9.5, color="#b2182b")

    # The cross term, pointed at the coupling blocks it is built from.
    ax.add_patch(
        FancyArrowPatch(
            (2.5, -0.35), (2.5, 0.55), arrowstyle="-|>", mutation_scale=12, color="#b2182b"
        )
    )
    ax.text(
        2.5,
        -0.75,
        r"$X = 2\,\mathrm{Im}\,\langle\Psi_d|\sum_n V_{dn}\varphi_n\rangle$",
        ha="center",
        va="top",
        fontsize=10,
        color="#b2182b",
    )
    ax.text(
        2.5,
        -1.75,
        "the DOORWAY talking to the ARMS:\nthe rate memory feeds amplitude back",
        ha="center",
        va="top",
        fontsize=8.5,
        style="italic",
    )
    ax.set_title(
        "$H_{ext}$ is an ARROW matrix — the arms couple to the doorway and not to\n"
        "each other. Eliminating them gives PRA 77 Eq. (55)'s kernel $F(E)$;\n"
        "keeping them makes the memory watchable.",
        fontsize=10,
        pad=14,
    )


def write_construction_figure(outdir: Path) -> Path:
    """What the observables are, and the two windows of a propagation."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = _load(outdir, "NO", None)
    t, s0 = d["time"], float(d["survival"][0])
    peak = float(t[int(np.argmax(d["arm_norm"]))])

    fig = plt.figure(figsize=(13.0, 9.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=0.32, wspace=0.22)

    _draw_arrow_hamiltonian(fig.add_subplot(gs[0, 0]))

    ax = fig.add_subplot(gs[0, 1])
    ax.semilogy(t, d["survival"] / s0, color="C0", lw=1.6, label=r"$S_d(t)/S_d(0)$  (doorway)")
    ax.semilogy(t, d["arm_norm"] / s0, color="C3", lw=1.6, label=r"$A(t)/S_d(0)$  (arms)")
    ax.axvline(peak, color="k", ls=":", lw=1.2)
    ax.axvspan(0, peak, color="0.85", zorder=0)
    ax.text(peak / 2, 2e-3, "arms\nfilling", ha="center", fontsize=9)
    ax.annotate(
        f"$t_{{peak}} = {peak:g}$",
        xy=(peak, 3e-1),
        xytext=(peak + 400, 3e-1),
        fontsize=9,
        arrowprops={"arrowstyle": "->"},
    )
    ax.set_xlabel("$t$ (a.u.)")
    ax.set_ylabel("normalized by $S_d(0)$")
    ax.set_title("NO: amplitude leaves the doorway and fills the arms", fontsize=10)
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.25)
    ax.set_ylim(1e-4, 3)

    ax = fig.add_subplot(gs[1, :])
    x, xl = d["exchange"] / s0, d["exchange_local"] / s0
    ax.plot(t, x, color="C0", lw=0.8, label=r"nonlocal  $X/S_d(0)$")
    ax.plot(t, xl, color="C1", lw=1.4, ls="--", label=r"Markovian  $X_{loc}/S_d(0)$")
    ax.fill_between(
        t, 0, x, where=x > 0, color="C2", alpha=0.65, label="$X>0$: amplitude RETURNING"
    )
    ax.axhline(0.0, color="k", lw=0.7)
    ax.axvline(peak, color="k", ls=":", lw=1.2)
    ax.axvspan(0, peak, color="0.85", zorder=0)
    ax.set_yscale("symlog", linthresh=1e-8)
    ax.set_xlabel("$t$ (a.u.)")
    ax.set_ylabel("exchange rate / $S_d(0)$  (a.u.$^{-1}$)")
    ax.set_title(
        r"$X_{loc} \leq 0$ wherever $\Gamma_{loc} \geq 0$ — a LOCAL model can only drain the "
        r"doorway. $X>0$ is what the LCP cannot represent.",
        fontsize=10,
    )
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(alpha=0.25)

    fig.suptitle(
        "The memory observables: where the amplitude is, and which way it is flowing",
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = outdir / "nrm-memory-construction.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def write_resolvability_figure(outdir: Path) -> Path:
    """Why the returning flux can be read on NO and on neither of the others."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.9))

    ax = axes[0]
    for molecule, suffixes in _REFINEMENT_FILES.items():
        base = _load(outdir, molecule, None)
        for k, suffix in enumerate(suffixes):
            other = _load(outdir, molecule, suffix)
            lifts = [
                coarse_grained_return(
                    base["exchange"],
                    other["exchange"],
                    float(base["dt"]),
                    float(other["dt"]),
                    w,
                )[0]
                for w in _BIN_WIDTHS
            ]
            ax.plot(
                _BIN_WIDTHS,
                lifts,
                marker="o",
                ms=4,
                color=_COLOUR[molecule],
                alpha=0.9,
                label=molecule if k == 0 else None,
            )
    ax.axhline(0.0, color="k", lw=1.0)
    ax.text(5.5, 0.02, "chance", fontsize=8)
    ax.set_xscale("log")
    ax.set_xticks(list(_BIN_WIDTHS))
    ax.set_xticklabels([f"{w:g}" for w in _BIN_WIDTHS])
    ax.minorticks_off()
    ax.set_xlabel("bin width (a.u.)")
    ax.set_ylabel("lift of returning bins over its own null")
    ax.set_title("does a finer run agree WHICH bins return?", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    ax = axes[1]
    for molecule, suffixes in _REFINEMENT_FILES.items():
        base = _load(outdir, molecule, None)
        for k, suffix in enumerate(suffixes):
            other = _load(outdir, molecule, suffix)
            mags = [
                coarse_grained_return(
                    base["exchange"],
                    other["exchange"],
                    float(base["dt"]),
                    float(other["dt"]),
                    w,
                )[3]
                for w in _BIN_WIDTHS
            ]
            ax.semilogy(
                _BIN_WIDTHS,
                mags,
                marker="s",
                ms=4,
                color=_COLOUR[molecule],
                alpha=0.9,
                label=molecule if k == 0 else None,
            )
    ax.axhline(1.0, color="k", lw=1.0, ls="--")
    ax.text(5.5, 1.15, "disagreement as large as the signal", fontsize=8)
    ax.set_xscale("log")
    ax.set_xticks(list(_BIN_WIDTHS))
    ax.set_xticklabels([f"{w:g}" for w in _BIN_WIDTHS])
    ax.minorticks_off()
    ax.set_xlabel("bin width (a.u.)")
    ax.set_ylabel(r"median $|\Delta X| / |X|$ on returning bins")
    ax.set_title("...and do they agree HOW MUCH returns?", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    ax = axes[2]
    width = 20.0
    labels, raws, nulls, colours = [], [], [], []
    for molecule, suffixes in _REFINEMENT_FILES.items():
        base = _load(outdir, molecule, None)
        other = _load(outdir, molecule, suffixes[0])
        lift, raw, _n, _m = coarse_grained_return(
            base["exchange"], other["exchange"], float(base["dt"]), float(other["dt"]), width
        )
        labels.append(molecule)
        raws.append(raw)
        nulls.append(raw - lift)
        colours.append(_COLOUR[molecule])
    xs = np.arange(len(labels))
    ax.bar(xs - 0.19, raws, 0.36, color=colours, label="raw concordance")
    ax.bar(xs + 0.19, nulls, 0.36, color=colours, alpha=0.4, hatch="//", label="its null")
    for i, (r, nl) in enumerate(zip(raws, nulls, strict=True)):
        ax.annotate(
            "",
            xy=(i - 0.19, r),
            xytext=(i - 0.19, nl),
            arrowprops={"arrowstyle": "<->", "color": "k", "lw": 1.4},
        )
        ax.text(
            i,
            max(r, nl) + 0.06,
            f"lift {r - nl:+.2f}",
            fontsize=9.5,
            ha="center",
            fontweight="bold" if r - nl > 0.3 else "normal",
        )
    ax.set_xticks(xs)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1.22)
    ax.set_ylabel("fraction of bins called returning")
    ax.set_title("why the RAW number misleads (20 a.u. bins)", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.25, axis="y")

    fig.suptitle(
        "Is the returning flux resolved?  Only on NO — and it takes both panels to see it",
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path = outdir / "nrm-memory-resolvability.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def write_ladder_figure(outdir: Path) -> Path:
    """The near-threshold inflation, the window that removes it, the ordering."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data: dict[str, list[tuple[float, float, float, float, bool]]] = {}
    for molecule, rungs in _LADDER_FILES.items():
        lo, hi = ENERGY_WINDOWS[molecule]
        rows = []
        for e_kin, suffix in rungs:
            z = _load(outdir, molecule, suffix)
            t = z["time"]
            num = np.abs(z["exchange"] - z["exchange_local"])
            den = np.abs(z["exchange_local"])
            rows.append(
                (
                    e_kin,
                    float(np.trapezoid(num, t) / np.trapezoid(den, t)),
                    nonlocality_post_peak(z["exchange"], z["exchange_local"], z["arm_norm"], t),
                    float(np.trapezoid(den, t)),
                    lo <= e_kin <= hi,
                )
            )
        data[molecule] = rows

    fig, axes = plt.subplots(1, 4, figsize=(18.5, 4.9))

    for ax, col, title, ylab in (
        (
            axes[0],
            1,
            "FULL-RUN integral: inflates near threshold",
            "nonlocality (full run)",
        ),
        (
            axes[1],
            2,
            "POST-PEAK integral: three clean bands",
            "nonlocality (from $t_{peak}$)",
        ),
    ):
        for molecule, rows in data.items():
            e = [r[0] for r in rows]
            y = [r[col] for r in rows]
            keep = [r[4] for r in rows]
            ax.plot(e, y, "-", color=_COLOUR[molecule], lw=1.2, alpha=0.5)
            ax.plot(
                [x for x, k in zip(e, keep, strict=True) if k],
                [v for v, k in zip(y, keep, strict=True) if k],
                "o",
                ms=6,
                color=_COLOUR[molecule],
                label=molecule,
            )
            ax.plot(
                [x for x, k in zip(e, keep, strict=True) if not k],
                [v for v, k in zip(y, keep, strict=True) if not k],
                "x",
                ms=8,
                mew=2,
                color=_COLOUR[molecule],
            )
        ax.set_xscale("log")
        ax.set_xlabel("$E_{kin}$ (Ha)")
        ax.set_ylabel(ylab)
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(alpha=0.25)
    axes[0].text(
        0.011,
        0.52,
        "x = outside that molecule's\ndeclared energy window\n(both were added by this ladder)",
        fontsize=8.5,
        va="top",
        bbox={"boxstyle": "round,pad=0.35", "fc": "w", "ec": "0.6"},
    )

    ax = axes[2]
    for molecule, rows in data.items():
        ax.loglog(
            [r[0] for r in rows],
            [r[3] for r in rows],
            "o-",
            ms=5,
            color=_COLOUR[molecule],
            label=molecule,
        )
    ax.set_xlabel("$E_{kin}$ (Ha)")
    ax.set_ylabel(r"$\int |X_{loc}|\,dt$   (the Markovian REFERENCE)")
    ax.set_title("the denominator collapses — that is the\ncause of the inflation", fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25, which="both")

    ax = axes[3]
    from matplotlib.patches import Rectangle

    for i, molecule in enumerate(("N2", "NO", "F2")):
        vals = [r[2] for r in data[molecule] if r[4]]
        lo_v, hi_v = min(vals), max(vals)
        # NO's band is 0.002 wide and would vanish; every band is drawn at a
        # visible minimum thickness, with the true range printed beside it.
        pad = max(0.0, (0.018 - (hi_v - lo_v)) / 2)
        ax.add_patch(
            Rectangle(
                (i - 0.26, lo_v - pad),
                0.52,
                (hi_v - lo_v) + 2 * pad,
                facecolor=_COLOUR[molecule],
                edgecolor="k",
                lw=0.8,
                alpha=0.85,
            )
        )
        ax.plot([i] * len(vals), vals, "k.", ms=5, zorder=3)
        ax.text(
            i,
            hi_v + pad + 0.03,
            f"{lo_v:.3f}-{hi_v:.3f}\n({len(vals)} rungs)",
            ha="center",
            fontsize=8.5,
        )
    n2 = [r[2] for r in data["N2"] if r[4]]
    no = [r[2] for r in data["NO"] if r[4]]
    f2 = [r[2] for r in data["F2"] if r[4]]
    for lo_v, hi_v, x0, lab in (
        (max(n2), min(no), 0.5, f"+{100 * (min(no) / max(n2) - 1):.1f}%"),
        (max(no), min(f2), 1.5, f"+{100 * (min(f2) / max(no) - 1):.1f}%"),
    ):
        ax.annotate(
            "", xy=(x0, hi_v), xytext=(x0, lo_v), arrowprops={"arrowstyle": "<->", "color": "k"}
        )
        ax.text(x0 + 0.06, (lo_v + hi_v) / 2, lab, fontsize=9, va="center")
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["N2", "NO", "F2"])
    ax.set_xlim(-0.55, 2.55)
    ax.set_ylim(0.12, 1.52)
    ax.set_ylabel("nonlocality (post-peak, in-window)")
    ax.set_title("the result: three bands, no overlap", fontsize=10)
    ax.grid(alpha=0.25, axis="y")

    fig.suptitle(
        "The energy ladder: what went wrong with the full-run integral, and the ordering that "
        "survives",
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path = outdir / "nrm-memory-energy-ladder.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main() -> None:
    """CLI: run one molecule, print its numbers, or assemble the figure."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "molecule", choices=[*sorted(DECKS), "figure", "report", "resolution", "explain"]
    )
    p.add_argument("--outdir", type=Path, default=_FIGURE_DIR)
    p.add_argument("--energy", type=float, default=None, help="override E_kin (hartree)")
    p.add_argument("--steps", type=int, default=None, help="override n_steps")
    p.add_argument("--dt", type=float, default=None, help="override dt")
    p.add_argument("--order", type=int, default=None, help="override the diagonal-Pade order")
    p.add_argument(
        "--against",
        default="N2",
        help="molecule whose campaign run the `resolution` mode compares its variants to",
    )
    p.add_argument("--smoke", action="store_true", help="20 steps; exercises the path only")
    p.add_argument(
        "--converge",
        action="store_true",
        help="print the electronic-box ladder instead of propagating (expensive; "
        "its results are recorded in E_BOX_LADDER)",
    )
    args = p.parse_args()

    if args.molecule == "explain":
        for build in (
            write_construction_figure,
            write_resolvability_figure,
            write_ladder_figure,
        ):
            print(f"wrote {build(args.outdir)}")
        return

    if args.molecule == "resolution":
        base = args.outdir / f"{args.against.lower()}-nrm-memory-observables.npz"
        variants = sorted(
            q
            for q in args.outdir.glob(f"{args.against.lower()}-nrm-memory-observables-*.npz")
            if "-dt" in q.name or "-order" in q.name
        )
        if not variants:
            raise SystemExit(f"no --dt / --order variant of {base.name} in {args.outdir}")
        for q in variants:
            compare_resolution(base, q)
        return

    if args.molecule in ("figure", "report"):
        paths = [
            args.outdir / f"{m.lower()}-nrm-memory-observables.npz"
            for m in ("N2", "F2", "NO")
            if (args.outdir / f"{m.lower()}-nrm-memory-observables.npz").exists()
        ]
        if not paths:
            raise SystemExit(f"no campaign .npz in {args.outdir} -- run the molecules first")
        loaded = [dict(np.load(q, allow_pickle=False)) for q in paths]
        summaries = [summarize(d) for d in loaded]
        for sm in summaries:
            print_summary(sm)
        _ordering_table(summaries, loaded)
        if args.molecule == "figure":
            print(f"wrote {write_figure(paths, args.outdir)}")
        return

    deck = DECKS[args.molecule]
    if args.energy is not None:
        deck = replace_energy(deck, args.energy)
    if args.steps is not None or args.dt is not None or args.order is not None:
        deck = replace_window(deck, args.steps, args.dt, args.order)
    if args.converge:
        converge(deck, tuple(E_BOX_LADDER[deck.name]))
        return
    path = run(deck, outdir=args.outdir, smoke=args.smoke)
    print_summary(summarize(dict(np.load(path, allow_pickle=False))))


def replace_energy(deck: CampaignDeck, e_kin: float) -> CampaignDeck:
    """`deck` at another incident energy.

    A separate propagation, never another column of the same one: `F(E)` and
    therefore `Gamma_loc` carry one total energy, and `MemorySpec` carries one
    `gamma_local`, so a multi-energy launch would read `exchange_local` at the
    wrong energy for every column but one.
    """
    return CampaignDeck(**{**deck.__dict__, "e_kin": e_kin})


def replace_window(
    deck: CampaignDeck, n_steps: int | None, dt: float | None, order: int | None = None
) -> CampaignDeck:
    """`deck` over another propagation window or propagator order.

    Both of the resolution checks the campaign runs go through here: halving
    `dt` at fixed order, and raising the order at fixed `dt`. They fail
    differently -- the first changes the number of solves and so the
    accumulated round-off as well as the truncation, the second changes only
    the truncation -- which is how the two can be told apart.
    """
    return CampaignDeck(
        **{
            **deck.__dict__,
            "n_steps": deck.n_steps if n_steps is None else n_steps,
            "dt": deck.dt if dt is None else dt,
            "order": deck.order if order is None else order,
        }
    )


if __name__ == "__main__":
    main()
