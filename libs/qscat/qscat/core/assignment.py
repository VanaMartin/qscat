"""Is an exact pole a resonance at all, and which quasi-bound state is it?

`qscat.core.exact_resonance_states` finds poles by two-angle ECS stability.
That test is **necessary and not sufficient** -- exactly as
`docs/physics/lcp-resonance-levels.md` warns for the 1-D case, and as measured
in 2-D on H2+, where 4 of 57 angle-stable poles turned out not to be resonances.
This module supplies the check that catches them, and the label that a bare
complex number does not carry.

## Why overlap and not energy

Energy proximity cannot distinguish a resonance from a rotated-continuum
eigenvalue that happens to land nearby: both are just numbers on the same axis.
Overlap can, because the two behave differently under it:

- A **genuine resonance** is, to the extent the Born-Oppenheimer picture holds
  at all, a product `phi_j(r; R) chi_v(R)` plus a correction. It has substantial
  overlap with ONE BO state and little with the rest, and the node and lobe
  counts agree -- a strong constraint, since neighbouring levels differ in
  curve index or vibrational quantum number and so in node count.
- A **discretized-continuum state** rotated onto the real axis by ECS has no
  such structure. It overlaps everything weakly and nothing strongly.

Measured on H2+ window 0 the two classes separate by three orders of magnitude
(genuine states 0.87-0.99, artefacts 6e-4 to 7e-3), so this is not a marginal
judgement. It is also how the published work assigns a cross-section feature to
a quasi-bound state (M. Vana, doctoral thesis, Charles University 2017,
Table 4.2, p. 73).

Energy then enters as a CHECK rather than as the assignment, because a large
overlap with an energetically distant level is real information: it means the
state is mixing across a gap, not that it is that level shifted.

## What the verdicts mean

`pair_by_overlap` returns one of seven, in priority order:

| verdict | meaning |
|---|---|
| `ok` | a clean identification; `shift_ev` is quotable |
| `spurious` | no partner, and the basis provably covers this energy -- not a resonance |
| `basis-limited` | no partner, but the basis does NOT cover this energy -- verdict withheld |
| `box-limited` | identified, but most of the state lies outside the unscaled region |
| `weak` | a best match too poor to call an identification |
| `mixed` | a near-equal blend of two levels; neither label describes it |
| `distant` | a clear identification whose partner lies further than a shift should |

The `spurious` / `basis-limited` split is the part that needs a physical
argument rather than a threshold, and it comes from
`qscat.core.bo.admissible_levels`: a Rydberg series is attached to a CLOSED
channel, so the set of levels that can exist at a given energy is finite and
computable. Without that split the two are indistinguishable -- and conflating
them once nearly discarded eight genuine states.

## The overlap is blind to a state leaving the box

The c-product cancels the rotated ECS tail by construction. That is what makes
it the right pairing -- and it means a state with 97 % of its probability
outside the unscaled region still pairs at 0.99 with the BO product it
genuinely is. The overlap answers "which state is this" correctly and says
nothing about whether the grid holds it.

On H2+ window 0 the two diverge sharply as the Rydberg series climbs toward its
threshold: `real_weight` falls 0.998 -> 0.68 -> 0.29 -> 0.12 -> 0.031 -> 0.008
-> 0.0003 while every overlap stays near 0.99. Those orbitals are larger than
the 300-bohr box. `real_weight` measures it and the `box-limited` verdict
reports it; supply `localization=` or the check does not run.

Across the three H2+ DR windows this reclassifies **18 of 57 poles**, and every
statistic computed over the survivors moves with it: the measured regime split
in the Born-Oppenheimer shift went from 0.264 / 3.375 meV over 40 rows to
0.457 / 3.702 meV over 24. N2's poles sit at `real_weight` 0.96 and are
untouched.

## Thresholds are calibrated, not derived

The module defaults were measured on the H2+ DR windows. They are exposed as
keyword arguments because a different system may separate differently; a run on
a new model should check that its genuine and artefact populations still
separate cleanly before trusting the verdicts.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from scipy.optimize import linear_sum_assignment

from qscat.dvr import TensorGrid
from qscat.linalg import c_product
from qscat.units import HARTREE_TO_EV

from .bo import BoBasis, basis_covers

__all__ = [
    "MAX_SHIFT_EV",
    "MIN_LOCALIZATION",
    "MIXED_RATIO",
    "NO_PARTNER",
    "WEAK",
    "OverlapPair",
    "PeakAlignment",
    "overlap",
    "pair_by_overlap",
    "pair_one_to_one",
    "peak_alignment",
    "peak_positions",
    "real_weight",
]

# An overlap below this is "no partner in this basis". On its own it does NOT
# distinguish a spurious state from a real one whose partner was never built;
# `basis_covers` is what separates them.
NO_PARTNER = 0.10

# Between NO_PARTNER and this, the best match exists but is not an
# identification: on H2+ a genuine assignment scores 0.87-0.99 and even the
# strongly mixed states reach 0.63-0.78, so anything under half is a state the
# basis does not really describe.
WEAK = 0.50

# `second/best` above this means the state is a blend, not an assignment: at the
# measured Ry_16 v=1 / omega_2^5 crossing the ratio reaches 0.95 and neither
# label fits.
MIXED_RATIO = 0.70

# Below this fraction of |psi|^2 inside the unscaled region, a state has left
# the box and nothing about it is quotable -- however well it pairs. Half is a
# deliberately generous floor: the H2+ states this catches sit at 0.29 and
# below, and the ones it clears at 0.68 and above.
MIN_LOCALIZATION = 0.5

# A partner further than this in energy is flagged even when the overlap is
# high. Overlap alone would happily pair a state with a level several meV away,
# and past a few meV that is a statement about mixing rather than a shift worth
# quoting. 0.020 eV = 20 meV: every clean H2+ pairing sits well inside it, and
# the rows that exceed it are the known curve crossings.
MAX_SHIFT_EV = 0.020


@dataclass(frozen=True)
class OverlapPair:
    """One exact pole paired to a BO level BY OVERLAP, with the checks."""

    pole_energy: float  # Re E, Hartree
    level: tuple[int, int] | None  # (curve, vib); None when nothing matched
    overlap: float
    second_level: tuple[int, int] | None
    second_overlap: float
    shift_ev: float  # Re E - level energy, in eV (NaN when level is None)
    # ok | spurious | basis-limited | box-limited | weak | mixed | distant
    verdict: str

    @property
    def shift_mev(self) -> float:
        return 1000.0 * self.shift_ev

    @property
    def label(self) -> str:
        r"""`omega^j_v`, or `-` when no level was assigned."""
        if self.level is None:
            return "-"
        return rf"$\omega^{{{self.level[0]}}}_{{{self.level[1]}}}$"

    @property
    def is_quotable(self) -> bool:
        """True only for `ok`. Every other verdict withholds something."""
        return self.verdict == "ok"

    @property
    def is_resonance(self) -> bool:
        """False only for `spurious` -- the one verdict that rejects a pole.

        `basis-limited` and `weak` withhold a label without rejecting the state,
        which is the honest reading: the evidence is about the basis, not about
        the pole.
        """
        return self.verdict != "spurious"


def overlap(a: npt.NDArray[np.complex128], b: npt.NDArray[np.complex128]) -> float:
    """`|<a|b>|` under the c-product, with both sides c-normalized.

    The c-product (bilinear, NOT conjugated) is the ECS-correct pairing -- the
    same convention `qscat.core`'s cross sections use. A conjugated dot would
    weight the exponentially growing ECS tail instead of cancelling it, and on a
    rotated grid that is not a small difference.

    **This can exceed 1, and that is not a bug.** The c-product is a bilinear
    form, so Cauchy-Schwarz does not bound it. With both states c-normalized the
    value is exactly `|c(a, b)|`, and the inflation over the Hermitian intuition
    is `1/sqrt(rho_a rho_b)`, where::

        rho = |c(psi, psi)| / ||psi||^2

    measures how close to REAL-VALUED (up to one global phase) a state is: 1 for
    a real vector, falling as the state acquires an internal phase profile.
    Measured on N2's broad anion resonances the six clean identifications score
    1.02 to 1.19, rising monotonically with `Gamma` as `rho` falls 0.66 to 0.42;
    H2+'s narrow Rydberg resonances stay at 0.87-0.99 because they are nearly
    real.

    `rho` is NOT `real_weight`, and substituting one for the other is wrong.
    N2's poles have `rho` of 0.42-0.66 while sitting 96 % inside the unscaled
    region -- well localized, and genuinely complex where they live, because a
    broad resonance is. H2+'s high Rydberg poles have both collapse together,
    for the different reason that their orbitals leave the box. `rho` explains
    the inflation; `real_weight` catches the escape.

    **Do NOT "fix" this by dividing by the Euclidean norms.** `|c(a,b)| /
    (||a|| ||b||)` is bounded by 1 and is WRONG here, which was measured rather
    than argued. Its denominator weights the exponentially growing ECS tail --
    reintroducing exactly the contamination the c-product's numerator cancels,
    which is the same error as using `vdot`. On H2+ window 0 it collapses to
    0.03, 0.008 and 0.006 for three states whose node counts identify them
    unambiguously as `Ry_14..Ry_16 v=1`, and it re-ranks all three onto the
    wrong (compact) partner. It penalizes diffuse states for being diffuse.

    The right response to a state whose `real_weight` is small is not a
    different overlap: it is to report the localization and treat the state as
    box-limited. `pair_by_overlap` does that.

    Returns `0.0` when either vector is (numerically) self-orthogonal, which is
    the honest answer: the normalization it would need does not exist.
    """
    na = complex(c_product(a, a))
    nb = complex(c_product(b, b))
    if na == 0 or nb == 0:
        return 0.0
    return float(abs(complex(c_product(a, b)) / np.sqrt(na * nb)))


def real_weight(psi: npt.NDArray[np.complex128], tgrid: TensorGrid) -> float:
    """Fraction of `|psi|^2` inside the UNSCALED region of `tgrid`.

    How much of a state the grid actually holds, as opposed to pushes into the
    rotated tail. `qscat.core.lcp.ResonanceLevels` reports exactly this for 1-D
    levels; the 2-D path had no equivalent, and the gap hid a real defect.

    Distinct from the `rho = |c(psi,psi)| / ||psi||^2` that inflates `overlap`
    above 1 (see there): `rho` measures how close to real-valued a state is,
    this measures where it lives. N2's poles sit at `real_weight` 0.96 with
    `rho` 0.42-0.66 -- localized, and genuinely complex where they live. The two
    move together only when the cause is a tail.

    **The overlap cannot see this and is not supposed to.** The c-product
    cancels the rotated tail by construction, so a state that is 97 % tail still
    pairs at 0.99 with the BO product it genuinely is -- a correct answer to the
    question asked, and a badly misleading summary of the state. Measured on H2+
    window 0, `real_weight` falls 0.998 -> 0.68 -> 0.29 -> 0.12 -> 0.031 ->
    0.008 -> 0.0003 as the Rydberg series climbs toward its threshold: those
    orbitals are simply larger than the 300-bohr box, and their overlaps stay at
    0.99 throughout.

    Returns a value in `[0, 1]`. Near 1 the grid holds the state; small values
    mean the answer is about the box, not the physics.
    """
    p2 = np.abs(np.asarray(psi, dtype=np.complex128)) ** 2
    total = float(p2.sum())
    if total == 0.0:
        return 0.0
    return float(p2[tgrid.real_mask()].sum() / total)


def pair_by_overlap(
    pole_energy: complex,
    pole_state: npt.NDArray[np.complex128],
    basis: BoBasis,
    thresholds: npt.ArrayLike | None = None,
    *,
    no_partner: float = NO_PARTNER,
    weak: float = WEAK,
    mixed_ratio: float = MIXED_RATIO,
    max_shift_ev: float = MAX_SHIFT_EV,
    n_eff_max: float | None = None,
    basis_complete: bool | None = None,
    localization: float | None = None,
    min_localization: float = MIN_LOCALIZATION,
) -> OverlapPair:
    """Pair one exact pole to the BO level it most resembles, with checks.

    Parameters
    ----------
    pole_energy : complex
        The pole, `E_r - i Gamma/2`. Only `Re E` enters the shift.
    pole_state : ndarray of complex
        Its eigenvector, flat over the same `(r, R)` layout as `basis`.
    basis : BoBasis
        Reference states, from `qscat.core.bo.bo_basis`. Must carry states.
    thresholds : array_like of float, optional
        Channel thresholds for the closed-channel admissibility check. **Omit at
        your peril**: without it a poor overlap cannot be attributed to the
        state rather than to the basis, and every such pole is reported
        `basis-limited` rather than `spurious`.
    no_partner, weak, mixed_ratio, max_shift_ev : float, optional
        Verdict thresholds; see the module docstring on calibration.
    n_eff_max : float, optional
        Passed to the admissibility check to cut the Rydberg series short of an
        accumulation region.
    basis_complete : bool, optional
        Assert directly that the basis holds every state the BO picture admits
        here, overriding the `thresholds` computation. This exists for the
        NEUTRAL case, where there is no Rydberg series to count: an anion
        resonance curve has one electronic state and the basis is complete once
        all its vibrational levels are built. **It is a claim, not a check** --
        asserting it wrongly converts every `basis-limited` verdict into
        `spurious`, which is the one verdict that rejects a pole.
    localization : float, optional
        This pole's `real_weight`. Supply it: the overlap is blind to a state
        that has left the box, by design, so without this a `box-limited` pole
        is reported as a clean identification.
    min_localization : float, optional
        Below this `real_weight`, the verdict is `box-limited`.

    Returns
    -------
    OverlapPair

    Raises
    ------
    ValueError
        If `basis` carries no states (built with `with_states=False`).
    """
    if not basis.has_states:
        raise ValueError(
            "pair_by_overlap needs a basis carrying product states: rebuild the "
            "ElectronicCurves with with_states=True"
        )

    ranked = sorted(
        ((overlap(s.psi, pole_state), key) for key, s in basis.items()), key=lambda t: -t[0]
    )
    best_v, best_k = ranked[0]
    second_v, second_k = ranked[1] if len(ranked) > 1 else (0.0, None)
    e_r = float(np.real(pole_energy))
    shift = (e_r - float(np.real(basis[best_k].energy))) * HARTREE_TO_EV

    # Energy decides which reading a poor overlap supports. If every level that
    # could exist at this energy is in the basis, a poor overlap is a fact about
    # the STATE; otherwise it is a fact about the BASIS. Checked for any weak
    # identification, not only a vanishing one: a state whose true partner is
    # absent still scores moderately against a wrong partner -- one measured H2+
    # pole reached 0.21 against a level it is not, while the Ry_18 it actually
    # is had never been built.
    if basis_complete is not None:
        covered = basis_complete
    elif thresholds is not None:
        covered = basis_covers(e_r, thresholds, basis, n_eff_max=n_eff_max)
    else:
        covered = False

    if best_v < weak and not covered:
        verdict, level = "basis-limited", None
    elif best_v < no_partner:
        verdict, level = "spurious", None
    elif localization is not None and localization < min_localization:
        # Ordered ahead of the match-quality verdicts deliberately: a state the
        # grid does not hold can still pair beautifully, and saying "ok" about
        # it is the failure this verdict exists to prevent. The label is kept
        # (the identification may well be right) but nothing is quotable.
        verdict, level = "box-limited", best_k
    elif best_v < weak:
        verdict, level = "weak", best_k
    elif best_v > 0 and second_v / best_v > mixed_ratio:
        verdict, level = "mixed", best_k
    elif abs(shift) > max_shift_ev:
        verdict, level = "distant", best_k
    else:
        verdict, level = "ok", best_k

    return OverlapPair(
        pole_energy=e_r,
        level=level,
        overlap=best_v,
        second_level=second_k,
        second_overlap=second_v,
        shift_ev=shift if level is not None else float("nan"),
        verdict=verdict,
    )


def pair_one_to_one(
    pole_energy: npt.ArrayLike,
    level_energy: npt.ArrayLike,
    *,
    max_distance: float = 1.0e-3,
) -> dict[int, int]:
    """Match poles to levels ONE-TO-ONE by POSITION, minimising total displacement.

    Only real parts are compared; widths play no part in the matching.

    Returns `{pole index: level index}`, omitting any pair further apart than
    `max_distance` (default 1 mHa).

    **A cross-check, not the assignment.** `pair_by_overlap` is the primary
    criterion; where the two disagree, overlap wins. This is here because an
    independent pairing that agrees is evidence, and one that disagrees marks a
    row worth looking at by hand.

    **Why not nearest-neighbour.** Assigning each pole to whichever level is
    closest treats the poles independently, so two poles can claim one level
    while another goes unclaimed -- and it fails worst exactly where the physics
    is most interesting. On H2+ two BO levels 20 uHa apart correspond to poles
    154 uHa apart; nearest-neighbour hands both poles to one level and calls the
    second pairing ambiguous, when it is fully determined once the matching is
    required to be a bijection. A minimum-total-cost assignment uses the
    constraint the physics supplies -- distinct states pair with distinct levels.

    It is not infallible: where a whole neighbourhood is denser than the shifts,
    the global optimum can still permute members within it, so a pairing is only
    as trustworthy as the local spacing. Report the per-pair distance so that
    stays visible.
    """
    # Real parts, taken explicitly: both sides are routinely complex here
    # (`E_r - i Gamma/2`), and an implicit float cast would discard the width
    # with a warning rather than saying that only positions are being matched.
    pe = np.real(np.asarray(pole_energy)).astype(np.float64)
    le = np.real(np.asarray(level_energy)).astype(np.float64)
    if pe.size == 0 or le.size == 0:
        return {}
    cost = np.abs(pe[:, None] - le[None, :])
    rows, cols = linear_sum_assignment(cost)
    return {int(r): int(c) for r, c in zip(rows, cols, strict=True) if cost[r, c] <= max_distance}


def peak_positions(
    energy: npt.NDArray[np.float64],
    sigma: npt.NDArray[np.float64],
    *,
    prominence: float = 10.0,
) -> npt.NDArray[np.float64]:
    """Energies of prominent local maxima in a cross-section sweep.

    A local maximum counts when it exceeds both neighbours AND is `prominence`
    times the smaller of the two samples two steps out -- a deliberately crude
    test, because the sweeps this is used on are resonance spectra where a
    genuine peak towers over its surroundings and a sophisticated peak finder
    would mostly add tunable knobs.

    Returns an empty array for a sweep shorter than five samples.
    """
    e = np.asarray(energy, dtype=np.float64)
    s = np.asarray(sigma, dtype=np.float64)
    if e.size != s.size:
        raise ValueError(f"energy and sigma must have equal length, got {e.size} and {s.size}")
    if e.size < 5:
        return np.empty(0, dtype=np.float64)
    out = [
        e[i]
        for i in range(2, e.size - 2)
        if s[i] > s[i - 1] and s[i] > s[i + 1] and s[i] > prominence * min(s[i - 2], s[i + 2])
    ]
    return np.asarray(out, dtype=np.float64)


@dataclass(frozen=True)
class PeakAlignment:
    """How well a set of computed positions lands on a set of observed peaks."""

    n_marks: int
    n_peaks: int
    median_widths: float  # median distance to the nearest peak, in widths
    within_one_width: int
    distances: npt.NDArray[np.float64]  # per mark, in widths

    def __str__(self) -> str:
        return (
            f"{self.n_marks} marks vs {self.n_peaks} peaks: median "
            f"{self.median_widths:.1f} widths, within 1 width "
            f"{self.within_one_width}/{self.n_marks}"
        )


def peak_alignment(
    marks: npt.ArrayLike,
    peaks: npt.ArrayLike,
    *,
    width: float,
) -> PeakAlignment:
    """Distance from each computed position to the nearest observed peak.

    **Distances are quoted in units of a resonance width, and that is the
    point.** The energy axis is not the scale on which "lands on the peak" means
    anything -- a 3 uHa disagreement is negligible against a 1 meV-wide feature
    and fatal against a 0.02 meV-wide one. `width` is the FWHM of a typical
    resonance in the set being compared.

    Marks with no peak on the same side are not special-cased: every mark takes
    its distance to the nearest peak anywhere in `peaks`, so a mark far outside
    the peak range reports a large distance rather than being dropped.

    Raises `ValueError` if `width <= 0`, or if either array is empty -- an
    alignment of nothing against something is not a number worth returning.
    """
    m = np.atleast_1d(np.asarray(marks, dtype=np.float64))
    p = np.atleast_1d(np.asarray(peaks, dtype=np.float64))
    if width <= 0:
        raise ValueError(f"width must be positive, got {width}")
    if m.size == 0 or p.size == 0:
        raise ValueError(f"peak_alignment needs both marks ({m.size}) and peaks ({p.size})")
    d = np.abs(m[:, None] - p[None, :]).min(axis=1) / width
    return PeakAlignment(
        n_marks=int(m.size),
        n_peaks=int(p.size),
        median_widths=float(np.median(d)),
        within_one_width=int(np.sum(d <= 1.0)),
        distances=d,
    )
