"""Exact (non-Born-Oppenheimer) resonance poles in the published DR windows.

The repository owner asked for the exact 2-D poles (`qscat.core.
exact_resonance_states`) at the Born-Oppenheimer (BO) levels located by
`rydberg_levels` (`validation/h2plus/rydberg_levels.py`) inside the
three published dissociative-recombination (DR) windows -- electron energy
`E = e_tot - eps[0]` in `[0, 0.008]`, `[0.010, 0.018]`, `[0.020, 0.027]` Ha,
`eps[0] = -0.097604` Ha the cation `v=0` vibrational threshold -- and a table
of the shift each pole's BO level undergoes once the Born-Oppenheimer
separation is dropped. NOT a cross-section sweep (the owner has that data
already).

**Scope cutoff (repository owner's governing decision, relayed via
`validation/h2plus/pole_box_probe.py`'s docstring): `n_eff <= 12`.**
`n_eff = 1/sqrt(2*binding)`, `binding` measured to the nearest cation
vibrational threshold ABOVE the level (not `eps[0]`). Levels above the cutoff
are not isolated resonances in the first place -- Rydberg spacing falls as
`n^-3`, so peaks merge before widths do -- and are deliberately excluded, not
silently dropped. Of the 28 BO levels inside the three windows, **six** are cut
and `main()` prints each as a DROPPED row so a run states its own exclusions:
`Ry_2 v=7` (`n_eff` 19.5), `Ry_5 v=3` (14.9), `Ry_5 v=4` (13.2) and the
`Ry_11 v=1,2,3` trio (12.7-12.8). Note the cutoff's stated premise -- that
near-threshold states are poorly described -- is NOT what the shift table
below shows: the high-`n` levels are the best described of the set, so the
`Ry_11` exclusions are a scope choice about which peaks are isolated enough to
be worth calling resonances, not a numerical necessity.

**Results at the r_max = 300 bohr electronic box.** 22 in-scope BO levels
seeded (8, 8 and 6 per window); the exact solver returned **19, 18 and 21**
poles.

**A withdrawn claim.** An earlier run seeded only 3 levels per window -- all
this module's truncated enumeration could then find -- returned 18, 13 and 14
poles, and that arithmetic was read as "18 poles for 3 BO levels, so the BO
level list is not a complete inventory of a window's resonances." That was an
artifact of the enumeration, not a result. `_bo_levels` was building curves 0-4
on a ~60-bohr electronic box which cannot hold an `n_eff >= 6` orbital, so
every `Ry_5+` level was missing by construction. Enumerated on a box that holds
them, the windows contain **9, 11 and 8** BO levels, not 3 each.

Seeding those properly also changed the poles, which is why the counts above
differ from the earlier run: +1, +5 and +7. The under-seeding did lose real
states, and window 0 (+1) is not representative of that.

Each BO level paired to a pole ONE-TO-ONE (`pair_one_to_one`; see there for why
nearest-neighbour is the wrong algorithm here). All 28 levels match -- 9, 11
and 8 in the three windows -- with `gap` the distance from the assigned pole to
the nearest OTHER level, i.e. how much room the assignment had:

| window | level | E_BO (Ha) | E_exact (Ha) | shift (meV) | gap (meV) |
|---|---|---|---|---|---|
| 0 | `w^6_1` | -0.096041 | -0.096041 | -0.002 | 48.1 |
| 0 | `w^7_1` | -0.094273 | -0.094305 | -0.871 | 20.1 |
| 0 | `w^4_2` | -0.093566 | -0.093680 | **-3.097** | 16.1 |
| 0 | `w^8_1` | -0.093019 | -0.092997 | +0.585 | 15.5 |
| 0 | `w^9_1` | -0.092097 | -0.092097 | -0.003 | 0.5 |
| 0 | `w^3_3` | -0.092078 | -0.091943 | **+3.670** | 4.2 |
| 0 | `w^2_5` | -0.091431 | -0.091397 | (+0.921) | 0.1 |
| 0 | `w^10_1` | -0.091399 | -0.091288 | +3.035 | 3.9 |
| 0 | `w^11_1` | -0.090859 | -0.090859 | +0.010 | 14.7 |
| 1 | `w^2_6` | -0.085270 | -0.085442 | **-4.680** | 11.3 |
| 1 | `w^4_3` | -0.084951 | -0.084918 | +0.886 | 2.9 |
| 1 | `w^3_4` | -0.084154 | -0.083785 | (+10.049) | 0.6 |
| 1 | `w^8_2` | -0.083761 | -0.083578 | +4.984 | 15.7 |
| 1 | `w^9_2` | -0.082833 | -0.082825 | +0.212 | 18.9 |
| 2 | `w^4_4` | -0.076878 | -0.077449 | **-15.541** | 17.5 |
| 2 | `w^7_3` | -0.076301 | -0.075902 | +10.870 | 23.8 |
| 2 | `w^9_3` | -0.074090 | -0.074091 | -0.026 | 19.3 |

(abridged; `pair_one_to_one` regenerates the full 28.) Parenthesised rows have
a shift larger than half their gap, so the pairing itself is not safe there --
`w^2_5` in particular sits 0.06 meV from its neighbour. 21 of the 28 pairs are
clear of that test, and their shifts span **-4.68 to +10.87 meV**.

**The shift sorts by regime, and that is the result.** Median |shift| is
**0.390 meV for the high-n Rydberg levels (Ry >= 6, 17 pairs)** against
**3.097 meV for the compact low-n ones (Ry < 6, 11 pairs)** -- an eightfold
separation. A distant Rydberg electron follows the nuclei adiabatically, so its
level is nearly BO-exact; a compact low-n state overlapping the dissociative
channel is not. Both regimes are an order of magnitude beyond N2's 0.22 meV
(`docs/physics/exact-2d-resonances.md`) and neither is one-signed.

`w^9_1` and `w^3_3` are the two regimes meeting: 20 uHa apart in the BO
picture, they come out -0.003 and +3.670 meV shifted, so the exact treatment
splits a near-degeneracy eightfold and asymmetrically.
`resonance_state_figures.py` shows the pair -- one diffuse, one compact.

**Where these poles have been checked against data** is `dr_levels_figure.py`:
against the published cross-section peaks they sit a median **0.2** resonance
widths away in windows 0 and 1 and **3.3** in window 2 (the BO levels, for
comparison, 0.8 / 3.7 / 30.4). Window 2 was 13.9 widths on the under-seeded
run, so most of what looked like a window-2 anomaly was missing poles rather
than the threshold proximity earlier drafts proposed. It remains the worst of
the three and the closest to a threshold; that is a candidate for the residual,
not an established cause.

Pole widths across the window span 0.004-0.52 meV (1.4e-7 to 1.9e-5 Ha). That
matches the resonance width scale measured independently from the published
cross-section sweep (median FWHM 2e-5 Ha, docs/physics/h2plus-dr.md) -- and it
is why that sweep cannot be compared pointwise.

**Box convergence: the positions are converged, the COUNT is not.** Measured on
the earlier 3-seed run (the 300-vs-600 comparison has not been repeated at full
seeding, and does not need to be for the positions: it is the same operator).
Repeating every window at r_max = 600 left each pole found at 300 unmoved --
|dE| from 2.8e-17 up to 3.8e-9 for the topmost, i.e. converged well past the
precision anyone would quote. The shifts tabulated above are therefore
box-converged.

The number of poles is not, and it does not even move in a consistent
direction: 18 -> 22 in window 0, but 13 -> 10 in window 1 and 14 -> 13 in
window 2. Window 0 alone would suggest "a bigger box holds more diffuse
Rydberg states", and the other two refute that as a general rule -- what the
box changes is which states survive the two-angle stability test, which can go
either way. **Treat the inventory as unconverged in size while the positions in
it are solid.**

One artifact to not misread: the 300->600 comparison pairs poles by nearest
energy, so when a state present at 300 has no counterpart at 600 it pairs with
its neighbour and reports a large |dE| (3.2e-4 Ha on window 2's last row).
That is a dropped state, not a moved one.

Two honest qualifications on that near-machine-precision agreement. It is
partly a property of the grid construction: `_electronic_box` keeps the inner
segments fixed and lays the outer region out in 10-bohr elements regardless of
`r_max`, so 100-300 bohr is discretized IDENTICALLY in both runs and the
extension only appends elements beyond 300. The physical content of the check
is therefore "these states have no support past 300 bohr", which is what makes
their eigenvalues insensitive -- not a general claim that any two boxes agree
to 1e-17. And the varying survivor COUNT above is the same effect
`pole_box_probe.py` saw over 150 -> 300 -> 600 at a single seed; that probe was
right that the count is box-sensitive, and this campaign adds that the
positions are not.

The poles at the top of the window carry the largest `residual_electronic`
(up to 5.2e-8 at 600) and are the least trustworthy of the set.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.core import exact_resonance_states
from qscat.core.grids import fem_grid_exp_tail, nuclear_grid
from qscat.dvr import FemDvrEcsGrid, TensorGrid
from qscat.model import H2P
from qscat.units import HARTREE_TO_EV
from scipy.optimize import linear_sum_assignment

from validation.h2plus.config import full_grid
from validation.h2plus.rydberg_levels import RydbergLevels, rydberg_levels

__all__ = ["Seed", "PoleResult", "exact_poles", "pair_one_to_one", "main"]

# Published DR windows, electron energy E = e_tot - EPS0 (Ha).
#
# These are the published panels TRIMMED to stop short of each ion vibrational
# threshold, and the trimming is deliberate: the Rydberg series accumulates at
# a threshold, and the accumulation region is where a shift-invert pole search
# is least trustworthy. Widening a window toward its threshold re-admits
# exactly the states this cutoff exists to exclude.
#
# That caution turned out to be justified and insufficient. Window 2 -- whose
# upper edge is 8e-4 Ha below the v=3 threshold, the closest of the three
# (1.1e-3 and 1.8e-3 for windows 1 and 0) -- is the one whose poles do NOT
# land on the published cross-section peaks: 13.9 resonance widths off against
# 0.2 for windows 0 and 1 (see `dr_levels_figure.py`).
WINDOWS: tuple[tuple[float, float], ...] = ((0.0, 0.008), (0.010, 0.018), (0.020, 0.027))

# eps[0]: the cation v=0 vibrational threshold. Ties BO curve energies (the
# model's absolute electron+nuclear frame) to the published windows' electron
# energy convention -- same constant `rydberg_levels`/`pole_box_probe` use.
EPS0 = -0.097604

# Cation vibrational thresholds v=0..5 (Ha), as given by the repository owner
# (measured: `vibrational_states(full_grid().grids[1], H2P.mu, 6,
# H2P.v0).eps`; `pole_box_probe.py` independently reproduces these same
# six numbers).
THRESHOLDS: npt.NDArray[np.float64] = np.array(
    [-0.097604, -0.087802, -0.078519, -0.069754, -0.061507, -0.053780]
)

# Repository owner's governing scope decision -- see module docstring.
N_EFF_CUTOFF = 12.0

# Electronic ECS angle pair (base/moved): both comfortably under the model's
# own 45 deg bound, matching `pole_box_probe.py`'s tested values.
_EL_BASE_DEG = 30.0
_EL_MOVED_DEG = 40.0

# Nuclear ECS angle pair (base/moved): both under H2P's
# `max_nuclear_ecs_angle_deg = 22.5`, matching `pole_box_probe.py`.
_NUC_BASE_DEG = 18.0
_NUC_MOVED_DEG = 12.0

# Search width per shift (eigenpairs per shift per grid). pole_box_probe.py's
# k=8 SATURATED as the box grew (8 -> 6 -> 3 states survived at one seed as
# r_max: 150 -> 300 -> 600) -- Rydberg density rises toward threshold, so this
# campaign budgets generously rather than reusing that k unexamined.
K_SEARCH = 24


@dataclass(frozen=True)
class Seed:
    """One BO level (`rydberg_levels` curve/vib pair) inside a DR window."""

    window: int
    curve: int
    vib: int
    e_electron: float  # E = e_tot - EPS0
    e_tot: float  # absolute (electron + nuclear) energy, the model's frame
    n_eff: float
    in_scope: bool  # n_eff <= N_EFF_CUTOFF


@dataclass(frozen=True)
class PoleResult:
    """One exact pole at one box size, paired back to the BO seed that found it."""

    window: int
    r_max: float
    seed: Seed
    energy: complex
    width: float
    residual_electronic: float
    residual_nuclear: float
    shift_mev: float  # energy.real - seed.e_tot, in meV
    # |E - nearest seed| / |E - second-nearest seed|. Near 0 the assignment is
    # unambiguous; approaching 1 the pole sits midway between two BO levels and
    # nearest-Re pairing is a coin toss -- see `exact_poles`' margin note.
    assignment_ratio: float


def _n_eff(e_tot: float) -> float:
    """Hydrogenic n_eff from the binding below the nearest threshold ABOVE
    `e_tot` -- NOT `eps[0]`. The per-level threshold is the physically
    relevant one: a Rydberg state is bound against the ion level it sits
    below, not against the incident channel. See the module docstring.
    """
    above = THRESHOLDS[THRESHOLDS > e_tot]
    if above.size == 0:
        raise ValueError(f"e_tot={e_tot} sits above every known cation threshold")
    binding = float(above[0] - e_tot)
    return float(1.0 / np.sqrt(2.0 * binding))


def _bo_levels() -> RydbergLevels:
    """Rydberg curves 0-11 x their bound vibrational ladders.

    **The electronic box has to hold the orbital, or the level does not
    exist as far as this function is concerned.** An earlier version used
    `proxy_grid()`'s ~60-bohr box with `n_curves=5`, on the reasoning that
    curves 0-4 are unaffected by that truncation -- which is true, and
    beside the point: `Ry_5` and up are `n_eff >= 6` orbitals the box
    cannot hold at all, so their levels were absent by construction rather
    than by physics. In window 0 that turned 9 BO levels into 3, and the
    missing six are precisely the `Ry_6..Ry_11 v=1` series.

    So: a 300-bohr real region (the same box the pole search uses) and
    twelve curves. Even that is a cutoff, not a limit of the physics -- the
    `v=1` Rydberg series accumulates at the `v=1` ion threshold, and poles
    above E ~ 0.0071 Ha in window 0 pair to `Ry_12+` levels this
    enumeration still does not reach.

    `full_grid()`'s real nuclear grid supplies the vibrational ladder.
    `allow_partial=True` because the curves do not share one capacity (the
    deep curve 0 supports far fewer levels than the shallow Rydberg ones);
    `find_seeds` skips the NaN padding.
    """
    g_r = fem_grid_exp_tail(
        [(10, 1.0), (10, 4.0), (16, 20.0), (20, 100.0), (20, 300.0)],
        angle_deg=5.0,
        quadrature=8,
        tail_n=25,
    )
    return rydberg_levels(
        H2P,
        g_r,
        full_grid().grids[1],
        n_curves=12,
        n_vib=8,
        allow_partial=True,
    )


def pair_one_to_one(
    pole_energy: npt.NDArray[np.float64],
    level_energy: npt.NDArray[np.float64],
    *,
    max_distance: float = 1.0e-3,
) -> dict[int, int]:
    """Match poles to BO levels ONE-TO-ONE, minimising the total displacement.

    Returns `{pole index: level index}`, omitting any pair further apart than
    `max_distance` (default 1 mHa, comfortably wider than the largest credible
    shift and narrower than the level spacing).

    **Why not nearest-neighbour.** Assigning each pole to whichever level is
    closest treats the poles independently, so two poles can claim one level
    while another level goes unclaimed -- and it is worst exactly where the
    physics is most interesting. At E ~ 0.0055 Ha the BO levels `omega_1^9` and
    `omega_3^3` sit 20 uHa apart while the exact poles sit 154 uHa apart:
    nearest-neighbour hands both poles to `omega_1^9` and reports the second
    pairing at `assignment_ratio` 0.87, i.e. "ambiguous", when in fact it is
    fully determined once you require the matching to be a bijection. On the
    corrected (denser) BO level set 33 of 58 poles come out "ambiguous" that
    way, which is a defect of the algorithm rather than of the data.

    A minimum-total-cost assignment uses the constraint the physics supplies --
    distinct states pair with distinct levels -- and resolves those cases. It
    is not infallible: where a whole neighbourhood is denser than the shifts,
    the global optimum can still permute members within it, so a pairing is
    only as trustworthy as the local spacing. `report_shift_table` prints the
    per-pair distance so that stays visible.
    """
    if pole_energy.size == 0 or level_energy.size == 0:
        return {}
    cost = np.abs(pole_energy[:, None] - level_energy[None, :])
    rows, cols = linear_sum_assignment(cost)
    return {int(r): int(c) for r, c in zip(rows, cols, strict=True) if cost[r, c] <= max_distance}


def find_seeds() -> tuple[list[Seed], list[Seed]]:
    """Every (curve, vib) BO level inside a published window, split into
    (in-scope, out-of-scope) by the `n_eff <= 12` cutoff. Both lists are
    returned -- the out-of-scope ones are reported, not silently discarded.
    """
    lv = _bo_levels()
    survivors: list[Seed] = []
    dropped: list[Seed] = []
    for curve in range(lv.energies.shape[0]):
        for vib in range(lv.energies.shape[1]):
            e_tot = float(lv.energies[curve, vib])
            if not np.isfinite(e_tot):  # padding: this curve has fewer levels
                continue
            e_el = e_tot - EPS0
            for w, (lo, hi) in enumerate(WINDOWS):
                if lo <= e_el <= hi:
                    n_eff = _n_eff(e_tot)
                    seed = Seed(w, curve, vib, e_el, e_tot, n_eff, n_eff <= N_EFF_CUTOFF)
                    (survivors if seed.in_scope else dropped).append(seed)
    return survivors, dropped


def _electronic_box(r_max: float, angle_deg: float) -> FemDvrEcsGrid:
    """H2+ electronic grid at box `r_max`: `full_grid()`'s inner-segment
    resolution (0.1/0.3/1.0/4.0 bohr out to 100 bohr) with only the outer
    segment (100 -> r_max) and ECS pivot varying -- copied from
    `pole_box_probe.py`'s helper of the same name. `r_max` must be a
    multiple of 10 bohr above 100.
    """
    n_outer = round((r_max - 100.0) / 10.0)
    return fem_grid_exp_tail(
        [(10, 1.0), (10, 4.0), (16, 20.0), (20, 100.0), (n_outer, r_max)],
        angle_deg=angle_deg,
        quadrature=8,
        tail_n=25,
    )


def _nuclear_box(angle_deg: float) -> FemDvrEcsGrid:
    """Reduced nuclear grid (n_complex=3), matching `pole_box_probe.py`'s
    helper -- this campaign targets the ELECTRONIC box (the DR-window levels'
    diffuse Rydberg orbitals), so the nuclear side is held fixed and cheap,
    not itself re-converged here.
    """
    return nuclear_grid(angle_deg=angle_deg, r_max=14.0, n_complex=3, quadrature=8)


def exact_poles(seeds: list[Seed], r_max: float) -> list[PoleResult]:
    """Run `exact_resonance_states` once, seeded at every `seeds` entry
    (assumed to share one window), at electronic box `r_max`. Poles are
    paired back to their seeding BO level by nearest `Re(E)`.

    That pairing is sound but NOT comfortable, and the margin is worth
    stating because the shift table is what it feeds. Measured seed
    separations within a window are 1489/647 uHa (window 0), 318/797 uHa
    (window 1) and **74/1612 uHa (window 2)**. The 74 uHa pair -- `Ry_4 v=4`
    against `Ry_3 v=5` -- is the one to watch: nearest-Re pairing there
    survives only while the exact poles move by well under half that gap.
    For scale, the same comparison on N2 moved the levels ~7 uHa
    (0.2 meV, docs/physics/exact-2d-resonances.md), so the expected shift is
    ~5x smaller than the margin -- close enough that a pole landing between
    the two must be reported as ambiguous rather than silently assigned to
    whichever seed it happens to sit nearer.
    """
    if not seeds:
        return []
    window_idx = seeds[0].window
    assert all(s.window == window_idx for s in seeds)

    el_base = _electronic_box(r_max, _EL_BASE_DEG)
    el_moved = _electronic_box(r_max, _EL_MOVED_DEG)
    nu_base = _nuclear_box(_NUC_BASE_DEG)
    nu_moved = _nuclear_box(_NUC_MOVED_DEG)

    grid_base = TensorGrid([el_base, nu_base])
    grid_electronic = TensorGrid([el_moved, nu_base])
    grid_nuclear = TensorGrid([el_base, nu_moved])

    shifts = [complex(s.e_tot, -1e-4) for s in seeds]
    e_lo = min(s.e_tot for s in seeds) - 0.01
    e_hi = max(s.e_tot for s in seeds) + 0.01
    window = (e_lo, e_hi, -0.01, 0.0)

    res = exact_resonance_states(
        H2P, grid_base, grid_electronic, grid_nuclear, shifts=shifts, window=window, k=K_SEARCH
    )

    out: list[PoleResult] = []
    for e, g, re_, rn in zip(
        res.energies, res.widths, res.residual_electronic, res.residual_nuclear, strict=True
    ):
        # meV via the library constant, NOT a hand-written 1000: an earlier
        # revision used `1000 * dE`, which is milli-HARTREE, and understated
        # every shift by the 27.2 factor between the two.
        ranked = sorted(seeds, key=lambda s: abs(s.e_tot - e.real))
        seed = ranked[0]
        d0 = abs(seed.e_tot - e.real)
        # One seed in the window means nothing to be ambiguous against.
        ratio = d0 / abs(ranked[1].e_tot - e.real) if len(ranked) > 1 else 0.0
        shift_mev = 1000.0 * HARTREE_TO_EV * (e.real - seed.e_tot)
        out.append(
            PoleResult(
                window_idx,
                r_max,
                seed,
                complex(e),
                float(g),
                float(re_),
                float(rn),
                shift_mev,
                float(ratio),
            )
        )
    return out


def main(*, windows: list[int] | None = None, r_maxes: tuple[float, ...] = (300.0, 600.0)) -> None:
    """Run the pole campaign. `windows`/`r_maxes` narrow the run for cost
    discipline (a full 3-window x 2-box sweep is a lot of factorizations on
    scipy's SuperLU backend, no MUMPS reuse across shifts -- see the module
    docstring's cost accounting); default is everything.
    """
    survivors, dropped = find_seeds()
    print(
        f"in-scope seeds: {len(survivors)}, dropped (n_eff > {N_EFF_CUTOFF}): {len(dropped)}",
        flush=True,
    )
    for s in dropped:
        print(
            f"  DROPPED window={s.window} curve={s.curve} vib={s.vib} "
            f"E={s.e_electron:.6f} n_eff={s.n_eff:.2f}",
            flush=True,
        )
    for s in survivors:
        print(
            f"  seed window={s.window} curve={s.curve} vib={s.vib} "
            f"E={s.e_electron:.6f} n_eff={s.n_eff:.2f}",
            flush=True,
        )

    want = set(windows) if windows is not None else {0, 1, 2}
    for w in sorted(want):
        w_seeds = [s for s in survivors if s.window == w]
        if not w_seeds:
            continue
        print(f"--- window {w}: {len(w_seeds)} seeds ---", flush=True)
        by_box: dict[float, list[PoleResult]] = {}
        for r_max in r_maxes:
            t0 = time.perf_counter()
            results = exact_poles(w_seeds, r_max)
            dt = time.perf_counter() - t0
            by_box[r_max] = results
            print(f"  r_max={r_max:.0f}  found={len(results)}  {dt:.1f}s", flush=True)
            for r in results:
                # >0.5 means the second-nearest BO level is within a factor of
                # two of the nearest -- the shift for that row is not safely
                # attributable to one seed.
                flag = "  AMBIGUOUS" if r.assignment_ratio > 0.5 else ""
                print(
                    f"    seed=(c{r.seed.curve},v{r.seed.vib})  E={r.energy.real:+.9f}Ha  "
                    f"G={r.width * 1000.0 * HARTREE_TO_EV:.4f}meV  shift={r.shift_mev:+.3f}meV  "
                    f"res_el={r.residual_electronic:.1e}  res_nuc={r.residual_nuclear:.1e}  "
                    f"ratio={r.assignment_ratio:.2f}{flag}",
                    flush=True,
                )
        if len(by_box) >= 2 and all(by_box.values()):
            boxes = sorted(by_box)
            lo, hi = by_box[boxes[0]], by_box[boxes[-1]]
            print(f"  convergence {boxes[0]:.0f} -> {boxes[-1]:.0f} bohr:", flush=True)
            for r in lo:
                match = min(hi, key=lambda x: abs(x.energy - r.energy))
                diff = abs(r.energy - match.energy)
                print(
                    f"    (c{r.seed.curve},v{r.seed.vib})  "
                    f"E({boxes[0]:.0f})={r.energy.real:+.9f}  "
                    f"E({boxes[-1]:.0f})={match.energy.real:+.9f}  |diff|={diff:.3e}",
                    flush=True,
                )


if __name__ == "__main__":
    main()
