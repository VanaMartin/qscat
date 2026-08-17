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
silently dropped. Of the 11 BO levels inside the three windows, exactly two
are cut, both sitting just under a threshold: `Ry_2 v=7` (E = 0.017768 Ha,
`n_eff` = 19.5) and `Ry_2 v=9` (E = 0.026313 Ha, `n_eff` = 18.0). `main()`
prints both as DROPPED rows so a run states its own exclusions.

**Results at the r_max = 300 bohr electronic box.** Three BO levels seeded per
window; the exact solver returned 18, 13 and 14 poles.

**A withdrawn claim.** An earlier reading of that arithmetic said "18 poles for
3 BO levels, so the BO level list is not a complete inventory of a window's
resonances." That was an artifact of THIS module's own enumeration, not a
result. `_bo_levels` was building curves 0-4 on a ~60-bohr electronic box which
cannot hold an `n_eff >= 6` orbital, so every `Ry_5+` level was missing by
construction. Enumerated on a box that holds them (see `_bo_levels`), the
windows contain **9, 11 and 8** BO levels, not 3 each -- and the poles that
still look unpaired sit where the `v=1` Rydberg series accumulates and pair to
`Ry_12+` levels beyond even the widened cutoff. Every exact pole appears to have
a BO partner once the series is enumerated far enough.

What survives is better than the withdrawn claim, because it is about physics
rather than counting: **the shift depends on which regime the level is in.**
High-`n` Rydberg levels are nearly BO-exact -- `Ry_6..Ry_11` at `v=1,2,3` shift
by <= 0.5 meV -- while compact low-`n`, high-`v` levels shift by several meV
(`Ry_3 v=3` +3.60, `Ry_4 v=2` -3.15, `Ry_2 v=6` -4.78, `Ry_4 v=4` -15.6). A
distant Rydberg electron follows the nuclei adiabatically; a compact one
overlapping the dissociative channel does not. `resonance_state_figures.py`
shows the two regimes at the same energy, repelling.

Each BO level paired to the pole nearest it:

| window | BO level | E_BO (Ha) | E_exact (Ha) | shift (meV) | Gamma (meV) | ratio |
|---|---|---|---|---|---|---|
| 0 | `Ry_4 v=2` | -0.093564 | -0.093680 | **-3.161** | 0.0152 | 0.07 |
| 0 | `Ry_3 v=3` | -0.092075 | -0.092097 | **-0.598** | 0.1822 | 0.03 |
| 0 | `Ry_2 v=5` | -0.091428 | -0.091397 | **+0.832** | 0.1513 | 0.05 |
| 1 | `Ry_2 v=6` | -0.085266 | -0.085237 | **+0.790** | 0.1900 | 0.10 |
| 1 | `Ry_4 v=3` | -0.084948 | -0.084918 | **+0.815** | 0.6188 | 0.09 |
| 1 | `Ry_3 v=4` | -0.084151 | -0.083785 | +9.966 | 0.4014 | 0.31 |
| 2 | `Ry_2 v=8` | -0.075189 | -0.075173 | **+0.446** | 0.1568 | 0.01 |
| 2 | `Ry_4 v=4` | -0.076875 | -0.076355 | (+14.155) | 0.9012 | **1.17** |
| 2 | `Ry_3 v=5` | -0.076801 | -0.076355 | (+12.141) | 0.9012 | **0.86** |

Bold shifts are the ones worth quoting. They run **-3.2 to +0.8 meV**, an order
of magnitude larger than the N2 result in `docs/physics/exact-2d-resonances.md`
(0.22 meV) and, unlike N2's, **not one-signed** -- `Ry_4 v=2` and `Ry_3 v=3`
move down while the rest move up, so the Born-Oppenheimer correction here is
not a uniform lowering.

The last two rows are **not individually assignable** and their shifts are
parenthesised for that reason. `Ry_4 v=4` and `Ry_3 v=5` are the 74 uHa pair
this module's `exact_poles` docstring flags in advance: both come out nearest
to the SAME pole (-0.076355), and `Ry_4 v=4`'s `assignment_ratio` of 1.17
exceeds 1, meaning that pole is closer to the OTHER BO level than to its own.
`Ry_3 v=4`'s 0.31 is marginal rather than broken. Everything else is
unambiguous at <= 0.10, with angle residuals at or below 1e-17.

**Where these poles have been checked against data**, and where that check
fails, is in `dr_levels_figure.py`: in windows 0 and 1 they sit a median 0.2
resonance widths from the published cross-section peaks, but in window 2 they
are 13.9 widths off -- and so are the BO levels (15.5), which is why that is
read as a problem with the window rather than with the approximation. Window
2's per-level numbers should not be quoted.

Pole widths across the window span 0.004-0.52 meV (1.4e-7 to 1.9e-5 Ha). That
matches the resonance width scale measured independently from the published
cross-section sweep (median FWHM 2e-5 Ha, docs/physics/h2plus-dr.md) -- and it
is why that sweep cannot be compared pointwise.

**Box convergence: the positions are converged, the COUNT is not.** Repeating
every window at r_max = 600 leaves each pole found at 300 unmoved -- |dE| from
2.8e-17 up to 3.8e-9 for the topmost, i.e. converged well past the precision
anyone would quote. The shifts tabulated above are therefore box-converged.

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

from validation.h2plus.config import full_grid
from validation.h2plus.rydberg_levels import RydbergLevels, rydberg_levels

__all__ = ["Seed", "PoleResult", "exact_poles", "main"]

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
