"""The published DR cross sections with computed resonance positions marked.

The figure this analysis exists to produce, in the form of the published
version (M. Vana, doctoral thesis, Charles University 2017, Fig. 4.7): sigma_DR
against incident electron energy across the three published windows, with the
Born-Oppenheimer levels `omega_i^j` drawn on it -- so the question "do the
computed resonance positions land on the cross-section peaks" is answered by
looking.

Two series are drawn over the published curves:

- **BO levels** (dashed), from `rydberg_levels` -- the quasi-bound levels of the
  Rydberg curves, which is what the published figure marks.
- **Exact 2-D poles** (solid), from `exact_poles` -- the approximation-free
  resonance positions, which the published figure does not have. Where the two
  differ is the Born-Oppenheimer error, and it is the point of the comparison.

**Measured result** (`report_peak_alignment`, median distance to the nearest
published peak in units of a resonance width -- FWHM 2e-5 Ha, the only scale on
which "lands on the peak" means anything -- with the count landing inside one
width):

| series | window 0 | window 1 | window 2 |
|---|---|---|---|
| exact poles | **0.2** (11/12) | **0.2** (14/18) | **0.3** (9/17) |
| BO levels | 0.8 (5/9) | 3.7 (2/11) | 30.4 (0/8) |

**All three windows agree**: the exact poles reproduce, to within a third of a
resonance width, the peak positions of a sweep they were never fitted to. The
BO levels degrade steadily across the three, which is the Born-Oppenheimer
error measured against data and is the point of the figure.

**Both rows reached those numbers by correction, and the route matters, because
each intermediate value was wrong in a way that looked like physics.**

The BO row was once 1.8 / 8.2 / 15.5, measured against a level set truncated to
3-4 marks per window by an electronic box too small to hold `Ry_5+` (see
`bo_levels`) -- and the missing levels were exactly the ones sitting on the
lower-window peaks.

The pole row's window 2 was once 13.9 widths, then 3.3, now 0.3. The first drop
came from seeding the campaign at every BO level instead of three per window;
the second from deleting four states that are not resonances at all
(`bo_overlap`, and the note above `EXACT_POLES`). **There is no window-2
anomaly.** Two explanations were offered for it along the way -- threshold
proximity, then residual under-seeding -- and neither survived; it was
contaminated input twice over.

**Reduced mass.** The levels here are computed at `REFERENCE_MU = 918.25`, the
value the published sweep was computed with, so that the marks and the curves
are on the same footing -- a level and a peak separated by the 2e-6 Ha the mass
correction is worth would be ~10% of a resonance width and could be misread as
a shift. This repository's shipped `H2P.mu` is the corrected `918.076` (Vana
2017 Table 1.2, Hvizdos 2016 Table 1.1, Hvizdos et al. 2018 Sec. II A), and the
full sweep is to be recomputed at that value later; every figure caption says
which mass produced it.

Run as::

    uv run python -m validation.h2plus.dr_levels_figure

Writes `docs/physics/figures/h2p-dr-levels.png`. The BO levels cost a few
minutes; exact poles are read from `exact_poles`' cached run if present and
simply omitted from the figure if not.
"""

from __future__ import annotations

import pathlib
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from qscat.core import peak_alignment, peak_positions
from qscat.core.grids import fem_grid_exp_tail

from validation.h2plus.config import full_grid
from validation.h2plus.exact_poles import EPS0, THRESHOLDS, WINDOWS
from validation.h2plus.reference_sweep import (
    REFERENCE_MU,
    DrSweep,
    load,
    mu_matched_model,
)
from validation.h2plus.rydberg_levels import rydberg_levels

FIGURES = pathlib.Path(__file__).resolve().parents[2] / "docs" / "physics" / "figures"
LEVEL_CACHE = pathlib.Path(__file__).with_name("dr_levels_figure.levels.npz")

# Curves at or above this index are the diffuse high-n Rydberg series whose
# levels are labelled along the TOP; the compact low-n curves label along the
# BOTTOM. Splitting on 6 reproduces the published figure's own arrangement and
# separates the two physical regimes -- see the annotate() call in `main`.
_HIGH_RYDBERG = 6


def _CURVE_COLOR(j: int) -> str:
    """One colour per Rydberg curve, as the published figure uses."""
    palette = ("#1f77b4", "#9467bd", "#d62728", "#ff7f0e", "#2ca02c", "#8c564b")
    return palette[j % len(palette)]


# Exact 2-D pole positions (absolute energy, Ha) from the r_max=300 campaign,
# with the four states `bo_overlap` identifies as NOT resonances removed:
# E = 0.004479, 0.021782, 0.022796 and 0.026065 have best BO overlaps of
# 6e-4, 7e-4, 7e-3 and 7e-3 where genuine states score 0.87-0.99, and every
# level energetically admissible at those energies is in the basis, so the
# poor overlap is a fact about the state rather than about the basis. They
# passed the two-angle ECS stability test that found them, which is the
# point: stability is necessary and not sufficient.
#
# Removing them is what makes window 2 agree with the others -- its median
# distance to a published peak falls from 3.3 resonance widths to 0.3. One
# of the four sat 0.4 widths from a peak by coincidence and was flattering
# the figure; the other three sat 6-32 widths away and were spoiling it.
#
# Keyed by window, and recorded rather than recomputed: each window is a
# ~10-30 min multi-shift 2-D solve, far too slow for a figure script.
EXACT_POLES: dict[int, tuple[float, ...]] = {
    0: (
        -0.097748890,
        -0.096041061,
        -0.094304949,
        -0.093680172,
        -0.092997199,
        -0.092096961,
        -0.091942982,
        -0.091397427,
        -0.091287909,
        -0.090858824,
        -0.090433412,
        -0.090091431,
        -0.089813206,
        -0.089586571,
        -0.089423775,
        -0.089344464,
        -0.089209371,
        -0.089074490,
    ),
    1: (
        -0.086823362,
        -0.085441899,
        -0.085236950,
        -0.084918035,
        -0.083784756,
        -0.083577996,
        -0.082825025,
        -0.082131676,
        -0.081593201,
        -0.081170577,
        -0.080850737,
        -0.080694597,
        -0.080508256,
        -0.080288534,
        -0.080097879,
        -0.079935442,
        -0.079877848,
        -0.079796340,
    ),
    2: (
        -0.078166690,
        -0.077449152,
        -0.076354825,
        -0.075901901,
        -0.075172625,
        -0.074999717,
        -0.074091044,
        -0.073397529,
        -0.072877765,
        -0.072593147,
        -0.072357456,
        -0.072040067,
        -0.071764636,
        -0.071534086,
        -0.071378514,
        -0.071340399,
        -0.071176496,
        -0.071036632,
    ),
}


@dataclass(frozen=True)
class BoLevels:
    """BO levels with the identity each one needs to be labelled `omega^j_i`."""

    energy: npt.NDArray[np.float64]  # (n,) absolute energy, Ha
    curve: npt.NDArray[np.int_]  # (n,) Rydberg curve j -- the SUPERscript
    vib: npt.NDArray[np.int_]  # (n,) vibrational level i -- the SUBscript

    def label(self, k: int) -> str:
        return rf"$\omega^{{{int(self.curve[k])}}}_{{{int(self.vib[k])}}}$"


def bo_levels(mu: float = REFERENCE_MU) -> BoLevels:
    """BO level energies (absolute, Ha) at the given reduced mass, with labels.

    Enumerated on a 300-bohr electronic box with twelve Rydberg curves, NOT
    the ~60-bohr proxy box with five that an earlier version used. That box
    cannot hold an `n_eff >= 6` orbital, so the `Ry_5+` levels were absent by
    construction -- it put 3 marks in window 0 where there are 9, and the
    missing ones are the `Ry_6..Ry_11 v=1` series that the published peaks in
    the lower half of that window belong to. See `exact_poles._bo_levels`.

    Cached to a git-ignored `.npz`: 818 nuclear points x a dense electronic
    eigensolve is ~15 minutes, and this figure gets restyled. Delete the file
    to recompute.
    """
    model = mu_matched_model() if mu == REFERENCE_MU else None
    if model is None:
        raise ValueError("only REFERENCE_MU is wired up; pass mu=REFERENCE_MU")
    if LEVEL_CACHE.exists():
        z = np.load(LEVEL_CACHE)
        return BoLevels(energy=z["energy"], curve=z["curve"], vib=z["vib"])
    g_r = fem_grid_exp_tail(
        [(10, 1.0), (10, 4.0), (16, 20.0), (20, 100.0), (20, 300.0)],
        angle_deg=5.0,
        quadrature=8,
        tail_n=25,
    )
    res = rydberg_levels(
        model,
        g_r,
        full_grid().grids[1],
        n_curves=12,
        n_vib=8,
        allow_partial=True,
    )
    j_idx, i_idx = np.meshgrid(
        np.arange(res.energies.shape[0]), np.arange(res.energies.shape[1]), indexing="ij"
    )
    keep = np.isfinite(res.energies)
    out = BoLevels(
        energy=np.asarray(res.energies[keep], dtype=np.float64),
        curve=np.asarray(j_idx[keep], dtype=np.int_),
        vib=np.asarray(i_idx[keep], dtype=np.int_),
    )
    np.savez(LEVEL_CACHE, energy=out.energy, curve=out.curve, vib=out.vib)
    return out


def main() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    sweep = load()
    levels = bo_levels()

    # Tall panels with generous gaps: each panel carries a row of rotated
    # `omega^j_i` labels above AND below it, and those need vertical room that
    # a default layout does not leave.
    fig, axes = plt.subplots(len(WINDOWS), 1, figsize=(10.0, 13.0))
    fig.subplots_adjust(hspace=0.55, top=0.94, bottom=0.10)
    for w, (lo, hi) in enumerate(WINDOWS):
        ax = axes[w]
        m = (sweep.energy >= lo) & (sweep.energy <= hi)
        ax.semilogy(sweep.energy[m], sweep.sigma[m, 0], lw=1.0, color="C0", label=r"$DR_0$")
        ax.semilogy(sweep.energy[m], sweep.sigma[m, 1], lw=1.0, color="C1", label=r"$DR_1$")

        # Exact poles go UNDERNEATH the BO levels and lighter, so that where the
        # two coincide the dashed BO line reads as sitting ON a pale solid one
        # and the overlap stays visible. Drawn first, low zorder.
        for e_tot in EXACT_POLES.get(w, ()):
            e = e_tot - EPS0
            if lo <= e <= hi:
                ax.axvline(e, color="C3", ls="-", lw=1.6, alpha=0.30, zorder=0)

        # BO levels on top, dashed, coloured by Rydberg curve and labelled
        # `omega^j_i` as the published figure does.
        ee = levels.energy - EPS0
        for k in np.flatnonzero((ee >= lo) & (ee <= hi)):
            e, j = float(ee[k]), int(levels.curve[k])
            ax.axvline(e, color=_CURVE_COLOR(j), ls="--", lw=0.9, alpha=0.9, zorder=3)
            # Top row for the high-n Rydberg series, bottom row for the compact
            # low-n levels. That is the published figure's own split, and it is
            # not merely cosmetic: those are the two regimes (adiabatic vs not),
            # and separating them is what keeps a near-degenerate pair such as
            # omega^9_1 / omega^3_3 legible at the same energy.
            # Upright, INSIDE the axes, just to the right of its own line.
            # High-n Rydberg labels ride high and compact low-n ones ride low --
            # the published figure's split, which keeps a near-degenerate pair
            # such as omega^9_1 / omega^3_3 legible at one energy, and which is
            # the physical division (adiabatic vs not) besides.
            top = j >= _HIGH_RYDBERG
            ax.annotate(
                levels.label(int(k)),
                xy=(e, 0.955 if top else 0.045),
                xycoords=("data", "axes fraction"),
                xytext=(3, 0),
                textcoords="offset points",
                ha="left",
                va="center",
                fontsize=7,
                color=_CURVE_COLOR(j),
                clip_on=True,
            )

        ax.set_xlim(lo, hi)
        # Headroom at top and bottom so the two label rows sit in empty space
        # rather than over the curves. Both ends, because the low-n labels ride
        # along the bottom.
        y0, y1 = ax.get_ylim()
        ax.set_ylim(y0 / 12.0, y1 * 25.0)
        ax.set_ylabel(r"$\sigma_{DR}$ (bohr$^2$)")

        # The window caption sits OUTSIDE the axes, as a title. The labels are
        # inside now, so nothing collides. The next ion vibrational threshold is
        # named rather than drawn: it lies beyond the right edge (the window
        # stops short of it deliberately), so a marker at the edge would claim
        # the threshold is AT the edge.
        thr = THRESHOLDS[w + 1] - EPS0
        ax.set_title(
            f"window {w}:  {lo}–{hi} Ha"
            rf"      (next threshold $v_{w + 1}$ at {thr:.4f} Ha, "
            rf"{(thr - hi) * 1e3:.1f} mHa beyond the right edge)",
            fontsize=9,
            color="0.25",
            pad=6,
        )
        if w == 0:
            handles, labels = ax.get_legend_handles_labels()
            handles += [
                Line2D([], [], color="0.45", ls="--", lw=0.9),
                Line2D([], [], color="C3", ls="-", lw=1.6, alpha=0.30),
            ]
            labels += [r"BO levels $\omega^j_i$", "exact 2-D poles"]
            # Top-left, but anchored just under the high-n label row so the two
            # do not overlap.
            ax.legend(
                handles,
                labels,
                fontsize=8,
                loc="upper left",
                bbox_to_anchor=(0.005, 0.90),
                framealpha=0.92,
            )
    axes[-1].set_xlabel("incident electron energy (Ha)")

    fig.suptitle(
        "H$_2^+$ dissociative recombination: published cross sections with computed "
        "resonance positions",
        fontsize=11,
    )
    fig.text(
        0.5,
        0.005,
        f"Cross sections: published sweep, rescaled by $1/2\\pi$ to this repository's "
        f"convention. BO levels computed at $\\mu$ = {REFERENCE_MU}, the reduced mass "
        f"the sweep itself used, so marks and curves are directly comparable. Exact "
        f"poles were computed at this repository's corrected $\\mu$ = 918.076 (Váňa "
        f"2017 Tab. 1.2; Hvizdoš 2016 Tab. 1.1; Hvizdoš 2018 §II A) and are NOT "
        f"mass-matched to the curves: the correction is worth ~2e-6 Ha, about 10% of "
        f"a resonance width, so a pole may sit that far off its peak for that reason "
        f"alone. The full sweep is to be recomputed at 918.076. Exact poles: "
        f"$r_{{max}}$ = 300 bohr, all three windows.",
        ha="center",
        va="bottom",
        fontsize=7.0,
        wrap=True,
    )
    # NOT tight_layout: it recomputes spacing from artist extents and collapses
    # the gaps the label rows need (the labels are annotations in axes-fraction
    # coordinates, which it does not account for). The explicit
    # `subplots_adjust` above is the layout.

    FIGURES.mkdir(parents=True, exist_ok=True)
    out = FIGURES / "h2p-dr-levels.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")

    n_in = sum(1 for e in levels.energy - EPS0 for lo, hi in WINDOWS if lo <= e <= hi)
    n_poles = sum(len(v) for v in EXACT_POLES.values())
    print(f"BO levels inside the three windows: {n_in}; exact poles drawn: {n_poles}")
    report_peak_alignment(sweep, levels.energy)


def report_peak_alignment(sweep: DrSweep, levels: npt.NDArray[np.float64]) -> None:
    """How close do the marks sit to the published peaks they should explain?

    The comparison the figure exists to make, stated numerically rather than
    left to the eye. Distances are quoted in units of a resonance width (median
    FWHM 2e-5 Ha), since that -- not the energy axis -- is the scale on which
    "lands on the peak" means anything.
    """
    fwhm = 2.0e-5
    peaks = peak_positions(sweep.energy, sweep.sigma[:, 1])
    print(f"\nprominent DR_1 peaks in the published sweep: {peaks.size}")
    for label, marks in (
        ("BO levels", levels - EPS0),
        ("exact poles", np.array([e - EPS0 for v in EXACT_POLES.values() for e in v])),
    ):
        for w, (lo, hi) in enumerate(WINDOWS):
            m = marks[(marks >= lo) & (marks <= hi)]
            p = peaks[(peaks >= lo) & (peaks <= hi)]
            if m.size == 0 or p.size == 0:
                continue
            print(f"  {label:<12} window {w}: {peak_alignment(m, p, width=fwhm)}")


if __name__ == "__main__":
    main()
