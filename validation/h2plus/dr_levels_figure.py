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

**Measured result** (`report_peak_alignment`, distances in units of a resonance
width, median FWHM 2e-5 Ha -- the only scale on which "lands on the peak" means
anything):

| series | window 0 | window 1 | window 2 |
|---|---|---|---|
| exact poles | **0.2** | **0.2** | 13.9 |
| BO levels | 1.8 | 8.2 | 15.5 |

In windows 0 and 1 the exact poles reproduce the peak positions of a sweep they
were never fitted to, while the BO levels miss by 2-8 widths -- the
Born-Oppenheimer error, measured against data.

**Window 2 is the exception and it is not explained.** There both series are
comparably wrong (13.9 and 15.5 widths), so it is not the approximation
failing. It is not peak sparsity either: window 2 holds 13-15 prominent peaks
at normal spacing, and the poles are still 8.7 widths off when allowed to match
a peak in EITHER channel. The one thing distinguishing it is threshold
proximity -- its upper edge sits 8e-4 Ha below the `v=3` vibrational threshold,
the closest of the three windows (1.1e-3 and 1.8e-3 for windows 1 and 0) -- and
this repository's sigma_DR is independently known to misbehave within ~1e-3 Ha
of a threshold (docs/physics/h2plus-dr.md). Suggestive, not established.
Window 2's per-level numbers should not be quoted.

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

import numpy as np
import numpy.typing as npt
from qscat.core.grids import fem_grid_exp_tail

from validation.h2plus.config import full_grid
from validation.h2plus.exact_poles import EPS0, WINDOWS
from validation.h2plus.reference_sweep import (
    REFERENCE_MU,
    DrSweep,
    load,
    mu_matched_model,
)
from validation.h2plus.rydberg_levels import rydberg_levels

FIGURES = pathlib.Path(__file__).resolve().parents[2] / "docs" / "physics" / "figures"
LEVEL_CACHE = pathlib.Path(__file__).with_name("dr_levels_figure.levels.npz")

# Exact 2-D pole positions (absolute energy, Ha) from the r_max=300 campaign,
# keyed by window. Recorded rather than recomputed: each window is a ~10-30 min
# multi-shift 2-D solve (see exact_poles.py), far too slow for a figure script.
EXACT_POLES: dict[int, tuple[float, ...]] = {
    0: (
        -0.096041061,
        -0.094304949,
        -0.093680172,
        -0.093125072,
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
    ),
    2: (
        -0.078166690,
        -0.077449152,
        -0.076354825,
        -0.075901901,
        -0.075822041,
        -0.075172625,
        -0.074999717,
        -0.074808161,
        -0.074091044,
        -0.073397529,
        -0.072877765,
        -0.072593147,
        -0.072357456,
        -0.072040067,
    ),
}


def bo_levels(mu: float = REFERENCE_MU) -> npt.NDArray[np.float64]:
    """BO level energies (absolute, Ha) at the given reduced mass, flattened.

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
        return np.asarray(np.load(LEVEL_CACHE)["levels"], dtype=np.float64)
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
    e = res.energies.ravel()
    out = np.asarray(e[np.isfinite(e)], dtype=np.float64)
    np.savez(LEVEL_CACHE, levels=out)
    return out


def main() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    sweep = load()
    levels = bo_levels()

    fig, axes = plt.subplots(len(WINDOWS), 1, figsize=(9.0, 10.5))
    for w, (lo, hi) in enumerate(WINDOWS):
        ax = axes[w]
        m = (sweep.energy >= lo) & (sweep.energy <= hi)
        ax.semilogy(sweep.energy[m], sweep.sigma[m, 0], lw=1.0, color="C0", label=r"$DR_0$")
        ax.semilogy(sweep.energy[m], sweep.sigma[m, 1], lw=1.0, color="C1", label=r"$DR_1$")

        for e in levels - EPS0:
            if lo <= e <= hi:
                ax.axvline(e, color="0.45", ls="--", lw=0.9, zorder=0)
        for e_tot in EXACT_POLES.get(w, ()):
            e = e_tot - EPS0
            if lo <= e <= hi:
                ax.axvline(e, color="C3", ls="-", lw=0.9, alpha=0.8, zorder=1)

        ax.set_xlim(lo, hi)
        ax.set_ylabel(r"$\sigma_{DR}$ (bohr$^2$)")
        ax.set_title(f"window {w}: {lo}–{hi} Ha", fontsize=10)
        if w == 0:
            handles, labels = ax.get_legend_handles_labels()
            handles += [
                Line2D([], [], color="0.45", ls="--", lw=0.9),
                Line2D([], [], color="C3", ls="-", lw=0.9),
            ]
            labels += ["BO levels $\\omega_i^j$", "exact 2-D poles"]
            ax.legend(handles, labels, fontsize=8, loc="best")
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
    fig.tight_layout(rect=(0, 0.045, 1, 0.97))

    FIGURES.mkdir(parents=True, exist_ok=True)
    out = FIGURES / "h2p-dr-levels.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")

    n_in = sum(1 for e in levels - EPS0 for lo, hi in WINDOWS if lo <= e <= hi)
    n_poles = sum(len(v) for v in EXACT_POLES.values())
    print(f"BO levels inside the three windows: {n_in}; exact poles drawn: {n_poles}")
    report_peak_alignment(sweep, levels)


def _peaks(
    energy: npt.NDArray[np.float64], sigma: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    """Energies of prominent local maxima in a channel of the published sweep."""
    out = []
    for i in range(2, len(energy) - 2):
        if (
            sigma[i] > sigma[i - 1]
            and sigma[i] > sigma[i + 1]
            and sigma[i] > 10.0 * min(sigma[i - 2], sigma[i + 2])
        ):
            out.append(energy[i])
    return np.asarray(out, dtype=np.float64)


def report_peak_alignment(sweep: DrSweep, levels: npt.NDArray[np.float64]) -> None:
    """How close do the marks sit to the published peaks they should explain?

    The comparison the figure exists to make, stated numerically rather than
    left to the eye. Distances are quoted in units of a resonance width (median
    FWHM 2e-5 Ha), since that -- not the energy axis -- is the scale on which
    "lands on the peak" means anything.
    """
    fwhm = 2.0e-5
    peaks = _peaks(sweep.energy, sweep.sigma[:, 1])
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
            d = np.array([np.min(np.abs(p - x)) for x in m])
            print(
                f"  {label:<12} window {w}: {m.size:>2} marks, "
                f"median distance to nearest peak {np.median(d) / fwhm:>6.1f} widths, "
                f"within 1 width: {int(np.sum(d <= fwhm))}/{m.size}"
            )


if __name__ == "__main__":
    main()
