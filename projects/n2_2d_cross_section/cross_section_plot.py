"""Generic sigma(E) plotting utility -- no physics, no Houfek, no N2.

`projects/` must not import `validation/` (see `CLAUDE.md`'s "Global
Constraints"), so this module never reads any golden-data file itself: it
takes whatever reference curve the caller wants to overlay as a plain
`(E_ref, sigma_ref)` array pair. The N2-specific driver that actually reads
Houfek's `CSVE.V00.J00` lives in `validation/n2/ti_curve.py` and calls this
function; this module would be equally happy plotting a different
molecule's channels, or no reference at all.

`Agg` is selected before `pyplot` is imported (mirrors
`projects/n2_2d_td_cross_section/observation.py`) so this module never
requires a display, works headlessly in CI/tests, and never leaks a figure
window across test runs (every call closes its `Figure`).
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")  # non-interactive backend, set before pyplot import

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt

__all__ = ["plot_cross_sections"]

_PathLike = str | os.PathLike[str]


def plot_cross_sections(
    E_grid: npt.NDArray[np.float64],
    sigma: npt.NDArray[np.float64],
    *,
    channels: list[int] | None = None,
    reference: tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]] | None = None,
    thresholds: dict[int, float] | None = None,
    title: str | None = None,
    path: _PathLike,
) -> None:
    """Plot one sigma(E) curve per channel, save to `path` as a PNG.

    `E_grid` is `(M,)` Hartree; `sigma` is `(M, C)` bohr^2, column `c` is
    channel `channels[c]` (default `range(C)` when `channels` is `None`).

    `reference`, if given, is `(E_ref, sigma_ref)` with `sigma_ref` shaped
    `(N, C)` -- the SAME column order/count as `sigma` -- drawn as dashed,
    marked curves in the same color as their matching computed channel, so
    the eye pairs "solid vs. dashed of the same color" rather than reading
    a legend. This is a plain overlay: this function does no interpolation,
    resampling, or alignment between `E_grid` and `E_ref` -- the two curves
    are simply drawn on the same axes.

    `thresholds`, if given, is `{channel: threshold_energy_Ha}`; channels
    present in both `channels` and this dict get a dotted vertical line at
    their threshold, colored to match that channel's curve.

    Non-positive cross sections (e.g. a closed channel returning exactly
    `sigma=0` below its own threshold) are masked to `NaN` before plotting
    on the log-scaled y-axis: matplotlib's log locator silently drops them
    from the visible curve rather than raising, but leaving raw zeros in
    would still emit a `RuntimeWarning` at `log10(0)` -- this masks that
    cleanly rather than suppressing the warning.

    No physics lives here: this function does not know what a channel, a
    threshold, or a molecule is -- it plots whatever arrays it is given.
    """
    e = np.asarray(E_grid, dtype=np.float64)
    s = np.asarray(sigma, dtype=np.float64)
    n_channels = s.shape[1]
    ch_labels = channels if channels is not None else list(range(n_channels))
    if len(ch_labels) != n_channels:
        raise ValueError(
            f"channels has {len(ch_labels)} entries but sigma has {n_channels} columns"
        )

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    s_masked = np.where(s > 0.0, s, np.nan)
    for j, ch in enumerate(ch_labels):
        color = colors[j % len(colors)]
        ax.plot(e, s_masked[:, j], "-", color=color, label=f"v'={ch} (computed)")

        if reference is not None:
            e_ref, sigma_ref = reference
            e_ref_arr = np.asarray(e_ref, dtype=np.float64)
            s_ref_arr = np.asarray(sigma_ref, dtype=np.float64)
            s_ref_masked = np.where(s_ref_arr[:, j] > 0.0, s_ref_arr[:, j], np.nan)
            ax.plot(
                e_ref_arr, s_ref_masked, "--", marker="o", markersize=2.5,
                linewidth=0.8, color=color, alpha=0.6, label=f"v'={ch} (reference)",
            )

        if thresholds is not None and ch in thresholds:
            ax.axvline(thresholds[ch], color=color, linestyle=":", alpha=0.5, linewidth=1.0)

    ax.set_xlabel("E (Hartree)")
    ax.set_ylabel(r"$\sigma$ (bohr$^2$)")
    ax.set_yscale("log")
    if title:
        ax.set_title(title)
    ax.legend(fontsize="small", ncol=2)
    ax.grid(True, which="both", alpha=0.2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
