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

matplotlib is imported LAZILY inside `plot_cross_sections`, not at module
scope: it is an OPTIONAL dependency (the `qscat[plot]` extra), so merely
importing `qscat.core` -- which re-exports this function -- must not require
matplotlib. See `libs/qscat/pyproject.toml`'s `[project.optional-dependencies]`
and `tests/test_no_matplotlib_at_import.py`.
"""

from __future__ import annotations

import os
from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

__all__ = ["plot_cross_sections", "plot_resonance_levels"]

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

    Requires the optional `qscat[plot]` extra (matplotlib); raises
    `ModuleNotFoundError` with an actionable hint if it is not installed.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")  # non-interactive backend, set before pyplot import
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:  # pragma: no cover - trivial guard
        raise ModuleNotFoundError(
            "qscat.core.plot_cross_sections requires matplotlib. "
            "Install the plotting extra: pip install 'qscat[plot]'."
        ) from exc

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

    # Parse the (E_ref, sigma_ref) overlay once, not per channel.
    e_ref_arr = s_ref_arr = None
    if reference is not None:
        e_ref, sigma_ref = reference
        e_ref_arr = np.asarray(e_ref, dtype=np.float64)
        s_ref_arr = np.where((tmp := np.asarray(sigma_ref, dtype=np.float64)) > 0.0, tmp, np.nan)

    s_masked = np.where(s > 0.0, s, np.nan)
    for j, ch in enumerate(ch_labels):
        color = colors[j % len(colors)]
        ax.plot(e, s_masked[:, j], "-", color=color, label=f"v'={ch} (computed)")

        if e_ref_arr is not None and s_ref_arr is not None:
            ax.plot(
                e_ref_arr,
                s_ref_arr[:, j],
                "--",
                marker="o",
                markersize=2.5,
                linewidth=0.8,
                color=color,
                alpha=0.6,
                label=f"v'={ch} (reference)",
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


def plot_resonance_levels(
    levels: dict[str, npt.NDArray[np.complex128]],
    *,
    path: _PathLike,
    curves: dict[str, tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]] | None = None,
    band: tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]
    | None = None,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    baseline: str | None = None,
    pairing: dict[str, Sequence[tuple[int, int]]] | None = None,
    title: str | None = None,
    xlabel: str = "R (bohr)",
    ylabel: str = "energy (Hartree)",
) -> None:
    """Plot complex resonance levels `omega_i = E_r - i*Gamma/2`, save to `path`.

    The left panel is the level diagram in the conventional form: any `curves`
    are drawn as solid lines against the internuclear coordinate, `band` shades
    a width envelope around one of them, and every level in `levels` becomes a
    horizontal dashed line at its `Re E`, one colour per series. The right panel
    plots each series' difference from `baseline`, in meV, against level index.

    Parameters
    ----------
    levels : dict of str to ndarray of complex
        One entry per series, e.g. `{"exact 2-D": ..., "BO / LCP": ...}`.
    curves : dict of str to (ndarray, ndarray), optional
        Background potential curves as `(x, y)` pairs, e.g. the neutral `V0(R)`
        and the resonance curve `E_res(R)`.
    band : tuple of (ndarray, ndarray, ndarray), optional
        `(x, centre, half_width)` shaded envelope — the `Gamma(R)/2` band around
        a resonance curve.
    baseline : str, optional
        Series the others are differenced against. Omit for a single panel.
    pairing : dict of str to sequence of (int, int), optional
        Explicit `(series index, baseline index)` pairs per non-baseline series.
        **Supply this whenever the difference panel is meant to be read as a
        physical shift.** Without it the series are paired level-by-level in
        ascending `Re E` and truncated to the shortest, which is only correct
        when both sets are complete and ordered alike -- and that assumption
        fails exactly where the physics is interesting. Measured on H2+: two BO
        levels 20 uHa apart correspond to exact poles 154 uHa apart, so index
        pairing crosses them and reports two shifts that belong to neither
        level. `qscat.core.assignment.pair_by_overlap` produces a pairing that
        is defensible; `pair_one_to_one` produces one that is at least a
        bijection. When `pairing` is given the series are used in the order
        supplied, NOT sorted.

    Notes
    -----
    No physics: it takes arrays and labels, so it serves any pair of level sets
    on any coordinate.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from qscat.units import HARTREE_TO_EV

    # A supplied pairing indexes the caller's own ordering; sorting would
    # silently re-map every pair it names.
    series = {
        k: (
            np.asarray(v, dtype=np.complex128)
            if pairing is not None
            else np.sort_complex(np.asarray(v, dtype=np.complex128))
        )
        for k, v in levels.items()
    }
    if not series:
        raise ValueError("plot_resonance_levels: no series given")
    if baseline is not None and baseline not in series:
        raise ValueError(f"baseline {baseline!r} is not one of {sorted(series)}")

    ncols = 1 if baseline is None else 2
    fig, axes = plt.subplots(1, ncols, figsize=(6.0 * ncols, 4.6), squeeze=False)
    ax = axes[0][0]

    if band is not None:
        x_b, centre, half = (np.asarray(a, dtype=np.float64) for a in band)
        ax.fill_between(x_b, centre - half, centre + half, color="gold", alpha=0.45, zorder=1)
    if curves is not None:
        for label, (x_c, y_c) in curves.items():
            ax.plot(np.asarray(x_c), np.asarray(y_c), linewidth=1.8, label=label, zorder=2)

    x0, x1 = xlim if xlim is not None else ax.get_xlim()
    for i, (label, w) in enumerate(series.items()):
        style = ("--", "-.")[i % 2]
        color = f"C{i + len(curves or {})}"
        for j, e in enumerate(w):
            ax.plot(
                [x0, x1],
                [e.real, e.real],
                style,
                color=color,
                linewidth=1.0,
                alpha=0.9,
                zorder=3,
                label=label if j == 0 else None,
            )
        # Index the levels once, off the first series. Labelling every series
        # would just overprint: series that differ by less than a pixel are the
        # normal case here, and saying so is the right panel's job.
        if i == 0:
            for j, e in enumerate(w):
                ax.annotate(
                    rf"$\omega_{{{j}}}$",
                    xy=(x0, e.real),
                    xytext=(3, 2),
                    textcoords="offset points",
                    fontsize="x-small",
                    color=color,
                )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.legend(fontsize="small", loc="lower right")
    ax.grid(True, alpha=0.2)
    if title:
        ax.set_title(title)

    if baseline is not None:
        ax2 = axes[0][1]
        base = series[baseline]
        for label, w in series.items():
            if label == baseline:
                continue
            if pairing is not None:
                pairs = list(pairing.get(label, ()))
                if not pairs:
                    continue
                left = np.array([w[i] for i, _ in pairs], dtype=np.complex128)
                right = np.array([base[j] for _, j in pairs], dtype=np.complex128)
                idx = np.arange(len(pairs))
            else:
                m = min(w.size, base.size)
                left, right = w[:m], base[:m]
                idx = np.arange(m)
            d_e = (left.real - right.real) * HARTREE_TO_EV * 1000.0
            d_g = (-2.0 * left.imag + 2.0 * right.imag) * HARTREE_TO_EV * 1000.0
            ax2.plot(idx, d_e, marker="o", label=rf"$\Delta E_r$ ({label} − {baseline})")
            ax2.plot(
                idx,
                d_g,
                marker="s",
                linestyle="--",
                label=rf"$\Delta\Gamma$ ({label} − {baseline})",
            )
        ax2.axhline(0.0, color="0.6", linewidth=0.8)
        ax2.set_xlabel("level index")
        ax2.set_ylabel("difference (meV)")
        ax2.legend(fontsize="small")
        ax2.grid(True, alpha=0.2)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
