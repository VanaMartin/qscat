"""Numeric-output + visualization layer for the 2-D TD N2 VE cross section
(sub-project #7).

The user's primary goal for this whole sub-project is *observing the
transient anion form and decay from the correlation functions* -- the
numeric arrays (`t`, `c_{v'}(t)`, `norm(t)`, the snapshot densities,
`sigma(E)`) are the deliverable; the PNGs below are a thin
`matplotlib.use("Agg")` layer drawn on top of them, not a replacement for
them. This module does no physics of its own: it only serializes and plots
objects Tasks 1-5 already produced --
`projects.n2_2d_td_cross_section.td_propagation.PropagationResult` (the
correlation functions and density/norm snapshots) and a `sigma(E)` curve
(`td_cross_section.td_ve_cross_section_2d` or
`convergence.sigma_curve`), optionally overlaid on #6's exact TI oracle
`projects.n2_2d_cross_section.cross_section_2d.ve_cross_section_2d`.

Functions:
  * `save_numeric_outputs` -- the primary deliverable: a documented `.npz`.
  * `plot_snapshots` -- rho(R,t) / rho(r,t) at the snapshot times (incoming
    packet -> transient anion at the molecule -> decay), plus ||Psi(t)||.
  * `plot_correlation` -- |c_v'(t)| (and Re/Im) vs t, per channel -- the
    correlation build-up.
  * `plot_sigma_vs_ti` -- sigma_TD(E) overlaid on sigma_TI(E), distinguishing
    the spectral window where the Tannor-Weeks deconvolution is trustworthy
    (solid) from outside it, where dividing by a small `|eta_incident(E)|`
    amplifies residual noise (faded, never presented as signal).

    This function used to draw a THIRD region, a "finite-T resolution floor"
    below which the boomerang sub-features were said to be narrower than the
    propagation's `2*pi/T` and therefore unresolvable. That region was an
    artefact of the order-1 Crank-Nicolson propagator this module used at the
    time, not a property of the method: it was removed once order-3 Pade
    (`qscat.evolution.make_pade_stepper`) made the same energies agree with
    the exact TI oracle. Nothing here should re-introduce a "regions" concept
    without a measurement showing the region is real.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # non-interactive backend, set before pyplot import

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from qscat.dvr import TensorGrid

from projects.n2_2d_td_cross_section.td_propagation import PropagationResult

__all__ = ["save_numeric_outputs", "plot_snapshots", "plot_correlation", "plot_sigma_vs_ti"]

_PathLike = str | os.PathLike[str]


def save_numeric_outputs(
    result: PropagationResult,
    sigma_E: npt.ArrayLike,
    E_grid: npt.ArrayLike,
    path: _PathLike,
    *,
    dt: float,
    wp_in: dict[str, float],
    wp_out: dict[str, float],
) -> None:
    """Write the sampled dynamics + cross-section arrays to a documented `.npz`.

    This file IS the primary numeric deliverable of sub-project #7 -- the
    user's stated goal is observing the transient-anion formation/decay in
    the raw correlation functions and densities, not just the final
    `sigma(E)` curve. It is written to be *self-contained*: the propagation
    time step `dt` and the incoming/outgoing wavepacket parameters
    (`wp_in`/`wp_out`) that the Tannor-Weeks transform needs are stored
    alongside `c(t)`, so a pure-`.npz` consumer can reload `t`/`c` and re-run
    `td_cross_section.sigma_from_correlations(..., dt=..., wp_in=...,
    wp_out=...)` to reproduce `sigma_E` without any externally-remembered
    metadata.

    Keys and shapes (`n_t` = number of propagation samples, `n_ch` = number
    of outgoing channels `len(vprimes)`, `n_snap` = number of density
    snapshots, `n_E` = number of energies in `E_grid`):

      * `t`        (n_t,)          -- sample times `n*dt` (`result.t`).
      * `c`         (n_t, n_ch)     -- correlation functions `c_{v'}(t_n)`
                                       (`result.c`), recorded every step --
                                       the raw material of both the
                                       Tannor-Weeks transform and the
                                       "formation from the correlation
                                       functions" observation goal.
      * `norm`      (n_t,)          -- `||Psi(t_n)||`, Hermitian L2
                                       (`result.norm`); decays as the ECS
                                       contour absorbs outgoing flux.
      * `times`     (n_snap,)       -- the (coarse-cadence) times at which a
                                       density snapshot was recorded
                                       (`s.time` for `s in result.snapshots`).
      * `rho_R`     (n_snap, n_R)   -- nuclear density `sum_r |Psi|^2` at
                                       each snapshot time, stacked in the
                                       same order as `times`
                                       (`s.rho_R` for each snapshot; `n_R` is
                                       the full nuclear-axis length, real
                                       region + ECS tail).
      * `rho_r`     (n_snap, n_r)   -- electronic density `sum_R |Psi|^2` at
                                       each snapshot time, stacked likewise
                                       (`s.rho_r`; `n_r` the full electronic
                                       axis length).
      * `E_grid`    (n_E,)          -- the collision energies (Hartree) the
                                       cross section was evaluated at.
      * `sigma_E`   whatever shape `sigma_E` was given in (typically
                                       `(n_E, n_ch)`, matching
                                       `td_ve_cross_section_2d`'s array-`E`
                                       convention) -- `sigma_{v_init->v'}(E)`
                                       in bohr^2.
      * `dt`        scalar float64   -- the propagation time step used to
                                       generate `c(t)`; needed (with `t`) to
                                       re-run the energy transform.
      * `wp_in`     0-d string       -- the incoming-wavepacket parameter
                                       dict (keys `r0`/`p0`/`sigma`), JSON-
                                       encoded; reload with
                                       `json.loads(str(data["wp_in"]))`.
      * `wp_out`    0-d string       -- the outgoing-channel wavepacket
                                       parameter dict (keys
                                       `r0_out`/`p0_out`/`sigma_out`), JSON-
                                       encoded; reload with
                                       `json.loads(str(data["wp_out"]))`.

    Together `dt`, `wp_in`, and `wp_out` make the file self-contained: they
    are exactly the keyword arguments `sigma_from_correlations` needs beyond
    the reloaded `t`/`c`.

    `rho_R`/`rho_r` are NOT masked to the real (unscaled) region here -- they
    are saved exactly as `td_propagation.propagate` computed them (see that
    module's `_densities`); a consumer that wants only the physically
    meaningful real-region density should mask by the grid's
    `real_points <= R0` at load time (see `plot_snapshots` below, which does
    this when given a `tgrid`).
    """
    snaps = result.snapshots
    times = np.array([s.time for s in snaps], dtype=np.float64)
    rho_R = np.stack([s.rho_R for s in snaps]) if snaps else np.zeros((0, 0))
    rho_r = np.stack([s.rho_r for s in snaps]) if snaps else np.zeros((0, 0))

    np.savez(
        Path(path),
        t=result.t,
        c=result.c,
        norm=result.norm,
        times=times,
        rho_R=rho_R,
        rho_r=rho_r,
        E_grid=np.asarray(E_grid, dtype=np.float64),
        sigma_E=np.asarray(sigma_E, dtype=np.float64),
        dt=np.float64(dt),
        wp_in=np.asarray(json.dumps(wp_in)),
        wp_out=np.asarray(json.dumps(wp_out)),
    )


def plot_snapshots(
    result: PropagationResult,
    path: _PathLike,
    tgrid: TensorGrid | None = None,
) -> None:
    """rho(R,t) and rho(r,t) at the snapshot times, plus ||Psi(t)|| vs t.

    Three stacked panels: nuclear density snapshots (incoming packet ->
    transient anion at the molecule -> decay), electronic density snapshots,
    and the norm-decay curve with the snapshot times marked.

    If `tgrid` is given, densities are masked to each axis's real (unscaled)
    region and plotted against the physical `R`/`r` (bohr) grid, with the
    ECS pivot `R0` marked -- the complex tail carries outgoing flux, not
    probability density, so plotting it as density would be physically
    misleading (same masking discipline as
    `projects.n2_2d_cross_section.nuclear_density`). Without `tgrid` (e.g. a
    quick test with no grid handy), the full stored arrays are plotted
    against a plain grid-point index.
    """
    snaps = result.snapshots
    if not snaps:
        raise ValueError("plot_snapshots: result has no density snapshots to plot")

    if tgrid is not None:
        R_full = tgrid.grids[1].real_points
        r_full = tgrid.grids[0].real_points
        R_mask = R_full <= tgrid.grids[1].R0
        r_mask = r_full <= tgrid.grids[0].R0
        R_axis = R_full[R_mask]
        r_axis = r_full[r_mask]
        R0_R: float | None = tgrid.grids[1].R0
        R0_r: float | None = tgrid.grids[0].R0
        xlabel_R, xlabel_r = "R (bohr)", "r (bohr)"
    else:
        R_mask = slice(None)
        r_mask = slice(None)
        R_axis = np.arange(snaps[0].rho_R.size)
        r_axis = np.arange(snaps[0].rho_r.size)
        R0_R = R0_r = None
        xlabel_R, xlabel_r = "nuclear grid index", "electronic grid index"

    n_snap = len(snaps)
    cmap = plt.get_cmap("viridis")

    fig, axes = plt.subplots(3, 1, figsize=(7.5, 11.5))

    for i, snap in enumerate(snaps):
        color = cmap(i / max(1, n_snap - 1))
        axes[0].plot(R_axis, snap.rho_R[R_mask], color=color, label=f"t={snap.time:g}")
        axes[1].plot(r_axis, snap.rho_r[r_mask], color=color, label=f"t={snap.time:g}")
    if R0_R is not None:
        axes[0].axvline(R0_R, color="k", ls=":", lw=1, label="ECS pivot R0")
    if R0_r is not None:
        axes[1].axvline(R0_r, color="k", ls=":", lw=1, label="ECS pivot R0")

    axes[0].set_xlabel(xlabel_R)
    axes[0].set_ylabel("rho(R, t)")
    axes[0].set_title("Nuclear density snapshots -- incoming -> transient anion -> decay")
    axes[0].legend(fontsize=7, ncol=2)

    axes[1].set_xlabel(xlabel_r)
    axes[1].set_ylabel("rho(r, t)")
    axes[1].set_title("Electronic density snapshots")
    axes[1].legend(fontsize=7, ncol=2)

    axes[2].plot(result.t, result.norm, color="tab:blue", lw=1.5)
    for snap in snaps:
        axes[2].axvline(snap.time, color="gray", ls=":", lw=0.6)
    axes[2].set_xlabel("t (a.u.)")
    axes[2].set_ylabel("||Psi(t)||")
    axes[2].set_title("Norm decay: resonance formation + decay (ECS-absorbed flux)")

    fig.tight_layout()
    fig.savefig(Path(path), dpi=150)
    plt.close(fig)


def plot_correlation(
    result: PropagationResult,
    path: _PathLike,
    labels: list[str] | None = None,
) -> None:
    """`|c_v'(t)|` (top) and `Re/Im c_v'(t)` (bottom) vs t, one line per channel.

    `labels` names each of `result.c`'s `n_ch` columns (default `"v'=0"`,
    `"v'=1"`, ...); pass the actual `vprimes` list (as strings) for
    physically meaningful legends.
    """
    n_ch = result.c.shape[1]
    labels = labels if labels is not None else [f"v'={j}" for j in range(n_ch)]
    if len(labels) != n_ch:
        raise ValueError(f"plot_correlation: {len(labels)} labels for {n_ch} channels")

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 8), sharex=True)
    for j in range(n_ch):
        axes[0].plot(result.t, np.abs(result.c[:, j]), label=labels[j])
        axes[1].plot(result.t, result.c[:, j].real, label=f"Re c_{labels[j]}")
        axes[1].plot(result.t, result.c[:, j].imag, "--", label=f"Im c_{labels[j]}")

    axes[0].set_ylabel("|c_v'(t)|")
    axes[0].set_title("Correlation function build-up (formation observed here)")
    axes[0].legend(fontsize=8)

    axes[1].set_xlabel("t (a.u.)")
    axes[1].set_ylabel("Re / Im c_v'(t)")
    axes[1].legend(fontsize=7, ncol=2)

    fig.tight_layout()
    fig.savefig(Path(path), dpi=150)
    plt.close(fig)


def _as_channel_columns(a: npt.ArrayLike, n_e: int) -> npt.NDArray[np.float64]:
    """Reshape a 1-D `(n_e,)` or 2-D `(n_e, n_ch)` array-like to 2-D `(n_e, n_ch)`."""
    arr = np.asarray(a, dtype=np.float64)
    if arr.ndim == 1:
        if arr.shape[0] != n_e:
            raise ValueError(f"expected length {n_e}, got shape {arr.shape}")
        return arr.reshape(n_e, 1)
    if arr.shape[0] != n_e:
        raise ValueError(f"expected first axis length {n_e}, got shape {arr.shape}")
    return arr


def plot_sigma_vs_ti(
    E_grid: npt.ArrayLike,
    sigma_td: npt.ArrayLike,
    sigma_ti: npt.ArrayLike,
    usable: tuple[float, float],
    path: _PathLike,
) -> None:
    """sigma_TD(E) overlaid on the exact sigma_TI(E), inside vs outside the window.

    `sigma_td`/`sigma_ti` are `(n_E,)` or `(n_E, n_ch)` array-likes aligned
    with `E_grid`; `usable = (E_lo, E_hi)` (e.g. from
    `convergence.usable_window`) is the sub-interval where
    `|eta_incident(E)|` is large enough that the Tannor-Weeks deconvolution
    is trustworthy (see `convergence.usable_window`'s docstring). Points
    inside it are drawn solid; points outside are faded and dashed, since
    dividing by a small `|eta_incident|` there amplifies residual
    truncation/discretization noise into something that looks like signal.

    **A third region used to be drawn here and was wrong.** It marked
    `E < 0.13` Ha as "finite-T unresolved", on the reasoning that the
    boomerang sub-features (~0.004 Ha wide) were narrower than the
    propagation's energy resolution `2*pi/T ~ 0.0042` Ha, and cited measured
    ratios (5.7 at E=0.09, 0.37 at E=0.11) as evidence of a limit that "a
    longer T would sharpen" but re-tuning could not fix. Those ratios were
    produced by the order-1 Crank-Nicolson propagator this module used at the
    time, which under-converges badly over a long propagation. With order-3
    Pade (`qscat.evolution.make_pade_stepper`, now the default) the same
    energies track the exact TI oracle, so the region described a solver
    artefact as if it were physics. It is gone, along with the
    `validated_anchors` mechanism that existed only to except two energies
    from it.

    Do not re-introduce a regions concept without a measurement showing the
    region is real at the propagator actually in use.

    No sigma value is altered by this function -- only how each point is
    drawn and labeled.
    """
    E = np.asarray(E_grid, dtype=np.float64)
    order = np.argsort(E)
    E = E[order]
    td = _as_channel_columns(sigma_td, len(E_grid))[order]
    ti = _as_channel_columns(sigma_ti, len(E_grid))[order]
    e_lo, e_hi = usable
    n_ch = ti.shape[1]

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.axvspan(e_lo, e_hi, color="tab:green", alpha=0.12, label="usable window", zorder=0)

    below = E < e_lo
    above = E > e_hi
    inside = (E >= e_lo) & (E <= e_hi)

    for j in range(n_ch):
        color = f"C{j}"
        ti_label = "sigma_TI (exact)" if n_ch == 1 else f"sigma_TI v'={j}"
        ax.plot(E, ti[:, j], "-", color=color, lw=2, alpha=0.9, label=ti_label)

        td_label = "sigma_TD" if n_ch == 1 else f"sigma_TD v'={j}"
        if np.any(inside):
            ax.plot(E[inside], td[inside, j], "o-", color=color, ms=4, label=td_label)

        noise_label = "sigma_TD outside window (eta-deconvolution noise)"
        first_outside = True
        for mask in (below, above):
            if np.any(mask):
                ax.plot(
                    E[mask],
                    td[mask, j],
                    "o--",
                    color=color,
                    alpha=0.3,
                    ms=3,
                    label=noise_label if (j == 0 and first_outside) else None,
                )
                first_outside = False

    ax.set_xlabel("E (Hartree)")
    ax.set_ylabel("sigma (bohr^2)")
    ax.set_title("TD vs TI vibrational-excitation cross section")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(Path(path), dpi=150)
    plt.close(fig)
