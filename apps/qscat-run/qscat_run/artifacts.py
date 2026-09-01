"""Artifact writers for a `qscat-run` experiment: the cross-section
CSV/NPZ/PNG trio (TI and/or TD -- `methods: [ti, td]` overlays both on the
SAME `cross_section.png`, since their keys share disjoint `"ti:"`/`"td:"`
prefixes), the TD-only moment-resolved `cross_section_vs_time` NPZ/PNG, the
TD-only opt-in `correlations.npz`, TI/TD wavefunction-density snapshots
(the SAME writer for both -- a `WavefunctionSnapshot`'s `kind`/`label` are
the only thing that differs), the resolved config (for reproducibility),
and the run manifest (provenance).

`matplotlib.use("Agg")` is set at import time (before `pyplot`), so this
module never needs a display and is safe to import in CI/tests. Nothing here
calls `datetime.now()`/`time.time()` at import time or otherwise generates
its own timestamp -- the CLI captures it once and passes it in via
`timestamp`, keeping this module a pure function of its arguments.
"""

from __future__ import annotations

import csv
import dataclasses
import json
import os
import platform
import re
import subprocess
import warnings
from importlib import metadata
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import yaml

from qscat_run.config import ExperimentConfig
from qscat_run.reference import config_base_dir, load_reference
from qscat_run.runner import (
    EigenStates,
    ExperimentResult,
    ResonanceLevelsRun,
    ResonanceState,
    WavefunctionSnapshot,
)

__all__ = ["write_artifacts"]


# `git rev-parse` walks up from `cwd` to find the repo root, so pinning the
# probe to this file's own location (rather than inheriting the caller's
# working directory) finds the SHA wherever `qscat-run` is invoked from.
# Overridden in tests to simulate "no repo here".
_REPO_PROBE_DIR = Path(__file__).resolve().parent

_SHA_RE = re.compile(r"\A[0-9a-f]{40}\Z")


def _git_sha() -> str:
    """The commit SHA the run's artifacts belong to, or `"unknown"`.

    A RECORD of what produced the numbers, not an address for them: published
    artifacts are addressed by content digest, so a missing SHA costs
    traceability but breaks nothing. It therefore warns rather than raising --
    a local figure should not fail over forty missing bytes of metadata.

    Publishing is where provenance actually matters, and that is where the
    strictness lives instead: the publisher refuses a manifest whose `git_sha`
    is not a commit, so nothing reaches the store unciteable.
    """
    # Baked in at image build time: the Docker build context excludes `.git`
    # (see .dockerignore), so `git rev-parse` inside a container has no repo
    # to read. `docker/build.sh` and `docker/run.sh` pass the host's SHA as a
    # build arg.
    #
    # Only a real SHA is honoured. `ARG GIT_SHA=unknown` in the Dockerfile
    # means an image built without `--build-arg GIT_SHA` carries the literal
    # string "unknown" here; treating any non-empty value as authoritative
    # let that sentinel outrank the working fallback below, silently. That is
    # how three O2 sweeps recorded "unknown" from a machine whose repo was
    # readable the whole time.
    baked = os.environ.get("QSCAT_GIT_SHA", "").strip().lower()
    if _SHA_RE.match(baked):
        return baked
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
            cwd=_REPO_PROBE_DIR,
        )
        probed = out.stdout.strip().lower()
        if _SHA_RE.match(probed):
            return probed
    except Exception:
        pass
    if not os.environ.get("QSCAT_ALLOW_UNKNOWN_SHA", "").strip():
        warnings.warn(
            "cannot determine the commit SHA for this run; its manifest will "
            'record "unknown" and the artifacts will not be traceable to the '
            "code that produced them. Inside a container pass the host SHA as "
            "QSCAT_GIT_SHA (docker/build.sh and docker/run.sh do this via "
            "--build-arg GIT_SHA). Set QSCAT_ALLOW_UNKNOWN_SHA=1 to silence "
            "this when there is genuinely no provenance to record.",
            RuntimeWarning,
            stacklevel=2,
        )
    return "unknown"


def _qscat_version() -> str:
    try:
        return metadata.version("qscat")
    except metadata.PackageNotFoundError:
        return "unknown"


def _config_to_dict(cfg: ExperimentConfig) -> dict[str, Any]:
    return dataclasses.asdict(cfg)


def _write_cross_section_csv(
    path: Path, energies: npt.NDArray[np.float64], series: dict[str, npt.NDArray[np.float64]]
) -> None:
    keys = list(series)
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["energy", *keys])
        for i, e in enumerate(energies):
            w.writerow([e, *(float(series[k][i]) for k in keys)])


def _write_cross_section_npz(
    path: Path, energies: npt.NDArray[np.float64], series: dict[str, npt.NDArray[np.float64]]
) -> None:
    # mypy false-positive: unpacking a `dict[str, ndarray]` here makes it
    # (over-cautiously) check the values against `savez`'s unrelated
    # `allow_pickle: bool` keyword too, and complain about the mismatch.
    np.savez(path, energy=energies, **series)  # type: ignore[arg-type]


def _disambiguated_labels(label: str | None, keys: list[str]) -> dict[str, str]:
    """Map one reference spec's series keys to their legend labels.

    `label` (from `ReferenceSpec.label`) is used verbatim when the spec
    produced exactly one series (e.g. a single-channel reference). With
    MULTIPLE channels from the same reference, the shared label alone would
    put several identical legend entries on the plot, so each key's own
    channel suffix (`"ch0"`, `"ch1"`, ...) is appended to keep them
    distinguishable. `label is None` returns `{}` -- `_write_cross_section_png`
    then falls back to the series key itself, same as a computed curve.
    """
    if label is None:
        return {}
    if len(keys) <= 1:
        return dict.fromkeys(keys, label)
    return {k: f"{label} ({k.rsplit(':', 1)[-1]})" for k in keys}


def _write_cross_section_png(
    path: Path,
    energies: npt.NDArray[np.float64],
    series: dict[str, npt.NDArray[np.float64]],
    reference: dict[str, tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]] | None = None,
    reference_labels: dict[str, str] | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    for key, sigma in series.items():
        masked = np.where(sigma > 0.0, sigma, np.nan)
        ax.plot(energies, masked, "-", label=key)
    # A reference dataset keeps ITS OWN energy axis (never interpolated onto
    # `energies`), so it is overlaid against that axis directly -- dashed and
    # thinner so it reads as "someone else's data", visually distinct from
    # the computed curves above. `reference_labels` (built from
    # `ReferenceSpec.label` by `_disambiguated_labels`) takes precedence over
    # the raw series key when the config named the reference.
    for key, (r_energy, r_sigma) in (reference or {}).items():
        masked_r = np.where(r_sigma > 0.0, r_sigma, np.nan)
        label = (reference_labels or {}).get(key, key)
        ax.plot(r_energy, masked_r, "--", linewidth=1.0, alpha=0.8, label=label)
    ax.set_xlabel("E (Hartree)")
    ax.set_ylabel(r"$\sigma$ (bohr$^2$)")
    ax.set_yscale("log")
    ax.legend(fontsize="small", ncol=1 if len(series) <= 8 else 2)
    ax.grid(True, which="both", alpha=0.2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _write_reference(
    out_dir: Path, series: dict[str, tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]]
) -> None:
    """Reference datasets keep their OWN energy axis, so they get their own
    files rather than columns in `cross_section.csv` (whose rows are the run's
    energies). Interpolating published data onto our grid would fabricate
    values and present them as the reference."""
    if not series:
        return
    with (out_dir / "reference.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["series", "energy", "sigma"])
        for key, (energy, sigma) in series.items():
            for e, s in zip(energy, sigma, strict=True):
                w.writerow([key, float(e), float(s)])
    flat: dict[str, npt.NDArray[np.float64]] = {}
    for key, (energy, sigma) in series.items():
        flat[f"{key}:energy"] = energy
        flat[f"{key}:sigma"] = sigma
    np.savez(out_dir / "reference.npz", **flat)  # type: ignore[arg-type]


def _write_correlations_npz(path: Path, correlations: dict[str, npt.NDArray[Any]]) -> None:
    """The raw per-step series behind each TD extractor's transform (opt-in,
    `cfg.artifacts.correlations`): `"{label}:t"`/`"{label}:c"`
    (`TannorWeeks`/`Dirac`) or `"{label}:t"`/`"{label}:b"`/`"{label}:d"`
    (`Flux`) -- see `runner.ExperimentResult.correlations`'s docstring for
    the exact keying. Mixed float (`t`) and complex (`c`/`b`/`d`) dtypes
    across keys are fine -- `np.savez` stores each named array independently.
    """
    np.savez(path, **correlations)  # type: ignore[arg-type]


def _write_wavefunction_snapshot(out_dir: Path, wf: WavefunctionSnapshot) -> None:
    # NOTE: build the full filename with an f-string rather than
    # `stem.with_suffix(...)` -- `wf.label` (e.g. "E0.05") contains a literal
    # dot, and `Path.with_suffix` REPLACES whatever follows the last dot
    # rather than appending, which would silently truncate "E0.05" to "E0".
    npz_path = out_dir / f"psi_{wf.label}.npz"
    png_path = out_dir / f"psi_{wf.label}.png"
    arrays: dict[str, npt.NDArray[Any]] = {
        "rho_r": wf.rho_r,
        "rho_R": wf.rho_R,
        "r": wf.r,
        "R": wf.R,
    }
    if wf.psi is not None:
        arrays["psi"] = wf.psi  # full complex field (n_r, n_R), for qscat.viz
    np.savez(npz_path, **arrays)  # type: ignore[arg-type]

    fig, (ax_r, ax_R) = plt.subplots(1, 2, figsize=(10, 4))
    ax_r.plot(wf.r, wf.rho_r)
    ax_r.set_xlabel("r (bohr)")
    ax_r.set_ylabel(r"$|\Psi^{(+)}|^2$ (summed over R)")
    ax_r.set_title("electronic")
    ax_R.plot(wf.R, wf.rho_R)
    ax_R.set_xlabel("R (bohr)")
    ax_R.set_title("nuclear")
    fig.suptitle(f"{wf.kind} {wf.label}")
    fig.tight_layout()
    fig.savefig(png_path)
    plt.close(fig)

    if wf.psi is not None:
        _write_wavefunction_field_png(out_dir / f"psi_{wf.label}_field.png", wf)


def _write_wavefunction_field_png(png_path: Path, wf: WavefunctionSnapshot) -> None:
    """Domain-coloured (phase->hue, magnitude->brightness) render of the full
    complex Psi field on the real-region r x R block, via `qscat.viz`'s pure-numpy
    `complex_to_rgb` -- the phase-carrying view the density marginals discard."""
    from qscat.viz import complex_to_rgb

    assert wf.psi is not None
    rgb = complex_to_rgb(wf.psi)  # (n_r, n_R, 3)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.imshow(
        rgb,
        origin="lower",
        aspect="auto",
        extent=(float(wf.R[0]), float(wf.R[-1]), float(wf.r[0]), float(wf.r[-1])),
    )
    ax.set_xlabel("R (bohr)")
    ax.set_ylabel("r (bohr)")
    ax.set_title(rf"{wf.kind} {wf.label}: $\Psi(r, R)$ (phase=hue, |·|=brightness)")
    fig.tight_layout()
    fig.savefig(png_path)
    plt.close(fig)


def _write_eigenstates(out_dir: Path, es: EigenStates) -> None:
    """`eigenstates_{label}.npz` (energies + fields + axis) plus a png: each
    field's |psi|^2 offset to its own energy -- the levels-and-their-state-
    wavefunctions view (vibrational) or the scattering states at each collision
    energy (LCP)."""
    stem = f"eigenstates_{es.label.replace(':', '_')}"
    np.savez(out_dir / f"{stem}.npz", energies=es.energies, states=es.states, axis=es.axis)

    is_vib = es.kind == "vibrational"
    fig, ax = plt.subplots(figsize=(7, 5))
    dens = np.abs(es.states) ** 2  # (n, len(axis))
    spread = float(es.energies[-1] - es.energies[0]) if len(es.energies) > 1 else 1.0
    scale = 0.6 * spread / max(1, len(es.energies))
    for i, e in enumerate(es.energies):
        d = dens[i]
        peak = d.max() or 1.0
        ax.axhline(float(e), color="0.8", lw=0.6)
        ax.plot(es.axis, float(e) + scale * d / peak, label=(f"v={i}" if is_vib else f"E={e:g}"))
    ax.set_xlabel("R (bohr)")
    ylabel = r"$|\chi_v|^2$" if is_vib else r"$|\psi_{sc}|^2$"
    ax.set_ylabel(f"energy (Hartree) + {ylabel} (offset)")
    ax.set_title(f"{es.kind} ({es.label})")
    ax.legend(fontsize="small", ncol=2)
    fig.tight_layout()
    fig.savefig(out_dir / f"{stem}.png")
    plt.close(fig)


def _write_resonance_state(out_dir: Path, rs: ResonanceState) -> None:
    """`resonance_{label}.npz` (complex pole energy, width, R, electronic
    eigenfunction + axis) plus a png of |phi_res(r)|^2 and Re/Im."""
    stem = f"resonance_{rs.label.replace(':', '_')}"
    np.savez(
        out_dir / f"{stem}.npz",
        energy=np.array(rs.energy),
        width=np.array(rs.width),
        R=np.array(rs.R),
        state=rs.state,
        axis=rs.axis,
    )
    real = rs.axis <= rs.axis.max()  # plot the real region only (axis is real_points)
    r = rs.axis[real]
    phi = rs.state[real]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(r, np.abs(phi) ** 2, "k-", label=r"$|\phi_{res}|^2$")
    ax.plot(r, phi.real, "C0-", lw=0.8, alpha=0.7, label=r"Re $\phi$")
    ax.plot(r, phi.imag, "C1-", lw=0.8, alpha=0.7, label=r"Im $\phi$")
    ax.set_xlabel("r (bohr)")
    ax.set_title(
        f"resonance state ({rs.label}) at R={rs.R:.3g} bohr: "
        rf"$E_r$={rs.energy.real:.4g}, $\Gamma$={rs.width:.3g} Ha"
    )
    ax.legend(fontsize="small")
    fig.tight_layout()
    fig.savefig(out_dir / f"{stem}.png")
    plt.close(fig)


def _resonance_levels_real_mask(
    R_axis: npt.NDArray[np.float64], R0: float
) -> npt.NDArray[np.bool_]:
    """Boolean mask selecting the physical real-region nodes of a
    `ResonanceLevelsRun`'s `R_axis` (`R_axis <= R0`).

    `R_axis`/`Vd`/`Gamma` all cover the grid's FULL node set (real nodes plus
    the complex-rotated ECS tail) -- despite its name, `R_axis` is not
    restricted to the real region, it is simply real-*valued* (the tail
    nodes' real pre-image coordinate, not their complex ECS position). Past
    `R0`, `Vd` is `model.v0` analytically continued at a complex argument (an
    absorbing-boundary artifact, not physics) plotted against that real
    pre-image -- which would stretch the figure's x-axis out to the tail's
    far edge and squash the physically relevant curve into a sliver. Factored
    out as its own function so the boundary logic is unit-testable without
    rendering a figure.
    """
    return R_axis <= R0


def _resonance_levels_ylim(
    Vd: npt.NDArray[np.float64], re_e: npt.NDArray[np.float64]
) -> tuple[float, float]:
    """Physical y-limits for the levels figure: the well, not the wall.

    `V_d(R)` runs from a well bottom around -0.15 Ha up a repulsive wall that
    reaches ~200 Ha as `R -> 0`. Autoscaling on that range collapses the well,
    every level bar and `Gamma(R)` onto a single line at `y ~ 0`. So the frame
    is set from the physics instead: bottom just under the well minimum, top
    just over the highest `Re E_v` (the wall then simply leaves the top of the
    frame, which is what it should do). With no levels at all, the curve's own
    outermost real value -- the anion dissociation limit, `Vd` at the largest
    real `R` -- stands in for the top; that is where any level would have to
    lie below, and unlike a percentile it does not depend on how the grid
    happens to distribute its nodes across the wall.

    Factored out so the choice is unit-testable without rendering a figure.
    """
    well = float(Vd.min())
    top = float(re_e.max()) if re_e.size else float(Vd[-1])
    span = top - well
    if not np.isfinite(span) or span <= 0.0:
        span = max(abs(well), 1.0)  # degenerate curve -- keep the frame finite
    return well - 0.10 * span, top + 0.20 * span


def _write_resonance_levels(out_dir: Path, run: ResonanceLevelsRun) -> None:
    """`resonance_levels_{label}.{csv,npz,png}` -- the quasi-bound level table.

    The csv/npz carry the full per-level data (energies, widths, residuals,
    real-weight fractions, and the golden-rule comparator where available --
    `golden_rule` legitimately may be all-`nan`, comparator unavailable or
    disabled, or `nan` on individual levels a distance guard rejected; those
    are written as literal `nan` text/values, not coerced to `0.0`). The png
    draws each level as a horizontal bar across the REAL-REGION `V_d(R)`
    curve at its `Re E_v`, with bar thickness proportional to `Gamma_v` -- so
    a broad, short-lived level reads as a thick smear and a long-lived one as
    a hairline. `golden_rule` never feeds the png (only `Re E_v`/`Gamma_v`
    do), so a `nan` comparator value cannot distort or silently drop a level
    bar there. The frame is set from the physical range
    (`_resonance_levels_ylim`) rather than autoscaled, and `Gamma(R)` gets its
    own right-hand axis -- otherwise `V_d`'s ~200 Ha repulsive wall at small
    `R` flattens everything worth seeing onto one line.
    """
    stem = f"resonance_levels_{run.label.replace(':', '_')}"
    lv = run.levels

    with (out_dir / f"{stem}.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        # Atomic units throughout (energies/Gamma in hartree, v the
        # vibrational quantum number); Re_E0/Gamma_v_1 are the golden-rule
        # comparator's pole energy/width (nan where unavailable/rejected).
        w.writerow(["v", "Re_E", "Gamma_v", "residual", "real_weight", "Re_E0", "Gamma_v_1"])
        for v in range(lv.energies.size):
            # `float(...)` before `!r` -- numpy>=2's scalar repr is
            # `"np.float64(...)"`, not a bare float literal (NEP 51); casting
            # to a plain Python float keeps the csv a clean, round-trippable
            # number (including a literal `"nan"`, never a numpy wrapper).
            w.writerow(
                [
                    v,
                    f"{float(lv.energies[v].real)!r}",
                    f"{float(lv.widths[v])!r}",
                    f"{float(lv.residuals[v])!r}",
                    f"{float(lv.real_weight[v])!r}",
                    f"{float(lv.golden_rule[v].real)!r}",
                    f"{float(-2.0 * lv.golden_rule[v].imag)!r}",
                ]
            )

    np.savez(
        out_dir / f"{stem}.npz",
        energies=lv.energies,
        widths=lv.widths,
        states=lv.states,
        residuals=lv.residuals,
        real_weight=lv.real_weight,
        golden_rule=lv.golden_rule,
        R_axis=run.R_axis,
        Vd=run.Vd,
        Gamma=run.Gamma,
        R0=np.array(run.R0),
    )

    real = _resonance_levels_real_mask(run.R_axis, run.R0)
    R = run.R_axis[real]
    Vd = run.Vd[real].real
    Gamma = run.Gamma[real]

    fig_h = 4.5
    fig, ax = plt.subplots(figsize=(7.0, fig_h))
    (line_vd,) = ax.plot(R, Vd, label=r"$V_d(R)$", color="tab:blue")

    lo, hi = _resonance_levels_ylim(Vd, lv.energies.real)
    ax.set_ylim(lo, hi)

    # `Gamma(R)` gets its OWN axis: it is a width (typically 1e-15 to 1e-5 Ha
    # here), not an energy on `V_d`'s scale, and sharing the left axis pins it
    # to a flat line at the bottom of the frame.
    ax_g = ax.twinx()
    (line_g,) = ax_g.plot(R, Gamma, label=r"$\Gamma(R)$", color="tab:red", ls="--", lw=1.0)
    ax_g.set_ylabel(r"$\Gamma(R)$ [hartree]", color="tab:red")
    ax_g.tick_params(axis="y", labelcolor="tab:red", labelsize=8)
    ax_g.set_ylim(bottom=0.0)

    # Annotate greedily from the bottom up, skipping any label that would land
    # within ~10 points of the previous one: with 40 levels in a 0.01 Ha window
    # the labels otherwise overprint into an unreadable smear.
    min_sep = (hi - lo) * 10.0 / (0.78 * fig_h * 72.0)
    last_label = -np.inf
    for v in range(lv.energies.size):
        e_v = float(lv.energies[v].real)
        ax.axhline(e_v, color="k", lw=max(0.6, 400.0 * float(lv.widths[v])), alpha=0.55)
        if e_v - last_label >= min_sep:
            ax.annotate(
                rf"$\omega_{{{v}}}$",
                xy=(R[-1], e_v),
                xytext=(-24, 2),
                textcoords="offset points",
                fontsize=8,
            )
            last_label = e_v

    ax.set_xlabel(r"$R$ [bohr]")
    ax.set_ylabel("energy [hartree]")
    ax.set_title(f"quasi-bound levels -- {run.label}")
    ax.legend(handles=[line_vd, line_g], loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / f"{stem}.png", dpi=150)
    plt.close(fig)


def write_artifacts(
    result: ExperimentResult,
    cfg: ExperimentConfig,
    output_dir: str | Path,
    *,
    timestamp: str,
) -> None:
    """Write every artifact `cfg.artifacts` requests to `output_dir`.

    `timestamp` is an ISO-8601 string captured by the caller (the CLI) --
    this function is otherwise pure (same inputs -> same files, module-load
    time excepted for the `matplotlib.use("Agg")` side effect above).
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Reference datasets keep their own energy axis (see `_write_reference`'s
    # docstring) -- resolved once here, then both overlaid on the PNG and
    # written to their own reference.{csv,npz}, whichever/however many
    # artifacts are requested below.
    ref_base = config_base_dir(cfg.config_dir)
    reference_series: dict[str, tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]] = {}
    reference_labels: dict[str, str] = {}
    for ref in cfg.reference:
        series_i = load_reference(ref, ref_base)
        reference_series.update(series_i)
        reference_labels.update(_disambiguated_labels(ref.label, list(series_i)))

    if cfg.artifacts.cross_section and result.cross_sections:
        e, series = result.energies, result.cross_sections
        _write_cross_section_csv(out_dir / "cross_section.csv", e, series)
        _write_cross_section_npz(out_dir / "cross_section.npz", e, series)
        _write_cross_section_png(
            out_dir / "cross_section.png", e, series, reference_series, reference_labels
        )

    # Unconditional (not gated on `cfg.artifacts.cross_section`), by design:
    # a `reference:` block is its own independent artifact request, the same
    # way `result.wavefunctions`/`result.resonance_states` below are written
    # whenever the result carries them, regardless of any OTHER artifact
    # flag. `_write_reference` is a no-op when no `reference:` is configured,
    # so this costs nothing on every run that doesn't use the feature.
    _write_reference(out_dir, reference_series)

    cvt_spec = cfg.artifacts.cross_section_vs_time
    if cvt_spec is not None and cvt_spec.moments and result.cross_section_vs_time:
        e, cvt_series = result.energies, result.cross_section_vs_time
        # Reuses the cross-section NPZ/PNG writers verbatim -- they are
        # already generic over (energies, series) and need no TD-specific
        # logic: `cvt_series`'s keys already carry the "@t{t_i}" moment
        # suffix, so the PNG legend reads one curve per (series, moment).
        _write_cross_section_npz(out_dir / "cross_section_vs_time.npz", e, cvt_series)
        _write_cross_section_png(out_dir / "cross_section_vs_time.png", e, cvt_series)

    if cfg.artifacts.correlations and result.correlations:
        _write_correlations_npz(out_dir / "correlations.npz", result.correlations)

    if result.wavefunctions:
        wf_dir = out_dir / "wavefunction"
        wf_dir.mkdir(parents=True, exist_ok=True)
        for wf in result.wavefunctions:
            _write_wavefunction_snapshot(wf_dir, wf)

    if cfg.artifacts.eigenstates and result.eigenstates:
        es_dir = out_dir / "eigenstates"
        es_dir.mkdir(parents=True, exist_ok=True)
        for es in result.eigenstates:
            _write_eigenstates(es_dir, es)

    if result.resonance_states:
        rs_dir = out_dir / "resonance"
        rs_dir.mkdir(parents=True, exist_ok=True)
        for rs in result.resonance_states:
            _write_resonance_state(rs_dir, rs)

    for run in result.resonance_levels:
        _write_resonance_levels(out_dir, run)

    (out_dir / "config.resolved.yaml").write_text(
        yaml.safe_dump(_config_to_dict(result.resolved_cfg), sort_keys=False)
    )

    manifest = {
        "qscat_version": _qscat_version(),
        "git_sha": _git_sha(),
        "timestamp": timestamp,
        "backend": cfg.backend,
        "timings": result.timings,
        "platform": platform.platform(),
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
