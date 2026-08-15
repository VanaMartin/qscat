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
import platform
import subprocess
from importlib import metadata
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")  # noqa: E402 -- must precede pyplot import

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import numpy.typing as npt  # noqa: E402
import yaml  # noqa: E402

from qscat_run.config import ExperimentConfig
from qscat_run.runner import EigenStates, ExperimentResult, ResonanceState, WavefunctionSnapshot

__all__ = ["write_artifacts"]


def _git_sha() -> str:
    """The current commit SHA, best-effort: `"unknown"` if `git` is absent,
    this isn't a repo, or anything else goes wrong.

    `cwd` is pinned to this file's own location (rather than inherited from
    the caller's working directory) so the SHA is found regardless of where
    `qscat-run` is invoked from -- `git rev-parse` walks up from `cwd` to
    find the repo root, so any path inside the repo works.
    """
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
            cwd=Path(__file__).resolve().parent,
        )
        return out.stdout.strip() or "unknown"
    except Exception:  # noqa: BLE001 -- provenance is best-effort, never fatal
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


def _write_cross_section_png(
    path: Path, energies: npt.NDArray[np.float64], series: dict[str, npt.NDArray[np.float64]]
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    for key, sigma in series.items():
        masked = np.where(sigma > 0.0, sigma, np.nan)
        ax.plot(energies, masked, "-", label=key)
    ax.set_xlabel("E (Hartree)")
    ax.set_ylabel(r"$\sigma$ (bohr$^2$)")
    ax.set_yscale("log")
    ax.legend(fontsize="small", ncol=1 if len(series) <= 8 else 2)
    ax.grid(True, which="both", alpha=0.2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


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

    if cfg.artifacts.cross_section and result.cross_sections:
        e, series = result.energies, result.cross_sections
        _write_cross_section_csv(out_dir / "cross_section.csv", e, series)
        _write_cross_section_npz(out_dir / "cross_section.npz", e, series)
        _write_cross_section_png(out_dir / "cross_section.png", e, series)

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
