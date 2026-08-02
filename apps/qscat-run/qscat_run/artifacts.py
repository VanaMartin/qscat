"""Artifact writers for a `qscat-run` experiment: the cross-section
CSV/NPZ/PNG trio, TI wavefunction-density snapshots, the resolved config
(for reproducibility), and the run manifest (provenance).

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
from qscat_run.runner import ExperimentResult, WavefunctionSnapshot

__all__ = ["write_artifacts"]


def _git_sha() -> str:
    """The current commit SHA, best-effort: `"unknown"` if `git` is absent,
    this isn't a repo, or anything else goes wrong."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
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


def _write_wavefunction_snapshot(out_dir: Path, wf: WavefunctionSnapshot) -> None:
    # NOTE: build the full filename with an f-string rather than
    # `stem.with_suffix(...)` -- `wf.label` (e.g. "E0.05") contains a literal
    # dot, and `Path.with_suffix` REPLACES whatever follows the last dot
    # rather than appending, which would silently truncate "E0.05" to "E0".
    npz_path = out_dir / f"psi_{wf.label}.npz"
    png_path = out_dir / f"psi_{wf.label}.png"
    np.savez(npz_path, rho_r=wf.rho_r, rho_R=wf.rho_R, r=wf.r, R=wf.R)

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

    if result.wavefunctions:
        wf_dir = out_dir / "wavefunction"
        wf_dir.mkdir(parents=True, exist_ok=True)
        for wf in result.wavefunctions:
            _write_wavefunction_snapshot(wf_dir, wf)

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
