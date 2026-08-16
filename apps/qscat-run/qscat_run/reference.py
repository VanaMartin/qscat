"""Published reference datasets overlaid on a computed cross section.

Reads a data file BY PATH, named in the config. This module deliberately does
NOT import `validation` -- `qscat_run` must stay independent of it (see
`tests/test_no_validation_import.py`), and naming the file in the config is
what keeps that true while still letting a run cite committed data.

A reference keeps its OWN energy axis. It is never interpolated onto the run's
energies: doing so would fabricate values and present them as someone else's
measurement, in the one figure whose purpose is honest external comparison.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt

__all__ = ["ReferenceSpec", "REFERENCE_FORMATS", "load_reference", "resolve_path"]

Series = dict[str, tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]]


@dataclass(frozen=True)
class ReferenceSpec:
    """One `reference:` entry: a data file, its format, and how to label it."""

    path: str
    format: str
    label: str | None = None
    channels: tuple[int, ...] | None = None


def resolve_path(spec: ReferenceSpec, base_dir: Path) -> Path:
    """Absolute path for `spec.path`, resolved against `base_dir` when relative."""
    p = Path(spec.path)
    return p if p.is_absolute() else (Path(base_dir) / p)


def _load_columns(path: Path, channels: tuple[int, ...] | None) -> Series:
    """Whitespace-delimited table: column 0 is energy (Hartree), columns 1..N
    are sigma (bohr^2) for successive final channels."""
    raw = np.loadtxt(path)
    if raw.ndim != 2 or raw.shape[1] < 2:
        raise ValueError(
            f"{path}: expected a 2-D table with an energy column and at least one "
            f"cross-section column, got shape {raw.shape}"
        )
    energy = raw[:, 0].astype(np.float64)
    n_channels = raw.shape[1] - 1
    wanted = tuple(range(n_channels)) if channels is None else channels
    bad = [c for c in wanted if c < 0 or c >= n_channels]
    if bad:
        raise ValueError(
            f"{path}: requested channel(s) {bad} but the file has {n_channels} "
            f"(valid 0..{n_channels - 1})"
        )
    return {f"ref:ve:ch{c}": (energy, raw[:, c + 1].astype(np.float64)) for c in wanted}


REFERENCE_FORMATS = {"houfek": _load_columns}


def load_reference(spec: ReferenceSpec, base_dir: Path) -> Series:
    """Load one reference dataset as `{series_key: (energies, sigma)}`."""
    loader = REFERENCE_FORMATS.get(spec.format)
    if loader is None:
        raise ValueError(
            f"unknown reference format {spec.format!r}; "
            f"choose one of {sorted(REFERENCE_FORMATS)}"
        )
    return loader(resolve_path(spec, base_dir), spec.channels)
