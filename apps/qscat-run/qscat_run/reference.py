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

__all__ = [
    "REFERENCE_FORMATS",
    "ReferenceSpec",
    "bad_channels",
    "config_base_dir",
    "load_reference",
    "peek_n_channels",
    "resolve_path",
]

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


def config_base_dir(config_dir: str | None) -> Path:
    """The directory a `reference.path` (or any other config-relative path)
    resolves against: `config_dir` (the directory the config YAML itself
    lives in) when set, else the process's current working directory.

    Shared by `config.py` (validate-time) and `artifacts.py` (write-time) so
    the two base-directory computations can never drift apart.
    """
    return Path(config_dir) if config_dir else Path.cwd()


def peek_n_channels(path: Path) -> int:
    """The number of data columns in a `"houfek"`-format reference file,
    cheaply: reads only the first non-blank line rather than parsing the
    whole table with `np.loadtxt`.

    `validate_config` calls this at validate-time, which must stay fast and
    solve-free -- it must not pay for a full parse of a potentially large
    (e.g. 400-row) file just to bounds-check `channels`.
    """
    with path.open() as f:
        for line in f:
            fields = line.split()
            if fields:
                return len(fields) - 1
    raise ValueError(f"{path}: file is empty, cannot determine its column count")


def bad_channels(channels: tuple[int, ...] | None, n_channels: int) -> list[int]:
    """The entries of `channels` that are out of bounds for a file with
    `n_channels` data columns (empty when `channels` is `None` or every
    index is valid). Shared by `_load_columns` (load-time) and
    `validate_config` (validate-time) so the two bounds checks can never
    drift apart.
    """
    if channels is None:
        return []
    return [c for c in channels if c < 0 or c >= n_channels]


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
    bad = bad_channels(channels, n_channels)
    if bad:
        raise ValueError(
            f"{path}: requested channel(s) {bad} but the file has {n_channels} "
            f"(valid 0..{n_channels - 1})"
        )
    wanted = tuple(range(n_channels)) if channels is None else channels
    return {f"ref:ve:ch{c}": (energy, raw[:, c + 1].astype(np.float64)) for c in wanted}


REFERENCE_FORMATS = {"houfek": _load_columns}


def load_reference(spec: ReferenceSpec, base_dir: Path) -> Series:
    """Load one reference dataset as `{series_key: (energies, sigma)}`."""
    loader = REFERENCE_FORMATS.get(spec.format)
    if loader is None:
        raise ValueError(
            f"unknown reference format {spec.format!r}; choose one of {sorted(REFERENCE_FORMATS)}"
        )
    return loader(resolve_path(spec, base_dir), spec.channels)
