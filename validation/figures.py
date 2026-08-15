"""The committed-figures directory (`docs/physics/figures/`), shared by every
validation driver that writes a figure — a single definition instead of the
per-driver `Path(__file__).resolve().parents[...] / "docs" / "physics" /
"figures"` copies."""

from __future__ import annotations

from pathlib import Path

__all__ = ["FIGURE_DIR"]

FIGURE_DIR = Path(__file__).resolve().parent.parent / "docs" / "physics" / "figures"
