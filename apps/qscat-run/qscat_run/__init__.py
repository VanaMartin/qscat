"""qscat-run — a config-driven CLI for 2-D electron-diatomic model experiments.

See `docs/superpowers/specs/2026-08-01-qscat-run-cli-design.md` for the
design. This package depends only on `qscat` (+ click/pyyaml/matplotlib/
numpy) — never on `validation.*` or `projects.*` (enforced by
`tests/test_no_validation_import.py`); per-molecule numerical decks in
`presets.py` are copied in as literal data, not imported.
"""

from __future__ import annotations

__all__: list[str] = []
