# qModeling

Quantum-mechanics research monorepo (QSCAT). Python-first, HPC-capable.

Development lifecycle: **design → Python toy model → validate → optimize hot paths
in Rust → promote to the `qscat` standard library.** Everything runs on CPU locally
and is containerizable.

## Layout
- `libs/qscat/` — the standard library (validated, reusable)
- `native/` — Rust kernels (PyO3/maturin)
- `projects/` — research toy models
- `validation/` — analytic benchmarks & golden data
- `reference/` — read-only old codebases (eMoScat, libXcuda)
- `docs/` — specs, physics notes, ADRs
- `docker/` — CPU test / runtime images

## Quickstart
```bash
uv sync
uv run pytest
```

See `CLAUDE.md` for the full operating manual.
