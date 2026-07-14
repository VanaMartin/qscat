---
name: qscat-conventions
description: Reference for qModeling/QSCAT shared conventions — atomic units, FEM-DVR-ECS notation, tolerance defaults, and standard-library layout. Consult when unsure how the project names or measures things.
---

# qscat-conventions

## Overview

Lookup reference for the vocabulary and defaults shared across qModeling —
not a process to follow, just facts to check.

## Units

- **Atomic units are the default** everywhere physics is computed: energy in
  Hartree, length in Bohr, unless a spec explicitly states otherwise.
- Conversions live in `libs/qscat/qscat/units.py` (e.g. `hartree_to_ev`,
  `ev_to_hartree`, using CODATA 2018 constants). Convert at the boundary
  (I/O, reporting) — keep internal computation in atomic units.

## FEM-DVR-ECS Notation

- **DVR** (Discrete Variable Representation) — grid-based basis where the
  potential is diagonal; lives in `libs/qscat/qscat/dvr/`.
- **ECS** (Exterior Complex Scaling) — contour deformation technique for
  handling continuum/scattering boundary conditions; lives in
  `libs/qscat/qscat/ecs/`.
- **FEM-DVR** — finite-element DVR, combining piecewise DVR grids across
  elements; the combination these two subpackages are meant to compose
  toward as methods are added.

## `qscat` Subpackage Map (`libs/qscat/qscat/`)

| Subpackage  | Purpose                                                |
|-------------|---------------------------------------------------------|
| `special`   | Special functions (analytic benchmarks live here too)   |
| `dvr`       | Discrete Variable Representation grids/bases            |
| `ecs`       | Exterior Complex Scaling                                |
| `evolution` | Time evolution / propagators                            |
| `linalg`    | Linear algebra helpers                                  |
| `units`     | Atomic-unit conversions (`units.py`, not a subpackage dir)|

Only validated, reusable code lives here (see `qm-method-lifecycle` step 5) —
`projects/<name>/` is for in-progress toy models.

## Tolerance Defaults

- Never compare floats with bare `==`; always state `rtol`/`atol` explicitly
  (see `numerical-validation` for the full technique set).
- Typical bands: `rtol=1e-8`–`1e-10` for analytic-benchmark comparisons;
  `rtol=1e-12` or tighter for differential tests between two implementations
  of the same deterministic arithmetic (e.g. Python vs. a Rust kernel on
  identical inputs, as in `native/qscat-kernels/tests/test_l2_norm.py`);
  looser (`1e-6` or method-dependent) for convergence-study error bounds,
  since the point there is a trend, not a fixed target.

## Naming

- Python: `snake_case` for functions/modules, matching `libs/qscat/qscat/`
  (e.g. `hartree_to_ev`, `l2_norm`).
- Rust kernel crates live under `native/<crate-name>/` and compile to a
  Python module conventionally named `<crate>_kernels` (e.g. the
  `qscat-kernels` crate builds the `qscat_kernels` Python module — see
  `[tool.maturin] module-name` in `native/qscat-kernels/pyproject.toml`).
- Reference oracles under `reference/` (`reference/eMoScat`,
  `reference/libXcuda`) are read-only — never edited, only read for
  algorithms/expected outputs.

## CPU-First

Everything in qModeling runs on CPU and is containerizable (see
`containerize-and-run`); there is no GPU runtime dependency even though
`reference/libXcuda` exists as a read-only algorithmic oracle.
