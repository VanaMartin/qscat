# Changelog

All notable changes to `qscat` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `qscat.viz` animation: `WavefunctionArtist` (draw/update one HSV panel in a
  caller-supplied axes — isolated from figure creation, so panels compose),
  `animate_wavefunction` (a decoupled sequence of `psi(t)` frames → `.mp4`/`.gif`,
  fixed mag + static potential overlay, per-frame `|psi|`), and `animate_artists`
  (several panels in one figure). `plot_wavefunction_2d` is now a thin wrapper
  over `WavefunctionArtist`. (TODO noted: an inverse-value/print colour mode.)
- `qscat.viz`: combined contour overlay — the solid-white `|psi|` contours plus a
  DEDICATED dotted potential overlay (`potential=` + `potential_levels=`) drawn
  together, with inline energy labels. `qscat.viz.energy_contour_levels` picks
  physically relevant levels (vibrational thresholds `eps_v` and total energies
  `eps[v_init]+E` = classical turning surfaces) via `potential_levels="auto"`.
- `qscat.viz.plot_wavefunction_2d` contour overlay: thin white lines at 0.6
  opacity (all overridable), tracing either the wavefunction `magnitude`
  (`|psi|`, levels derived from `mag`) or a `potential` — supplied as nodal
  values on the same tensor grid (projected via the new
  `EquidistantProjector.project_values`; the robust path for a numerically-
  evaluated potential) or as a callable `V(r, R)` for an analytic one.
- `qscat.viz`: 2-D wavefunction visualisation (ported from eMoScat's
  `EquidistantProjector2d` + `display_wf.py`). `EquidistantProjector` caches a
  sparse projection of a 2-D `TensorGrid` state onto a uniform grid (build once,
  apply per frame); `complex_to_rgb`/`complex_to_hsv` domain-colour the field
  (phase→hue, magnitude→brightness); `plot_wavefunction_2d` renders it
  (matplotlib, lazy). The building block for time-dependent animations.
- `qscat.dvr.dvr_interpolation_matrix`: the sparse FEM-DVR-ECS field-value
  operator (Lagrange interpolation at arbitrary points, ECS-mapped, with the
  dropped-Dirichlet-endpoint bridge factors) — the projector's kernel.
- Docs: numpydoc-format docstrings for the getting-started public API
  (`electronic_grid`, `nuclear_grid`, `vibrational_states`, `IncidentSpec`,
  `GridError`) so the API reference renders parameter/return/raises tables; the
  Sphinx build is now strict (`-W`, 0 warnings) with a doctest step in CI.
- Typed exception hierarchy `qscat.exceptions` (`QscatError` base +
  `GridError`, `ModelError`, `BackendError`, `ConvergenceError`); each also
  subclasses the built-in it replaces, so `except ValueError`/`RuntimeError`
  still catches them. Grid-construction and sparse-backend raises now use it.
- Lazy top-level namespace (PEP 562): `qscat.dvr`, `qscat.core`, … import on
  first access; the exception classes are exposed as `qscat.QscatError` etc.
- Public API stability policy (ADR 0004); the wide solver signatures are marked
  provisional pending the context-object refactor.
- Property-based tests (hypothesis) for the D-general primitives (`ecs_map`,
  `gll_nodes_weights`, `kron_sum`).
- Packaging for PyPI: complete `[project]` metadata (readme, authors,
  classifiers, keywords, URLs), BSD-3-Clause `LICENSE`, `CITATION.cff`, and a
  `py.typed` marker so downstream type-checkers see qscat's types.
- `qscat[plot]` optional dependency for the matplotlib-based figure helpers.
- Single-sourced version (`qscat.__version__`, read dynamically by the build).

### Changed
- `qscat.core.plot_cross_sections` now imports matplotlib lazily, so importing
  `qscat.core` / `qscat.tuning` no longer requires matplotlib (it is the
  optional `plot` extra). Previously a clean install without matplotlib crashed
  on `import qscat.core`.
- `qscat.core.time_dependent.free_hamiltonian` and `qscat.core.Flux.series` are
  now public (were `_free_hamiltonian` / `Flux._arrays`).

[Unreleased]: https://github.com/VanaMartin/qscat/compare/main...HEAD
