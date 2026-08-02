# Changelog

All notable changes to `qscat` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
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
