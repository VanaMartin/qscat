# Changelog

All notable changes to `qscat` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `qscat.core.lcp` resonance/scattering observables: `lcp_da_cross_section(...,
  return_wavefunction=True)` returns the 1-D nuclear resolvent `psi_sc(R)` per
  energy (the DA scattering state; same overload convention as
  `ve`/`da`/`dr_cross_section`); `resonance_eigenstate` and
  `resonance_eigenstate_at_peak_width` return the resonance pole energy
  (`E_r - iGamma/2`) + its c-product-normalized electronic eigenfunction (the
  eigenstate counterpart of `local_complex_potential`, which keeps only the
  pole energy), the width-peak variant robustly skipping the frozen small-R
  continuation tail.
- `qscat-run` LCP artifacts: with `artifacts.eigenstates`, an LCP run now also
  emits the **resonance electronic eigenstate** (`ResonanceState`, `resonance/`
  npz+png) at the width peak; with `wavefunction_snapshots.full_field` +
  `ti_energies`, the **LCP nuclear scattering states** `psi_sc(R)` at those
  collision energies (an `EigenStates` of `kind="lcp_scattering"`). Closes the
  last "all observables" gaps for the LCP method (resonance states + LCP
  wavefunctions).
- Docker base image now leverages OpenMP-threaded numerics: the **OpenMP
  OpenBLAS** variant is selected (`update-alternatives`) and the hard
  `OPENBLAS_NUM_THREADS=1` pin is dropped, and an **OpenMP-enabled MUMPS 5.9.1
  is built from source** (scivision/mumps CMake) to replace Debian's non-OpenMP
  `libmumps-seq`, so BLAS/LAPACK and MUMPS share one libgomp thread pool. The
  MUMPS backend stays exact vs SuperLU (`max|dx|~6e-16`). (Measured: on the
  2-D ECS decks this halves the single-thread factor via the newer MUMPS but
  does not scale with threads — small independent fronts; see GitHub #3/#4.)
- Docker images render `qscat.viz` animations: the base image ships **ffmpeg**
  (matplotlib's `FFMpegWriter` backend) and the `build` stage installs the
  `plot` extra (matplotlib), so the `runtime` image can write `.mp4`/`.gif`
  out of the box and the `test` image exercises a new ffmpeg-gated `.mp4` viz
  test (`test_animate_wavefunction_writes_mp4`) that skips on a bare Mac.
- `qscat.viz` phase reference: `animate_wavefunction(..., phase_reference=E_ref)`
  colours each frame after a global phase `e^{+i E_ref * t}`, showing the phase
  RELATIVE to a channel base energy (removes the fast base-energy hue spin);
  `WavefunctionArtist.update(state, phase=...)` / `plot_wavefunction_2d(..., phase=)`
  are the static knob. A global phase, so `|psi|`/contours are unchanged; applied
  on the projected field (projector stays phase-agnostic), ~zero overhead.
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
