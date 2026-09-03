# Changelog

All notable changes to `qscat` are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

This file covers the whole qModeling monorepo, which is what the published
`qscat` package links to. Entries naming `apps/qscat-run`, `docker/`,
`projects/` or `validation/` describe the surrounding repository rather than the
installed package.

## [Unreleased]

### Added
- **`energies: {ranges: [...]}` — a sweep written as `np.arange` segments.**
  A level-aware mesh is a coarse background sweep plus a dense window around
  each resonance level: a union of uniform segments, and nothing more. Written
  out point by point, O₂'s VE mesh was 3343 lines generated from 27 numbers —
  a file long enough that nobody reads it, and in which the background step,
  the window width and the level a point belongs to are all invisible. The
  segments record all three, and the three O₂ configs drop from 42 kB to
  2.9 kB. Segments may overlap; the mesh is their union, so an energy two
  segments both reach is solved once. `stop` is exclusive, as in numpy.
  Rewriting the O₂ mesh this way also repaired it: those configs had been
  rounded to 6 decimals, a 1.0e-6 Ha grid against their own finest step of
  1.1e-6 Ha — a ratio of 1.10, so the rounding was nearly as coarse as the
  mesh it was quantising. That collapsed 12 energies onto duplicates of their
  neighbours. Expanded from segments the mesh carries no duplicates and sits
  within 0.014 meV of the published axis everywhere. The axis a run actually
  solved is stored in its own `cross_section.npz` regardless, so the config
  records the recipe and the artifact records the result.

- **Computed artifacts move out of git, and `qscat-run fetch` brings them
  back.** 37% of every byte this repository had ever stored was computed
  output — 23.7 MB of blobs under `docs/physics/figures` and `validation`
  against a 39 MB packed clone — none of it source, none of it read by a
  human. Converged sweeps are now published to public object storage and the
  run directory keeps a KB-sized `artifacts.json` naming a URL prefix and a
  sha256 per file; `qscat-run fetch DIR` downloads what it names and verifies
  every byte, leaving nothing on disk that fails its digest. Reads are
  anonymous HTTPS, so a URL works from a notebook or a `curl` with no account
  and no client library. Objects are addressed by CONTENT — the sha256 goes
  into the filename — so a cited URL cannot come to mean something else, and
  a re-run that reproduces its numbers republishes to the same key and costs
  nothing. What stays in git is what a test or a note needs in
  order to stand alone: golden inputs, fit reports, cited figures. The
  classification and the measurements behind it are in
  `docs/adr/0008-computed-artifacts-live-in-public-object-storage.md`; the
  reader's side is `docs/artifacts.md`. Nothing is migrated yet — this adds
  the machinery and the rule.
- **Rendered mathematics in the documentation.** `docs/conf.py` gained
  `sphinx-copybutton`, `sphinx-design`, `sphinxext-opengraph` and
  `sphinx.ext.githubpages`, plus equation numbering (`numfig`,
  `math_eqref_format`) and Furo source links. `sphinx.ext.mathjax` and MyST
  `dollarmath` had been loaded and unused since the site was created — every
  physics note contained zero `$` characters. Every note is now typeset: three
  converted as pilots (`femdvr-ecs`, `nonlocal-resonance-model`,
  `n2-2d-cross-section`), and the remaining twenty followed in the same branch
  rather than the staged "as they are next touched" rollout first planned —
  group by group along the sidebar's sections, each file checked with a
  numeric-token multiset diff so no measured value could change under cover of
  a notation change. The checklist is in
  `docs/superpowers/plans/2026-08-19-docs-latex-and-theory-ia.md`.
- **A portability rule the notes are held to.** Notes under `docs/physics/`
  are read both as files in a clone and as pages on the site, so they use only
  what GitHub's renderer and MathJax both understand: plain `$...$`/`$$...$$`,
  `\begin{aligned}` inside `$$`, `\tag{n}` for published equation numbers —
  no MathJax macros, no `sphinx-design` directives. `tests/test_docs_portability.py`
  enforces it; the convention and a canonical symbol table live in the
  `qscat-conventions` skill.
- **A navigable Theory section.** The sidebar's flat 22-entry list became seven
  collapsible groups behind new index pages, plus `docs/physics/validation-harnesses.md`
  as a flat entry recording what each harness under `validation/` gates, against
  which oracle, at what tolerance. No note moved — 443 lines across 165 tracked
  files cite those paths.
- **Per-molecule guides** (`docs/molecules/{n2,no-f2,h2plus}.md`): what is known
  about each molecule, the headline numbers with their caveats, and links into
  the method notes. Nothing on the site previously answered "what do we know
  about F₂".
- `qscat.core.Verdict`: the `Literal` of the seven legal `assignment` verdict
  strings is now exported from `qscat.core`, so callers can type-check an
  `OverlapPair.verdict` without reaching into `qscat.core.assignment`.
- `qscat.core.exact_resonance_states` (+ `ExactResonanceStates`): poles of the FULL
  2-D S-matrix — eigenvalues `E_r − iΓ/2` of the complex-scaled electronic × nuclear
  Hamiltonian, with no Born-Oppenheimer separation, no discrete state and no local
  approximation. These are what `resonance_levels` approximates, so the pair makes
  the non-adiabatic error directly measurable. Identification generalizes ECS angle
  stability to TWO angles: three spectra (base, electronic-angle moved,
  nuclear-angle moved) with a state accepted only if it survives both comparisons,
  and both residuals reported. Validated against the exact separable-limit oracle
  (pole to 8.6e-15, width to the electronic width alone, eigenvector to the product
  state). On N₂ it returns the anion vibrational ladder and quantifies the
  non-adiabatic error: with the electronic grid converged (order converges at 8,
  `r_max` swept 24→72 bohr), the `v=0` exact pole lies 0.22 meV BELOW the BO/LCP
  level in position and 0.30 meV below in width. The exact pole is below BO in both
  at every level and every box, but only `v=0`'s difference is converged — widths
  converge more slowly than positions, and higher levels later than lower ones. Also
  `qscat.core.plot_resonance_levels` (generic complex-level plotting) and committed
  N₂ figures in the published convention. See docs/physics/exact-2d-resonances.md.
- `qscat.linalg.ShiftInvertEigs`: the `k` eigenpairs of a sparse complex-symmetric
  matrix nearest a complex shift (shift-invert Arnoldi with `SparseLU` as the inner
  solve, so it inherits the MUMPS backend). Resonances are INTERIOR eigenvalues,
  which a plain Krylov iteration cannot reach. Because `A - sigma*I` keeps its
  sparsity pattern for every shift, repeated `near(sigma)` calls reuse the symbolic
  analysis via `SparseLU.refactor` — the eigenvalue analogue of the time-independent
  energy sweep. Validated in 1-D against dense `np.linalg.eig`, eigenvalues (rtol
  1e-9) and eigenvectors (c-norm overlap), on synthetic complex-symmetric matrices
  and on the real N₂ electronic FEM-DVR-ECS Hamiltonian, where the recovered pole
  gives `E_res = 2.441 eV` / `Gamma = 0.4535 eV` against the published 2.445 /
  0.455 eV. See docs/physics/shift-invert-eigensolver.md for the conventions (the
  `A - sigma*I` sign, nearest-shift ordering) and the measured working range.
- `qscat.core.lcp` BO/LCP resonance levels: `lcp_resonance_levels` diagonalizes the
  nuclear Hamiltonian `H_N = T(mu) + V_d(R) - i*Gamma(R)/2` in the LCP complex curve
  (`ResonanceLevels`: complex `energies`/`widths`, c-product-normalized `states`, a
  two-angle stability `residuals` diagnostic, `real_weight`, and a `golden_rule`
  perturbative comparator), the Born-Oppenheimer approximation to the 2-D model's
  resonance energies — the thesis's `omega_j`, promoted from eMoScat's real-part-only
  levels to genuine complex quasi-bound states. `resonance_levels(model, ...)` runs
  the electronic pole walk once and adds one extra nuclear diagonalization on a
  second, differently-angled grid for the two-angle selection
  (`qscat.ecs.match_angle_stable`, the multi-state sibling of `find_resonance_pole`).
  `IonicResonanceModel.max_nuclear_ecs_angle_deg = 22.5` records the H₂⁺ nuclear ECS
  divergence bound (Hvizdoš et al., Phys. Rev. A 97, 022704 (2018), Sec. II) that
  `resonance_levels` now enforces. `qscat-run` exposes this as the `resonance_levels`
  observable kind (no `energies:` block needed) and as an opt-in
  `artifacts.resonance_levels` flag on an existing LCP run, writing
  `resonance_levels_{label}.{csv,npz,png}`. These are complex-scaled (ECS) resonance
  eigenstates, explicitly **not** Siegert pseudostates — see
  docs/physics/lcp-resonance-levels.md for the terminology, normalization, and the
  first F2 numbers.
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
- `docs/related-work.md`: a survey of what this functionality already is and is
  not published as code (PyPI, GitHub, Zenodo, Linux packages), the overlap with
  `quantumgrid` (the one other released FEM-DVR-ECS package), how the *ab initio*
  electron–molecule suites differ, what has no released counterpart, and the
  limits of that claim. Summarized in a "Relation to existing work" section of
  the package README.

### Changed
- **`qscat.org` is the project's canonical identity, and the install route is a
  clone.** The landing page went live, which turned two long-standing
  inaccuracies load-bearing. The package metadata and `CITATION.cff` named
  GitHub as the homepage because there was nothing else to name; `Homepage` is
  now `https://qscat.org`, a `Documentation` URL points at the GitHub Pages API
  reference for the first time, and `Repository`/`Issues`/`Changelog` stay on
  GitHub — four public locations, named distinctly, with `data.qscat.org`
  serving the published artifacts. `docs/related-work.md` called `qscat.com` the
  project's domain, in the very section teaching a reader to tell this `qscat`
  apart from the QGIS shoreline tool and NASA's QuikSCAT, where a wrong domain
  is worse than none. The maintainer's `martin@qscat.com` mailbox is
  deliberately unchanged — a website moving is not a mailbox moving — and both
  metadata files now say so where the address is written, so the mismatch does
  not read as a defect. Separately, the install instructions told readers to
  `pip install` a distribution that does not exist and will not until the
  citation article is out, contradicting `libs/qscat/README.md`; that had six
  sites, not the two first noticed — `docs/getting-started.md`, `docs/api/viz.md`,
  and four error messages and a docstring in the shipped library, which reach a
  user at the moment their import fails. All six now name the
  `uv sync --all-packages [--extra plot]` route. Nothing asserted on any of
  those strings; `tests/test_project_identity.py` now does, over the shipped
  sources, the published docs and both metadata files.
- Two public names renamed to end collisions with an identically-named
  sibling elsewhere in `qscat.core`/`qscat.tuning` (api-surface-pass,
  2026-08-25): `qscat.core.nrm.scattering.free_hamiltonian` →
  `electronic_free_hamiltonian` (was shadowed by
  `qscat.core.time_dependent.free_hamiltonian`), and
  `qscat.tuning.resonance_curve` → `resonance_curve_arrays` (was shadowed by
  `qscat.core.bo.resonance_curve`). The old names were kept briefly as
  deprecated aliases and have since been removed — see Removed, below.
- `qscat.core.plot_cross_sections` now imports matplotlib lazily, so importing
  `qscat.core` / `qscat.tuning` no longer requires matplotlib (it is the
  optional `plot` extra). Previously a clean install without matplotlib crashed
  on `import qscat.core`.
- `qscat.core.time_dependent.free_hamiltonian` and `qscat.core.Flux.series` are
  now public (were `_free_hamiltonian` / `Flux._arrays`).

### Removed
- All deprecation and legacy-compatibility machinery, per ADR 0004 point 2
  (pre-1.0 minor releases may break the public API) — the package is still
  `0.1.0.dev0`, so the one-release-cycle grace period these existed to serve
  never applies. The last commit where each removed name/behaviour existed is
  `ccb3fa2`.
  - `qscat._deprecation` (the shared `deprecated_alias_getattr` `__getattr__`
    factory) — deleted outright; nothing replaces it, its only callers were
    the aliases below.
  - `qscat.core.nrm.scattering.free_hamiltonian` (deprecated alias) — use
    `electronic_free_hamiltonian`.
  - `qscat.tuning.resonance_curve` / `qscat.tuning.resonance.resonance_curve`
    (deprecated aliases) — use `resonance_curve_arrays`.
  - `qscat.core.time_dependent._propagate` (deprecated alias) — use
    `propagate_wavepacket`.
  - `qscat.core.time_dependent._s_vector_one_energy` /
    `_sigma_one_energy` (deprecated aliases) — use `s_vector_one_energy` /
    `sigma_one_energy`.
  - `qscat.core.correlation.hankel_point_value`'s legacy grid-first calling
    form (`hankel_point_value(grid, z_position, k, l, charge=0, *, mass=1.0)`)
    — call `hankel_point_value(z_position, k, l, charge=0, *, mass=1.0)`.
  - `qscat.core.correlation.outgoing_surface_wave`'s legacy grid-first calling
    form — call `outgoing_surface_wave(z_surface, k, l, charge=0.0, *,
    mass=1.0)`.
  - `qscat.core.dissociation.dr_cross_section`'s `return_wavefunction` /
    `return_amplitude` flag-shaped tuple return — call `dr_solve(...,
    store_wavefunction=..., store_amplitude=...)` and read the `DrResult`
    fields.
  - `qscat.tuning.propose_grid`'s `rtol` parameter — it was never consumed;
    simply omit it.

### Fixed
- **The coupled-partial-wave summaries no longer assert the withdrawn 58 %
  width result.** The molecule guide (`docs/molecules/no-f2.md`), the
  resonance index (`docs/physics/resonances.md`) and
  `validation/coupled/cross_section.py`'s module docstring still said the
  fixed-wave width and cross section miss the coupled result by tens of
  percent — the headline the physics note itself had already withdrawn. That
  58 % was a POSITION artifact: the two truncations were compared at points
  5-10 mHa apart on `E_res(R)`, where `Gamma` falls steeply with `E_res`, so
  the position difference manufactured a width difference with no angular
  physics in it; pinned to the same `E_res`, the median difference over
  R in [1.6, 2.2] is 0.56 %. The companion 11.8x cross-section factor was
  measured on a model whose anion is unbound at every R, because splitting the
  well hands the deeper centre only `(1+kappa)/2` of `lam`. The three pages
  now carry the surviving result: only `l = 1` hosts a resonance at all (O⁻
  has one bound orbital, 2p), which EXPLAINS the single pole, and the
  truncation costs 2-7 % on the angle-integrated VE cross section against a
  reference converged to 0.3-0.5 % — resolved, but small, because a
  low-energy electron cannot resolve the anisotropy. No physics code, result
  data or analysis section changed; the full account was already in
  `docs/physics/coupled-partial-waves.md`.
- **A published run's inputs stay in git, so an offline clone can reproduce
  it.** The artifact store exists to move expensive *output* out of the
  repository, but the three published O₂ sweeps had also pushed out their
  `config.resolved.yaml` — 43 kB of resolved input, ignored by `.gitignore`
  and listed as fetch-only — so a clone with no network could see that a sweep
  existed and not what it was a run of, while
  `docs/adr/0008-computed-artifacts-live-in-public-object-storage.md` and the
  `artifact_store` docstring both promised otherwise. The three resolved
  configs are now tracked, and `config.resolved.yaml` and `manifest.json` are
  no longer named in any pointer; the remaining entries and their digests are
  untouched, so every published URL still resolves. Nothing failed before,
  because nothing can: a directory whose config arrived by fetch looks exactly
  like one whose config was cloned. The invariant is therefore held by a guard
  over the committed pointers (`apps/qscat-run/tests/test_artifact_store.py`),
  which fails if a published directory lacks either record, if a pointer lists
  one of them, or if the allow-list lets one exist untracked.
- **`qscat.core`'s installed API inventory describes the `dr_cross_section`
  that ships.** Removing the flag-shaped tuple returns (see Removed, above)
  left the module docstring still calling `dr_cross_section` "a thin
  flag-shaped-tuple wrapper kept for one deprecation cycle" — prose promising
  a grace period the package had already dropped, in the one place a reader
  with only the installed package looks before reading a signature. The
  inventory now describes the sigma-only callable as it is, including its
  `n_channels` and `ordering` keywords and the shape and channel ordering of
  what it returns, and points at `dr_solve`/`DrResult` for a wavefunction, a
  T-matrix amplitude, or any other detailed result — noting that the object-API
  method `ScatteringProblem.dr_cross_section`, a separate callable, does still
  carry flags. Nothing executes a docstring, which is how it fell behind; a new
  `libs/qscat/tests/test_public_api_prose.py` ties the inventory to
  `inspect.signature`, so a keyword the prose omits, a keyword it invents, and
  a promise of deprecation machinery the package does not carry each fail a
  fast test rather than a user's call.
- **`manifest.json` records a real commit, or the run fails.** Three committed
  O₂ sweeps carried `"git_sha": "unknown"`, so the artifacts behind the
  spin–orbit VE figure could not be tied to the code that produced them.
  Three independent causes, each sufficient alone: `ARG GIT_SHA=unknown` put
  the literal string into the environment, where it outranked the
  `git rev-parse` fallback written for exactly that case; the Docker
  `runtime` stage never received the variable at all (`ENV` does not cross a
  `FROM`, and `COPY --from=` copies files, not environment); and nothing
  tested any of it. Only a real 40-hex SHA is now honoured, an undeterminable
  one warns rather than writing a plausible-looking manifest in silence
  (`QSCAT_ALLOW_UNKNOWN_SHA=1` silences it where there is genuinely no
  provenance), and a Dockerfile guard fails if a stage starting from a fresh
  base omits the `ARG`/`ENV` pair. Publishing is where the strictness sits:
  the publisher refuses a manifest that cannot name its commit. The
  three sweeps were re-run: cross sections bit-identical across all 3343 × 7
  values, so the repair carries no physics change.
- `qscat.viz` contour colours now follow `inverse`. The `|psi|` contours and the
  dotted potential overlay had fixed defaults (white, `0.75` grey), so the
  `inverse` (light-ground) render drew white on white and the overlay was barely
  visible. They now default to white / `0.75` grey on the dark render and black /
  `0.25` grey under `inverse`; an explicit `contour_color` / `potential_color`
  still wins, so no existing call changes behaviour.
- `H2P.mu` (`qscat.model.library`): `918.25` → `918.076` (`m_p/2` for the modern
  proton mass, `1836.15267/2`). The old value was inherited from eMoScat's JSON
  deck; Vana 2017 Table 1.2 and Hvizdoš et al., Phys. Rev. A 97, 022704 (2018)
  Sec. II A both give `918.076`. The 0.019% error shifted H₂⁺ vibrational
  spacings by ~1e-4 relative — harmless qualitatively, but wrong for reproducing
  published numbers.

[Unreleased]: https://github.com/VanaMartin/qscat/compare/main...HEAD
