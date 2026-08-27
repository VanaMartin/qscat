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

## Mathematics in Documentation

The notes under `docs/physics/` are read in two renderers: as files in a
clone, and as pages on the published site. So they use only what both
handle.

- **Backticks are for code identifiers only** — `ve_cross_section`,
  `SparseLU`, `backend="mumps"`. A backticked `sigma` that means σ is a
  defect; so is a backticked Greek letter — `` `Γ(R)` `` and
  `` `μ = 918.25` `` are maths and are written `$\Gamma(R)$`,
  `$\mu = 918.25$`. A symbol is either code or maths, never Greek dressed
  in backticks. (Enforced: `tests/test_docs_portability.py` flags any
  backtick span containing a Greek letter — zero false positives, because
  this repository's identifiers are ASCII.)
- **Inline mathematics:** `$\sigma_{v\to v'}$`.
- **Display mathematics:** `$$...$$` on its own lines. Multi-line goes
  inside, as `$$ \begin{aligned} ... \end{aligned} $$` — a *bare*
  `\begin{align}` does not render on github.com.
- **Published equation numbers** use `\tag{15}`, with the page locator in
  the surrounding prose (see `mastering-references`).
- **No `{eq}` cross-references in a note.** The `$$...$$ (eq-name)` label and
  the `{eq}` role both show as literal text on github.com, and the portability
  test flags the role. Where a note genuinely writes "Eq. (3)", use `\tag{3}`
  on the equation and refer to it in prose. Labelled equations are available on
  site-first pages, where nobody reads the source.
- **No custom macros.** Spell out `\Psi^{(+)}`, `\sigma_\mathrm{DA}`.
- Atomic units are stated once per note, not per equation.
- **Headings stay plain unicode.** A heading also renders in the sidebar and in
  `toctree` entries, where MathJax does not run — `$^2\Pi_g$` in a title shows
  as literal source there. Write `²Π_g` in the heading and use maths in the body.
- **Inside a Markdown table cell, never write a bare `|` in maths** — it is
  read as a column separator, which splits the row and silently drops a cell
  from the built page without any build warning. Use `\vert` (or `\lvert` /
  `\rvert`).
- **Short level labels stay backticked** — `v = 0`, `v'=1`, `Ry₄`. This is a
  deliberate exception to the backticks-are-for-code rule and the dominant
  practice in the notes (`exact-2d-resonances.md` alone carries ~20 of them).
  A backticked label needs no MathJax, so it survives a table cell without a
  `\vert` hazard and reads identically in both renderers. Use maths for a level
  only when it appears inside a larger expression.
- **MyST directives and roles are for site-first pages only** — `docs/molecules/`
  and the section index pages. That covers `{toctree}`, `{figure}` and the
  `sphinx-design` set (`{dropdown}`, `{grid}`, `{tab-set}`) as well as inline
  roles like `{doc}`. In a note, embed a figure as a plain markdown image with
  the caption as ordinary prose, and link to a sibling note with a relative
  markdown link — ``[`other-note.md`](other-note.md)`` — which resolves in a
  clone and on the site alike.

### What is enforced, and what is convention-by-review

`tests/test_docs_portability.py` scans every `docs/physics/*.md` outside its
`SITE_FIRST_PAGES` allowlist and **mechanically enforces six** of the rules
above:

| rule | detector |
|---|---|
| no MyST directives (any name, ` ``` ` or `:::` form) | `find_myst_directives` |
| no MyST inline roles (`{doc}`, `{ref}`, …) | `find_myst_roles` |
| no `$` in a heading | `find_math_in_headings` |
| no bare `\begin{align\|gather\|equation}` outside `$$` | `find_bare_math_environments` |
| no LaTeX macro definitions, in a note **or** in `docs/conf.py` | `find_macro_definitions`, `test_conf_py_defines_no_mathjax_macros` |
| no Greek letters inside backtick spans | `find_greek_in_backticks` |

Everything else on this page is **convention, checked by review**: backticks
reserved for code identifiers, `$...$` vs `$$...$$` placement, `\tag{}` for
published equation numbers, bare `|` inside a table cell, level labels staying
backticked, and stating atomic units once per note. Backticked maths is detected for the **Greek-letter subset** only
(`find_greek_in_backticks`): a Greek letter inside a backtick span is
always maths here, since the repository's identifiers are ASCII. The
broader family — superscript/subscript unicode in a span — stays
convention-by-review, because it collides with legitimate spans: molecule
names in real paths (`reference/eMoScat/input/{NO,F₂}/grids.txt`),
reaction equations, and the level-label exception above (`Ry₄`).

### Canonical symbols

| Symbol | Meaning |
|---|---|
| $\theta$ | ECS rotation angle |
| $R_0$ | ECS pivot radius (always on an element boundary). **Reserved for this sense only** — the per-molecule equilibrium bond length shares the bare name `R0` as a `qscat.model.N2`-style model attribute (alongside `mu`, `ell`, `D0`, `alpha0`); it is a CODE IDENTIFIER, stays in backticks as `R0`, and must never be typeset as $R_0$. If a display symbol for the bond length is ever wanted, use $R_\mathrm{e}$, not $R_0$. |
| $z(x)$ | the ECS coordinate map |
| $x$ | the unscaled radial coordinate |
| $r$ | electronic coordinate |
| $R$ | internuclear coordinate |
| $\mu$ | reduced mass |
| $E$ | total energy |
| $k$ | wavenumber |
| $\Psi^{(+)}$ | the outgoing-wave driven solution |
| $\chi_v$ | vibrational state $v$ of the neutral |
| $\phi_d$ | the discrete (resonant) electronic state |
| $V_d(R)$ | the resonance curve |
| $\Gamma(R)$ | the resonance width |
| $\sigma_{v\to v'}$ | vibrational-excitation cross section |
| $\sigma_\mathrm{DA}$, $\sigma_\mathrm{DR}$ | dissociative attachment / recombination |
| $\mathbb{1}$ | identity operator/matrix (written as $\mathbb{1}$ rather than $I$ to avoid ambiguity with indices or current) |
| $T$ | **overloaded three ways — disambiguate by CONTEXT, not by shape.** (1) the kinetic-energy operator/matrix: bare $T$, or coordinate-subscripted $T_R$; (2) a T-matrix: bare $T$, channel-subscripted $T_{v\to v'}$, or superscripted by contribution $T^\mathrm{bg}$, $T^\mathrm{res}$ (as in $T = T^\mathrm{bg} + T^\mathrm{res}$); (3) a matrix transpose, as the superscript in $H^{T}$. A bare $T$ is therefore ambiguous on its own: a note that uses more than one sense must say which it means in prose at first use. In practice each note uses one sense throughout — kinetic in `femdvr-ecs`, T-matrix in `n2-2d-cross-section` and `nonlocal-resonance-model` |
| $S$ | S-matrix elements |

## Physics-note structure

A note under `docs/physics/` opens with the standard header block
(**Location** / **Origin** or **Source** / **Units**), then follows this
top-level spine:

1. **Key result** — 5–10 lines, the measured numbers and the one-sentence
   claim first. A reader who stops here leaves with the result and where
   it was measured.
2. **Physical picture** — what the method is and why it exists here.
3. **Method** — the equations, conventions, and implementation choices.
4. **Validation** — what was measured against what, with the evidence
   (and the negative results; the notes record limitations too).

This is the rule for NEW notes. Existing notes migrate opportunistically —
when a note is being edited anyway — never as a big-bang re-heading.
Closest to the spine today (results-first content early, though not yet
under these exact headings): `n2-resonance`, `n2-cross-section`,
`n2-td-cross-section`, `nd-tensor-hamiltonian`, `femdvr-ecs`,
`exact-2d-resonances`; results-first retrofits (a Key result opening
prepended or moved up, body prose otherwise unchanged):
`discretisation-tuning`, `mumps-sparse-backend`,
`h2plus-dr`, `optimization-targets`, `potential-factory`, and the split-out
`nrm-vibrational-excitation`.

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
| `core`      | Model-independent electron–diatomic scattering engine (driven/TI, TD, DA/DR, LCP, NRM, resonance verification) |
| `model`     | The `ResonanceModel` protocol + the per-molecule registry (`N2`, `NO`, `F2`, `H2P`, and the fitted `O2`/`O2_SO12`/`O2_SO32` via `qscat.model.flexible`) |
| `tuning`    | Automatic FEM-DVR-ECS discretisation tuner                |
| `viz`       | Wavefunction rendering and animation (plot extra)         |

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
- **Cross-arch BLAS floor:** a golden value produced *through a sparse
  solve* is compared at `rtol=1e-9`, never `1e-12`. The identical solve on
  a different BLAS (CI's Linux OpenBLAS vs a Mac's Accelerate) differs at
  ~1e-12 (2–3 ULPs), and a `1e-12` gate on such a golden has already failed in
  CI for exactly this reason. The `1e-12` band above is for two
  implementations of the *same deterministic arithmetic* on the same
  machine, not for cross-architecture goldens.

## Naming

- Python: `snake_case` for functions/modules, matching `libs/qscat/qscat/`
  (e.g. `hartree_to_ev`, `l2_norm`).
- **Solver argument order** (docs/adr/0007): observable solvers take the
  discretisation first — `f(grid(s), model, eps, chi, v_init, <per-call
  physics>, *, options)` — while model-derived builders (curves, level
  pipelines, tuning) take the model first. `lcp_da_cross_section` is the one
  documented exception (curve-first, no model — see the ADR).
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
