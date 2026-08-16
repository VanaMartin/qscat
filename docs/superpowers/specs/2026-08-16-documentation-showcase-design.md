# Documentation: structured API reference + honest showcase gallery

**Date:** 2026-08-16
**Status:** approved (design); execution split into two plans

## Goal

Make the qscat documentation carry two things it does not carry today:

1. **A navigable API reference** — one page per subpackage instead of a single
   119-name scroll, with every public name reachable.
2. **A showcase gallery** — eight worked examples, each with its model, its
   parametrisation, a figure produced by a committed config on real hardware,
   and a reading of what the figure shows *including where the physics
   disagrees*.

The organising discipline for Part B: **no figure ships that a reader cannot
regenerate.** Every gallery figure is produced by a committed `qscat-run` YAML,
run under Docker on `sadaharu`, and lands next to that run's `manifest.json`.

## Positioning: what the gallery is allowed to claim

`docs/related-work.md` (added to `main` on 2026-08-16) surveyed whether this
functionality already exists as released code. Its finding governs how every
gallery page is framed:

- **The numerical core is not distinctive.** `quantumgrid` (McCurdy, Streeter &
  Barbalinardo; PyPI 2.0.30) publishes the same FEM-DVR grid, bridge functions,
  ECS contour, and Crank-Nicolson propagator. A page that leads with "we have
  FEM-DVR-ECS grids" is advertising the one layer that has a published
  counterpart.
- **The physics layer above it has no released counterpart that turned up** —
  the electronic × nuclear driven solve for VE/DA/DR, the LCP approximation,
  the TD extractors, the discretisation tuner, and the complex-symmetric MUMPS
  sweep. That is what the gallery advertises.

Concretely: every showcase page leads with an **observable and a comparison**
(a cross section against published data, an approximation against its oracle,
three routes against each other), never with a discretisation technique. The
grid is stated in the Cost section as a parameter, not sold as a feature.

The eight entries chosen below already satisfy this; the rule exists so that
page copy written later does not drift back toward selling the numerics.

**A follow-up this survey implies, deliberately out of scope here:**
`quantumgrid` is an independent implementation of the same Rescigno & McCurdy
(2000) method and is therefore a legitimate external cross-check of
`qscat.dvr` / `qscat.ecs` that this repository does not currently use. Adding
it is a validation task, not a documentation task.

## Execution split

Two plans. **Plan A merges before Plan B is written.**

- **Plan A** — Part A below. Mechanical, no compute, ~12 files. Unblocks
  nothing else, so it should not wait behind a day of solver runs.
- **Plan B** — Parts B and C below. Two library features, a day of production
  runs, eight pages, figure dedup, README rework.

---

# Part A — structured API reference

## Current state

`docs/api.md` is a single page: eight `automodule` blocks, 119 public names in
one scroll, `qscat.core` alone contributing 37. Two submodules are **absent
entirely** — `qscat.viz` (9 public names) and `qscat.units` (4) — so 13 public
names have no rendered documentation at all.

## Target structure

```
docs/api/index.md     layer map + the ADR 0004 contract; toctree over the pages below
docs/api/core.md      37 names, grouped under H2 sections (see below)
docs/api/model.md     ResonanceModel protocol, Diatomic/Ionic forms, N2/NO/F2/H2P registry
docs/api/dvr.md       15 names
docs/api/ecs.md        3 names
docs/api/linalg.md     6 names
docs/api/evolution.md  4 names
docs/api/special.md    7 names
docs/api/tuning.md    22 names
docs/api/viz.md        9 names   -- NEW, currently undocumented
docs/api/base.md      units (4) + exceptions (5)  -- units currently undocumented
```

`docs/api.md` is deleted; `docs/index.md`'s toctree points at `api/index`.

**`qscat.core` sectioning.** 37 names in source order is the problem the split
is meant to fix, so this page groups them by job:

| section | names | n |
|---|---|---|
| The problem object | `ScatteringProblem` | 1 |
| Cross sections (time-independent) | `ve_cross_section`, `da_cross_section`, `dr_cross_section` | 3 |
| Cross sections (time-dependent) | `td_ve_cross_section`, `td_ve_cross_sections_all`, `td_da_cross_section`, `td_da_cross_sections_all` | 4 |
| Grids | `electronic_grid`, `nuclear_grid`, `fem_grid_exp_tail`, `segmented_grid` | 4 |
| Channels | `channel_vector`, `anion_electronic_states`, `v_dr_diag`, `outgoing_channel`, `outgoing_channel_nuclear` | 5 |
| Vibrational structure | `vibrational_states`, `VibrationalBasis` | 2 |
| Wavepacket & correlation | `gaussian_coeffs`, `initial_state`, `eta_incident`, `eta_outgoing`, `hankel_point_value`, `outgoing_surface_wave`, `propagate`, `sigma_from_correlations` | 8 |
| TD energy extractors | `Extractor`, `TannorWeeks`, `Dirac`, `Flux` | 4 |
| LCP approximation | `local_complex_potential`, `lcp_da_cross_section`, `resonance_levels`, `lcp_resonance_levels`, `ResonanceLevels` | 5 |
| Plotting helpers | `plot_cross_sections` | 1 |

That is all 37 names in `qscat.core.__all__`, verified against the installed
package. If `__all__` changes before implementation, the coverage test below is
what catches the drift.

## Page template

Each subpackage page:

```markdown
# qscat.<name>

One paragraph: what this layer is and when a user reaches for it.

```{eval-rst}
.. currentmodule:: qscat.<name>

.. autosummary::
   :nosignatures:

   <every name in __all__>
```

## <Section>

```{eval-rst}
.. autofunction:: qscat.<name>.<fn>
```
```

The `autosummary` table at the top is the scannable index the current page has
no equivalent of; the detail follows underneath.

## The coverage test

`libs/qscat/tests/test_api_docs_coverage.py` (or the docs-adjacent equivalent):

> For every submodule in `qscat._SUBMODULES` plus `qscat.exceptions`, every name
> in that module's `__all__` appears somewhere under `docs/api/`.

This is the load-bearing part of Part A. The layout is cosmetic; the test is
what stops the `viz`/`units` gap from silently recurring. It must fail if a new
public name is added without a docs entry.

## Also in Part A: publish the theory notes

`docs/conf.py` currently lists `physics` in `exclude_patterns`, so 18
substantial theory notes are invisible on the published site. Part A
un-excludes them and adds a **Theory** toctree section to `docs/index.md`.

The work this actually implies is link-fixing: the notes contain relative links
to repo paths (`validation/n2/...`, `libs/qscat/qscat/...`) that resolve in a
clone but 404 from a rendered site. Each such link is either

- rewritten to an absolute GitHub URL (`https://github.com/VanaMartin/qscat/blob/main/...`), or
- converted to a plain-text path reference (backticks, no link) where the point
  is "this lives at X", not "click here".

The notes also cross-reference `docs/superpowers/` plans and specs, which stay
excluded from the site — those links become plain-text references, consistent
with the `mastering-github` self-sufficiency rule that main must stand alone.

## Gate

- `uv run sphinx-build -b html -W docs docs/_build/html` clean. CI already runs
  with `-W`, so any unresolved reference or orphaned page fails the build.
- The coverage test passes.
- No page under `docs/api/` is orphaned from a toctree.

---

# Part B — the showcase gallery

## B1. Prerequisite feature: `reference:` in the qscat-run config

The flagship figure overlays this solver's N₂ cross section on Karel Houfek's
independent published data. That data **is** committed
(`validation/n2/data/CSVE.V00.J00`), so the comparison is reproducible from a
bare clone — but `qscat-run` has no reference-overlay capability, so today the
overlay exists only inside `validation/n2/ti_curve.py`.

Add to the config schema:

```yaml
reference:
  - path: validation/n2/data/CSVE.V00.J00   # relative to the config, or absolute
    format: houfek                           # named loader
    label: "Houfek (2006), CSVE.V00.J00"
```

- Loader lives in a new `apps/qscat-run/qscat_run/reference.py`.
- **It reads a file path. It does not import `validation`** — the layering rule
  (`tests/test_no_validation_import.py`) stays intact, because the data file is
  named by the config, not imported as a module.
- Reference series are drawn on `cross_section.png` (dashed/scatter,
  distinguishable from solver curves) and written as extra columns in
  `cross_section.csv` under a `ref:` key prefix, disjoint from the existing
  `ti:`/`td:`/`lcp:` prefixes.
- `validate_config` errors actionably when `path` is missing or `format` is
  unknown.
- A missing reference file is a config error at validate time, not a silent
  skip at plot time.

## B2. Prerequisite feature: region-split domain colouring in `qscat.viz`

`complex_to_hsv(z, mag)` normalizes the whole field by a single scalar `mag`.
In a time-dependent field the incident packet's amplitude dwarfs the resonant
and outgoing amplitude by orders of magnitude, so a global normalization
renders everything interesting as black. This is a genuine limitation, not a
tuning preference.

Change:

- `complex_to_hsv` / `complex_to_rgb` accept an **array-valued** `mag`
  broadcastable to `z.shape`. Scalar `mag` keeps working unchanged (existing
  callers and tests are unaffected).
- A new helper builds a piecewise scale along one axis — separate normalization
  for the interaction region and the asymptotic region, split at a configurable
  coordinate, each region scaled to a percentile of its own magnitude rather
  than its max (so a single hot pixel doesn't re-flatten the region).
- `plot_wavefunction_2d` / `animate_wavefunction` gain the parameter needed to
  pass it through.

**This will need visual iteration.** The published TD panels (Váňa & Houfek,
PRA 95 (2017), Figs. 1-25; Hvizdoš 2016 thesis Figs. 2.6-2.12) are the target.
The split radius and the per-region percentile are expected to take two or
three passes before the figure reads well; the plan budgets for that rather
than assuming one shot. The deliverable is a figure where the resonant and
outgoing structure is visible at the same time as the incident packet.

Note: `coloring.py` already carries a deferred TODO for a print-mode brightness
inversion. That stays deferred — it is a separate concern from region scaling.

## B3. Configs

Gallery configs live in `apps/qscat-run/examples/gallery/*.yaml`.

`apps/qscat-run/tests/test_examples.py` uses `EXAMPLES_DIR.glob("*.yaml")`,
which is **not recursive** — it must become `rglob("*.yaml")` (or
`glob("**/*.yaml")`) so gallery configs are schema-gated by the existing
`test_example_validates_clean` parametrisation for free. The existing
"directory is not empty" guard already protects against a typo'd glob.

Each gallery config carries a header comment naming the gallery page it feeds
and stating that it is a production deck (not a fast local example, unlike the
existing `examples/*.yaml`).

## B4. Running on sadaharu

`sadaharu` is the x86_64 Docker host: 32 cores, 123 GB RAM, `qmodeling-base`
and `qmodeling:runtime-cpu` already built. Production decks need the `test`
image (it carries MUMPS; `runtime` deliberately does not).

A committed driver — `docs/gallery/run-gallery.sh` — iterates the gallery
configs, invokes `docker/run.sh <config> <out>` for each, and collects
`cross_section.png` (and the entry's other artifacts) plus `manifest.json` into
`docs/gallery/figures/<entry>/`.

The script is the reproduction recipe, not a convenience: each gallery page
quotes the single `docker/run.sh` line for its own entry so a reader can run
exactly that one.

**Compute budget (estimate, to be replaced by measured values from the
manifests):** N₂ dense sweep tens of minutes; each TD entry roughly an hour;
H₂⁺ DR at ~1.15M unknowns is the long pole and may run for hours. Roughly a day
of mostly-unattended wall-clock across all eight entries. Actual per-entry
wall-clock is read from each `manifest.json` and quoted on the page.

## B5. Provenance layout

```
docs/gallery/figures/<entry>/
    <figure>.png            the committed figure
    manifest.json           qscat version, git SHA, timestamp, backend,
                            per-stage timings, host platform
```

The manifest is small JSON and is the honesty mechanism: it ties the figure to
a version, a commit, and a backend. A figure without its manifest does not ship.

## B6. The eight pages

`docs/gallery/index.md` is a contact sheet — thumbnail plus one line per entry.

**`models.md` — what you are solving.** Comes first, before any cross section.
Per molecule (N₂, NO, F₂, H₂⁺): the neutral curve V₀(R), the resonance position
V_d(R), the width Γ(R), the vibrational ladder, a parameter table, and the
published source (Houfek *et al.* PRA 73 (2006) for N₂/NO; PRA 77 (2008) for
F₂; Váňa & Houfek PRA 95 (2017) and the 2017 thesis for H₂⁺). Cheap to compute
— V₀ is analytic and V_d/Γ come from `qscat.core.lcp.local_complex_potential`.

| page | what it advertises |
|---|---|
| `n2-ve.md` **(hero)** | Exact 2-D TI σ_{0→v'}(E) overlaid on Houfek's independent data, boomerang oscillations resolved point-by-point. The credibility anchor: someone else's data, reproducible from a bare clone. |
| `f2-da-lcp.md` | F₂ σ_DA(E), exact 2-D vs the LCP approximation on one axis. ~11% agreement away from threshold and a **visible failure near threshold**. The repo's actual thesis: the exact solution is the oracle, the approximation is under test. |
| `three-routes.md` | The same cross section from the TI driven solve and from TD propagation with three energy extractors (Tannor-Weeks / Dirac / flux), agreeing to ~1-3%. Independent numerics converging. |
| `h2p-dr.md` | H₂⁺ dissociative recombination: Coulomb incident, Rydberg exit series, ~1.15M unknowns under Docker/MUMPS. Breadth and scale. |
| `td-wavefunction.md` | A small committed **GIF** of Ψ(r,R,t) with region-split domain colouring (B2). The visual hook; GIF so it renders on GitHub, PyPI, and the docs site alike. |
| `tuning.md` | The discretisation tuner deriving a grid a priori from the potential — 37% fewer points than the hand-tuned F₂ deck, with the honest caveat that the 1-D probes are necessary but not sufficient (the 2-D spot-check finding). |
| `mumps.md` | MUMPS vs SuperLU on the real N₂ matrices: 72× factor time, 9× peak RSS at the 143k deck. A table plus the `benchmarks/mumps_vs_superlu` invocation. |

## B7. Page template

Every showcase page, in this order:

1. **What it shows** — one paragraph.
2. **The model** — one line, linking to `models.md`.
3. **The config** — `literalinclude` of the committed YAML. Never a
   paraphrase; the file itself.
4. **The figure.**
5. **How to read it** — the physics, and **explicitly where it disagrees or
   is limited**. A gallery entry that only reports agreement is advertising,
   not documentation. Each entry names its known limitation (LCP's
   near-threshold departure; TD's near-threshold behaviour; the tuner's 1-D
   probe insufficiency; the N₂/NO deck-size finding).
6. **Cost** — grid sizes, unknowns, backend, wall-clock, read from the manifest.
7. **Reproduce** — the exact `docker/run.sh` line.
8. **Sources** — links to the `reference/literature/*.md` notes.

## B8. Figure dedup

Gallery figures supersede several existing `docs/physics/figures/*.png` that
were produced by the now-retired validation drivers. **The dedup rule is
narrower than "delete what's superseded", because of a constraint the figure
map exposed:** most figures are also referenced by frozen plans and specs under
`docs/superpowers/`, which are historical records of what was done and must not
be rewritten.

> **Rule.** Delete a superseded figure only when *every* reference to it lives
> in a live document — a `docs/physics/` note or the README — which is
> repointed at the gallery figure in the same change. If any frozen
> `docs/superpowers/` plan or spec references it, **keep the file** and leave
> the frozen document untouched.

Supersession also requires the gallery figure to show the same thing. A physics
note may analyse a detail the gallery figure does not reproduce (e.g. the
correlation-function panel, the nuclear-density panel); those are not
superseded and stay. Each candidate is checked individually — no blanket pass.

---

# Part C — README and site landing

## Root README

- **Hero figure** — the N₂-vs-Houfek overlay, near the top.
- **Highlights block** — three or four thumbnails linking to their gallery
  pages. Short; the site carries the depth.
- **Fix the stale reference at `README.md:104`**, which cites
  `validation/diatomic/da_curves.py` and `curves.py`. Both were **deleted** in
  the qscat-run consolidation, so the README currently points a reader at files
  that do not exist in a clone. Repoint at the gallery entry and its config.
- The two remaining `docs/superpowers/specs/...` references (lines 154, 178)
  are re-checked against the self-sufficiency rule while the file is open.

## `docs/index.md`

Toctree becomes: `getting-started` → `gallery/index` → `api/index` →
`related-work` → Theory (`physics/*`). The gallery sits between the tutorial
and the reference, which is where a reader deciding whether to adopt the
library will look; `related-work` (already in the toctree as of `main`
c327865) keeps its place after the reference.

## Package README (`libs/qscat/README.md`)

The PyPI long-description. Gets the hero figure as an absolute raw-GitHub URL
(relative image paths do not render on PyPI) and a link to the gallery.

---

# Non-goals

- Replicating every published result. The gallery showcases highlights; the
  validation suite remains the exhaustive check.
- New physics. Every number in the gallery comes from capabilities that already
  exist and are already validated.
- The print-mode brightness inversion deferred in `coloring.py`.
- Rewriting frozen `docs/superpowers/` plans and specs.
- Publishing `docs/adr/` or `docs/superpowers/` to the site.

# Success criteria

**Plan A**
- `docs/api/` renders, every page reachable from a toctree, `sphinx-build -W` clean.
- The `__all__` coverage test passes and fails when a public name is undocumented.
- `qscat.viz` and `qscat.units` are documented for the first time.
- The 18 theory notes render on the published site with no broken links.

**Plan B**
- Eight gallery pages, each with a figure produced by a committed config and
  accompanied by that run's `manifest.json`.
- Every gallery config validates clean under the recursive `test_examples.py` glob.
- The `reference:` overlay works and `test_no_validation_import.py` still passes.
- The TD GIF shows resonant and outgoing structure simultaneously with the
  incident packet.
- Every showcase page names a limitation as well as a result.
- `README.md:104`'s dead reference is gone.
