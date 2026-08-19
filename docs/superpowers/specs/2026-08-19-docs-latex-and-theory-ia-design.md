# Rendered mathematics and a navigable Theory section

Date: 2026-08-19

## Purpose

The published documentation at <https://vanamartin.github.io/qscat> does not look
like scientific documentation. Two concrete defects:

1. **No mathematics is rendered.** Across all 22 notes in `docs/physics/` the
   character `$` appears zero times. Every equation is ASCII inside a fenced
   block or backticks — `sigma_{v->v'} = 4 pi^3 |T|^2 / k^2` where the reader
   expects $\sigma_{v\to v'} = 4\pi^3 |T|^2 / k^2$. This is not a broken build:
   `docs/conf.py` already loads `sphinx.ext.mathjax` and enables the MyST
   `dollarmath` extension. The capability is present and unused.
2. **The Theory navigation is a flat list of 22 entries.** `docs/index.md`
   carries one `toctree` with a `:caption: Theory`, which Furo renders as an
   always-expanded block. `docs/physics/README.md` already groups the same notes
   into five themed sections; that grouping is invisible in the sidebar. There is
   no per-molecule entry point at all — a reader who wants "what do we know about
   F₂" has no page to open.

This spec fixes both, and writes down the convention that keeps them fixed.

## Scope

**In scope.** The documentation build configuration; the Theory information
architecture; three new molecule guide pages; seven new section index pages; one
new validation note; conversion of three pilot notes to rendered mathematics; the
written convention and a test that enforces it.

**Out of scope.** Converting the other nineteen notes. That is deliberate — see
"Rollout" below. No library code changes. No changes to `docs/api/`, which is
autodoc-generated. No move of any existing `docs/physics/*.md` file.

## Constraints that shaped the design

**C1 — 443 lines across 165 tracked files cite `docs/physics/<name>.md` paths.** They appear in
`CLAUDE.md`, in library source (`projects/n2_resonance/potential.py` and others),
in tests, in `docs/api/*.md`, in `.claude/skills/`, and in the notes themselves.
Physically reorganising the notes into subdirectories would mean rewriting all of
them, with a rename-shaped diff obscuring the real change. **Therefore no note
moves; only navigation changes.**

**C2 — Sphinx permits a document in exactly one toctree.** Sphinx 9.1 emits
`toc.multiple_toc_parents` when a document is reachable from two, and the docs
workflow builds with `-W`, so that warning is fatal. A note therefore has exactly
one home in the sidebar, and any second reference to it must be a prose `{doc}`
link. This is what rules out "the method notes appear under both the technical
tree and their molecule".

**C3 — The notes are read in two renderers.** They are the working notes of the
repository, reached as files in a clone (`mastering-github`'s rule that main must
stand alone) and as pages on the site. GitHub's Markdown renderer handles
`$...$` and `$$...$$` natively; it does not handle `conf.py`-defined MathJax
macros, and it renders MyST directives such as `{dropdown}` as visible junk.
**Therefore the physics notes use only constructs both renderers understand.**

**C4 — Furo collapses toctree parents, not captions.** A `:caption:` renders as a
permanently expanded heading — today's problem. Collapsibility requires each
group to be a real page owning its own `toctree`.

## Decisions

### D1 — Portable mathematics in the notes; directives only in site-first pages

Files under `docs/physics/` may use: standard Markdown, `$inline$`, `$$display$$`,
`\begin{aligned}` *inside* `$$`, and `\tag{n}`. They may not use `sphinx-design`
directives or custom macros.

Files under `docs/molecules/` and the section index pages are site-first — nobody
reads them raw — and may use the full `sphinx-design` toolkit.

Rejected: defining MathJax macros (`\Psip`, `\cprod`) in `conf.py`. They were
selected during brainstorming and then dropped, because C3 confines them to the
handful of site-first pages, which is a great deal of configuration for almost no
use and creates two dialects of mathematics in one tree.

### D2 — No MyST `amsmath`; multi-line math uses `aligned` inside `$$`

The MyST `amsmath` extension parses a bare `\begin{align}` block, which GitHub
does not render. `$$ \begin{aligned} ... \end{aligned} $$` produces the same
output under MathJax and renders in both places. One form, no dialect to police.

### D3 — Literature equation numbers become `\tag{}`

`nonlocal-resonance-model.md` alone carries 91 `Eq.` references, currently written
as right-aligned comments inside fenced blocks
(`p. 012710-3, Eq. (15)`). The equation number moves into the mathematics as
`\tag{15}`; the page locator moves into the surrounding prose, which is where the
`mastering-references` skill wants it. `\tag` is core MathJax and renders on
github.com, so the number survives both views.

### D4 — Molecule guides are self-contained pages, not toctree parents

Every N₂-named note (`n2-2d-td-cross-section` and the rest) is a *method* note
that happens to report its results on N₂; the order-3 Padé finding is not an N₂
fact. Under C2, letting a molecule own those notes would remove them from the
technical tree, which is where a reader looking for "how does the TD extractor
work" will go. So a molecule owns no notes.

The consequence, stated plainly because it departs from the original request: the
molecule entries are **not** collapsible parents in the sidebar. They are pages
whose own section headings expand beneath them when active, and whose long
parameter tables and result sets collapse in-page via `{dropdown}`. If true
`▸ N₂` parents are wanted later, the alternative is to re-home the method notes
under molecules, at the cost above.

### D5 — Three pilot notes, then note-by-note conversion

A transcription error in an equation is invisible to every test in this
repository: nothing executes the mathematics in a note. A 380 KB single-pass
sweep would therefore trade a visible defect for an invisible one. Three notes
are converted now, chosen to stress the convention; the rest are converted as
they are next touched.

### D6 — `validation-harnesses` is a flat entry, not a section of one

The harness note was originally given its own "Validation & tooling" group. A
collapsible section holding a single page is nav overhead with no navigation in
it, and the obvious ways to pad it are worse: `optimization-targets` is an open
direction, not a validation harness, and moving it would leave "Open directions"
a singleton instead.

The note is cross-cutting reference material — it describes what every harness in
`validation/` gates, across all molecules and methods — so it sits at the top
level of the technical toctree beside the Theory landing page, where a
cross-cutting document belongs. Seven collapsible groups, two flat entries.

## The information architecture

`docs/index.md` gains two captioned toctrees in place of the current single one:

```
Contents
    Getting started
    API reference
    Related work

Theory — technical
    Theory notes                    docs/physics/README.md (landing)
    Validation harnesses            validation-harnesses (new, flat — see D6)
  ▸ Discretisation                  femdvr-ecs · nd-tensor-hamiltonian ·
                                    discretisation-tuning
  ▸ Linear algebra & solvers        mumps-sparse-backend · ti-energy-sweep-reuse ·
                                    shift-invert-eigensolver
  ▸ The scattering engine           qscat-core-scattering · n2-resonance ·
                                    n2-cross-section · n2-2d-cross-section
  ▸ Dissociation & approximations   diatomic-ve-cross-sections · h2plus-dr ·
                                    nonlocal-resonance-model
  ▸ Time-dependent routes           n2-td-cross-section · n2-2d-td-cross-section ·
                                    td-extractors · td-da
  ▸ Resonances & levels             lcp-resonance-levels · exact-2d-resonances ·
                                    h2plus-resonance-states
  ▸ Open directions                 angular-coupled-channels · optimization-targets

Theory — molecules
    N₂ — the benchmark target
    NO / F₂ — dissociative attachment
    H₂⁺ — the ionic target
```

Every one of the 22 existing notes keeps its path. Each `▸` entry is a new index
page under `docs/physics/` owning a `toctree` of its notes; its prose is the
one-line description that already exists in `docs/physics/README.md`, so the
seven index pages are largely a split of that file rather than new writing.
`docs/physics/README.md` remains the Theory landing page and keeps portable
Markdown (it is read raw as often as any note).

## New files

| Path | Kind | Content |
|---|---|---|
| `docs/physics/discretisation.md` | index | orientation + `toctree` of 3 notes |
| `docs/physics/solvers.md` | index | orientation + `toctree` of 3 notes |
| `docs/physics/engine.md` | index | orientation + `toctree` of 4 notes |
| `docs/physics/dissociation.md` | index | orientation + `toctree` of 3 notes |
| `docs/physics/time-dependent.md` | index | orientation + `toctree` of 4 notes |
| `docs/physics/resonances.md` | index | orientation + `toctree` of 3 notes |
| `docs/physics/open-directions.md` | index | orientation + `toctree` of 2 notes |
| `docs/physics/validation-harnesses.md` | note | what each harness in `validation/` gates, and how to run it |
| `docs/molecules/n2.md` | guide | site-first N₂ guide |
| `docs/molecules/no-f2.md` | guide | site-first NO/F₂ guide |
| `docs/molecules/h2plus.md` | guide | site-first H₂⁺ guide |
| `tests/test_docs_portability.py` | test | enforces D1 |

### The validation note

`docs/physics/validation-harnesses.md` answers the "code/repo validations"
half of the request. Today this material exists only scattered through
`CLAUDE.md`. It records, per harness: what it gates, against which oracle, at
what tolerance, and the command that runs it —
`validation/n2` (Houfek `CSVE.V00.J00` anchors, groups C5/D1/E1/F1),
`validation/h2plus` (the published ω_i^j table, the σ_DR sweep),
`validation/diatomic` (no independent golden data; the exact solver is the
oracle), and `validation/tuning` (the `C = 0.10` calibration and the eMoScat deck
gates). It also states where the `@slow` boundary falls and why the F1 group
reports recorded NOTE rows rather than running a ~210 s propagation in-harness.

### The molecule guides

Each guide is one page carrying: the model form and parameters from
`qscat.model` (inside a `{dropdown}`, since a parameter table is reference
material and not something to read top to bottom); what has actually been
computed on that molecule and the headline numbers with their measured
tolerances; the committed figures from `docs/physics/figures/`; and a `{grid}` of
cards linking into the method notes that report on it.

The NRM results that arrived on main are part of this content: the N₂ guide
carries `n2-ve-nrm-vs-exact.png`, and the NO/F₂ guide carries
`f2-da-nrm-vs-lcp-vs-exact.png` together with the honest pair of findings — NRM
reproduces the exact oracle on F₂ to 0.06–1.9 %, and collapses on NO by 5–8
orders in dissociative attachment, unresolved.

## Build configuration

`pyproject.toml`, `[dependency-groups].docs` gains `sphinx-copybutton`,
`sphinx-design`, `sphinxext-opengraph`. All are pure-Python wheels; the Pages
workflow needs no change beyond its existing `uv sync --package qscat --group docs`.

`docs/conf.py`:

| Change | Reason |
|---|---|
| `extensions += sphinx_copybutton, sphinx_design, sphinxext.opengraph, sphinx.ext.githubpages` | copy buttons; cards/dropdowns; link previews; `.nojekyll` |
| `copybutton_prompt_text` matching `>>> `, `... `, `$ ` | copying a shell or doctest block yields runnable text |
| `myst_enable_extensions += "deflist"` (keep `dollarmath`, `colon_fence`; **not** `amsmath`) | D2 |
| `numfig = True`, `math_number_all = False`, `math_eqref_format = "Eq. {number}"` | number only equations a note actually references |
| `html_baseurl`, `ogp_site_url`, `ogp_image` | shared links preview |
| `html_theme_options`: `source_repository` / `source_branch` / `source_directory`; light and dark brand colours | "Edit on GitHub" link; theming |
| `-W --keep-going` unchanged | C2 depends on it |

## The convention

Written into `.claude/skills/qscat-conventions/`, referenced from
`CONTRIBUTING.md`:

- **Backticks are for code identifiers only** — `ve_cross_section`, `SparseLU`,
  `backend="mumps"`. A backticked `sigma` that means σ is a defect.
- **Inline mathematics** is `$\sigma_{v\to v'}$`. **Display mathematics** is
  `$$...$$` on its own lines; multi-line via `\begin{aligned}`.
- **Labelled equations** (`$$...$$ (eq-driven-ls)` plus the `{eq}` role) only
  where a note genuinely writes "Eq. (3)". The trailing label is literal text on
  GitHub, so it is a deliberate cost rather than a default.
- **No custom macros.** `\Psi^{(+)}`, `\sigma_\mathrm{DA}` spelled out.
- **A canonical symbol table** lives in the conventions skill, so that θ is the
  ECS angle, μ the reduced mass, `R0` the ECS pivot, and so on, consistently
  across 22 notes.
- Atomic units are stated once per note, not per equation.

### Fence audit

The 32 unlabelled fenced blocks are triaged into three buckets: genuine equations
become `$$` display mathematics; terminal transcripts become `console`; algorithm
sketches and schematic layouts become `text`. The 9 `python`, 6 `bash` and 2 `yaml` fences
are already correct and gain only copy buttons. A fence retagged to `python`
becomes visible to the doctest builder, so the audit must not retag prose.

## Pilot notes

1. **`femdvr-ecs.md`** — foundational. Its symbols (θ, `R0`, `z(x)`, the GLL
   nodes, the weight Jacobian `hz`) seed the canonical symbol table that every
   other note inherits, so it is converted first and the table is derived from it.
2. **`nonlocal-resonance-model.md`** — the hardest case at 68 KB, 13 fences and 91
   equation references, with bra-kets, integrals and a nonlocal kernel. It is the
   test of D2 and D3 under maximum stress. If the convention survives here it
   survives everywhere.
3. **`n2-2d-cross-section.md`** — the driven Lippmann–Schwinger block and the
   canonical $\sigma = 4\pi^3|T|^2/k^2$, plus result tables and a committed
   figure, so the pilot covers a mixed prose/mathematics/table/figure page.

## Rollout

The remaining nineteen notes are converted as each is next touched. The spec's
companion implementation plan carries the checklist; a converted note is one where
`test_docs_portability.py` passes *and* no ASCII equation remains outside a code
fence. No deadline is attached: an unconverted note renders exactly as it does
today, so the tree is correct at every intermediate state.

## Verification

Three gates, all runnable, all in CI:

1. `uv run sphinx-build -b html -W --keep-going docs docs/_build/html` stays
   clean. Under C2 this is also what proves every note has exactly one home: a
   note left out of the new index pages, or listed in two of them, fails the
   build rather than silently vanishing from the sidebar.
2. `uv run sphinx-build -b doctest docs docs/_build/doctest` unchanged — the
   fence audit must not break a runnable example.
3. `tests/test_docs_portability.py` — makes D1 testable rather than aspirational.
   It checks the **notes**, defined as every `docs/physics/*.md` except
   `README.md` and the seven section index pages, which the test names in an
   explicit `SITE_FIRST_PAGES` allowlist alongside `docs/molecules/`. Keeping the
   exemption an explicit list rather than a path rule is deliberate: adding an
   eighth site-first page is then a visible edit to the test, not a silent opt-out.

   Three assertions per note:

   - no `sphinx-design` directive — `{dropdown}`, `{grid}`, `{grid-item-card}`,
     `{tab-set}`, `{tab-item}`;
   - no `\begin{align}`, `\begin{gather}` or `\begin{equation}` outside a `$$`
     block (D2 — GitHub does not render these bare);
   - no `\newcommand`, `\def` or `\renewcommand` anywhere, and
     `conf.py` defines no `mathjax3_config` macros (D1 — this is the whole of
     "no custom macros", stated as something a test can check, rather than the
     unbounded "not in MathJax's core set").

   Without this test the portability rule decays the first time someone reaches
   for a better-looking directive.

Manual check, once, at the end: build locally and read the three pilot notes on
the rendered site and on github.com, confirming the mathematics is identical in
both.

## Risks

| Risk | Mitigation |
|---|---|
| An equation is mistranscribed during conversion; no test can catch it | Only 3 notes convert here, each reviewed against its source equations; the literature-anchored ones (`nonlocal-resonance-model`) are checked against `reference/literature/houfek-2008-pra77-012710.md` under `mastering-references` |
| The eight index pages become stale stubs nobody maintains | Their prose is moved from `docs/physics/README.md`, not newly invented, so there is one description per note rather than two |
| `sphinx-design` output looks wrong in Furo's dark theme | The theming change is in the same work; the manual check covers both themes |
| MathJax is CDN-loaded, so an offline build shows raw LaTeX | Accepted. The published site is the deliverable; offline readers get the same portable source GitHub renders |

## Out of scope, recorded for later

- Converting the remaining twenty notes (tracked by the plan's checklist).
- True collapsible per-molecule sidebar parents, which requires re-homing method
  notes under molecules — see D4.
- `sphinx-build -b linkcheck` in CI. The tree cites many external papers and a
  flaky publisher redirect would fail unrelated builds.
