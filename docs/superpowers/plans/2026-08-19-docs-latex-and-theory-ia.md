# Rendered Mathematics and Theory IA Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the qscat documentation site render its mathematics and give the Theory section a navigable, collapsible structure with per-molecule entry points.

**Architecture:** No existing note moves — 443 lines across 165 tracked files cite `docs/physics/<name>.md` paths, so navigation changes only. Seven new section index pages own collapsible `toctree`s; three new site-first molecule guides use `sphinx-design`; the notes themselves stay portable to GitHub's Markdown renderer (plain `$...$`, `aligned` inside `$$`, `\tag{}`, no macros, no directives), enforced by a new test. Three pilot notes convert to rendered mathematics now; nineteen convert as they are next touched.

**Tech Stack:** Sphinx 9.1, MyST-Parser 5.1, Furo, `sphinx-copybutton` 0.5.2, `sphinx-design` 0.7.0, `sphinxext-opengraph` 0.13.0, uv dependency groups, pytest.

**Spec:** `docs/superpowers/specs/2026-08-19-docs-latex-and-theory-ia-design.md`

## Global Constraints

- **No file under `docs/physics/` moves or is renamed.** (Spec C1.)
- **A document appears in exactly one `toctree`.** Sphinx 9.1 emits `toc.multiple_toc_parents` otherwise, and the docs workflow builds with `-W`, making it fatal. (Spec C2.)
- **Notes are portable to two renderers.** Files under `docs/physics/` may use only: standard Markdown, `$inline$`, `$$display$$`, `\begin{aligned}` *inside* `$$`, and `\tag{n}`. They may **not** use `sphinx-design` directives or any custom macro. (Spec C3, D1, D2.)
- **Site-first pages are exempt** from the rule above: `docs/molecules/*.md` and the seven section index pages named in `SITE_FIRST_PAGES` (Task 2).
- **No MathJax macros anywhere.** `conf.py` defines no `mathjax3_config` macros. (Spec D1.)
- **Backticks are for code identifiers only.** A backticked `sigma` that means σ is a defect. (Spec, "The convention".)
- **Atomic units throughout**, stated once per note, not per equation.
- The docs build command is `uv run sphinx-build -b html -W --keep-going docs docs/_build/html`. It must be clean at the end of every task **except Task 3**, which introduces a link to the note Task 4 creates. Tasks 3 and 4 are therefore a pair: do not stop between them, and do not push with only Task 3 landed.
- Verified before planning: `sphinx-design` 0.7.0 declares `sphinx>=7,<10`, `sphinx-copybutton` 0.5.2 `sphinx>=1.8`, `sphinxext-opengraph` 0.13.0 `Sphinx>=6.0` — all satisfy the pinned Sphinx 9.1. No downgrade is needed or permitted.

---

### Task 1: Build configuration and dependencies

**Files:**
- Modify: `pyproject.toml` (the `[dependency-groups].docs` list, lines 29-33)
- Modify: `docs/conf.py`

**Interfaces:**
- Consumes: nothing.
- Produces: the `sphinx_design`, `sphinx_copybutton`, `sphinxext.opengraph` extensions, available to every later task. `numfig = True`, `math_number_all = False`, `math_eqref_format = "Eq. {number}"`.

- [ ] **Step 1: Add the three packages to the docs dependency group**

In `pyproject.toml`, replace the `docs` group:

```toml
docs = [
    "sphinx>=8",
    "myst-parser>=4",
    "furo>=2024.8",
    # Copy button on code blocks; card/dropdown/tab directives for the
    # site-first molecule guides; OpenGraph tags so shared links preview.
    "sphinx-copybutton>=0.5.2",
    "sphinx-design>=0.7",
    "sphinxext-opengraph>=0.13",
]
```

- [ ] **Step 2: Sync and confirm the three import**

Run: `uv sync --package qscat --group docs`
Then: `uv run --no-sync python -c "import sphinx_design, sphinx_copybutton, sphinxext.opengraph; print('ok')"`
Expected: `ok`

- [ ] **Step 3: Update `docs/conf.py`**

Replace the `extensions` list and the settings below it. Note `sphinx.ext.mathjax` and MyST `dollarmath` are **already present** — this adds to them, it does not introduce math support.

```python
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.githubpages",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinxext.opengraph",
]
```

Add after the existing `myst_enable_extensions` line:

```python
# dollarmath: $...$ and $$...$$. colon_fence: ::: directives. deflist:
# definition lists for the symbol table.
#
# NOT amsmath: it parses a BARE \begin{align} block, which github.com does
# not render. The notes are read in two renderers (spec C3), so multi-line
# mathematics is written as $$ \begin{aligned} ... \end{aligned} $$, which
# both handle. One form, no dialect to police.
myst_enable_extensions = ["dollarmath", "colon_fence", "deflist"]

# Number only the equations a note actually cross-references with {eq};
# numbering all of them adds visual noise to notes that never refer back.
numfig = True
math_number_all = False
math_eqref_format = "Eq. {number}"

# No mathjax3_config macros are defined, deliberately: a macro renders on
# the site but not on github.com, and the notes must render in both
# (spec D1). Equations spell \Psi^{(+)} and \sigma_\mathrm{DA} out longhand.

# Strip interpreter and shell prompts so a copied block is runnable.
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True

html_baseurl = "https://vanamartin.github.io/qscat/"
ogp_site_url = html_baseurl
ogp_description_length = 200
```

Replace the `html_theme` block at the end:

```python
html_theme = "furo"
html_title = f"qscat {release}"
html_theme_options = {
    "source_repository": "https://github.com/VanaMartin/qscat/",
    "source_branch": "main",
    "source_directory": "docs/",
    "light_css_variables": {
        "color-brand-primary": "#1a5490",
        "color-brand-content": "#1a5490",
    },
    "dark_css_variables": {
        "color-brand-primary": "#7cb3e8",
        "color-brand-content": "#7cb3e8",
    },
}
```

- [ ] **Step 4: Build and confirm the tree is still warning-clean**

Run: `uv run sphinx-build -b html -W --keep-going docs docs/_build/html`
Expected: exit 0, `build succeeded`. If it fails with `Extension error`, one of the three packages is incompatible with Sphinx 9.1 — stop and report; do **not** downgrade Sphinx.

- [ ] **Step 5: Confirm mathematics now renders end to end**

Create a scratch file `docs/_mathcheck.md` containing exactly:

```markdown
# Math check

Inline $\sigma_{v\to v'} = 4\pi^3 |T|^2 / k^2$ and display:

$$P = 1 - Q, \qquad Q = |\phi_d\rangle\langle\phi_d| \tag{15}$$
```

Run: `uv run sphinx-build -b html -q docs docs/_build/html 2>&1 | head`
Then: `grep -c 'mathjax' docs/_build/html/_mathcheck.html`
Expected: a non-zero count, confirming MathJax is wired to the page.
Then delete the scratch file: `rm docs/_mathcheck.md docs/_build/html/_mathcheck.html`

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml uv.lock docs/conf.py
git commit -m "docs(build): copy buttons, design directives and equation numbering

sphinx.ext.mathjax and MyST dollarmath were already loaded and unused. This
adds what was actually missing: sphinx-copybutton, sphinx-design for the
site-first pages, OpenGraph tags, Furo source links and brand colours, and
equation numbering that only numbers what a note references.

Deliberately NOT the MyST amsmath extension, and deliberately no
mathjax3_config macros -- both render on the site but not on github.com,
and the notes are read in both."
```

---

### Task 2: The portability test

**Files:**
- Create: `tests/test_docs_portability.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `SITE_FIRST_PAGES: frozenset[str]` (basenames exempt from the note rules), and three detector functions `find_design_directives(text) -> list[str]`, `find_bare_math_environments(text) -> list[str]`, `find_macro_definitions(text) -> list[str]`, each returning the offending snippets. Task 3 adds the seven index page names to `SITE_FIRST_PAGES`.

This is the test that makes the portability rule real rather than aspirational (spec D1, verification gate 3). TDD applies to the **detectors**: they are unit-tested against known-bad and known-good strings, which genuinely fail before the functions exist. The tree scan is then a thin loop over those detectors.

- [ ] **Step 1: Write the failing test**

Create `tests/test_docs_portability.py`:

```python
"""The docs/physics notes must render in two places.

They are the working notes of the repository -- read as files in a clone --
and pages on the published site. GitHub's Markdown renderer handles
`$...$` and `$$...$$`, but NOT conf.py-defined MathJax macros, and it shows
MyST directives such as `{dropdown}` as visible junk.

So the notes are restricted to constructs both renderers understand. This
test enforces that; see docs/superpowers/specs/2026-08-19-docs-latex-and-
theory-ia-design.md, decisions D1 and D2.

Site-first pages -- the molecule guides and the section index pages, which
nobody reads raw -- are exempt, by an explicit allowlist rather than a path
rule, so adding an eighth is a visible edit here and not a silent opt-out.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PHYSICS_DIR = REPO_ROOT / "docs" / "physics"
CONF_PY = REPO_ROOT / "docs" / "conf.py"

# Exempt from the note rules: read on the site, never raw.
SITE_FIRST_PAGES = frozenset({"README.md"})

_DESIGN_DIRECTIVES = (
    "dropdown",
    "grid",
    "grid-item-card",
    "tab-set",
    "tab-item",
    "card",
)

# ``` fenced blocks: their contents are not prose and are skipped.
_FENCE_RE = re.compile(r"^```.*?^```", re.DOTALL | re.MULTILINE)
_DISPLAY_MATH_RE = re.compile(r"\$\$.*?\$\$", re.DOTALL)


def _strip_fences(text: str) -> str:
    return _FENCE_RE.sub("", text)


def find_design_directives(text: str) -> list[str]:
    """Return every sphinx-design directive opener found in `text`."""
    found = []
    for name in _DESIGN_DIRECTIVES:
        # MyST spells these ```{dropdown} or :::{dropdown}
        pattern = re.compile(r"^(?:```+|:::+)\{" + re.escape(name) + r"\}", re.MULTILINE)
        found.extend(m.group(0) for m in pattern.finditer(text))
    return found


def find_bare_math_environments(text: str) -> list[str]:
    """Return `\\begin{align|gather|equation}` used OUTSIDE a `$$` block.

    GitHub renders these only when wrapped in `$$`; bare, they show as
    literal LaTeX source.
    """
    outside = _DISPLAY_MATH_RE.sub("", _strip_fences(text))
    pattern = re.compile(r"\\begin\{(align\*?|gather\*?|equation\*?)\}")
    return [m.group(0) for m in pattern.finditer(outside)]


def find_macro_definitions(text: str) -> list[str]:
    """Return every LaTeX macro definition; none may appear in a note."""
    pattern = re.compile(r"\\(newcommand|renewcommand|def)\b")
    return [m.group(0) for m in pattern.finditer(text)]


def _notes() -> list[Path]:
    return sorted(
        p for p in PHYSICS_DIR.glob("*.md") if p.name not in SITE_FIRST_PAGES
    )


# --- detector unit tests -------------------------------------------------


def test_find_design_directives_flags_a_dropdown():
    text = "Some prose.\n\n```{dropdown} Parameters\ncontent\n```\n"
    assert find_design_directives(text) == ["```{dropdown}"]


def test_find_design_directives_flags_colon_fence_form():
    text = ":::{grid} 2\n:::\n"
    assert find_design_directives(text) == [":::{grid}"]


def test_find_design_directives_passes_a_plain_note():
    text = "Prose with $\\sigma$ and a table.\n\n| a | b |\n|---|---|\n"
    assert find_design_directives(text) == []


def test_find_bare_math_environments_flags_a_bare_align():
    text = "\\begin{align} a &= b \\end{align}\n"
    assert find_bare_math_environments(text) == ["\\begin{align}"]


def test_find_bare_math_environments_allows_aligned_inside_display_math():
    text = "$$\n\\begin{aligned} a &= b \\end{aligned}\n$$\n"
    assert find_bare_math_environments(text) == []


def test_find_bare_math_environments_allows_align_inside_display_math():
    text = "$$\n\\begin{align} a &= b \\end{align}\n$$\n"
    assert find_bare_math_environments(text) == []


def test_find_macro_definitions_flags_newcommand():
    text = "\\newcommand{\\Psip}{\\Psi^{(+)}}\n"
    assert find_macro_definitions(text) == ["\\newcommand"]


def test_find_macro_definitions_passes_plain_latex():
    text = "$\\Psi^{(+)}$ and $\\sigma_\\mathrm{DA}$\n"
    assert find_macro_definitions(text) == []


# --- tree scan -----------------------------------------------------------


@pytest.mark.parametrize("note", _notes(), ids=lambda p: p.name)
def test_note_has_no_design_directives(note: Path):
    found = find_design_directives(note.read_text())
    assert not found, (
        f"{note.name} uses sphinx-design {found}, which github.com renders as "
        f"visible junk. Directives belong in docs/molecules/ or a section "
        f"index page (see SITE_FIRST_PAGES)."
    )


@pytest.mark.parametrize("note", _notes(), ids=lambda p: p.name)
def test_note_has_no_bare_math_environments(note: Path):
    found = find_bare_math_environments(note.read_text())
    assert not found, (
        f"{note.name} uses {found} outside a $$ block; github.com will show it "
        f"as literal LaTeX. Wrap it: $$ \\begin{{aligned}} ... \\end{{aligned}} $$"
    )


@pytest.mark.parametrize("note", _notes(), ids=lambda p: p.name)
def test_note_defines_no_macros(note: Path):
    found = find_macro_definitions(note.read_text())
    assert not found, f"{note.name} defines a LaTeX macro {found}; macros do not render on github.com."


def test_conf_py_defines_no_mathjax_macros():
    """The other half of 'no custom macros' (spec D1)."""
    assert "mathjax3_config" not in CONF_PY.read_text(), (
        "conf.py defines mathjax3_config. Macros defined there render on the "
        "site but not on github.com, and the notes are read in both."
    )


def test_the_scan_actually_covers_the_notes():
    """Guard against the glob silently matching nothing."""
    assert len(_notes()) >= 20
```

- [ ] **Step 2: Run the detector tests to verify they fail**

Run: `uv run pytest tests/test_docs_portability.py -q`
Expected: collection succeeds and the whole file passes **except** — confirm first that it genuinely exercises the detectors by temporarily breaking one: change `find_macro_definitions` to `return []` and re-run.
Expected: `test_find_macro_definitions_flags_newcommand` FAILS. Restore the function.

- [ ] **Step 3: Run the full file**

Run: `uv run pytest tests/test_docs_portability.py -q`
Expected: all pass. The tree scan passes today because no note uses a directive or a macro yet — that is the point: the test locks in a property that currently holds and would otherwise silently decay.

- [ ] **Step 4: Verify the scan has teeth against the real tree**

Run: `printf '\n```{dropdown} x\ny\n```\n' >> docs/physics/femdvr-ecs.md && uv run pytest tests/test_docs_portability.py -q; git checkout docs/physics/femdvr-ecs.md`
Expected: `test_note_has_no_design_directives[femdvr-ecs.md]` FAILS, then the checkout restores the file. Confirm `git status` is clean afterwards.

- [ ] **Step 5: Commit**

```bash
git add tests/test_docs_portability.py
git commit -m "test(docs): keep the physics notes renderable on github.com

The notes are read both as files in a clone and as pages on the site.
GitHub renders \$...\$ but not conf.py MathJax macros, and shows MyST
directives as junk -- so the notes are restricted to what both renderers
handle, and this test is what stops that rule from decaying the first time
someone reaches for a better-looking directive.

Site-first pages are exempt through an explicit allowlist, so adding one is
a visible edit here rather than a silent opt-out."
```

---

### Task 3: Section index pages and the navigation rewiring

**Files:**
- Create: `docs/physics/discretisation.md`, `docs/physics/solvers.md`, `docs/physics/engine.md`, `docs/physics/dissociation.md`, `docs/physics/time-dependent.md`, `docs/physics/resonances.md`, `docs/physics/open-directions.md`
- Modify: `docs/index.md` (the `Theory` toctree, lines 25-48)
- Modify: `docs/physics/README.md` (becomes the Theory landing page)
- Modify: `tests/test_docs_portability.py` (`SITE_FIRST_PAGES`)

**Interfaces:**
- Consumes: `SITE_FIRST_PAGES` from Task 2.
- Produces: seven index documents, each owning a `toctree`. Task 4 adds `physics/validation-harnesses` as a flat entry to `docs/index.md`; Task 5-6 add `docs/molecules/*` under a second caption.

The one-line descriptions below are **moved** from `docs/physics/README.md`, not newly written, so each note keeps exactly one description in the tree. `README.md` is left as a short landing page pointing at the seven sections.

- [ ] **Step 1: Add the seven index pages to the test allowlist**

In `tests/test_docs_portability.py`, replace `SITE_FIRST_PAGES`:

```python
# Exempt from the note rules: read on the site, never raw.
SITE_FIRST_PAGES = frozenset(
    {
        "README.md",
        "discretisation.md",
        "solvers.md",
        "engine.md",
        "dissociation.md",
        "time-dependent.md",
        "resonances.md",
        "open-directions.md",
    }
)
```

- [ ] **Step 2: Create the seven index pages**

`docs/physics/discretisation.md`:

```markdown
# Discretisation

How a continuous radial coordinate becomes a finite matrix, and how the grid
is chosen rather than guessed.

```{toctree}
:maxdepth: 1

femdvr-ecs
nd-tensor-hamiltonian
discretisation-tuning
```

- {doc}`femdvr-ecs` — the FEM-DVR grid with an exterior-complex-scaled tail,
  and the four analytic benchmarks that pin it down.
- {doc}`nd-tensor-hamiltonian` — the N-dimensional sparse tensor Hamiltonian.
- {doc}`discretisation-tuning` — deriving a grid from the potential instead of
  hand-picking element lengths, and where the 1-D probes are not sufficient.
```

`docs/physics/solvers.md`:

```markdown
# Linear algebra and solvers

The sparse factorizations and eigensolves the two-dimensional problems run
on, and the reuse tricks that make a sweep affordable.

```{toctree}
:maxdepth: 1

mumps-sparse-backend
ti-energy-sweep-reuse
shift-invert-eigensolver
```

- {doc}`mumps-sparse-backend` — the complex-symmetric MUMPS backend.
- {doc}`ti-energy-sweep-reuse` — reusing the symbolic factorization across an
  energy sweep.
- {doc}`shift-invert-eigensolver` — the eigenpairs nearest a complex shift, for
  resonances that sit in the interior of the spectrum. Validated in 1-D only.
```

`docs/physics/engine.md`:

```markdown
# The scattering engine

The model-independent solver and the time-independent routes to a cross
section.

```{toctree}
:maxdepth: 1

qscat-core-scattering
n2-resonance
n2-cross-section
n2-2d-cross-section
```

- {doc}`qscat-core-scattering` — the model-independent engine and the
  model/engine split it enforces.
- {doc}`n2-resonance` — locating the resonance pole.
- {doc}`n2-cross-section` — the one-dimensional time-independent route.
- {doc}`n2-2d-cross-section` — the exact two-dimensional driven solve, gated
  against independent published data.
```

`docs/physics/dissociation.md`:

```markdown
# Dissociation and approximations

The channels where the molecule comes apart, and the two reduced models —
local and nonlocal — measured against the exact solver.

```{toctree}
:maxdepth: 1

diatomic-ve-cross-sections
nonlocal-resonance-model
h2plus-dr
```

- {doc}`diatomic-ve-cross-sections` — NO and F₂, and the local-complex-potential
  approximation measured against the exact oracle.
- {doc}`nonlocal-resonance-model` — the rung above the LCP: a nonlocal,
  energy-dependent kernel for dissociative attachment. Reproduces the oracle on
  F₂; collapses on NO for reasons not established.
- {doc}`h2plus-dr` — dissociative recombination for an ionic target.
```

`docs/physics/time-dependent.md`:

```markdown
# Time-dependent routes

Wavepacket propagation as the second, independent route to the same cross
sections.

```{toctree}
:maxdepth: 1

n2-td-cross-section
n2-2d-td-cross-section
td-extractors
td-da
```

- {doc}`n2-td-cross-section` — wavepacket propagation in one dimension.
- {doc}`n2-2d-td-cross-section` — the exact two-dimensional time-dependent
  route, including why order-1 Crank–Nicolson was not enough.
- {doc}`td-extractors` — three energy extractors sharing one propagation.
- {doc}`td-da` — the dissociative-attachment generalization.
```

`docs/physics/resonances.md`:

```markdown
# Resonances and levels

Quasi-bound states: the Born–Oppenheimer approximation to them, the exact
two-dimensional poles, and how to tell a resonance from an artefact.

```{toctree}
:maxdepth: 1

lcp-resonance-levels
exact-2d-resonances
h2plus-resonance-states
```

- {doc}`lcp-resonance-levels` — Born–Oppenheimer quasi-bound levels in the
  complex curve.
- {doc}`exact-2d-resonances` — the same levels without the approximation: poles
  of the full 2-D S-matrix, and what the Born–Oppenheimer error actually
  measures on N₂.
- {doc}`h2plus-resonance-states` — the same comparison on H₂⁺, against a σ_DR
  sweep: the Born–Oppenheimer error sorted by regime, and the four "resonances"
  that turned out not to be.
```

`docs/physics/open-directions.md`:

```markdown
# Open directions

Work that is designed but parked, and where the remaining cost is.

```{toctree}
:maxdepth: 1

angular-coupled-channels
optimization-targets
```

- {doc}`angular-coupled-channels` — the parked angular extension.
- {doc}`optimization-targets` — where the remaining hot paths are.
```

- [ ] **Step 3: Rewrite `docs/physics/README.md` as the landing page**

Replace the whole file:

```markdown
# Theory notes

One note per method: the derivation, the equations, the unit conventions
(atomic units throughout), the validation evidence, and the literature it
comes from. These are the working notes behind the implementation, so they
record limitations and negative results as well as what works.

The notes are grouped into seven sections in the sidebar:

- {doc}`discretisation` — turning a coordinate into a matrix, and choosing the
  grid.
- {doc}`solvers` — the sparse factorizations and eigensolves underneath.
- {doc}`engine` — the model-independent solver and the time-independent routes.
- {doc}`time-dependent` — the wavepacket route to the same answers.
- {doc}`dissociation` — the channels where the molecule comes apart, and the
  reduced models measured against the exact solver.
- {doc}`resonances` — quasi-bound states, exactly and approximately.
- {doc}`open-directions` — designed but parked.

{doc}`validation-harnesses` cuts across all of them: what each harness in
`validation/` gates, and how to run it.

For a per-molecule view — what has been computed on N₂, NO, F₂ and H₂⁺, and
which notes report it — start from the Molecules section of the sidebar.
```

Note: `validation-harnesses` does not exist until Task 4, so this link will
break the `-W` build until then. That is expected and is why Steps 4-6 below
defer the build check to Task 4. Run the *test* now, the *build* after Task 4.

- [ ] **Step 4: Rewire `docs/index.md`**

Replace the second `toctree` (the one captioned `Theory`) with:

```markdown
```{toctree}
:maxdepth: 2
:caption: Theory — technical

physics/README
physics/discretisation
physics/solvers
physics/engine
physics/time-dependent
physics/dissociation
physics/resonances
physics/open-directions
```
```

Leave the `Contents` toctree above it untouched.

- [ ] **Step 5: Run the portability test**

Run: `uv run pytest tests/test_docs_portability.py -q`
Expected: all pass. The seven new index pages are in `SITE_FIRST_PAGES`, so their `toctree` directives do not trip the directive scan.

- [ ] **Step 6: Confirm every note is claimed exactly once**

Run:

```bash
uv run --no-sync python - <<'PY'
import pathlib, re, collections
idx = pathlib.Path("docs/physics")
notes = {p.stem for p in idx.glob("*.md")} - {"README"}
owners = collections.defaultdict(list)
for page in list(idx.glob("*.md")) + [pathlib.Path("docs/index.md")]:
    for block in re.findall(r"```\{toctree\}(.*?)```", page.read_text(), re.S):
        for line in block.splitlines():
            line = line.strip().removeprefix("physics/")
            if line and not line.startswith(":"):
                owners[line].append(page.name)
missing = sorted(n for n in notes if n not in owners)
dupes = {k: v for k, v in owners.items() if len(v) > 1}
print("unclaimed:", missing)
print("claimed twice:", dupes)
PY
```

Expected exactly:

```
unclaimed: []
claimed twice: {}
```

`validation-harnesses` does not exist yet, so it cannot appear in either list. A name under `unclaimed` was dropped from an index page — add it. A name under `claimed twice` will fail the `-W` build in Task 4 — remove one of the two entries.

- [ ] **Step 7: Commit**

```bash
git add docs/index.md docs/physics/README.md docs/physics/discretisation.md \
        docs/physics/solvers.md docs/physics/engine.md \
        docs/physics/dissociation.md docs/physics/time-dependent.md \
        docs/physics/resonances.md docs/physics/open-directions.md \
        tests/test_docs_portability.py
git commit -m "docs(theory): group the notes into collapsible sections

The Theory sidebar was one flat list of 22 entries. Furo collapses toctree
parents but not captions, so each group becomes a short index page owning
its own toctree -- the prose is moved from physics/README.md, which keeps
one description per note rather than creating a second.

No note moves: 443 lines across 165 tracked files cite these paths.

The build is intentionally red until the next commit adds
validation-harnesses, which README.md now links."
```

---

### Task 4: The validation-harnesses note

**Files:**
- Create: `docs/physics/validation-harnesses.md`
- Modify: `docs/index.md` (add the flat entry)

**Interfaces:**
- Consumes: the `Theory — technical` toctree from Task 3.
- Produces: `physics/validation-harnesses`, linked from `docs/physics/README.md`. This is the last change needed to make the `-W` build green again.

This note is a **note**, not a site-first page — it is subject to the portability rules. Spec D6: it is a flat top-level entry, not a section of one.

- [ ] **Step 1: Gather the facts from the source, not from memory**

Read each harness's module docstring and gate constants before writing:

```bash
sed -n '1,40p' validation/n2/experiment.py
sed -n '1,30p' validation/n2/exact2d.py
grep -rn "GATED_RTOL\|ANCHOR_FACTOR" validation/ | head
sed -n '1,30p' validation/tuning/calibrate.py
sed -n '1,30p' validation/h2plus/reference_levels.py
sed -n '1,30p' validation/diatomic/config.py
```

Every tolerance, group name and command in the note must come from what these
print. Do not carry a number over from `CLAUDE.md` without confirming it here —
`CLAUDE.md` is a summary and can lag.

- [ ] **Step 2: Write the note**

Create `docs/physics/validation-harnesses.md` with this structure, filling the
table from Step 1's output:

```markdown
# Validation harnesses

Every method in qscat is checked against an analytic benchmark, a conservation
law, a convergence study, or an independent reference. This note is the index
of those checks: what each harness in `validation/` gates, against which
oracle, at what tolerance, and the command that runs it.

Atomic units throughout.

## The harnesses

| Harness | Oracle | Gate | Run it |
|---|---|---|---|
| `validation/n2` | Houfek's independent `CSVE.V00.J00` data | `GATED_RTOL = 1e-3` for the exact 2-D solver (`validation/n2/exact2d.py:82`); the LCP is held only to the cross-model `ANCHOR_FACTOR = 3.0` band (`validation/n2/reference.py:41`) | `uv run python -m validation.n2.experiment` |
| `validation/h2plus` | the published ω_i^j level table (`reference_levels`) | overlap verdicts against the Born–Oppenheimer basis | `uv run python -m validation.h2plus.bo_overlap` |
| `validation/diatomic` | none published — the exact 2-D solver **is** the oracle | well-posedness (finite, non-negative) plus a per-molecule resonance floor `_V1_FLOOR = {"F2": 0.15, "NO": 10.0}` (`test_diatomic.py:101`, asserted at :133). NOT `ANCHOR_FACTOR` — that name appears in this package only as prose at `test_ve_nrm.py:374` | `uv run pytest validation/diatomic -q` |
| `validation/tuning` | the eMoScat per-molecule decks | `_2D_SPOT_CHECK_RTOL = 0.02` (`test_emoscat_decks.py:257`), `_2D_CONVERGENCE_RTOL = 0.15` (`test_resonance_aware.py:107`) | `uv run python -m validation.tuning.calibrate` |

Confirm each constant still holds before committing — Step 1's grep is the
source, this table is a transcription of it.

## What each one actually establishes

*(One short subsection per harness: the claim it supports, and the claim it
does NOT support. For `validation/diatomic`, state plainly that no independent
golden data exists for NO or F₂ — only N₂ has Houfek's — so agreement there is
self-consistency, not external validation.)*

## Where the `@slow` boundary falls

*(Explain that the harness has a per-group time budget, so the N₂ F1 group
reports the time-dependent 2-D agreement as recorded, cited NOTE rows rather
than running a ~210-250 s propagation in-harness; the live PASS/FAIL gate for
that claim lives in `projects/n2_2d_td_cross_section/test_td_cross_section.py`
under `@pytest.mark.slow`. Confirm the timing figure from the source before
quoting it.)*

## Reading a NOTE row

*(Explain the PASS / NOTE distinction: a NOTE row is a documented, reported
limitation carried in the harness output rather than a silent omission.)*
```

- [ ] **Step 3: Add the flat entry to `docs/index.md`**

In the `Theory — technical` toctree, insert `physics/validation-harnesses`
immediately after `physics/README`:

```markdown
physics/README
physics/validation-harnesses
physics/discretisation
```

- [ ] **Step 4: Run the portability test**

Run: `uv run pytest tests/test_docs_portability.py -q`
Expected: all pass, now including `[validation-harnesses.md]` in the scan (it is a note, not exempt).

- [ ] **Step 5: Build — the tree must be green again**

Run: `uv run sphinx-build -b html -W --keep-going docs docs/_build/html`
Expected: exit 0, `build succeeded`. This is the first green build since Task 3, and it proves both that every note has exactly one toctree home and that `README.md`'s link resolves.

- [ ] **Step 6: Commit**

```bash
git add docs/physics/validation-harnesses.md docs/index.md
git commit -m "docs(theory): index what every validation harness gates

The claims each harness supports lived only in CLAUDE.md, scattered. This
collects them: oracle, tolerance, command, and -- for validation/diatomic --
the fact that no independent golden data exists for NO or F2 at all, so
agreement there is self-consistency rather than external validation.

Flat top-level entry rather than a section of one (spec D6)."
```

---

### Task 5: The N₂ molecule guide (establishes the template)

**Files:**
- Create: `docs/molecules/n2.md`
- Modify: `docs/index.md` (new `Theory — molecules` caption)

**Interfaces:**
- Consumes: `sphinx_design` from Task 1.
- Produces: the guide template — a `{grid}` of `{grid-item-card}`s linking to method notes, a `{dropdown}` holding the parameter table — that Task 6 replicates for NO/F₂ and H₂⁺.

Site-first page: the full `sphinx-design` toolkit is allowed, and `docs/molecules/` is outside the portability scan by construction (the scan globs `docs/physics/*.md`).

- [ ] **Step 1: Gather the model parameters and headline numbers from source**

```bash
uv run --no-sync python -c "
from qscat.model import N2
m = N2 if not callable(N2) else N2()
print(type(m).__name__)
print({k: v for k, v in vars(m).items() if not k.startswith('_')})
"
grep -rn "0.973\|0.988\|1.0097\|GATED_RTOL" docs/physics/n2-2d-td-cross-section.md docs/physics/n2-2d-cross-section.md | head
ls docs/physics/figures/ | grep -i n2
```

Every number that appears in the guide must be traceable to a note or to
source. Where a note gives a range or a caveat, carry the caveat — this is a
summary page, not a highlight reel.

- [ ] **Step 2: Write the guide**

Create `docs/molecules/n2.md`:

```markdown
# N₂ — the benchmark target

N₂ is the only molecule in this repository with **independent published data**
to check against: Karel Houfek's `CSVE.V00.J00` vibrational-excitation cross
sections. Everything else — NO, F₂, H₂⁺ — is validated against qscat's own
exact solver. That makes N₂ the anchor: the solver earns the right to be an
oracle elsewhere by matching Houfek here.

Atomic units throughout (Hartree, bohr).

:::{dropdown} The model — form and parameters
:icon: table

`qscat.model.N2` is a `DiatomicResonanceModel` with `charge = 0`: the shared
neutral-diatomic form, a Morse neutral curve plus a sigmoid-and-Gaussian
resonant interaction. Reproduce the exact parameter values with

```python
from qscat.model import N2
print({k: v for k, v in vars(N2).items() if not k.startswith("_")})
```

and paste them into a table here, then link the thesis reference note in
`reference/literature/` that publishes them.
:::

## What has been computed

:::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} Vibrational excitation, exactly
:link: ../physics/n2-2d-cross-section
:link-type: doc

The exact 2-D driven Lippmann–Schwinger solve, gated against Houfek at
`GATED_RTOL = 1e-3`. Figure: `../physics/figures/n2-2d-ti-cross-section.png`.
:::

:::{grid-item-card} The same answer, time-dependent
:link: ../physics/n2-2d-td-cross-section
:link-type: doc

Wavepacket propagation reaching the same cross section — and why order-1
Crank–Nicolson was not good enough to get there. Figure:
`../physics/figures/n2-2d-td-vs-ti-vs-houfek.png`.
:::

:::{grid-item-card} The local-complex-potential approximation
:link: ../physics/n2-cross-section
:link-type: doc

The 1-D reduction, and where it departs from the exact result.
:::

:::{grid-item-card} Resonance states
:link: ../physics/exact-2d-resonances
:link-type: doc

Poles of the full 2-D S-matrix, and what the Born–Oppenheimer error actually
measures. Figure: `../physics/figures/n2-exact-2d-resonance-levels.png`.
:::

:::

## The nonlocal resonance model

*(Summarise from `docs/physics/nonlocal-resonance-model.md` §8: the NRM
reproduces the exact 2-D solver on N₂ in vibrational excitation to better than
0.7 % in both the elastic and first inelastic channel. Embed
`../physics/figures/n2-ve-nrm-vs-exact.png` with a caption naming what is
plotted against what.)*

## Where to read more

`{doc}` links to every note above, plus {doc}`../physics/n2-resonance` and
{doc}`../physics/td-extractors`. Prose links, **not** a `toctree` — these notes
are owned by the technical sections, and a second parent fails the `-W` build.
```

- [ ] **Step 3: Add the molecules caption to `docs/index.md`**

After the `Theory — technical` toctree, add:

```markdown
```{toctree}
:maxdepth: 1
:caption: Theory — molecules

molecules/n2
```
```

- [ ] **Step 4: Build**

Run: `uv run sphinx-build -b html -W --keep-going docs docs/_build/html`
Expected: exit 0. A failure naming `n2-2d-cross-section` as having two parents means a `toctree` was used in the guide instead of prose links — replace it with `{doc}` links.

- [ ] **Step 5: Check the cards render in both themes**

Run: `uv run python -m http.server -d docs/_build/html 8000` and open
`http://localhost:8000/molecules/n2.html`. Toggle Furo's light/dark switch.
Expected: cards and the dropdown are legible in both; no unstyled raw text.
Stop the server.

- [ ] **Step 6: Commit**

```bash
git add docs/molecules/n2.md docs/index.md
git commit -m "docs(molecules): an N2 entry point for the site

Nothing on the site answered 'what do we know about N2' -- the method notes
are organised by method, which is right for them and useless for that
question. This is the per-molecule view: the model, what has been computed,
the headline numbers with their caveats, and links into the notes.

Prose links rather than a toctree: the notes are owned by the technical
sections, and a second parent fails the -W build."
```

---

### Task 6: The NO/F₂ and H₂⁺ guides

**Files:**
- Create: `docs/molecules/no-f2.md`, `docs/molecules/h2plus.md`
- Modify: `docs/index.md`

**Interfaces:**
- Consumes: the card/dropdown template from Task 5.
- Produces: the complete `Theory — molecules` section.

- [ ] **Step 1: Gather each molecule's facts from source**

```bash
uv run --no-sync python -c "
from qscat.model import NO, F2, H2P
for m in (NO, F2, H2P):
    m = m() if callable(m) else m
    print(type(m).__name__, getattr(m, 'charge', None))
"
grep -rn "0.263\|1.736\|1.8e9\|0.06\|1.9 %" docs/physics/diatomic-ve-cross-sections.md docs/physics/nonlocal-resonance-model.md | head -20
ls docs/physics/figures/ | grep -Ei "f2|no-|h2p"
```

- [ ] **Step 2: Write `docs/molecules/no-f2.md`**

Follow Task 5's structure exactly — `{dropdown}` for parameters, `{grid}` of
`{grid-item-card}`s, prose links at the end. The content that must appear,
each verified against the note in Step 1:

- These two have **no independent published data**. The exact 2-D solver is the
  oracle; agreement is self-consistency, not external validation.
- The LCP's error is **systematic and energy-dependent, not a fixed
  percentage**: on F₂ the ratio LCP/exact sweeps 0.263 → 1.736 across
  0.010–0.050 Ha, crossing unity near E≈0.032. On NO it fails outright away
  from threshold — the exact σ_DA decays 13 orders of magnitude while the LCP
  stays flat, reaching a ratio of 1.8e9.
- The NRM result and its unresolved half: on F₂ it lands within 0.06–1.9 % of
  the exact oracle where the LCP is off by 11–74 %; on NO it collapses by 5–8
  orders in dissociative attachment, and eight hypotheses and two mechanisms
  were killed by measurement without resolving it. Say this plainly — it is the
  most interesting thing on the page.
- Embed `../physics/figures/f2-da-nrm-vs-lcp-vs-exact.png` and the committed
  `f2-2d-da-lcp-vs-exact.png` / `no-2d-da-lcp-vs-exact.png` figures.

- [ ] **Step 3: Write `docs/molecules/h2plus.md`**

Same structure. Content, each verified in Step 1:

- H₂⁺ is the **first ionic target** (`charge = -1`), which changes the incident
  channel from a free wave to a Coulomb one and adds a Rydberg exit series.
- It is the **first non-laptop model** — ~1.15M unknowns at full size, needing
  Docker and the MUMPS backend.
- The resonance-verification result: angle stability is necessary and not
  sufficient. Four of 57 angle-stable poles scored overlaps of 6e-4..7e-3
  against a Born–Oppenheimer basis where genuine states score 0.87–0.99, and a
  further 18 of 57 were `box-limited` — Rydberg orbitals larger than the
  300-bohr box, which the c-product overlap is blind to by construction.
- Embed the committed H₂⁺ figures listed by Step 1.

- [ ] **Step 4: Extend the molecules toctree**

```markdown
```{toctree}
:maxdepth: 1
:caption: Theory — molecules

molecules/n2
molecules/no-f2
molecules/h2plus
```
```

- [ ] **Step 5: Build and test**

Run: `uv run sphinx-build -b html -W --keep-going docs docs/_build/html && uv run pytest tests/test_docs_portability.py -q`
Expected: build succeeds, tests pass.

- [ ] **Step 6: Commit**

```bash
git add docs/molecules/no-f2.md docs/molecules/h2plus.md docs/index.md
git commit -m "docs(molecules): NO/F2 and H2+ entry points

Both pages lead with the limitation rather than burying it: NO and F2 have
no independent published data, so the exact solver is the oracle and
agreement is self-consistency. The NO NRM collapse and the four H2+ poles
that turned out not to be resonances are on the page, not in a footnote."
```

---

### Task 7: The convention, the symbol table, and the femdvr-ecs pilot

**Files:**
- Modify: `.claude/skills/qscat-conventions/SKILL.md`
- Modify: `CONTRIBUTING.md`
- Modify: `docs/physics/femdvr-ecs.md`

**Interfaces:**
- Consumes: the portability test from Task 2.
- Produces: the written convention and the canonical symbol table that Tasks 8-9 follow.

`femdvr-ecs.md` is converted **first** because its symbols seed the table every
other note inherits (spec, Pilot notes 1).

- [ ] **Step 1: Add the mathematics convention to the conventions skill**

Append to `.claude/skills/qscat-conventions/SKILL.md`, after the "Units"
section:

```markdown
## Mathematics in Documentation

The notes under `docs/physics/` are read in two renderers: as files in a
clone, and as pages on the published site. So they use only what both
handle.

- **Backticks are for code identifiers only** — `ve_cross_section`,
  `SparseLU`, `backend="mumps"`. A backticked `sigma` that means σ is a
  defect; write `$\sigma$`.
- **Inline mathematics:** `$\sigma_{v\to v'}$`.
- **Display mathematics:** `$$...$$` on its own lines. Multi-line goes
  inside, as `$$ \begin{aligned} ... \end{aligned} $$` — a *bare*
  `\begin{align}` does not render on github.com.
- **Published equation numbers** use `\tag{15}`, with the page locator in
  the surrounding prose (see `mastering-references`).
- **Labelled equations** (`$$...$$ (eq-name)` plus the `{eq}` role) only
  where the note genuinely writes "Eq. (3)" — the label is literal text on
  github.com, so it is a deliberate cost.
- **No custom macros.** Spell out `\Psi^{(+)}`, `\sigma_\mathrm{DA}`.
- Atomic units are stated once per note, not per equation.
- `sphinx-design` directives (`{dropdown}`, `{grid}`, `{tab-set}`) are for
  site-first pages only — `docs/molecules/` and the section index pages.

`tests/test_docs_portability.py` enforces all of this.

### Canonical symbols

| Symbol | Meaning |
|---|---|
| $\theta$ | ECS rotation angle |
| $R_0$ | ECS pivot radius (always on an element boundary) |
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
| $T$, $S$ | T-matrix and S-matrix elements |
```

- [ ] **Step 2: Point CONTRIBUTING.md at it**

Add to `CONTRIBUTING.md`, after the "Checks" section:

```markdown
## Writing documentation

Physics notes live in `docs/physics/`, one per method. They are read both as
files in a clone and as pages on <https://vanamartin.github.io/qscat>, so
they are restricted to Markdown and LaTeX that renders in both — see the
"Mathematics in Documentation" section of the `qscat-conventions` skill for
the rules and the canonical symbol table.

```bash
uv run pytest tests/test_docs_portability.py    # enforces those rules
uv run sphinx-build -b html -W --keep-going docs docs/_build/html
```
```

- [ ] **Step 3: Convert `docs/physics/femdvr-ecs.md`**

Work through the file top to bottom. The transformation, by example — this
passage:

> the two outermost grid points are dropped to enforce a Dirichlet (psi=0)
> boundary condition [...] rotates the true continuum spectrum by
> `arg(E) ~ -2*theta`

becomes:

> the two outermost grid points are dropped to enforce a Dirichlet
> ($\psi = 0$) boundary condition [...] rotates the true continuum spectrum by
> $\arg(E) \sim -2\theta$

And the ECS map, currently prose-with-backticks, becomes display mathematics:

```markdown
$$
z(x) = \begin{cases}
x & x \le R_0 \\
R_0 + (x - R_0)\,e^{i\theta} & x > R_0
\end{cases}
$$
```

Rules for this pass:
- Every symbol that denotes a *quantity* moves from backticks to `$...$`.
- Every identifier that denotes *code* (`ecs_map`, `GridSpec`, `nq`,
  `gll_nodes_weights`, `angle_deg`) **stays** in backticks. `nq` is a
  constructor argument, not a symbol — leave it.
- Unlabelled fences holding equations become `$$` blocks; unlabelled fences
  holding directory trees or transcripts get a `text` tag.
- Do not reword the physics. This is a typesetting pass; if a sentence seems
  wrong, leave it and note it separately.

- [ ] **Step 4: Verify the conversion did not change meaning**

Run: `git diff docs/physics/femdvr-ecs.md`
Read every hunk. Confirm each is a notation change, not a content change.
Then: `uv run pytest tests/test_docs_portability.py -q`
Expected: pass.

- [ ] **Step 5: Build and read the rendered page**

Run: `uv run sphinx-build -b html -W --keep-going docs docs/_build/html`
Then serve and open `http://localhost:8000/physics/femdvr-ecs.html`.
Expected: every equation typeset; no stray `$`; no raw `\theta`.

- [ ] **Step 6: Confirm it still renders on GitHub**

Run: `gh markdown --help >/dev/null 2>&1 || echo "no gh markdown; use the web preview instead"`
Push the branch and view `docs/physics/femdvr-ecs.md` on github.com, or paste
the file into a draft comment preview. Expected: the same equations render.
This is the check that D1 is actually being honoured.

- [ ] **Step 7: Commit**

```bash
git add .claude/skills/qscat-conventions/SKILL.md CONTRIBUTING.md docs/physics/femdvr-ecs.md
git commit -m "docs(physics): typeset the FEM-DVR-ECS note, and write down the rule

First of three pilot conversions. femdvr-ecs goes first because its symbols
-- theta, R0, z(x), the weight Jacobian -- seed the canonical symbol table
the other notes inherit, which now lives in the qscat-conventions skill
alongside the rule that backticks are for code identifiers and dollar signs
are for quantities.

Typesetting only: no physics was reworded."
```

---

### Task 8: The nonlocal-resonance-model pilot

**Files:**
- Modify: `docs/physics/nonlocal-resonance-model.md`

**Interfaces:**
- Consumes: the symbol table from Task 7.
- Produces: the worked pattern for `\tag{}`-carrying literature equations, which the remaining nineteen notes follow during rollout.

The hardest case: 68 KB, 13 fenced blocks, 91 `Eq.` references, bra-kets,
integrals and a nonlocal kernel. If the convention survives here it survives
everywhere (spec, Pilot notes 2).

- [ ] **Step 1: Re-read the reference note before touching an equation**

Run: `sed -n '1,60p' reference/literature/houfek-2008-pra77-012710.md`

Every equation in this note carries a page locator from that paper. Under the
`mastering-references` skill the locators must survive the conversion — they
move from inside the fence into the prose.

- [ ] **Step 2: Convert the Feshbach-split block**

This block:

````
```
Q = |φ_d><φ_d|,   P = 1 − Q                                p. 012710-3, Eq. (15)
V_d(R)    = V_0(R) + <φ_d|H_el|φ_d>                        p. 012710-3, Eq. (20)
V_dk^+(R) = <φ_d|H_el|φ_k^+>                               p. 012710-3, Eq. (21)
lim_{R→∞} φ_d(r;R) = φ_b(r)                                p. 012710-7, Eq. (67)
```
````

becomes prose carrying the locator plus tagged display mathematics:

```markdown
Eq. (15), (20) and (21) are on p. 012710-3; the asymptotic requirement
Eq. (67) is on p. 012710-7.

$$Q = |\phi_d\rangle\langle\phi_d|, \qquad P = 1 - Q \tag{15}$$

$$V_d(R) = V_0(R) + \langle\phi_d|H_\mathrm{el}|\phi_d\rangle \tag{20}$$

$$V_{dk}^{+}(R) = \langle\phi_d|H_\mathrm{el}|\phi_k^{+}\rangle \tag{21}$$

$$\lim_{R\to\infty} \phi_d(r;R) = \phi_b(r) \tag{67}$$
```

- [ ] **Step 3: Convert the nuclear-equation block**

```markdown
Eq. (52)-(54), p. 012710-5.

$$
\left[E - T_R - V_d(R)\right]\Psi_d^{+}(R)
  - \int F(E,R,R')\,\Psi_d^{+}(R')\,\mathrm{d}R'
  = V_{dk_i}^{+}(R)\,\chi_{v_i}(R)
\tag{52}
$$

$$
F(E,R,R') = \int V_{dk}^{+}(R)\,G_0^{+}(E,R,R')\,V_{dk}^{+}(R')^{*}\,k\,\mathrm{d}k
\tag{53}
$$

$$
\sigma_\mathrm{DA}(E) = \frac{2\pi^2}{k_i^2}\,\frac{K_\mathrm{DA}}{\mu}\,
  \lim_{R\to\infty}\left|\Psi_d^{+}(R)\right|^2
\tag{54}
$$
```

- [ ] **Step 4: Convert the remaining equation fences the same way**

Work through the rest of the file. Keep every `\tag{n}` matching the paper's
equation number exactly — a wrong tag is a citation error, not a typo. Fences
holding shell transcripts get a `console` tag; fences holding tabular output
get `text`.

- [ ] **Step 5: Verify every tag against the reference note**

Run:

```bash
uv run --no-sync python - <<'PY'
import pathlib, re
note = pathlib.Path("docs/physics/nonlocal-resonance-model.md").read_text()
tags = sorted({int(m) for m in re.findall(r"\\tag\{(\d+)\}", note)})
print("tagged equations:", tags)
ref = pathlib.Path("reference/literature/houfek-2008-pra77-012710.md").read_text()
anchored = sorted({int(m) for m in re.findall(r"Eq\.\s*\((\d+)\)", ref)})
print("anchored in reference note:", anchored)
print("tagged but NOT anchored:", [t for t in tags if t not in anchored])
PY
```

Expected: `tagged but NOT anchored: []`. Any equation number tagged here but
absent from the reference note must be anchored there first, under the
`mastering-references` skill — that skill's discipline has already caught three
real discrepancies in this repository.

- [ ] **Step 6: Test, build, review the diff**

Run: `uv run pytest tests/test_docs_portability.py -q && uv run sphinx-build -b html -W --keep-going docs docs/_build/html`
Expected: both clean.
Then `git diff docs/physics/nonlocal-resonance-model.md` and read every hunk:
notation only, no physics reworded, no locator dropped.

- [ ] **Step 7: Commit**

```bash
git add docs/physics/nonlocal-resonance-model.md
git commit -m "docs(nrm): typeset the nonlocal resonance model note

The hardest note in the tree and the real test of the convention: 91
equation references with page locators, bra-kets, and a nonlocal kernel.
Published equation numbers become \\tag{n} -- core MathJax, so they render
on github.com too -- and the page locators move into the prose, which is
where mastering-references wants them.

Every tag was checked against the anchored equation numbers in
reference/literature/houfek-2008-pra77-012710.md."
```

---

### Task 9: The n2-2d-cross-section pilot

**Files:**
- Modify: `docs/physics/n2-2d-cross-section.md`

**Interfaces:**
- Consumes: the symbol table from Task 7 and the `\tag{}` pattern from Task 8.
- Produces: the third and final pilot; the mixed prose/mathematics/table/figure case.

- [ ] **Step 1: Convert the driven Lippmann–Schwinger block**

This block, at roughly line 32:

````
```
Psi_i    = F_{E,l}(r) chi_v(R)                       [masked to the unscaled region]
Psi_sc   = (E_tot*I - H_2D)^-1 V_int Psi_i           [one sparse LU per energy]
Psi^(+)  = Psi_i + Psi_sc
T_{v->v'} = <chi_v' F_{E',l} | V_int | Psi^(+)>       [c-product, masked]
sigma_{v->v'} = 4 pi^3 |T|^2 / k^2                    [bohr^2]
```
````

becomes:

```markdown
$$
\begin{aligned}
\Psi_i        &= F_{E,l}(r)\,\chi_v(R) &&\text{masked to the unscaled region}\\
\Psi_\mathrm{sc} &= \left(E_\mathrm{tot}\mathbb{1} - H_\mathrm{2D}\right)^{-1} V_\mathrm{int}\,\Psi_i
                 &&\text{one sparse LU per energy}\\
\Psi^{(+)}    &= \Psi_i + \Psi_\mathrm{sc} \\
T_{v\to v'}   &= \langle \chi_{v'} F_{E',l} \,|\, V_\mathrm{int} \,|\, \Psi^{(+)} \rangle
                 &&\text{c-product, masked}\\
\sigma_{v\to v'} &= \frac{4\pi^3 |T|^2}{k^2} &&\text{bohr}^2
\end{aligned}
$$
```

Note `\begin{aligned}` **inside** `$$` — a bare `\begin{align}` fails the
portability test and does not render on github.com.

- [ ] **Step 2: Convert the remaining mathematics in the file**

Inline quantities move from backticks to `$...$`; code identifiers
(`ve_cross_section_2d`, `WORKING_GRID`, `SparseLU.refactor`, `BASELINE`,
`GATED_RTOL`) stay in backticks. Leave the results tables as Markdown tables —
they are data, not mathematics — but convert their *headers* (`σ_exact
(bohr²)`) only if they currently use ASCII.

- [ ] **Step 3: Leave the numbers alone**

This note carries measured values (`2.368e-06`, `GATED_RTOL = 1e-3`, the
six-anchor table). Do not reformat, round, or "tidy" any of them. Confirm with:

```bash
git diff --word-diff=porcelain docs/physics/n2-2d-cross-section.md | grep -E '^\-' | grep -E '[0-9]e-?[0-9]' | head
```

Expected: only lines where a number moved inside `$...$` unchanged, never a
line where the digits differ.

- [ ] **Step 4: Test, build, review**

Run: `uv run pytest tests/test_docs_portability.py -q && uv run sphinx-build -b html -W --keep-going docs docs/_build/html`
Expected: both clean. Read the rendered page and confirm the six-anchor table
and the committed figure still render.

- [ ] **Step 5: Commit**

```bash
git add docs/physics/n2-2d-cross-section.md
git commit -m "docs(n2): typeset the exact 2-D cross-section note

Third and last pilot -- the mixed case, with a derivation, a results table
and a figure on one page. The driven Lippmann-Schwinger chain becomes one
aligned block with its side conditions as text annotations rather than
bracketed comments.

Measured values are untouched; only their notation moved."
```

---

### Task 10: Fence audit for the unconverted notes

**Files:**
- Modify: the nineteen unconverted notes under `docs/physics/`, where a fence lacks a language

**Interfaces:**
- Consumes: nothing.
- Produces: correctly highlighted, copy-buttoned code blocks across the whole tree.

**Retagging only.** An equation sitting in an unconverted note gets a `text`
tag now and becomes `$$` mathematics when that note is converted during
rollout. This task must not convert any mathematics — that is what keeps it
mechanical and safe to review.

- [ ] **Step 1: List every unlabelled fence with its first line**

```bash
uv run --no-sync python - <<'PY'
import pathlib
pilots = {"femdvr-ecs.md", "nonlocal-resonance-model.md", "n2-2d-cross-section.md"}
for p in sorted(pathlib.Path("docs/physics").glob("*.md")):
    if p.name in pilots:
        continue
    lines = p.read_text().splitlines()
    inside = False
    for i, line in enumerate(lines):
        if line.startswith("```"):
            if not inside and line.strip() == "```":
                nxt = lines[i + 1] if i + 1 < len(lines) else ""
                print(f"{p.name}:{i+1}: {nxt[:70]}")
            inside = not inside
PY
```

- [ ] **Step 2: Tag each fence by what it holds**

For each line the script printed, open the file at that line and add exactly
one of:

- `console` — a shell transcript, or a command with its output.
- `bash` — a command with no output shown.
- `text` — everything else: equations awaiting conversion, tabular output,
  directory trees, schematic layouts.

Do **not** tag anything `python` in this task. A `python` fence becomes visible
to the doctest builder, and a non-runnable snippet tagged `python` would break
`sphinx-build -b doctest`.

- [ ] **Step 3: Confirm no unlabelled fences remain**

Re-run Step 1's script.
Expected: no output.

- [ ] **Step 4: Confirm the doctest builder is unaffected**

Run: `uv run sphinx-build -b doctest docs docs/_build/doctest`
Expected: exit 0, and the same number of tests as before this task — no fence
was accidentally tagged `python`.

- [ ] **Step 5: Build and test**

Run: `uv run sphinx-build -b html -W --keep-going docs docs/_build/html && uv run pytest tests/test_docs_portability.py -q`
Expected: both clean.

- [ ] **Step 6: Commit**

```bash
git add docs/physics/
git commit -m "docs(physics): tag the untagged fenced blocks

32 fences carried no language, so Pygments rendered them flat grey and the
copy button had nothing to attach to. Retagging only -- equations still
sitting in fences get 'text' now and become display mathematics when their
note is converted. Nothing is tagged python: that would hand a
non-runnable snippet to the doctest builder."
```

---

### Task 11: Rollout checklist and final verification

**Files:**
- Modify: `docs/superpowers/plans/2026-08-19-docs-latex-and-theory-ia.md` (this file — the checklist below)
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes: everything above.
- Produces: a green build, a green suite, and the tracked rollout list.

- [ ] **Step 1: Run every gate**

```bash
uv run sphinx-build -b html -W --keep-going docs docs/_build/html
uv run sphinx-build -b doctest docs docs/_build/doctest
uv run pytest tests/ -q
uv run ruff check .
```

Expected: all four clean. Report the actual output; do not summarise a failure
as a pass.

- [ ] **Step 2: Read the site in both themes**

Serve `docs/_build/html` and walk: the landing page, one section index, one
converted note, one unconverted note, and all three molecule guides. Toggle
light/dark on each. Confirm the sidebar sections collapse and expand.

- [ ] **Step 3: Confirm the dual-render promise on the real thing**

Push the branch and open all three converted notes on github.com. Every
equation must render there too. If one does not, that construct violates D1 and
must be rewritten — this is the check the whole portability rule exists for.

- [ ] **Step 4: Record the rollout checklist**

Append to this plan file:

```markdown
## Rollout: notes awaiting conversion

Converted when next touched. A note is done when the portability test passes
and no ASCII equation remains outside a code fence.

- [ ] nd-tensor-hamiltonian
- [ ] discretisation-tuning
- [ ] mumps-sparse-backend
- [ ] ti-energy-sweep-reuse
- [ ] shift-invert-eigensolver
- [ ] qscat-core-scattering
- [ ] n2-resonance
- [ ] n2-cross-section
- [ ] n2-td-cross-section
- [ ] n2-2d-td-cross-section
- [ ] td-extractors
- [ ] td-da
- [ ] diatomic-ve-cross-sections
- [ ] h2plus-dr
- [ ] lcp-resonance-levels
- [ ] exact-2d-resonances
- [ ] h2plus-resonance-states
- [ ] angular-coupled-channels
- [ ] optimization-targets
- [ ] validation-harnesses
```

- [ ] **Step 5: Add a CHANGELOG entry**

Follow the existing format in `CHANGELOG.md`. Record: rendered mathematics with
the portability constraint, the Theory regrouping, the molecule guides, and the
new validation-harnesses note.

- [ ] **Step 6: Commit**

```bash
git add docs/superpowers/plans/2026-08-19-docs-latex-and-theory-ia.md CHANGELOG.md
git commit -m "docs: record the conversion rollout and changelog

Nineteen notes still carry ASCII mathematics. They convert as they are next
touched rather than in one sweep: nothing in this repository executes the
mathematics in a note, so a transcription error is invisible to every test
here, and a 380 KB single pass would trade a visible defect for an
invisible one."
```

---

## Rollout: notes awaiting conversion

All twenty are converted. The staging this section described -- convert as
each note is next touched -- was overtaken: the whole rollout was carried out
in one branch after the pilot landed, group by group along the sidebar's own
sections, with a numeric-token multiset diff run per file so that no measured
value could change under cover of a notation change.

- [x] nd-tensor-hamiltonian
- [x] discretisation-tuning
- [x] mumps-sparse-backend
- [x] ti-energy-sweep-reuse
- [x] shift-invert-eigensolver
- [x] qscat-core-scattering
- [x] n2-resonance
- [x] n2-cross-section
- [x] n2-td-cross-section
- [x] n2-2d-td-cross-section
- [x] td-extractors
- [x] td-da
- [x] diatomic-ve-cross-sections
- [x] h2plus-dr
- [x] lcp-resonance-levels
- [x] exact-2d-resonances
- [x] h2plus-resonance-states
- [x] angular-coupled-channels
- [x] optimization-targets (no mathematics to convert; its fences are
      profiler output and stay `text`)
- [x] validation-harnesses

Two conventions were settled during the rollout and belong with the others in
the `qscat-conventions` skill:

- **Headings stay plain unicode.** A heading also renders in the sidebar and
  in `toctree` entries, where MathJax does not run, so `$^2\Pi_g$` in a title
  would show as literal source there.
- **Level labels in table cells stay backticked.** `v = 0`, `v = 1` read as
  labels rather than equations, and a table cell is where an unescaped pipe
  inside maths silently drops content -- the defect the NRM note had to fix.

One gap in the earlier fence audit was found and closed here: it matched
fences at line start, so three fences indented inside list items
(`angular-coupled-channels`, `n2-2d-td-cross-section`, `n2-cross-section`)
were never counted. All three held equations and are now display maths.
