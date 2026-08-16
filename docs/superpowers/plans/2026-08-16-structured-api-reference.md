# Structured API Reference Implementation Plan (Plan A)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the single 119-name `docs/api.md` with a navigable per-subpackage API reference, add a test that fails when a public name is undocumented, and publish the 18 `docs/physics/` theory notes that are currently excluded from the site.

**Architecture:** `docs/api/` gets one Markdown page per `qscat` subpackage. Every page lists each public name explicitly — via an `autosummary` index for the small modules, or sectioned `autoclass`/`autofunction` directives for `qscat.core` — so a repo-root pytest can assert coverage by reading the Markdown. Publishing the theory notes is a `conf.py` and toctree change only: a trial build proved it needs no link-fixing.

**Tech Stack:** Sphinx 8 + MyST-Parser + autodoc/autosummary/napoleon, furo theme, pytest.

**Source spec:** `docs/superpowers/specs/2026-08-16-documentation-showcase-design.md` (Part A). Plan B, the showcase gallery, is written separately after this merges.

## Global Constraints

- The docs build must pass the CI command exactly: `uv run sphinx-build -b html -W --keep-going docs docs/_build/html`. `-W` turns every warning into an error, so an unresolvable autodoc target or an orphaned page fails the build.
- The public API is exactly each submodule's `__all__` (ADR 0004). Nothing outside `__all__` gets documented; nothing inside it may be omitted.
- Never edit anything under `docs/superpowers/plans/` or `docs/superpowers/specs/` other than this plan file — those are frozen historical records.
- `docs/adr/` and `docs/superpowers/` stay excluded from the rendered site.
- Atomic units throughout; do not restate physics in API prose — link to `docs/physics/` instead.
- Preserve the repo's commit trailers on every commit:
  `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_01NxwNdBLUXBampDLrusGrAS`.
- Run tests with `uv run --no-sync pytest` when no dependency changed, to skip a redundant resolve.

## Verified facts this plan is built on

These were measured on 2026-08-16, not assumed. Do not re-derive them; do verify anything you change.

1. `qscat.core.__all__` has **37** names; `dvr` 15, `tuning` 22, `viz` 9, `model` 7, `special` 7, `linalg` 6, `exceptions` 5, `units` 4, `evolution` 4, `ecs` 3. Total **119**.
2. `qscat.viz` and `qscat.units` appear **nowhere** in the current `docs/api.md` — 13 public names with no rendered docs.
3. A trial build with `docs/physics/` un-excluded and all 18 notes in a toctree **succeeded under `-W --keep-going` with zero warnings**. All 15 image references (`figures/*.png`) are relative and were copied correctly. There are **no** repo-path hyperlinks in the notes to fix.
4. There is **no** repo-root `tests/` directory yet. `pyproject.toml` sets `pythonpath = ["."]` and `norecursedirs` does not exclude `tests`, so a root `tests/` package is collected automatically.
5. In `qscat.core`, exactly seven names are classes — `ScatteringProblem`, `VibrationalBasis`, `ResonanceLevels`, `Extractor`, `TannorWeeks`, `Dirac`, `Flux` — and the other 30 are functions. Using the wrong directive is a `-W` build failure.

---

## File Structure

**Create:**
- `tests/__init__.py` — makes the root test dir a package, matching repo convention.
- `tests/test_api_docs_coverage.py` — the `__all__`-coverage gate.
- `docs/api/index.md` — layer map, the ADR 0004 contract, toctree over the pages.
- `docs/api/core.md` — 37 names in eight sections.
- `docs/api/model.md`, `dvr.md`, `ecs.md`, `linalg.md`, `evolution.md`, `special.md`, `tuning.md`, `viz.md` — one page each.
- `docs/api/base.md` — `qscat.units` + `qscat.exceptions`.

**Modify:**
- `docs/index.md` — toctree `api` → `api/index`; add the Theory section.
- `docs/conf.py` — drop `"physics"` from `exclude_patterns`; update the module docstring.
- `docs/physics/README.md` — replace the three-line stub with a real section index.

**Delete:**
- `docs/api.md` — superseded by `docs/api/`.

---

### Task 1: The `__all__` coverage gate

This task deliberately ends **RED**. The test is the specification for Task 2; it cannot pass until `docs/api/` exists. Do not create any docs file in this task, and do not weaken the test to make it pass.

**Files:**
- Create: `tests/__init__.py`
- Create: `tests/test_api_docs_coverage.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: the contract Task 2 must satisfy — every name in every submodule's `__all__` must appear in some `docs/api/**/*.md` file, either as a bare name on its own line (an `autosummary` entry) or as the trailing component of an `.. auto<something>:: qscat.<mod>.<Name>` directive line.

- [ ] **Step 1: Create the test package marker**

Create `tests/__init__.py` containing exactly:

```python
"""Repository-wide tests (not part of any distributable package)."""
```

- [ ] **Step 2: Write the failing test**

Create `tests/test_api_docs_coverage.py`:

```python
"""Every public name in qscat's ``__all__`` must appear in the rendered API reference.

Guards the gap that motivated the split: the single-page ``docs/api.md``
documented eight submodules and silently omitted ``qscat.viz`` (9 public
names) and ``qscat.units`` (4), so 13 public names had no rendered
documentation at all. A layout can be re-prettified at any time; this test is
what stops the omission from recurring.

Lives at the repository root rather than under ``libs/qscat/tests/`` because it
reads ``docs/``, which ships in the repository but not in the qscat sdist.
"""

from __future__ import annotations

import importlib
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
API_DOCS = REPO_ROOT / "docs" / "api"

# qscat._SUBMODULES plus `exceptions`, which is a module rather than a lazily
# exposed submodule but is equally part of the public surface (ADR 0004).
MODULES = (
    "core",
    "model",
    "dvr",
    "ecs",
    "linalg",
    "evolution",
    "special",
    "tuning",
    "viz",
    "units",
    "exceptions",
)

# Two shapes count as "documented", and both are whole-line matches:
#   * an autosummary entry -- a bare indented name on its own line
#   * an autodoc directive -- `.. autofunction:: qscat.core.ve_cross_section`
# Matching whole lines only is what keeps an incidental prose mention of a name
# from counting as documentation for it.
_DOCUMENTED = re.compile(r"^\s*(?:\.\.\s+auto\w+::\s+[\w.]*?\.)?(\w+)\s*$")


def _documented_names() -> set[str]:
    names: set[str] = set()
    for path in sorted(API_DOCS.rglob("*.md")):
        for line in path.read_text(encoding="utf-8").splitlines():
            match = _DOCUMENTED.match(line)
            if match:
                names.add(match.group(1))
    return names


def test_api_docs_directory_is_populated() -> None:
    # Without this, a mistyped path would make `_documented_names()` return an
    # empty set and every check below would fail loudly -- but a mistyped
    # *glob* would make it return nothing while the directory exists, so assert
    # both the directory and at least one page.
    assert API_DOCS.is_dir(), f"no API reference directory at {API_DOCS}"
    assert list(API_DOCS.rglob("*.md")), f"no .md pages under {API_DOCS}"


@pytest.mark.parametrize("module", MODULES)
def test_every_public_name_is_documented(module: str) -> None:
    mod = importlib.import_module(f"qscat.{module}")
    public = set(mod.__all__)
    assert public, f"qscat.{module} exports no public names -- is __all__ missing?"
    missing = sorted(public - _documented_names())
    assert not missing, (
        f"qscat.{module} exports {len(missing)} public name(s) with no entry "
        f"under docs/api/: {missing}"
    )


def test_documented_names_do_not_include_private_names() -> None:
    # The reference documents the public surface only (ADR 0004). A leading
    # underscore in an autodoc directive means a private name leaked into it.
    leaked = sorted(n for n in _documented_names() if n.startswith("_"))
    assert not leaked, f"private names documented under docs/api/: {leaked}"
```

- [ ] **Step 3: Run the test and confirm it fails for the right reason**

Run: `uv run --no-sync pytest tests/test_api_docs_coverage.py -v`

Expected: 12 of the 13 tests FAIL — `test_api_docs_directory_is_populated`
with `no API reference directory at .../docs/api`, plus all 11 parametrized
cases. `test_documented_names_do_not_include_private_names` PASSES vacuously
(nothing is documented, so nothing private is). That is the intended RED state.

Confirm the failure is a missing directory, **not** an import error or a regex bug. If `qscat` fails to import, stop and fix the environment (`uv sync --all-packages`) before continuing.

- [ ] **Step 4: Commit the RED test**

```bash
git add tests/__init__.py tests/test_api_docs_coverage.py
git commit -m "test: require every public qscat name to appear in the API reference

Currently RED: docs/api/ does not exist yet. The test is the specification
for the per-subpackage split -- it fails today for all 119 public names, and
would have failed for the 13 names in qscat.viz and qscat.units that the
single-page docs/api.md silently omitted.

Lives at the repository root because it reads docs/, which ships in the repo
but not in the qscat sdist.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01NxwNdBLUXBampDLrusGrAS"
```

---

### Task 2: Build `docs/api/` and turn the gate GREEN

**Files:**
- Create: `docs/api/index.md`, `core.md`, `model.md`, `dvr.md`, `ecs.md`, `linalg.md`, `evolution.md`, `special.md`, `tuning.md`, `viz.md`, `base.md`
- Modify: `docs/index.md` (toctree entry `api` → `api/index`; the `{doc}` quick link)
- Delete: `docs/api.md`
- Test: `tests/test_api_docs_coverage.py` (from Task 1, unchanged)

**Interfaces:**
- Consumes: the coverage contract from Task 1.
- Produces: `docs/api/index.md` as the toctree root for the reference; Task 3 adds a sibling Theory section to `docs/index.md` and must not disturb this entry.

**Two page patterns. Use pattern B only for `qscat.core`.**

*Pattern A (every module except `core` and `base`)* — an `autosummary` index listing every public name, then `automodule` for the detail. The explicit index is not decoration: `automodule` alone puts no individual name in the Markdown, so the coverage test would fail. The index is what makes the page both scannable and testable.

*Pattern B (`core` only)* — 37 names grouped into eight H2 sections with explicit `autoclass`/`autofunction` directives. `core` is too large for one flat list, and the sections are its index.

- [ ] **Step 1: Write the reference landing page**

Create `docs/api/index.md`:

````markdown
# API reference

Generated from the source by autodoc, one page per subpackage.

## What is public

The public API is exactly the names exported in each submodule's `__all__`.
Anything with a leading underscore, and any name not in an `__all__`, is
private and may change without notice. While the version is `0.y.z`, minor
releases may contain breaking changes to the public API; they are called out in
`CHANGELOG.md`. See ADR 0004 (`docs/adr/0004-public-api-stability-policy.md`)
for the full policy, including the *provisional* marker and the post-1.0
deprecation rule.

Every recoverable error is a subclass of `qscat.exceptions.QscatError`, so
`except QscatError` is a valid catch-all.

## The layers

| layer | modules | what lives here |
|---|---|---|
| Engine | {doc}`core`, {doc}`model` | The model-independent scattering engine, and the per-molecule potentials it consumes. |
| Numerics | {doc}`dvr`, {doc}`ecs`, {doc}`linalg`, {doc}`evolution`, {doc}`special` | FEM-DVR grids, the complex-scaling map, sparse linear algebra, time propagators, radial/Coulomb special functions. |
| Tooling | {doc}`tuning`, {doc}`viz` | Deriving a grid from a potential; rendering wavefunctions. |
| Base | {doc}`base` | Unit conversions and the exception hierarchy. |

`qscat.core` never imports `qscat.model` at runtime — it depends only on the
`ResonanceModel` protocol. Adding a molecule is a registry entry, never solver
code.

```{toctree}
:maxdepth: 1

core
model
dvr
ecs
linalg
evolution
special
tuning
viz
base
```
````

- [ ] **Step 2: Write the `qscat.core` page (Pattern B)**

Create `docs/api/core.md`. Use `.. autoclass::` for exactly these seven —
`ScatteringProblem`, `VibrationalBasis`, `ResonanceLevels`, `Extractor`,
`TannorWeeks`, `Dirac`, `Flux` — and `.. autofunction::` for the other 30.
Classes take `:members:`.

````markdown
# qscat.core

The model-independent electron–diatomic scattering engine. `ScatteringProblem`
is the recommended entry point: it bundles the grid, model, and vibrational
basis once and exposes every observable as a method. The functional solvers
below are the low-level layer those methods call.

This module never imports `qscat.model` at runtime — it depends only on the
`ResonanceModel` protocol, so a new molecule needs no change here.

## The problem object

```{eval-rst}
.. autoclass:: qscat.core.ScatteringProblem
   :members:
```

## Cross sections (time-independent)

```{eval-rst}
.. autofunction:: qscat.core.ve_cross_section
.. autofunction:: qscat.core.da_cross_section
.. autofunction:: qscat.core.dr_cross_section
```

## Cross sections (time-dependent)

```{eval-rst}
.. autofunction:: qscat.core.td_ve_cross_section
.. autofunction:: qscat.core.td_ve_cross_sections_all
.. autofunction:: qscat.core.td_da_cross_section
.. autofunction:: qscat.core.td_da_cross_sections_all
```

## Grids

```{eval-rst}
.. autofunction:: qscat.core.electronic_grid
.. autofunction:: qscat.core.nuclear_grid
.. autofunction:: qscat.core.fem_grid_exp_tail
.. autofunction:: qscat.core.segmented_grid
```

## Channels

```{eval-rst}
.. autofunction:: qscat.core.channel_vector
.. autofunction:: qscat.core.anion_electronic_states
.. autofunction:: qscat.core.v_dr_diag
.. autofunction:: qscat.core.outgoing_channel
.. autofunction:: qscat.core.outgoing_channel_nuclear
```

## Vibrational structure

```{eval-rst}
.. autofunction:: qscat.core.vibrational_states
.. autoclass:: qscat.core.VibrationalBasis
   :members:
```

## Wavepacket and correlation

```{eval-rst}
.. autofunction:: qscat.core.gaussian_coeffs
.. autofunction:: qscat.core.initial_state
.. autofunction:: qscat.core.eta_incident
.. autofunction:: qscat.core.eta_outgoing
.. autofunction:: qscat.core.hankel_point_value
.. autofunction:: qscat.core.outgoing_surface_wave
.. autofunction:: qscat.core.propagate
.. autofunction:: qscat.core.sigma_from_correlations
```

## Time-dependent energy extractors

All three share one propagate-once protocol, so a single propagation can drive
every extractor. See `docs/physics/td-extractors.md`.

```{eval-rst}
.. autoclass:: qscat.core.Extractor
   :members:
.. autoclass:: qscat.core.TannorWeeks
   :members:
.. autoclass:: qscat.core.Dirac
   :members:
.. autoclass:: qscat.core.Flux
   :members:
```

## The LCP approximation

The local-complex-potential reduction and the Born–Oppenheimer resonance
levels built on it. These are the *approximation* under test against the exact
solvers above; see `docs/physics/diatomic-ve-cross-sections.md` and
`docs/physics/lcp-resonance-levels.md`.

```{eval-rst}
.. autofunction:: qscat.core.local_complex_potential
.. autofunction:: qscat.core.lcp_da_cross_section
.. autofunction:: qscat.core.resonance_levels
.. autofunction:: qscat.core.lcp_resonance_levels
.. autoclass:: qscat.core.ResonanceLevels
   :members:
```

## Plotting

```{eval-rst}
.. autofunction:: qscat.core.plot_cross_sections
```
````

- [ ] **Step 3: Write the eight Pattern-A pages**

Each file is the same shape. Substitute the module name, the intro paragraph,
and the name list from the table below:

````markdown
# qscat.<MODULE>

<INTRO>

```{eval-rst}
.. currentmodule:: qscat.<MODULE>

.. autosummary::
   :nosignatures:

   <one public name per line, indented three spaces, in __all__ order>
```

```{eval-rst}
.. automodule:: qscat.<MODULE>
   :members:
   :imported-members:
```
````

| file | module | intro | names, in `__all__` order |
|---|---|---|---|
| `model.md` | `qscat.model` | Everything tied to a specific molecule: the `ResonanceModel` protocol that `qscat.core` depends on, the shared neutral and ionic potential forms, and the molecule registry. Adding a molecule is a registry entry plus validation, never solver code. | `ResonanceModel`, `DiatomicResonanceModel`, `IonicResonanceModel`, `N2`, `NO`, `F2`, `H2P` |
| `dvr.md` | `qscat.dvr` | The FEM-DVR grid with an exterior-complex-scaled tail: grid construction, kinetic-energy assembly (dense and sparse), and the N-dimensional tensor layer. See `docs/physics/femdvr-ecs.md` and `docs/physics/nd-tensor-hamiltonian.md`. | `ElementSpec`, `GridSpec`, `FemDvrEcsGrid`, `kinetic`, `kinetic_sparse`, `dvr_first_derivative_at_node`, `dvr_interpolation_matrix`, `hamiltonian`, `eigen`, `gll_nodes_weights`, `diff_matrix`, `TensorGrid`, `kinetic_nd`, `potential_nd`, `hamiltonian_nd` |
| `ecs.md` | `qscat.ecs` | The exterior-complex-scaling coordinate map — the single source of the `z(x)` transform used by every complex tail — and the resonance-pole matchers built on angle stability. | `ecs_map`, `find_resonance_pole`, `match_angle_stable` |
| `linalg.md` | `qscat.linalg` | Dimension-general sparse linear algebra: Kronecker sums over arbitrary D, a cached sparse LU with a SuperLU and a complex-symmetric MUMPS backend, and the bilinear (non-conjugated) ECS inner product. See `docs/physics/mumps-sparse-backend.md`. | `kron_sum`, `c_product`, `SparseLU`, `default_backend`, `set_default_backend`, `get_default_backend` |
| `evolution.md` | `qscat.evolution` | Time propagators for `d/dt psi = -i H psi` with complex, possibly non-Hermitian `H`. The order-N diagonal-Padé stepper generalizes Crank–Nicolson (order 1); order 3 is what makes the time-dependent cross sections converge. | `make_cn_stepper`, `make_pade_stepper`, `make_sparse_cn_stepper`, `pade_roots` |
| `special.md` | `qscat.special` | Riccati–Bessel and Riccati–Hankel functions, their reduced-mass generalizations, and the Coulomb functions used for ionic targets. | `riccati_bessel_en`, `riccati_hankel_en`, `riccati_bessel_en_mass`, `riccati_hankel_en_mass`, `coulomb_f_en`, `coulomb_g_en`, `coulomb_h1_en` |
| `tuning.md` | `qscat.tuning` | The automatic discretisation tuner: derive a minimal FEM-DVR-ECS grid at a target precision from the potential and energy range, instead of hand-picking element lengths. See `docs/physics/discretisation-tuning.md`, including the documented limits of the 1-D convergence probes. | `IncidentSpec`, `PotentialProfile`, `ProbeResult`, `analyze_potential`, `equidistribution_elements`, `grid_cost`, `interaction_region`, `max_stable_angle`, `optimal_real_mesh`, `order_for_wavenumber`, `probe_channel_representation`, `probe_electronic`, `probe_nuclear`, `propose_grid`, `refine`, `refine_elements_in_window`, `refine_to_2d_convergence`, `required_extent`, `resonance_curve`, `tensor_cost`, `tune_ecs_tail`, `tw_analysis` |
| `viz.md` | `qscat.viz` | Rendering wavefunctions: a cached sparse projector onto an equidistant grid, complex-plane domain colouring, and static or animated 2-D field plots. Needs the `plot` extra (`pip install "qscat[plot]"`). | `EquidistantProjector`, `WavefunctionArtist`, `complex_to_hsv`, `complex_to_rgb`, `hsv_to_rgb`, `energy_contour_levels`, `plot_wavefunction_2d`, `animate_wavefunction`, `animate_artists` |

- [ ] **Step 4: Write the base page**

Create `docs/api/base.md`. Use explicit directives, **not** `automodule` — the
two `units` constants are plain module data and the exception classes need
`autoexception`, so `automodule` risks an unresolved-target warning under `-W`.

````markdown
# qscat.units and qscat.exceptions

## qscat.units

Atomic units are used throughout qscat: energies in Hartree, lengths in Bohr.
These conversions exist so unit handling is never scattered through method
code — there is exactly one place that knows the Hartree/eV factor.

```{eval-rst}
.. currentmodule:: qscat.units

.. autosummary::
   :nosignatures:

   HARTREE_TO_EV
   EV_TO_HARTREE
   hartree_to_ev
   ev_to_hartree

.. autodata:: qscat.units.HARTREE_TO_EV
.. autodata:: qscat.units.EV_TO_HARTREE
.. autofunction:: qscat.units.hartree_to_ev
.. autofunction:: qscat.units.ev_to_hartree
```

## qscat.exceptions

Every recoverable error qscat raises is a subclass of `QscatError`. Each also
subclasses the built-in it replaces, so catching the built-in remains valid.
Generic argument validation may still raise a plain `ValueError`/`TypeError`.

```{eval-rst}
.. currentmodule:: qscat.exceptions

.. autosummary::
   :nosignatures:

   QscatError
   GridError
   ModelError
   BackendError
   ConvergenceError

.. autoexception:: qscat.exceptions.QscatError
   :members:
.. autoexception:: qscat.exceptions.GridError
   :members:
.. autoexception:: qscat.exceptions.ModelError
   :members:
.. autoexception:: qscat.exceptions.BackendError
   :members:
.. autoexception:: qscat.exceptions.ConvergenceError
   :members:
```
````

- [ ] **Step 5: Run the coverage test — expect GREEN**

Run: `uv run --no-sync pytest tests/test_api_docs_coverage.py -v`

Expected: all 13 tests PASS.

If a module still reports missing names, add them to that page's `autosummary`
block — do not relax the test. If `test_documented_names_do_not_include_private_names`
fails, you documented something outside `__all__`; remove it.

- [ ] **Step 6: Delete the old page and repoint the site toctree**

Delete `docs/api.md`.

In `docs/index.md`, change the toctree entry `api` to `api/index`, and change
the quick link `- **API reference:** {doc}`api`` to `{doc}`api/index``.

- [ ] **Step 7: Build the docs with the CI command**

Run: `uv run sphinx-build -b html -W --keep-going docs docs/_build/html`

Expected: `build succeeded.` with no warnings. `-W` makes any unresolved
autodoc target fatal, so this is the real check that every directive names a
real object with the right kind (`autoclass` vs `autofunction` vs
`autoexception` vs `autodata`).

If a warning names a specific object, fix that directive. Do not add the
warning to `suppress_warnings`.

**One contingency, pre-checked.** CI builds docs with
`uv sync --package qscat --group docs`, which does **not** install the `plot`
extra, so matplotlib is absent there. `qscat.viz` was verified to import and
to have every public signature and type hint resolve with matplotlib blocked,
so `automodule:: qscat.viz` is expected to build clean. If it nevertheless
warns about an unresolved matplotlib annotation, the fix is to add
`--extra plot` to the `uv sync` line in `.github/workflows/docs.yml` — not to
suppress the warning and not to drop names from `viz.md`.

- [ ] **Step 8: Confirm nothing else referenced the deleted page**

Run: `grep -rn "docs/api\.md\|{doc}\`api\`" --include="*.md" --include="*.py" --include="*.yml" . | grep -v _build | grep -v "docs/superpowers"`

Expected: no output. Fix any live document that still points at `docs/api.md`.
Frozen files under `docs/superpowers/` are excluded from the grep on purpose —
leave them alone.

- [ ] **Step 9: Commit**

```bash
git add docs/api docs/index.md
git rm docs/api.md
git commit -m "docs: split the API reference into one page per subpackage

Replaces a single 119-name page with docs/api/, a landing page plus ten
subpackage pages. qscat.viz (9 names) and qscat.units (4) are documented for
the first time -- the single page omitted them entirely.

Small modules carry an autosummary index then automodule detail; qscat.core's
37 names are grouped into eight sections with explicit directives. The index
is load-bearing rather than decorative: automodule alone puts no individual
name in the Markdown, so it is what the coverage test reads.

The gate added in the previous commit is now green.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01NxwNdBLUXBampDLrusGrAS"
```

---

### Task 3: Publish the theory notes

18 substantial notes in `docs/physics/` are excluded from the site by
`conf.py`, so none of the project's depth is visible on GitHub Pages.

A trial build with them included **already succeeded under `-W --keep-going`
with zero warnings**, and all 15 image references resolved. Expect no
link-fixing. If the build does warn, fix the specific file it names rather
than re-excluding the directory.

**Files:**
- Modify: `docs/conf.py` (`exclude_patterns`, module docstring)
- Modify: `docs/index.md` (Theory toctree section, quick link)
- Modify: `docs/physics/README.md` (stub → section index)

**Interfaces:**
- Consumes: `docs/index.md` as left by Task 2 — the `api/index` toctree entry must survive unchanged.
- Produces: nothing later tasks depend on. This is the last task in Plan A.

- [ ] **Step 1: Un-exclude the physics directory**

In `docs/conf.py`, remove the `"physics",` entry from `exclude_patterns`,
leaving `"_build"`, `"Thumbs.db"`, `".DS_Store"`, `"adr"`, and `"superpowers"`.

Update the comment above it, which currently claims the site is
"index/getting-started/api only", to say that `physics/` is now published as
the Theory section while `adr/` and `superpowers/` remain repository-only.

In the module docstring, replace the sentence about numpydoc conversion being
"roadmap Part 2" only if it is inaccurate after this change — otherwise leave
the docstring's build instructions as they are.

- [ ] **Step 2: Rewrite the physics section index**

`docs/physics/README.md` is a three-line stub that tells authors to "link from
the method's spec" — fine as a directory note, useless as a published landing
page. Replace its entire contents with:

```markdown
# Theory notes

One note per method: the derivation, the equations, the unit conventions
(atomic units throughout), the validation evidence, and the literature it
comes from. These are the working notes behind the implementation, so they
record limitations and negative results as well as what works.

## Discretisation

- {doc}`femdvr-ecs` — the FEM-DVR grid with an exterior-complex-scaled tail,
  and the four analytic benchmarks that pin it down.
- {doc}`nd-tensor-hamiltonian` — the N-dimensional sparse tensor Hamiltonian.
- {doc}`discretisation-tuning` — deriving a grid from the potential instead of
  hand-picking element lengths, and where the 1-D probes are not sufficient.
- {doc}`mumps-sparse-backend` — the complex-symmetric MUMPS backend.
- {doc}`ti-energy-sweep-reuse` — reusing the symbolic factorization across an
  energy sweep.

## The scattering engine

- {doc}`qscat-core-scattering` — the model-independent engine and the
  model/engine split it enforces.
- {doc}`n2-resonance` — locating the resonance pole.
- {doc}`n2-cross-section` — the one-dimensional time-independent route.
- {doc}`n2-2d-cross-section` — the exact two-dimensional driven solve, gated
  against independent published data.

## Time-dependent routes

- {doc}`n2-td-cross-section` — wavepacket propagation in one dimension.
- {doc}`n2-2d-td-cross-section` — the exact two-dimensional time-dependent
  route, including why order-1 Crank–Nicolson was not enough.
- {doc}`td-extractors` — three energy extractors sharing one propagation.
- {doc}`td-da` — the dissociative-attachment generalization.

## Molecules and approximations

- {doc}`diatomic-ve-cross-sections` — NO and F₂, and the local-complex-potential
  approximation measured against the exact oracle.
- {doc}`h2plus-dr` — dissociative recombination for an ionic target.
- {doc}`lcp-resonance-levels` — Born–Oppenheimer quasi-bound levels in the
  complex curve.

## Open directions

- {doc}`angular-coupled-channels` — the parked angular extension.
- {doc}`optimization-targets` — where the remaining hot paths are.
```

- [ ] **Step 3: Add the Theory section to the site toctree**

In `docs/index.md`, add a second toctree block after the existing one:

````markdown
```{toctree}
:maxdepth: 1
:caption: Theory

physics/README
physics/femdvr-ecs
physics/nd-tensor-hamiltonian
physics/discretisation-tuning
physics/mumps-sparse-backend
physics/ti-energy-sweep-reuse
physics/qscat-core-scattering
physics/n2-resonance
physics/n2-cross-section
physics/n2-2d-cross-section
physics/n2-td-cross-section
physics/n2-2d-td-cross-section
physics/td-extractors
physics/td-da
physics/diatomic-ve-cross-sections
physics/h2plus-dr
physics/lcp-resonance-levels
physics/angular-coupled-channels
physics/optimization-targets
```
````

That is all 19 files in `docs/physics/` (18 notes plus the README index).
Every one must be listed, or `-W` fails the build on a document not included
in any toctree.

Also update the quick link in `docs/index.md`: replace
`- **Theory notes:** the `docs/physics/` directory in the repository` with
`- **Theory notes:** {doc}`physics/README``.

- [ ] **Step 4: Build the docs with the CI command**

Run: `uv run sphinx-build -b html -W --keep-going docs docs/_build/html`

Expected: `build succeeded.` with no warnings, and the log showing
`copying images...` for all 15 `physics/figures/*.png`.

- [ ] **Step 5: Run the doctest builder, as CI does**

Run: `uv run sphinx-build -b doctest docs docs/_build/doctest`

Expected: `3 passed and 0 failed.` and `build succeeded.`

CI runs this second builder after the HTML one, so a docs change that passes
`-W` can still break the pipeline here. Publishing the theory notes brings
their code blocks into doctest's scope; none currently contain `>>>` blocks,
so the count should stay at the 3 existing doctests. If a new failure appears,
it is a real docs/API drift — fix the example, do not add `:skipif:`.

- [ ] **Step 6: Verify the rendered output is actually there**

Run: `ls docs/_build/html/physics/ | head -25` and
`ls docs/_build/html/api/`

Expected: an `.html` for each of the 19 physics documents, and for each of the
11 API pages. A build can succeed while a page silently fails to be written,
so check the files exist rather than trusting the exit code alone.

- [ ] **Step 7: Run the full non-slow suite**

Run: `uv run --no-sync pytest -n 8 -m "not slow" -q`

Expected: no failures. Nothing in this plan touches library code, so any
failure is either pre-existing or caused by the new root `tests/` package
being collected — investigate rather than deselecting it.

- [ ] **Step 8: Commit**

```bash
git add docs/conf.py docs/index.md docs/physics/README.md
git commit -m "docs: publish the theory notes as the site's Theory section

conf.py excluded docs/physics/ from the build, so 18 notes -- the derivations,
validation evidence, and documented limitations behind every method -- were
invisible on the published site while being the deepest documentation the
project has.

Needed no link-fixing: every markdown link in the notes is a relative image
reference, and all 15 resolve and copy correctly. Verified with the CI
command (-W --keep-going): build succeeded, zero warnings.

docs/physics/README.md was a three-line note to authors; it becomes the
section's landing page, grouping the notes by topic.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01NxwNdBLUXBampDLrusGrAS"
```

---

## Definition of done

- `docs/api/` has 11 pages; `docs/api.md` is gone; nothing live references it.
- `tests/test_api_docs_coverage.py` passes, and fails if a public name is added
  without a docs entry. Verify this by hand once: add a throwaway name to some
  `__all__`, watch the test fail, then revert.
- `qscat.viz` and `qscat.units` are documented.
- All 19 `docs/physics/` documents render on the site with their figures.
- `uv run sphinx-build -b html -W --keep-going docs docs/_build/html` succeeds
  with zero warnings.
- `uv run --no-sync pytest -n 8 -m "not slow" -q` is green.

## Explicitly not in this plan

These belong to Plan B and must not be started here: the showcase gallery, the
`reference:` config block, region-split domain colouring in `qscat.viz`, any
run on sadaharu, figure dedup, and the README rework — including the stale
`validation/diatomic/da_curves.py` reference at `README.md:104`, which stays
broken until Plan B repoints it at a gallery entry.
