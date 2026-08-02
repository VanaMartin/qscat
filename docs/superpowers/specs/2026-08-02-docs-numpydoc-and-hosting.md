# Docs completion: numpydoc docstrings + hosting — Design Spec

**Date:** 2026-08-02
**Status:** Approved design (roadmap Part 2, the follow-on after the scaffold)
**Lifecycle:** docs/tooling on top of validated code.

## Goal

Turn the Sphinx scaffold (already committed: `docs/conf.py` + index/getting-started/
api, autodoc + furo + intersphinx, `docs.yml` CI) into a genuinely top-tier API
reference by (1) converting public docstrings to **numpydoc** so autodoc renders
parameter/return/raises tables, and (2) going live on GitHub Pages.

## Why

The site builds today, but the docstrings are prose (only ~8 structured-section
markers across all of qscat), so the API pages render descriptions without the
parameter/return tables users expect from scipy/numpy/astropy. numpydoc is the
scientific-Python standard and `napoleon` is already enabled.

## Non-goals

- No prose rewriting for its own sake — convert the *structure* (Parameters/
  Returns/Raises/Examples), preserve the existing (good) explanatory prose.
- No tutorials expansion beyond getting-started here (a later docs task).
- No change to code behavior.

## Design

- Convert public-API docstrings (the names in each module's `__all__`) to numpydoc
  sections: `Parameters`, `Returns`, `Raises`, `Examples`, `Notes`. Prioritize the
  user-facing surface: `ScatteringProblem`, the functional solvers, `qscat.dvr`
  grid builders, `qscat.model`, `qscat.linalg.SparseLU`, `qscat.exceptions`.
- Enforce with `ruff`'s `D` (pydocstyle, numpy convention) on the library package
  and the `numpydoc` validation checks; add both to CI incrementally
  (per-module `select`, not a big-bang repo-wide flip).
- Add runnable `Examples` where cheap; enable `doctest` in CI (`sphinx-build -b
  doctest`) for those examples so the docs cannot drift from the API.
- Tighten `docs.yml` to `sphinx-build -W` once the tree is warning-clean (it is at
  3 warnings now); enable **GitHub Pages** (repo Settings → Pages → GitHub
  Actions) so the existing `deploy` job publishes.
- A conventions page (atomic units, FEM-DVR-ECS notation) linked from the API ref,
  defined once.

## Validation

- `uv run sphinx-build -b html docs docs/_build/html` builds with **0 warnings**
  after conversion (currently 3).
- `sphinx-build -b doctest` passes for every `Examples` block.
- Spot-check rendered pages: `ScatteringProblem` and `ve_cross_section` show
  parameter/return tables.

## Deliverables

- numpydoc docstrings across the public API; `ruff`/`numpydoc` checks wired into
  CI per module; doctest job in `docs.yml`; `-W` enabled; Pages live; a
  conventions page.

## Verification

Docs CI green (build + doctest, `-W`); the published Pages URL renders the API
reference with tables; `ruff check` (with `D`) clean on the converted modules.
