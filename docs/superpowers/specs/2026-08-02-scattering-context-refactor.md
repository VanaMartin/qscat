# Full context-object refactor of the functional solvers — Design Spec

**Date:** 2026-08-02
**Status:** Approved design (roadmap Part 1, the deferred piece)
**Lifecycle:** cleanup — no new capability; a signature change across the repo.
**Execution:** best run as Subagent-Driven Development (fresh implementer +
review per package tree) — it is ~243 call sites; the full test suite (differential
oracles) is the safety net.

## Goal

Collapse the wide `(tgrid, model, eps, chi, v_init, …)` argument group of the
functional solvers into the already-shipped context. `ScatteringProblem` (the
facade) exists and is stable; this refactor makes the **functional layer** take
the bundle too, so there is one argument style, not a bundle-facade over a wide
core.

## Why

The facade delivered the clean user API safely, but the functional signatures are
still 7+ loose positional args threaded through libs/projects/validation/apps
(~243 call sites). Doing the actual reduction removes the duplication the facade
only hides, and is the "on par with top packages" cleanup. Do it BEFORE the
functional API is treated as frozen post-1.0 (the solvers are marked provisional
per ADR 0004 precisely to keep this door open).

## Non-goals

- No numerical change — every cross section is identical to round-off (the
  differential oracles enforce this).
- No change to `ScatteringProblem`'s public surface.
- Not the module splits (separate spec).

## Design

- Introduce `qscat.core.ScatteringContext` (frozen) bundling `(tgrid, model,
  eps, chi, v_init)` — the exact group `ScatteringProblem` already holds; the
  facade becomes a thin `ScatteringContext` + solved-basis wrapper.
- Change the functional solvers to `fn(ctx, <per-call args>, *, <options>)`:
  `ve_cross_section(ctx, vprimes, E, *, ...)`, `da_cross_section(ctx, E, *, ...)`,
  etc. Keep keyword-only options unchanged.
- Migrate all call sites, one package tree per SDD task (libs → projects →
  validation → apps), running that tree's tests between tasks.
- Remove the "provisional" notes from the docstrings once landed (ADR 0004).
- Decide: keep thin deprecated positional shims for one minor release, or hard
  cut (pre-1.0, a hard cut is defensible and cleaner — recommend hard cut with a
  CHANGELOG breaking-change note).

## Validation

- The existing differential/analytic/convergence suite is the oracle: it must
  stay green unchanged (no test's expected numbers change; only call syntax).
- Add one test that `ScatteringProblem.ve_cross_section(...)` and
  `ve_cross_section(ctx, ...)` return identical arrays (facade == functional).
- `ruff`/`mypy --strict` clean; the `test_core_no_model_import` guard still holds.

## Deliverables

- `ScatteringContext` + refactored functional signatures + all call sites migrated.
- Facade reimplemented on `ScatteringContext`; provisional notes removed;
  CHANGELOG breaking-change entry.

## Verification

Full `uv run pytest` green (incl. @slow in Docker) with zero expected-value
changes; `ruff` + `mypy libs/qscat/qscat` clean; grep shows no residual
`eps, chi, v_init` positional threading through the public solvers.
