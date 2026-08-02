# 4. Public API stability policy

Date: 2026-08-02

## Status

Accepted

## Context

`qscat` is heading for a PyPI release (see the hardening roadmap,
`docs/superpowers/plans/2026-08-02-qscat-hardening-roadmap.md`). At first publish
every importable name risks being treated as a stable contract. `qscat.core`
alone re-exports ~30 names, and several solver signatures are still expected to
change (Part 1 of the roadmap: bundling the `(tgrid, model, eps, chi, v_init)`
arguments into a context object). We need a written policy so users know what
they may depend on and we retain freedom to refactor pre-1.0.

## Decision

1. **Public surface.** The public API is exactly the names exported in each
   submodule's `__all__` (and re-exported by `qscat.core.__init__` etc.).
   Anything with a leading underscore, and any name not in an `__all__`, is
   private and may change without notice.

2. **SemVer, with a pre-1.0 caveat.** The project follows Semantic Versioning.
   While the version is `0.y.z`, minor releases MAY contain breaking changes to
   the public API (this is standard SemVer for pre-1.0). Breaking changes will
   be called out in `CHANGELOG.md`. From `1.0.0` on, breaking changes to the
   public API require a major-version bump.

3. **Provisional APIs.** An API still expected to change (currently: the wide
   solver signatures targeted for the context-object refactor) is marked
   *provisional* in its docstring. Provisional APIs may change in a minor
   release even after 1.0, until the "provisional" note is removed.

4. **Deprecation.** Post-1.0, a public API is deprecated for at least one minor
   release (emitting `DeprecationWarning` and noted in `CHANGELOG.md`) before
   removal.

5. **Exceptions.** qscat raises subclasses of `qscat.exceptions.QscatError` for
   the major recoverable categories (`GridError`, `ModelError`, `BackendError`,
   `ConvergenceError`). Each also subclasses the built-in it replaces, so
   catching the built-in remains valid. Generic argument validation may still
   raise a plain built-in `ValueError`/`TypeError`.

## Consequences

- Users can rely on `__all__` as the contract and on `except QscatError`.
- We keep refactoring freedom pre-1.0 and via the provisional mechanism, without
  surprising users, because breakage is disclosed in the changelog.
- Before tagging `1.0.0` we must audit every `__all__`, resolve provisional
  markers, and confirm the deprecation policy is in force.
