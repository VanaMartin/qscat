# Contributing to qModeling / qscat

Thanks for your interest in contributing. qModeling is a Python-first
quantum-scattering research monorepo; `libs/qscat` is the reusable library
published to PyPI. This file is the short version; `CLAUDE.md` is the full
operating manual.

## Development setup

Requires Python 3.12 and [uv](https://docs.astral.sh/uv/).

```bash
git clone https://github.com/VanaMartin/qscat
cd qscat
uv sync --all-packages   # NOT plain `uv sync` — that prunes the workspace members
```

`--all-packages` installs `qscat`, the `qscat-run` CLI, and builds the Rust
`qscat_kernels` extension.

## Checks (run before opening a PR)

```bash
uv run ruff check .            # lint
uv run ruff format .           # format
uv run mypy libs/qscat/qscat   # types (strict, clean over the library)

# The fast tier -- toy-scale, and the one CI gates on.
uv run pytest -m "not slow" -n auto --dist loadfile
```

CI runs the same on 3.12/3.13, plus a clean-venv import check and `twine check`.

### The two test tiers

The suite is tiered by COST, not by importance — see
`docs/adr/0005-test-tiers-fast-and-slow.md`. The default tier is toy-scale and
runs on every push; `@pytest.mark.slow` is production-scale physics (real
molecule decks, converged grids, multi-thousand-step propagations) and runs
only locally, in the Docker `test` image, or on demand in CI.

```bash
uv run pytest -m slow          # the production tier: SERIAL, minutes and gigabytes
docker/build.sh test           # both tiers, with MUMPS
```

Two rules worth knowing before you add a test:

- **A test needing more than a few seconds or ~0.5 GB belongs in `slow` — or
  wants a smaller deck.** Ask which applies. If what is under test is plumbing
  (shapes, cadences, dispatch), assert it on a toy deck and keep it in the gate.
- **Never pair `-n` with the full suite.** Parallelising the slow tier
  schedules multi-GB decks against each other with nothing bounding the total;
  one example alone peaks near 19 GB. Parallelise the fast tier or nothing.

### Asking CI to run the production tier

CI does not run the slow tier by default, so a change that moves a number can
merge green. When reviewing, apply the label covering what changed and the
**Validation** workflow runs that suite:

`validate:core` · `validate:n2` · `validate:diatomic` · `validate:factory` ·
`validate:h2plus` · `validate:tuning` · `validate:run` · `validate:all`

A path filter cannot make this judgement — a one-line kernel edit can
invalidate every cross section while a large refactor changes nothing — so the
workflow only nudges: when calculation-bearing source changes with no label, it
writes an advisory note to the run summary and does not fail the build.

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

## The method lifecycle

Every numerical capability moves through five stages (enforced by the
`qm-method-lifecycle` skill): **design** (`docs/physics/`) → **toy model**
(`projects/`) → **validate** (`validation/`, against analytic benchmarks /
conservation laws / convergence / a reference) → **optimize in Rust**
(`native/`, keeping the Python version as the differential oracle) → **promote
to `qscat`** (`libs/qscat`). Code is CPU-runnable on a laptop and containerizable
at every stage.

## Conventions

- **Atomic units** throughout (`qscat.units`); no ad-hoc conversions in method code.
- **Errors**: raise `qscat.exceptions` types (`GridError`, `ModelError`,
  `BackendError`, `ConvergenceError`) for recoverable categories; generic
  argument validation may raise built-in `ValueError`/`TypeError`.
- **Public API**: exported names in each module's `__all__`; see ADR 0004 for the
  stability policy. New public API needs a docstring and a test.
- **Testing**: numerical code is validated by differential/analytic/convergence
  tests, not exact equality (see the `numerical-validation` skill). Every Rust
  kernel keeps a Python differential oracle.
- **Commits**: conventional-commit style (`feat(scope):`, `fix(scope):`, …).
  Update `CHANGELOG.md` under `[Unreleased]` for user-facing changes.

## Adding a molecule

A new molecule is a registry entry in `qscat.model` + validation — never new
solver code. See `docs/physics/qscat-core-scattering.md`.

## Reporting issues

Open an issue at https://github.com/VanaMartin/qscat/issues with a minimal
reproducer (a small grid config is ideal).
