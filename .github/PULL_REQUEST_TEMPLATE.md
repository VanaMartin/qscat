## What

## Why

## Checks

- [ ] `uv run pytest -m "not slow" -n auto --dist loadfile` green
- [ ] `uv run ruff check .` and `uv run mypy libs/qscat/qscat apps/qscat-run/qscat_run` clean

## Validation tier

Does this change touch calculation-bearing source (`libs/qscat/qscat/`,
`native/`, `projects/`, `apps/qscat-run/qscat_run/`)?

- **No** — say so; no label needed.
- **Yes** — name the `validate:*` label that covers what changed and apply
  it, or state why the change cannot move a number (a pure refactor, a
  docstring, plumbing asserted on toy decks). The Validation workflow's
  advisory note will call out the paths either way; the judgement is yours,
  not the path filter's.
