---
name: mastering-github
description: Use when preparing a qModeling branch for review, flipping a draft PR to ready, cleaning up a fix-on-fix commit history, or deciding whether a file may cite a design spec, plan, scratch note, issue number, or PR number.
allowed-tools: Read, Grep, Glob, Edit, Write, Bash
---

# Mastering GitHub (qModeling)

Getting qModeling work from "the code passes" to "a stranger who cloned only
this repository can read it, review it, and build on it."

Two procedures live here, each with its own command:

| Command | Procedure | Reference |
|---|---|---|
| `/review-ready` | Full pre-flip pass: dissolve, de-reference, self-audit, tidy, flip draft → ready | `references/review-ready.md` |
| `/tidy-history` | Rewrite a fix-on-fix branch into readable commits, re-homed on current `main` | `references/tidy-history.md` |

## The core principle: main must stand alone

**A reader who has cloned the repository and nothing else must be able to
understand every shipped file.**

Not the session that produced it. Not the design spec. Not the plan, the
progress ledger, the PR thread, or the issue. A fork gets the tree and the git
history — nothing else. Anything a shipped file *needs* in order to make sense
must be *in* the tree.

This is not a style preference. `libs/qscat` publishes to PyPI; a docstring
that says "see docs/superpowers/specs/2026-08-15-…-design.md" is a dead pointer
for every user who installed the wheel.

Two behaviours follow, and they are the ones agents get wrong:

1. **Dissolve** durable content out of working files into permanent homes,
   *before* review — see `references/self-sufficiency.md`.
2. **Cite working files only when inlining would pollute** the code or the
   note — and then only as an optional pointer, never as a load-bearing one.
   The predicate and the test are in `references/self-sufficiency.md`.

## What counts as a working file

Anything produced to *make* the change rather than to *be* the change:

- `docs/superpowers/specs/*`, `docs/superpowers/plans/*` — design and plan docs
- `.superpowers/sdd/*` — ledgers, briefs, task reports, review packages
- session notes, scratch scripts, `/tmp` output
- GitHub issue numbers, PR numbers, PR URLs, review-comment links
- branch names, agent names, commit SHAs used as narrative ("as decided in
  a1b2c3d")

Git history and CHANGELOG entries are exempt: a commit message or changelog
line *may* reference a PR, because both travel with the clone and neither is
load-bearing for understanding the code.

## Red flags — stop and apply the skill

- "Link the spec from the PR so the reviewer can find it" — the baseline
  failure this skill exists to prevent. The reviewer may be a stranger; put the
  content where they already are.
- "The design rationale is in the plan, no need to repeat it"
- "I'll reference the issue for context"
- "History tidying is the branch owner's call" — it is part of review-ready
- "The spec explains why, so the docstring can be short"

All of these mean: dissolve the content first, then decide whether a pointer
still earns its place.

## When NOT to use this

- Mid-implementation. Dissolving and de-referencing belong at the end, once
  the content has settled; doing it early means doing it twice.
- On `main` directly. Both procedures assume a branch.
- For the design work itself — that is `superpowers:brainstorming` and
  `superpowers:writing-plans`. This skill is what happens *after* them.

## Repository facts these procedures rely on

- Default branch `main`; remote `origin` is `VanaMartin/qscat`.
- Tests: `uv run --no-sync pytest`. `libs/qscat/tests -m "not slow"` ≈ 5 min;
  the slow group ≈ 60 min; `apps/qscat-run/tests` ≈ 20 min with slow.
- Gates: `uv run ruff check .`; `uv run mypy libs/qscat/qscat` (clean).
  **`mypy libs/qscat` including tests has ~205 pre-existing findings, and the
  repo is NOT `ruff format`-clean at HEAD** — a repo-wide `ruff format` rewrites
  44 untouched files. Format only files you edited.
- Backgrounded `pytest` in this environment returns exit 0 with an empty output
  file. Run verification in the foreground.
- Commit trailers are required — see `references/tidy-history.md` before any
  rewrite, because rewriting drops them silently.
