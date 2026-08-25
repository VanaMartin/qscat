---
name: release-review
description: Use when asked for a harsh whole-repo review of qModeling, a release-readiness or pre-release audit, or a consistency review after a burst of recent commits or a multi-commit campaign (docs sweep, test re-tiering, feature landing).
---

# Release Review

## Overview

A whole-repo review run as parallel read-only agents, one per track, synthesized
by theme. Core principle, measured on this repo: **the dominant defect class is
not broken code — it is unverified prose.** The 2026-08-25 review found the
numerics mypy-strict-clean and ruff-clean while a dozen docstrings stated the
opposite of the adjacent code, nine file pointers were dead, and one policy was
stated four incompatible ways. Review the words against the code as hard as the
code itself.

## Track split

Dispatch all tracks in one message (they are independent and read-only). A
generic plan covers tracks 1–4; track 5 is the one that gets forgotten and the
one that catches what single-commit review structurally cannot.

| track | scope | judge it as |
|---|---|---|
| library | `libs/qscat` + its tests | a candidate public package |
| experiments | `projects/`, `validation/`, `apps/`, `benchmarks/` | conformance to the lifecycle + layering rules |
| documentation | `docs/`, READMEs, CLAUDE.md, skills, `reference/literature` | parity with code, not prose quality |
| tooling | packaging, typing, lint, CI, Docker, Rust | 2026 scientific-OSS baseline; run the tools, paste output |
| recent changes | `git log --since=<2 weeks>` grouped into campaigns | did each campaign apply its own rules everywhere it applies |

## The checks generic plans miss

These found the critical issues; put them verbatim in the agent briefs.

1. **Prose-drift sweep.** Every `path/to/file.py` named in a docstring or doc
   must resolve (nine didn't; two were the sole provenance of hard-coded
   constants). Every "enforced/gated by test X" claim is verified against the
   test body (a skill claimed 11 rules enforced; the test checked 4). Comments
   are checked against the lines they annotate — grid dicts, tolerances,
   anchor counts, "unchanged from" claims.
2. **N-places sweep.** For each policy or measured number, grep every file
   stating it and demand agreement (the slow-tier CI policy: 6 statements, 4
   stale). A campaign that updated fewer than all N is a finding even when
   each file is individually plausible.
3. **All-branches fix search.** `git log --all --grep` for commits titled as
   fixes to results the docs present as valid or open. An unmerged
   `fix(core): sigma_DA from the outgoing flux...` on a backup branch was the
   single most serious finding: a published negative result on a possibly
   defective oracle.
4. **Campaign self-consistency.** For each recent campaign, test its own rule
   against its own output: are `@slow` markers a function of *measured* cost
   (a 17 s test was unmarked, a 0.11 s test marked); did the docs convention
   land on every note it ticked converted; does the ADR's Consequences section
   still describe shipped CI.
5. **Measure, don't trust.** Re-run mypy/ruff/clippy and paste counts; time
   the tests you make claims about; rebuild the wheel; diff numeric-token
   multisets across "no numbers changed" commits.

## Orchestration mechanics

- Brief every agent: read-only; severity-ranked (Critical/Major/Minor); every
  finding cites `file:line` and quotes **both sides** of any contradiction; no
  nitpicks a linter already gates; and a required **"what held up"** section —
  verified strengths calibrate the harshness and catch one-sided reviewers.
- Agents idle without delivering. Tell them to `SendMessage` the report to
  `"main"`; on a bare idle notification, request the report explicitly.
- Synthesize by theme, not concatenation: name the one pattern behind the
  findings, lead with the verdict, end with a prioritized fix order sized by
  effort. Publish the report as an artifact and give the link.

## Common mistakes

| mistake | reality |
|---|---|
| Docs track does proofreading | Its job is parity: claim ↔ code, number ↔ gate, path ↔ file |
| Trusting a reviewer's claim | No evidence rule → confident wrong findings; require reproduction |
| Skipping the recent-changes track | Cross-commit inconsistency is invisible to per-commit review |
| Only judging HEAD | Unmerged branches can hold fixes to what HEAD publishes |
| Harsh = one-sided | Without "what held up", severity has no calibration and trust is lost |
