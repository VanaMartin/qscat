---
name: release-review
description: Use when asked for a harsh whole-repo review of qModeling, a release-readiness or pre-release audit, a consistency review after a burst of recent commits or a multi-commit campaign (audit mode), or a structural pass over code quality, duplication, dead code and placement (structure mode).
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

## Two modes

| mode | question | output |
|---|---|---|
| **audit** (below) | does the repository still say true things about itself | a severity-ranked findings report |
| **structure** | is the code compact and readable, and is it in the right place | a remediation plan, then a measured delta |

Audit mode is the five read-only tracks in the rest of this file. Structure mode is
the seven-phase pipeline below. They are independent; run either alone.

## Structure mode

Seven phases. Phases 0 and 6 run a script, not an agent, so the closing claim is a
diff of two measurements rather than a summary of what agents believed they did.

| phase | executor | rubric lives in |
|---|---|---|
| 0 · map | `.claude/skills/code-mapping/scripts/inventory.py` | `code-mapping` |
| 1 · classify | `symbol-classifier` | `code-mapping` |
| 2 · judge | `file-quality-judge`, one per review unit, in parallel | `code-quality-judging` |
| 3 · consolidate | `consolidation-analyst` | `code-consolidation` |
| 4 · plan | you | `superpowers:writing-plans` |
| 5 · execute | `superpowers:subagent-driven-development` | this file |
| 6 · verify | `.claude/skills/code-mapping/scripts/inventory.py` again | this file |

Workspace: `.superpowers/audit/<timestamp>/`, git-ignored, never committed.

### Intent isolation is the load-bearing rule

A dispatch for phases 1, 2 or 3 states a scope, the input artifact paths, and an
output path. **Nothing else.** No goal, no campaign context, no earlier phase's
conclusions.

A judge told "we are trying to make this compact" will find things to compact,
because that is what it was asked to help with. The same judge told "grade these
files against this rubric" reports what is there. The analyst never sees the
judges' verdicts, so it rules on similarity evidence rather than ratifying an
opinion.

Each agent loads its own rubric from its own skill. Rubrics are never pasted into a
dispatch.

### Review units

Partition the scope into units of about 2500 lines so each judge holds its unit at
once, and dispatch every judge in one message. `libs/qscat/qscat/core` is 11345
lines and splits six ways (measured 2026-08-27):

| unit | contents | lines |
|---|---|---|
| `nrm-ti` | `nrm/` minus its three time-dependent modules | 1788 |
| `nrm-td` | `nrm/{extended,propagation,td_cross_section}.py` | 1462 |
| `lcp` | `lcp/` plus `dissociation.py` | 1857 |
| `td` | `time_dependent.py` plus `td_extractors.py` | 2179 |
| `states` | `bo.py`, `assignment.py`, `resonance.py` | 1564 |
| `rest` | the remainder of `core/` | 2495 |

`nrm` splits rather than riding at 3250 lines because its two halves are separate
subsystems — a time-independent resolvent route and a time-dependent propagation
route — so the boundary is real rather than a size convenience.

### Task order in the plan

**Delete, then merge, then move, then rewrite.** Every deletion shrinks the surface
the later steps must reason about, and deletions are the cheapest changes to
verify. A plan that rewrites before deleting spends its effort on code that was
about to disappear.

### The two-lane gate

Every task declares a lane.

**Lane A — structural.** Dead code, comments, file splits, holder unification,
renames, moves.

```bash
uv run --no-sync pytest -m "not slow" -n auto --dist loadfile
uv run --no-sync mypy libs/qscat/qscat apps/qscat-run/qscat_run
uv run ruff check . && uv run ruff format --check .
```

The mypy target is wider than CLAUDE.md's `libs/qscat/qscat`, deliberately: a
structural pass moves code across the app boundary too. Measured 2026-08-27, the
two-target form is clean over 85 files.

**Lane B — a numerics path is touched.** Lane A's gate, plus: run
`.claude/skills/code-mapping/scripts/capture_observables.py` before the edit, make
the change, run it again, and compare **bitwise**.

Bitwise is the correct bar, not a harsh one — a behaviour-preserving refactor
reorders no floating-point operation. A mismatch has found exactly what the audit
exists to catch: two things that looked the same and were not. Explain it or revert
it; never loosen the comparison.

**An empty diff only covers the modules some case actually exercises.** The capture
script ships a small fixed case set, and each case declares the modules it reaches.
Before claiming a lane-B pass, check that the modules you edited appear in some case's
`modules` list — nothing enforces that intersection, so a refactor outside the covered
set produces an empty diff that means *nothing was measured*, not *nothing changed*.
Measured on the 2026-08-28 run: both refactors touched `qscat.core` and `qscat.tuning`,
which no case covers, so their empty diffs were vacuous and the real evidence was the
suite plus behaviour tests written for the change. When you edit outside the case set,
either add a case or say plainly in the report that the gate did not cover it.

### Phase 6

1. Re-map into `after/` and emit `delta.md`: symbols removed, clone clusters
   closed, holders unified, files still over length, symbols with no static
   consumer, total lines.
2. Confirm every defect the plan promised to fix is absent from the new map. One
   still present is reported, not dropped.
3. Re-run the prose-drift check (audit mode, check 1) over every touched file.
   Moving code is the standard way to strand a `path/to/file.py` reference.
4. Delete the workspace.

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

This section and Common mistakes below describe audit mode's five-track dispatch;
structure mode's orchestration is in the Structure mode section above.

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
