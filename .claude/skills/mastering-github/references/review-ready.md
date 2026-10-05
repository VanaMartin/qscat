# /review-ready

The pass that takes a branch from "the code works" to "a stranger can review
it." Ends with the draft → ready flip, which stays a human decision.

## Preconditions

Refuse and say why if any hold:

- On `main` — this operates on a branch.
- Uncommitted changes — settle the tree first; you cannot audit a moving target.
- Relevant checks not yet passing — settle failures before the cleanup pass.
  Use `superpowers:systematic-debugging` when available; otherwise reproduce and
  diagnose the failure directly.

## Step 1 — Establish state

```bash
git fetch origin
git branch --show-current && git status --short
git log --oneline origin/main..HEAD
gh pr view --json number,state,isDraft,title,baseRefName,headRefName
```

Report the branch, intended PR base, commit scope, and PR state. Confirm inherited
commits belong in that scope before publication. If the PR targets a different
base, inspect its delta against that base. Distinguish a missing PR from a failed
GitHub/authentication request. A non-draft PR makes the flip a no-op.

## Step 2 — Dissolve and prune

**REQUIRED: follow `self-sufficiency.md` in this directory.** Dissolve first,
prune second. This is the step agents skip; it is the reason this command
exists.

Report concretely: which durable blocks moved where, and every reference
removed or kept with the reason it survived the test.

## Step 3 — Remove what should not ship

Delete, do not tidy:

- scratch scripts, `/tmp` outputs, debug prints, commented-out code
- tests that assert nothing, or that only restate the implementation
- files whose lack of use is confirmed by the five dynamic-reference checks in
  `code-mapping`, not merely an import search
- `docs/` files that duplicate content now living elsewhere

Before deleting code or prose, apply `code-quality-judging`'s provenance check:
preserve unique rationale, assumptions, and measurements in their durable home.
An unresolved orphan or uncertain purpose is a finding to investigate, not
permission to delete. For a supposedly redundant test, identify the existing
test that covers the same realistic regression before removing it.

## Step 4 — Comment and docstring durability

For each comment/docstring the branch added, ask: **does this describe the code,
or the process that produced it?**

Process comments rot and mean nothing to a stranger. Rewrite them as statements
about the code, or delete them.

| Rot | Durable |
|---|---|
| "changed per review feedback" | (delete) |
| "the plan says to use 0.01 here" | "0.01 resolves the pole walk through the crossing; 0.05 leaves a spurious Γ ~2e-5" |
| "TODO: task 7 will wire this" | (delete, or a real issue outside the tree) |
| "workaround for now" | "…because `Vd` is complex in the ECS tail; a real cast would break the continuation" |

A number in a comment must say what it establishes, not where it came from.

## Step 5 — Verify

Choose checks from the actual change and current CI contracts. Await completion,
retain output, and report real exit status, counts, skips, and missing capabilities;
a launched process or an empty capture is not verification evidence. Run numerical
jobs centrally so parallel reviewers do not contend for memory or output paths.

```bash
uv run --no-sync pytest -q -m "not slow" -n auto --dist loadfile
uv run ruff check .
uv run ruff format --check .
uv run mypy libs/qscat/qscat apps/qscat-run/qscat_run
```

Pin Linux BLAS threads to one per worker. Select relevant production suites when
the change could move a number and run them serially. For guidance-only changes,
verify affected paths, commands, and instruction contracts instead of running
unrelated numerical calculations. State which checks were selected and why.

**Report actual numbers.** "Tests pass" is not a result; "374 passed, 9 skipped"
is. If you did not run something, say which and why — never imply coverage you
do not have.

Type-check shipped package paths, not `libs/qscat` including tests. Format only
edited files; the read-only repository format check matches CI. Do not assume
historical finding counts or formatting exceptions describe the current baseline.

## Step 6 — Decision point

Stop. Present to your human partner:

- what was dissolved, and where it went
- what was deleted
- every reference removed, and every one kept with its justification
- verification results, with numbers
- whether history needs tidying (Step 7) — and your recommendation

Wait for a decision. Steps 7-9 rewrite history and change PR visibility; neither
is yours to take unprompted.

## Step 7 — Tidy history

If the branch accreted "fix review comment" / "oops typo" / "revert that"
commits, follow `tidy-history.md` in this directory.

Skip on a branch whose commits already read as logical units — say so rather
than rewriting for its own sake.

## Step 8 — Make the PR body self-sufficient

The PR body is read by people who never open the tree. Apply the same rule:
**it must stand alone**, and it may not send the reviewer to a working file for
anything load-bearing.

It should carry:

- what changed and why, in the reader's terms
- what the change does **not** establish — limits, caveats, known rough edges
- how it was validated, with numbers
- anything the reviewer must decide

It should not carry: a milestone table, a task list, agent names, or "see the
plan for details."

If the branch's provenance genuinely matters (a port, a correction to published
values), state the finding itself, not a pointer to where it was discussed.

Reconcile the PR body with the final diff and verification results after the last
edit or history integration. Remove descriptions of implementations or tests that
were subsequently removed.

## Step 9 — Flip

Only after Step 6's approval:

```bash
gh pr ready
```

If no PR exists, create one with `gh pr create --base main` and the Step 8 body.

The flip is the human's call. If they have not said yes, stop at Step 6 and say
what remains.

## What this does NOT do

- It does not merge.
- It does not fix failing tests or review findings — those come first.
- It does not delete `docs/superpowers/` specs and plans. They stay as the
  record; Step 2 removes the *dependency* on them, not the files.

## Output shape

Report in this order: branch state → dissolved → deleted → references
(removed / kept-with-reason) → verification numbers → history recommendation →
what you did not do.
