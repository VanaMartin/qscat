# /review-ready

The pass that takes a branch from "the code works" to "a stranger can review
it." Ends with the draft → ready flip, which stays a human decision.

## Preconditions

Refuse and say why if any hold:

- On `main` — this operates on a branch.
- Uncommitted changes — settle the tree first; you cannot audit a moving target.
- Tests not yet passing — this is a cleanup pass, not a debugging pass. Use
  `superpowers:systematic-debugging`.

## Step 1 — Establish state

```bash
git branch --show-current && git status --short
git log --oneline $(git merge-base main HEAD)..HEAD
gh pr view --json number,state,isDraft,title 2>/dev/null || echo "no PR yet"
```

Report: branch, commit count, PR state. If a PR already exists and is *not* a
draft, say so — the flip step is then a no-op and the rest still applies.

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
- files added "to be safe" that nothing imports — check with `grep -rn`
- `docs/` files that duplicate content now living elsewhere

Deleting a whole unit is usually right where trimming it is wrong. If you are
unsure whether something ships, that uncertainty is the answer: it does not.

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

Foreground only — backgrounded `pytest` here returns exit 0 with empty output.

```bash
uv run --no-sync pytest libs/qscat/tests -q -m "not slow"   # ~5 min
uv run --no-sync pytest apps/qscat-run/tests -q             # ~20 min with slow
uv run ruff check .
uv run mypy libs/qscat/qscat
```

Run the slow library group (~60 min) when the branch touched `libs/qscat`.

**Report actual numbers.** "Tests pass" is not a result; "374 passed, 9 skipped"
is. If you did not run something, say which and why — never imply coverage you
do not have.

Known-clean baselines, so you do not chase them: `mypy libs/qscat` *including*
tests has ~205 pre-existing findings. The repo IS `ruff format`-clean —
a repo-wide run is a no-op, so a file it wants to rewrite is one of yours.

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
