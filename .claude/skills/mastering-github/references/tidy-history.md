# /tidy-history

Rewrite a fix-on-fix branch into a short chain of logical, readable commits
re-homed on current `origin/main` — without changing the settled end state.

For branches that accreted "fix review comment", "oops typo", "actually revert
that" while iterating. The tree is right; the history is not.

## Preconditions — hard refusals

Refuse, and say which one, if:

- **On `main`.** Never rewrite the default branch.
- **Uncommitted changes.** Settle the tree first — Step 2 depends on the end
  state being final.
- **Tests failing.** Rewriting a broken branch produces broken commits.
- **The branch is already clean.** If every commit is already a logical unit,
  say so and stop. Rewriting for its own sake destroys review anchors and
  invalidates any review comments already attached to those SHAs.
- **Someone else has commits on the branch,** unless your human partner
  explicitly says to proceed. Force-pushing over a collaborator is not yours to
  decide.

## Step 1 — Back up first, always

```bash
git branch backup/$(git branch --show-current)-$(git rev-parse --short HEAD)
git rev-parse HEAD          # record this; it is the recovery point
```

Do this even when confident. Especially then.

## Step 2 — Freeze the end state

```bash
MB=$(git merge-base main HEAD)
git rev-parse HEAD > /tmp/tidy-end-state
git diff $MB..HEAD > /tmp/tidy-end-state.diff
```

The final tree is the specification. Whatever you build in Step 4, it must
reproduce this diff exactly. Nothing in the working tree may change.

## Step 3 — Collapse onto the ORIGINAL merge-base

```bash
git reset --soft $MB
```

One staged snapshot on the base the branch actually started from. Collapse
*before* re-homing — mixing the two turns every unrelated upstream change into
a conflict inside your own rewrite.

## Step 4 — Rebuild as logical commits

Stage and commit in units a reviewer can hold in their head. Aim for units that
are individually reviewable and individually revertible; a good branch is
usually 1–6 of them, but let the work decide, not the number.

Split by **what changed and why**, never by file type or by task number:

| Good unit | Bad unit |
|---|---|
| "add the two-angle pole matcher" | "changes to pole.py" |
| "wire the CLI surface" | "task 6" |
| "fix the golden-rule crash on an empty comparator window" | "review fixes" |
| "correct the H2+ reduced mass to m_p/2" | "misc fixes" |

Rules:

- A pure refactor is its own commit, never mixed with behaviour.
- A fix that only repairs a defect introduced earlier *on this branch* belongs
  folded into the commit that introduced it — that defect never reached `main`,
  so its history has no audit value.
- A fix that corrects something already on `main` stays its own commit — that
  one is real history.
- Docs may ride with the code they document, or stand alone if substantial.

**Preserve the commit trailers.** This repo's commits carry
`Co-Authored-By:` and `Claude-Session:` lines. Rewriting drops them silently
unless you re-add them to every new message. Check with
`git log -1 --format=%b` on a pre-rewrite commit and reproduce the same
trailers.

## Step 5 — Re-home onto current main, by rebase

```bash
git fetch origin
git rebase origin/main
```

Rebase — never `reset --hard origin/main`, which discards your work while
looking like it worked. Resolve conflicts against *upstream's* intent, not by
reflexively keeping your side.

## Step 6 — Verify before claiming done

```bash
git diff $(cat /tmp/tidy-end-state) HEAD    # MUST be empty
```

**An empty diff is the whole gate.** Non-empty means the rewrite changed the
result: stop, restore from the backup branch, and start over. Do not
"reconcile" the difference — that is how a rewrite silently ships something
nobody wrote.

Then re-run verification, foreground (backgrounded `pytest` here returns exit 0
with empty output):

```bash
uv run --no-sync pytest libs/qscat/tests -q -m "not slow"
uv run ruff check . && uv run mypy libs/qscat/qscat
```

Rebasing onto a moved `main` can break a passing branch. The pre-rewrite pass
does not carry over.

## Step 7 — Push

```bash
git push --force-with-lease
```

`--force-with-lease`, never `--force`: it refuses if the remote moved under
you. If it refuses, someone else pushed — stop and ask, do not override.

## Anti-patterns

| Miss | Consequence |
|---|---|
| Rebasing onto new `main` *before* collapsing | Every upstream change becomes a conflict inside your rewrite |
| `reset --hard origin/main` to "re-home" | Silently discards the branch |
| Skipping the end-state diff check | Ships a tree nobody authored |
| Squashing everything to one commit | Throws away the review structure this command exists to create |
| Rewriting an already-clean branch | Destroys review anchors, orphans existing review comments |
| Dropping the commit trailers | Loses attribution on every commit |
| `--force` instead of `--force-with-lease` | Overwrites a collaborator |

## What this does NOT do

- It does not change the tree. If the end-state diff is non-empty, you broke it.
- It does not merge or flip a PR — see `review-ready.md`.
- It does not fix tests or review findings.
- It does not rewrite commits already on `main`.

## Output shape

Report: backup branch name → before/after commit counts → the new commit
subjects in order → end-state diff empty (yes/no) → verification numbers →
pushed (yes/no).
