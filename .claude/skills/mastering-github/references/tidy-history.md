# /tidy-history

Rebuild a settled fix-on-fix branch as logical commits on its original base,
prove that reconstruction preserves its tree, then integrate current `origin/main`.
History reconstruction and upstream integration have different verification gates.

## Preconditions

Stop if on `main`, the working tree is dirty, relevant verification is failing,
or the commits are already logical units. Obtain explicit approval before rewriting
published history or a branch containing another contributor's work. Never rewrite
commits already on the target base.

## 1. Establish scope and recovery points

Run the command blocks in one shell so their variables remain available. Record
the values and commands in an isolated temporary directory for recovery; do not
use shared fixed filenames across worktrees.

```bash
git fetch origin
BRANCH=$(git branch --show-current)
OLD_HEAD=$(git rev-parse HEAD)
NEW_MAIN=$(git rev-parse origin/main)
MB=$(git merge-base "$NEW_MAIN" "$OLD_HEAD")
BACKUP="backup/${BRANCH}-$(git rev-parse --short HEAD)"
git branch "$BACKUP" "$OLD_HEAD"
STATE_DIR=$(mktemp -d "${TMPDIR:-/tmp}/qscat-tidy.XXXXXX")
git diff --binary "$MB" "$OLD_HEAD" > "$STATE_DIR/original.diff"
printf '%s\n' "$BRANCH" "$OLD_HEAD" "$NEW_MAIN" "$MB" "$BACKUP" > "$STATE_DIR/refs"
REMOTE_HEAD=$(git ls-remote origin "refs/heads/$BRANCH" | cut -f1)
printf '%s\n' "$REMOTE_HEAD" > "$STATE_DIR/remote-head"
git log --oneline "$MB..$OLD_HEAD"
```

Check each command succeeded before proceeding. Confirm this commit set is the
intended PR scope; inherited unrelated work needs a base/scope decision first.
If the backup name already exists, choose a fresh name rather than overwriting it.
The recorded remote head is the lease for publication; fetching later must not
silently replace it.

Inspect persistent references to branch SHAs, including manifests and published
artifacts. A digest identifies content; a commit records provenance. If rewriting
would make required provenance unreachable, resolve its retention/publication with
the owner before proceeding. Do not change scientific records merely to hide an
orphaned SHA.

## 2. Rebuild on the original base

```bash
git reset --soft "$MB"
git restore --staged -- .
```

The snapshot is now in the working tree, unstaged. Stage and commit logical units,
using explicit paths or hunks, not an indiscriminate `git add -A`. Include intended
new files and deletions when rebuilding. Split by behavior and rationale rather
than file type or task number. Pure refactors are separate from behavior changes;
fixes to already-shipped defects retain their own rationale. Preserve applicable
`Co-Authored-By:` and `Claude-Session:` trailers from the commits being regrouped.

If staged/worktree state is unexpected, stop and inspect the backup. Recovery may
require discarding reconstruction changes; obtain approval before a destructive
reset and keep the backup and state directory until verification is complete.

## 3. Prove reconstruction identity before integration

```bash
git status --short
git diff --exit-code "$OLD_HEAD" HEAD
git rev-parse "$OLD_HEAD^{tree}" "HEAD^{tree}"
REBUILT_HEAD=$(git rev-parse HEAD)
printf '%s\n' "$REBUILT_HEAD" > "$STATE_DIR/rebuilt-head"
```

Require a clean working tree, an empty diff, and identical tree IDs. This compares
all tracked content, including file modes and deletions. A mismatch means the
reconstruction is incomplete or altered the result: stop and recover from the
backup, rather than explaining away the difference.

## 4. Integrate the recorded current main

```bash
git rebase --onto "$NEW_MAIN" "$MB"
```

Resolve conflicts against both the branch contract and upstream intent, rather
than reflexively keeping one side. Record every conflict resolution. On an
unexpected conflict or scope change, stop; `git rebase --abort` returns to the
verified reconstruction while a rebase is in progress.

## 5. Verify the integration separately

```bash
git range-diff "$MB..$REBUILT_HEAD" "$NEW_MAIN..HEAD"
git diff --stat "$REBUILT_HEAD" HEAD
git diff "$REBUILT_HEAD" HEAD
git diff --check "$NEW_MAIN...HEAD"
git status --short
```

Explain range-diff changes and inspect the integrated tree delta against upstream
changes and recorded conflict resolutions. Whole-tree equality with `OLD_HEAD`
is not required here: legitimate upstream changes are now present. An unexpected
branch behavior change is a failure, not an accepted consequence of rebasing.

Re-run verification affected by integration. Use the current CI targets:

```bash
uv run --no-sync pytest -m "not slow" -n auto --dist loadfile
uv run ruff check .
uv run ruff format --check .
uv run mypy libs/qscat/qscat apps/qscat-run/qscat_run
```

Pin Linux BLAS threads to one per worker. Run relevant production tests serially
when numerical behavior can change; record actual counts, skips, exit status,
and any justified exclusions. Await completion and retain output. A pre-rebase
pass does not establish correctness of the integrated tree.

## 6. Publish with the recorded lease

After approval and passing verification, push to the explicit remote branch:

```bash
if [ -n "$REMOTE_HEAD" ]; then
  git push "--force-with-lease=refs/heads/$BRANCH:$REMOTE_HEAD" origin "HEAD:refs/heads/$BRANCH"
else
  git push -u origin "$BRANCH"
fi
```

Never substitute `--force`. If the lease rejects, stop for a decision; do not
refresh it and retry automatically. Keep recovery points until publication is
confirmed. This procedure does not merge, change PR visibility, or fix review
findings as part of history reconstruction.

## Report

Report the backup/state location, original/base/upstream refs, before/after commit
counts and subjects, reconstruction identity result, integration/conflict review,
verification results, and publication status. Distinguish preserved branch content
from changes incorporated from upstream.
