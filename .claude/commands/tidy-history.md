---
description: Rewrite a fix-on-fix qModeling branch into a short chain of logical, reviewable commits re-homed on current origin/main, without changing the settled end state.
allowed-tools: Read, Grep, Glob, Bash
---

Tidy the current branch's commit history.

**Read `.claude/skills/mastering-github/references/tidy-history.md` and follow it
step by step.** It is the procedure; this file only launches it.

Three things that make the difference between a clean rewrite and a lost branch:

- **Check the preconditions and refuse if one holds** — especially "the branch
  is already clean." Rewriting readable history for its own sake destroys review
  anchors and orphans review comments already attached to those SHAs.
- **Back up before touching anything** (Step 1), and **collapse onto the
  ORIGINAL merge-base before re-homing** (Steps 3 and 5, in that order).
  Reversing them turns every upstream change into a conflict inside your rewrite.
- **The end-state diff must be empty** (Step 6). That check is the whole gate. If
  it is not empty, restore from the backup and start over — never reconcile the
  difference by hand.

Preserve this repo's `Co-Authored-By:` and `Claude-Session:` commit trailers;
rewriting drops them silently.

$ARGUMENTS
