---
description: Tidy a fix-on-fix qModeling branch using the shared history-rewriting procedure.
---

Tidy the current branch's commit history.

Read `.claude/skills/mastering-github/references/tidy-history.md` and follow it
step by step. It is the procedure; this file only launches it.

Check the preconditions and refuse if one holds, especially if the branch is
already clean. Back up before touching anything. Collapse onto the original
merge-base before re-homing, in that order.

The end-state diff must be empty. If it is not, restore from the backup and
start over; never reconcile the difference by hand.

Preserve existing `Co-Authored-By:` and `Claude-Session:` commit trailers.

$ARGUMENTS
