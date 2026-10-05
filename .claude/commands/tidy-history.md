---
description: Rebuild logical qModeling commits, prove tree identity on the original base, then separately verify integration with current origin/main.
allowed-tools: Read, Grep, Glob, Bash
---

Tidy the current branch's commit history.

**Read `.claude/skills/mastering-github/references/tidy-history.md` and follow it
step by step.** It is the procedure; this file only launches it.

Three things that make the difference between a clean rewrite and a lost branch:

- **Check the preconditions and refuse if one holds** — especially "the branch
  is already clean." Rewriting readable history for its own sake destroys review
  anchors and orphans review comments already attached to those SHAs.
- **Back up before touching anything** (Step 1), and **rebuild on the original
  merge-base before re-homing** (Steps 2 and 4, in that order). Keep upstream
  integration separate from the tree-preserving history reconstruction.
- **Prove tree identity before upstream integration** (Step 3 of the procedure).
  Then review conflict resolutions, the range-diff, and affected verification
  separately. Upstream changes can legitimately alter the integrated tree.

Preserve this repo's `Co-Authored-By:` and `Claude-Session:` commit trailers;
rewriting drops them silently.

$ARGUMENTS
