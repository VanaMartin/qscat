# Self-sufficiency: dissolving working files and pruning references

Two operations, in this order. Dissolve first — you cannot judge whether a
reference earns its place until the content it points at has a permanent home.

---

## 1. Dissolve

**Dissolving means: move the durable content out of the working file into the
permanent home where a reader would look for it, in that home's own voice.**

It is not copying, and it is not summarising the spec. A spec says "we will do
X because Y"; the physics note says "X, because Y." Design docs are written
forward-looking and provisional; permanent docs are written as settled fact.
Rewrite accordingly.

### What is durable

Ask of each block: **would this still be true and useful to someone six months
from now who never saw the branch?**

| Durable — dissolve it | Ephemeral — leave it | 
|---|---|
| The physics: what is computed, the operator, the conventions | Milestone tables, task lists, step ordering |
| Measured numbers and what they establish | "Task 4 will call this once" |
| Limits, caveats, failure modes, noise floors | Review findings and how they were fixed |
| Why an approach was chosen over a live alternative | Which agent did what, what the plan estimated |
| Citations to literature | Citations to the plan itself |
| Naming/terminology decisions a reader must share | Test counts at a point in time |

### Where it goes

| Content | Permanent home |
|---|---|
| Method, derivation, validation, limits | `docs/physics/<topic>.md` |
| What a function does, its contract, its units | the docstring |
| Where a capability lives, one-paragraph orientation | `CLAUDE.md` repo map |
| User-facing behaviour change | `CHANGELOG.md` |
| How to run it | `apps/qscat-run/README.md` + an `examples/*.yaml` |
| Architectural decision with lasting consequences | `docs/adr/` |

If content has no permanent home, that is the finding: either it is ephemeral
(leave it) or the repo is missing a doc (create it).

### Procedure

1. List the working files this branch produced or leaned on.
2. Read each; classify every block durable/ephemeral by the test above.
3. For each durable block, find its home in the table and write it there **in
   that home's voice**. Do not paste.
4. Re-read the permanent doc alone, as a stranger. Does it stand up with the
   spec deleted? If not, something durable is still stranded.
5. Only now prune references (below).

The working files stay in the repo. They are the record of how the work was
made. Dissolving is about the *shipped* files no longer depending on them.

---

## 2. Prune references

The rule is a conditional, not a ban:

> **Cite a working file only when inlining the content would pollute the file
> doing the citing. When you do cite, the reference must be optional — the
> reader must not need it.**

### The test — apply to each reference, in order

1. **Delete the reference mentally. Is the file still complete?**
   - No → the reference is load-bearing. It fails. Dissolve the content it
     carries, then re-apply the test.
   - Yes → continue.
2. **Would inlining what it points at pollute this file?** "Pollute" is
   observable, not a feeling: a docstring growing past its function, a physics
   note acquiring a milestone table, a README carrying a derivation.
   - Yes → the citation earns its place. Keep it, phrased as an optional
     pointer: "the design rationale is recorded in `<path>`", not "see `<path>`
     for the method."
   - No → inline the content and delete the reference.

A reference that survives both steps is one a reader may follow for *more*, never
one they must follow to understand what they are reading.

### Never citable from a shipped file

These fail step 1 by construction — they do not travel with a clone:

- issue numbers, PR numbers, PR/review URLs
- `.superpowers/sdd/*` — ledgers, briefs, task reports, review packages
- session notes, agent names, `/tmp` paths
- commit SHAs used as narrative ("as decided in a1b2c3d"). A SHA cited as
  *provenance* in a CHANGELOG or commit message is fine.

Exempt, because they travel with the clone and are not load-bearing:
**commit messages** and **CHANGELOG entries** may cite a PR number.

### Where to look

```bash
# working-file references in shipped content
grep -rn "docs/superpowers\|\.superpowers/sdd" --include="*.py" --include="*.md" \
  libs apps projects validation docs/physics CLAUDE.md README.md

# issue / PR references in shipped content
grep -rn "PR #[0-9]\|pull/[0-9]\|issue #[0-9]" --include="*.py" --include="*.md" \
  libs apps projects validation docs/physics CLAUDE.md README.md
```

Check the PR body too — it is read by people who never open the tree.

---

## Worked example

`qscat.core.lcp.lcp_resonance_levels`'s docstring ended with:

> See `docs/superpowers/specs/2026-08-15-bo-lcp-resonance-levels-design.md`.

Step 1: delete it — is the docstring complete? It was not. The docstring
described the call but not the two-angle selection criterion the user must
understand to trust the output. **Load-bearing → fails.**

Dissolve: the selection criterion, the normalization choice and the noise floor
went into `docs/physics/lcp-resonance-levels.md`; the contract and units went
into the docstring itself.

Re-apply: the docstring now cites `docs/physics/lcp-resonance-levels.md` — a
path inside the clone, carrying content a user needs. The spec is no longer
cited from shipped code at all, and nothing was lost: it remains in the repo as
the record of how the design was reached.

Contrast, a citation that legitimately survives: a physics note may say "the
convergence study behind this tolerance is recorded in
`docs/superpowers/plans/…`" — inlining a twenty-row table of grid refinements
would pollute the note, and a reader who skips the pointer still understands
the tolerance.

---

## Rationalizations

| Excuse | Reality |
|---|---|
| "Link the spec from the PR so the reviewer can find it" | The reviewer may have only the clone. Put the content where they already are. This is the observed baseline failure. |
| "The rationale is in the plan, no need to repeat it" | A fork has the plan only if it is committed — and even then, shipped code that *needs* it is incomplete. Dissolve, then decide. |
| "I'll reference the issue for context" | Issues do not travel with a clone. Never citable from shipped files. |
| "The spec explains why, so the docstring can be short" | Short is good; incomplete is not. Shorten by cutting words, not by outsourcing meaning. |
| "It's only a doc, not code" | `docs/physics/` is shipped. Same rule. |
| "Deleting the reference loses the audit trail" | The working file stays in the repo and in git history. Only the *dependency* is removed. |
