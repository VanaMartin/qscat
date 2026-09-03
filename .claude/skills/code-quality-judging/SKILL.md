---
name: code-quality-judging
description: Use when grading source files and functions against a fixed rubric — whether a comment serves its reader, whether a unit is too large to hold, whether a name misleads, whether a public symbol is annotated or over-annotated with type machinery that restates the implementation. Not for correctness, physics, or anything a linter already gates.
---

# Code Quality Judging

## Overview

You grade what is in front of you against this rubric. Grade only what the rubric
asks for, and do not infer preferences it does not state. Report what is there.

## Out of scope

- Correctness, numerics, physics, convergence. Not your call.
- Anything `ruff` or `ruff format` already gates — line length, import order,
  quote style, trailing commas. Reporting these is noise.
- Prose quality in documentation. You judge code and the comments inside it.
- **Anything the project has already decided.** Before filing a structural defect —
  a wide signature, a large class, a repeated shape — check whether the code says it
  is deliberate: a `provisional` marker in the docstring, a `docs/adr/` entry naming
  it, or a module docstring explaining the design. A decision already taken and
  recorded is not your finding, and reporting it as one buries the findings that are.
  If the decision exists but the code does not point at it, THAT is the defect, and
  its kind is `docstring-contradicts-code`.

  This has happened: a structural pattern that looked like the loudest finding
  available turned out to be named in an ADR as planned work, with a documented facade
  already standing over it. Check first; report either way on what you find.

## The comment rubric: whose chair is it written from

A comment is judged by **whom it serves**, never by its length. A long comment that
tells the reader something true and necessary is good. Verbosity is not a defect.

On a clear-cut comment this test agrees with ordinary good judgement, and it is not
claimed to beat it. What it adds is consistency across many files: the same closed
vocabulary, the same relocation precondition, and the same record shape, so that
independent reports merge into one account instead of twenty opinions. Apply it as
stated even where your own judgement would have arrived unaided — that sameness is
the point.

**Defects, at any length:**

| kind | what it looks like |
|---|---|
| `restates-code-comment` | `# loop over energies` above a loop over energies |
| `implementer-perspective-comment` | why we chose this, what we tried first, what it used to be, when it changed, an apology |

**Keep, at whatever length it needs:**

- what this thing **is**, and what a value **means**
- the invariant that holds here; what the reader must not assume
- the trap — the plausible edit that would silently break it
- the identifying source: the equation, paper section, or reference deck that says
  *this formula is that formula*

**When a comment is both.** A real comment is often both a trap and a history
note at once. Ask whether the comment prevents a specific, plausible edit. If a
competent reader who deleted it could make a silently-wrong change, it is a
trap: keep it. If it only explains why the current state was chosen, and no
plausible edit is prevented by knowing that, it is history: file it.

A comment that passes the trap test is kept, but trimmed to the part that does
the preventing. Measured numbers inside it follow the relocation rule
unchanged — secured in documentation first if they exist nowhere else.

Worked example: `# do not revert to CN here — it under-converges ~100%` is
kept, because reverting to CN is a plausible edit that the comment prevents.
The "~100%" figure is subject to the relocation rule below if it appears
nowhere else.

The line runs through the reader's chair. `# PRA 77 Eq. (37): the bra carries the
final channel energy` tells a reader what the line **is** — keep it. `# we use
MUMPS because it measured 72x faster on 2026-08-24` tells them how the code got
here — that belongs in documentation, not beside the code.

## Relocation precedes deletion

Before filing any rationale or measured number for removal, search `docs/physics/`,
`docs/adr/` and `CLAUDE.md` for that fact.

**Search for the bare numeral, not the number-plus-unit.** This repository writes
physics in Unicode — `µHa`, `Γ`, `σ`, `ω`, `×`, and U+2212 minus rather than hyphen —
so an ASCII-transliterated search term misses facts that are present. A judge searching
`154 uHa` found nothing while the document said `154 µHa`, and nearly filed a recorded
measurement as unrecorded. Search `154`, then read the hits. Search the whole of
`docs/physics/`, not only the file the comment names: a fact is often recorded under a
different note than the one the code points at.

- Recorded elsewhere → file the defect, and cite where the fact survives in `fix`.
- **Recorded nowhere else** → the kind is `provenance-at-risk`, not
  `implementer-perspective-comment`. The fix is to write the fact to documentation
  first; the deletion depends on it.

A measured number or citation that exists in exactly one comment is load-bearing.
Filing it for deletion without first securing it elsewhere is the worst outcome
this rubric permits, and the one it is written to prevent.

## Type machinery that costs more than it states

`untyped-public` catches a surface that says too little. This catches the
opposite: machinery that restates what the implementation already declares.

**The rule the project holds** (see `qscat-conventions`): one signature per
function, returning the `|` union of the shapes it can produce. `@overload`
stubs that differ only in a `Literal` flag are `redundant-overload` — each
restates a parameter list the implementation already states, so a signature
change has to be made in every copy and can drift in all but one. Twenty such
stubs across six methods once cost 221 lines to say what six implementation
signatures already said.

Two things this kind does NOT cover, so do not stretch it:

- **A union that is genuinely hard to read is not fixed by overloading it.**
  It is fixed by naming it — a type alias, or a dataclass when the shapes
  differ in more than arity. Report the unreadable union as `unclear-name`
  if it has no name and needs one.
- **Overloads that select on argument TYPE rather than on a literal flag**
  are not redundant; they express something a single signature cannot.

And its companion defect: with no overloads, a flag that changes the returned
shape is documented only in the docstring. A flag whose docstring does not say
which shape it selects is `docstring-contradicts-code` — the signature admits
several shapes and the prose names none.

## Verdicts and defects

**The report file is ONE object with exactly four top-level keys**, however many files
the unit contains:

```json
{"unit": "<the unit name you were given>",
 "files":   [{"file": "path/to/a.py", "verdict": "clean"}, ...],
 "defects": [ ... ],
 "held_up": [ ... ]}
```

`defects` and `held_up` are FLAT lists across the whole unit, each record carrying its
own `file` field — not nested per file, and the top level is never a list of per-file
objects. Six unit reports get merged programmatically, so a different shape breaks the
merge rather than degrading it. A judge has emitted a list of per-file objects and had
to be sent back with no fault in its judgements — the shape is not a formality.

Per file, one verdict from a closed set:

`clean` · `tighten` · `restructure` · `delete` · `promote` · `demote`

Per defect, one record:

```json
{"kind": "implementer-perspective-comment", "file": "path/to/file.py", "line": 204,
 "evidence": "# originally we tried the dense stepper here, too slow",
 "fix": "delete; the measurement is in docs/physics/nrm-time-dependent.md",
 "effort": "trivial"}
```

`kind` is closed: `restates-code-comment`, `implementer-perspective-comment`,
`provenance-at-risk`, `duplicate-logic`, `dead`, `commented-out-code`,
`misplaced-layer`, `oversized-unit`, `param-soup`, `unclear-name`,
`stringly-typed`, `docstring-contradicts-code`, `untyped-public`,
`redundant-overload`, `speculative-generality`.

`dead` is for a symbol nothing reaches; `commented-out-code` is for lines
disabled by commenting rather than deleted. There is deliberately no kind for a
magic number: in a repository whose modules carry physical constants, that would
report noise instead of findings.

`effort` is `trivial`, `contained`, or `invasive`.

`evidence` quotes the actual source line. A defect without a quote is not a defect.

## What held up

Every report ends with `held_up`: the things in this unit that are good, and why.
This is not politeness — without it, severity has no scale, and a report that finds
only faults is indistinguishable from a reviewer who was determined to find them.
