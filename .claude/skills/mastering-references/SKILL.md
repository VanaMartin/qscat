---
name: mastering-references
description: Use when adding a paper or thesis to reference/literature/, writing or updating a reference note for one, or citing a published result, equation, table or figure from qModeling code or docs.
allowed-tools: Read, Grep, Glob, Edit, Write, Bash, WebSearch
---

# Mastering references (qModeling)

The source PDFs in `reference/literature/` are **gitignored** — they do not
travel with a clone. The tracked `*.md` reference notes beside them are what the
repository actually cites, so a note has one job:

**Carry every published fact this repository relies on, each with a locator
precise enough to find it in the source in one jump.**

Not a summary of the paper. Not a review. The subset the repo depends on,
anchored.

## The locator rule

**Every extracted fact carries a locator, and every locator carries a page.**

A section name is not a locator — the reader still has to hunt. Write
`p. 032721-6, Eq. (32)`, not `Sec. III`. Journal articles with per-article
pagination (`032721-1` … `032721-17`) use that form; theses use plain page
numbers as printed.

| Fact | Locator |
|---|---|
| an equation | `p. 032721-4, Eq. (12)` |
| a parameter value | `p. 032721-9, Table I` |
| a numeric result or curve | `p. 032721-11, Fig. 7` |
| a stated conclusion | `p. 032721-13` |
| a definition or convention | `p. 5, Eq. (1.63)` (thesis) |

If the PDF's printed page differs from the extractor's page index, use the
**printed** page — that is what the reader sees — and say so once at the top if
the offset is confusing.

## Required structure

Write the note to `reference/literature/<same-stem-as-pdf>.md`. Every section
below is required; write "None." rather than dropping one.

```markdown
# <Author(s)>, <Journal> <vol>, <page> (<year>) — <short title>

**Source:** `reference/literature/<stem>.pdf` (gitignored) · <DOI or URL>
**Pagination:** <e.g. per-article, 032721-1 … 032721-17>

## Why this repository cares
One paragraph: what qModeling takes from this source, and what would be
unfounded without it.

## What this repository uses
A table — the load-bearing rows. One line per fact, each with its locator.

| Fact | Locator | Used by |
|---|---|---|

## Equations
The equations the repo implements, transcribed, each with its locator. Use the
paper's own numbering.

## Parameters and numeric values
Tables of published constants, with locators. State explicitly whether each was
checked against the repo, and the result.

## Findings and limits
The paper's own conclusions the repo relies on, and the limits it states on
them.

## Terminology map
Paper symbol → qModeling name, where they differ.

## Not used here
What the source contains that the repo deliberately does not take, so a later
reader does not assume coverage.
```

## What a note must not contain

- **Process narrative.** No "judgment calls I made", no TODOs, no "to verify".
  A note is a settled artifact; if something is unverified, say so as a fact in
  the relevant section ("not checked against the repo") rather than as a task.
- **Bulk verbatim text.** These are copyrighted. Transcribe equations, tables of
  constants, and short defining phrases; put everything else in your own words.
  A quoted sentence is fine where the exact wording is load-bearing — mark it as
  a quote and give its locator.
- **Facts without locators.** An unanchored claim is the thing this skill
  exists to prevent.

## Verify before you assert parity

When a note says a published value matches the repo, **check it** and say what
you checked:

```bash
grep -rn "918.076\|0.75102" libs/qscat/qscat/model/library.py
```

Write "matches `qscat.model.N2` exactly (verified)" or "differs: paper 918.076,
`H2P.mu` 918.25" — never "presumably matches". A disagreement between a paper
and the code is a finding, and has twice been a real one in this repo.

## Adding a new source

1. Copy the PDF to `reference/literature/<firstauthor>-<year>-<journal><vol>-<page>.pdf`
   (theses: `<firstauthor>-<year>-thesis.pdf`).
2. Extract text — see `references/extracting.md` for the command and the
   page-offset check.
3. Write the note per the structure above.
4. Add a row to `reference/literature/README.md`'s source table, and a short
   "What <source> covers" section pointing at the note.
5. `*.pdf` and `*.txt` are gitignored; the `*.md` note is tracked. Commit the
   note.

## When NOT to use this

- For qModeling's own physics write-ups — those are `docs/physics/`. A reference
  note describes a *source*; a physics note describes *our method*.
- For the eMoScat/libXcuda code snapshots in `reference/` — those are read by
  the `port-scout` agent, which extracts algorithms rather than citations.
