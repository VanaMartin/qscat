---
name: code-mapping
description: Use when a repository's structure must be known as measured fact — before a structure audit, a consolidation pass, a promotion decision, or any question of the form "what is stored where, who uses it, and what is duplicated".
---

# Code Mapping

## Overview

A deterministic AST scan that answers structural questions with measurements
instead of impressions. Agents read its tables; they never re-derive facts by
reading the tree.

## Running it

```bash
python .claude/skills/code-mapping/scripts/inventory.py \
  --root libs/qscat --root apps --root projects --root validation \
  --root benchmarks --root tests \
  --out .superpowers/audit/<timestamp>/before \
  --repo-root .
```

Stdlib only, under a minute on this repository, byte-identical across runs. The
output directory is git-ignored and never committed.

**The root is `libs/qscat`, not `libs/qscat/qscat`.** That package's tests live at
`libs/qscat/tests/`, beside the package rather than inside it, so rooting at the
package drops 65 test files from the caller table and inflates the orphan count.
Measured 2026-08-27: 75 apparent orphans in `qscat.core` against 52 with the root
correct. Before trusting any orphan list, check that every package's tests are
actually under one of the roots you passed.

## The tables

| file | one row per | the question it answers |
|---|---|---|
| `symbols.json` | function, method, class | where is it, how big, does it have a docstring |
| `imports.json` | import edge | which package depends on which |
| `callers.json` | symbol | how many sites use it, and in which packages |
| `duplicates.json` | clone cluster | what is written twice, exactly or nearly |
| `homonyms.json` | repeated name | what shares a name without sharing a body |
| `holders.json` | dataclass and friends | which result holders overlap, which methods echo |
| `hotspots.json` | symbol and file | what to look at first |

`hotspots.json` scores **review priority, not quality**. A high score means "look
here first", never "this is bad" — a long, branchy function may be exactly right.
Each row names the signals that produced its score.

## Two tables are lower bounds

**Near-clone groups.** `duplicates.json`'s near-clone pass is greedy and
anchor-relative: each group is built from one anchor compared against later
candidates. A chain of gradual drift — A close to B, B close to C, A not close to C —
is split across groups rather than gathered into one. A near-clone group is a lower
bound on its family, so treat two groups in the same file or subject area as possibly
one family, and read the source before ruling.

**Fan-in.** `callers.json` is a static scan. This repository dispatches by strings read from
YAML (`apps/qscat-run`), maps names to classes in registries, and uses `getattr`.
A symbol reachable only that way reports zero sites. **Zero sites is a candidate,
never a verdict.**

## Classification: the mismatch matrix

Reach is measured from `callers.json`: `shared` when consumers span two or more
packages, `local` when exactly one, `orphan` when none outside its own file.

Home is where the symbol lives, read against this repository's lifecycle rules:
`libs/qscat` is validated reusable code, `projects/` is toy and experiment stages,
`validation/` is harnesses, `apps/` is the execution surface, `benchmarks/` is
measurement.

Findings are the off-diagonal cells:

| reach | home | verdict | meaning |
|---|---|---|---|
| shared | `projects/` | `promote` | lifecycle stage 5 candidate |
| shared | `validation/` | `promote` | or `layering` if a consumer is `projects/` |
| local | `libs/qscat` | `demote` | speculative generality |
| orphan | exported | `dead-public` | only after the search below |
| orphan | private | `dead-private` | only after the search below |

## Before calling anything dead

For every candidate orphan, search and record the result:

1. `apps/qscat-run/examples/**/*.yaml` and any config for the bare name as a string.
2. Registries and dispatch dictionaries — grep the bare name in string literals.
3. `__all__` lists and package `__init__.py` re-exports.
4. `getattr(` calls whose attribute argument is a variable.
5. **Returned but never named.** A result-holder class or an exception type is
   constructed inside the module that defines it and handed back to callers, who
   then read its attributes without ever naming the type. The scan skips
   same-file references, so such a type reports zero consumers while being on
   every call path through its module. Check what the defining module's public
   functions RETURN before ruling on any class.

Expect this to be the common case, not the exception. Measured once on this repository:
89 of 421 exported non-test symbols reported zero static consumers, about a quarter of
them result-holder and exception classes of exactly this kind, and **none** proved dead.
A pass that reports most of its zero-consumer candidates as dead has not done this
search.

A candidate that survives all five is `dead-public` or `dead-private`. A candidate
you cannot resolve is `unresolved` — report it, never delete it.
