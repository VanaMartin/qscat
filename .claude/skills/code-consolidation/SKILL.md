---
name: code-consolidation
description: Use when ruling on whether code that looks alike is the same thing twice — clone clusters, functions sharing a name, result-holder classes with overlapping fields — and what should replace them.
---

# Code Consolidation

## Overview

You rule on similarity clusters. Your evidence is the cluster tables and the source
they point at. You do not receive, and must not ask for, anyone else's opinion
about the code.

## The three rulings

`unify` — these are the same thing and one should replace the rest. Requires a
target name, a target home, and **a statement of the behavioural difference between
the members and its explanation**.

`keep-separate` — these look alike and are not the same thing. Requires the reason.
This ruling is permanent: the cluster will resurface on every future run, and the
reason is what stops it being re-litigated.

**`keep-separate` is not the ruling that needs less work.** Every ruling carries a
non-empty `difference`, this one included — a `null` there is a schema violation, not
a lesser verdict. This has been got wrong at scale: four of thirty `keep-separate`
rulings once shipped with `difference: null`, and every one was a cluster the next
analyst would have had to re-derive from scratch.

`investigate` — the members look alike and you cannot explain the difference.

## The hard rule

**No cluster is marked `unify` until the difference between its members is stated
and explained.**

If two functions are nearly identical and you cannot say why they differ, the
ruling is `investigate` — never `unify`. An unexplained difference between
near-identical functions is worth more than the merge: it is usually a bug in one
of them, or an undocumented special case that nobody has written down.

## Before ruling: two checks that change the answer

**A cluster table is a lower bound, so ask whether several clusters are one family.**
Near-clone grouping is anchor-relative: each group is built from one anchor compared
against later candidates, so a chain of gradual drift splits across groups, and a
member whose body is a little shorter falls below the threshold and vanishes from the
family entirely. Two or more clusters in the same file, or over the same subject, are
therefore often fragments of one pattern.

This is not hypothetical: four separate near-clone clusters in one module were once
one nine-member facade, and four *further* methods with the identical shape were
missing from all four clusters because their bodies were shorter. Ruling on the
fragments would have produced four findings about a pattern that deserved one. Read
the file, not only the rows.

**Check whether the repository already knows.** Before ruling that a shape is a
defect, look for an existing decision covering it: a `provisional` marker in the
docstring, a `docs/adr/` entry naming it, or a module docstring explaining the design.
A pattern the project has already decided on is `keep-separate` with that decision
cited — not a finding. This has happened: a repeated shape and the facade standing over
it were both already recorded as deliberate, one of them in an ADR. Look before ruling.

## Same name is not same thing

A homonym cluster has four common explanations, and they get different rulings:

| pattern | ruling |
|---|---|
| a protocol method and its implementations | `keep-separate` — that is what a protocol is |
| a thin re-export or shim over one real definition | `unify` — collapse to the real one |
| the same logic written twice | `unify` |
| unrelated code that happens to share a word | `keep-separate` |

## Holders

For result-holding classes, propose `unify` only when **three or more** repeat the
same method body. Two is coincidence.

The preferred remedy is a shared mixin or base carrying the repeated method, not a
merged class. Holders differ in their fields for real reasons, and merging them to
share a `save` trades a small duplication for a large false abstraction.

## Output

```json
{"cluster": 7, "ruling": "unify",
 "members": ["a.py:212", "b.py:88"],
 "difference": "b.py:88 clamps the lower bound to zero; a.py:212 does not",
 "proposal": "keep a.py:212 as qscat.core.x.f, add an optional clamp argument"}
```

`difference` is required on `unify` and `investigate`. On `keep-separate` it holds
the reason they stay apart.
