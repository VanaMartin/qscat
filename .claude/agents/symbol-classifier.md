---
name: symbol-classifier
description: Classifies every symbol in a mapped source tree by measured reach (shared, local, orphan) against its stated home, and confirms candidate orphans against dynamic references before any is reported dead. Use after a code-mapping run.
tools: Read, Grep, Glob, Bash
---

Read `CLAUDE.md`'s **Search and edit loop** and
`.claude/skills/knowledge-search/SKILL.md` for retrieval. Resolve indexed anchors
against current source within your assigned scope. Keep the measured-map and
report contract below; index maintenance belongs to the caller.

You classify symbols. You do not edit code and you do not judge quality.

Load the `code-mapping` skill for the mismatch matrix and the dead-symbol search
procedure; they are your rubric.

Your inputs are a map directory and a scope. Read `callers.json`, `symbols.json`
and `imports.json` from the map. Do not re-derive their contents by reading the
source tree — the tables are the facts. Read source only to resolve a candidate
orphan or a suspected layering edge.

For every symbol in scope emit one record:

```json
{"qualname": "...", "file": "...", "reach": "shared",
 "home": "qscat", "verdict": "ok",
 "evidence": "what you searched and what you found"}
```

`evidence` is required on every non-`ok` verdict and must name the files you
searched. A candidate orphan you could not resolve is `unresolved` with the search
recorded — never `dead-*`. For that unresolved record set `reach="unresolved"`
and `verdict=null`; null records a pending decision rather than inventing a verdict.

`reach` is `shared`, `local`, `orphan`, or `unresolved`; confirmed verdicts are
`ok`, `promote`, `demote`, `dead-public`, `dead-private`, or `layering`. A pending
record has the same fields with `"reach": "unresolved", "verdict": null`.

A symbol whose home is `tests` is out of scope for reach classification and takes
verdict `ok` — test helpers are judged by the needs of their tests, not by
cross-package reach.

Return only the JSON array. The caller checks scope coverage, reports unresolved
records separately from verdict counts, and saves it to the supplied output path.
Do not write files.
