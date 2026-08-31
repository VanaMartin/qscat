# 8. Computed artifacts live in public object storage, not in git

Date: 2026-08-31

## Status

Accepted.

## Context

The repository had begun carrying computed results as if they were source.
Measured at `main`:

| | size |
|---|---|
| `.git`, packed | 39 MB |
| all history, uncompressed blobs | 64.4 MB over 3080 blobs |
| blobs under `docs/physics/figures` + `validation` | **23.7 MB over 499 blobs** — 37% of every byte the repository has ever stored |
| tracked at HEAD | ~18 MB, of which PNG 6.2 MB, CSV 1.7 MB, JSON 0.7 MB |

Small in absolute terms; the trajectory is the problem. Git cannot diff,
dedupe or partially clone any of it, every re-render of a figure adds its full
size to history forever, and a single spin-orbit study added three CSVs of
3343 rows each. One campaign output, `validation/coupled/results/screen.json`,
was 38 551 lines committed alongside the 365-line script that wrote it.

These files are not source and they are not data a human reads. They are a
function of `(config, commit, container image)` whose only value is the time
they save: hours of MUMPS solves, or a baseline to compare a new run against.

Some of them are nonetheless load-bearing. Tests read `o2-fit-report.json` and
Houfek's `CSVE.V00.J00`; physics notes cite committed figures. Removing them
indiscriminately breaks the gate.

## Decision

Classify by **what reads the file**, not by size or file type.

| class | example | home |
|---|---|---|
| golden / oracle inputs | `validation/n2/data/CSVE.V00.J00` | **git** — hand-verified, tests depend on them, changing one is a physics decision that must be reviewed |
| fit reports, locked constants, gate files | `results/*-fit-report.json` | **git** — KB-sized, lock tests read them |
| figures a note discusses | `docs/physics/figures/*.png` | **git** — the notes read like papers and are unreadable without them |
| sweep results, run outputs | `results/o2-ve/cross_section.csv` | **object storage** — reproducible from a committed config; nothing in the fast tier reads them |

The line: *the repository keeps what a test or a note needs in order to stand
alone, at KB scale; anything that is a function of (config, commit, image) at
MB scale moves out, with the config and the commit recorded.*

Storage is a **Cloudflare R2 bucket published read-only at
`https://data.qscat.org`**. Reads are anonymous HTTPS — no account, no
credentials, no client library, no vendor SDK. Writes are impossible over that
hostname: publishing goes through the S3 API with a per-maintainer token from
the private `qscat-infra` repository. **Public read-only is a property of the
architecture, not of a secret we keep.**

Addressing is by path, and the path carries the binding:

```
https://data.qscat.org/<scope>/<sha7>/<experiment>/<file>
                       main | branch | tag
```

Two rules the publisher enforces rather than documents:

1. **Published paths are immutable**; re-publishing a key is refused and there
   is no `--force`. Corrected numbers get a new commit and a new path, and the
   old URL keeps the old values — which is what a citation needs, and what
   makes a one-year `immutable` cache header safe.
2. **Every object's sha256 is recorded** in the committed `artifacts.json`
   pointer before upload, so a reader verifies what they downloaded.

`branch/` expires after 90 days; `main/` and `tag/` are permanent, because an
artifact a note cites must outlive the branch that produced it.

## Consequences

`qscat-run fetch DIR` downloads and verifies. A file already present and
correct is skipped, so re-running is free and an interrupted fetch resumes.

CI must never depend on the store: the fast tier reads nothing that is
fetched. A network outage may not turn into a red build.

The repository still stands alone in the sense that matters. The *inputs* a
claim depends on stay in git, and every published directory carries its
`config.resolved.yaml`, so a clone with no network can regenerate any fetched
artifact — slowly. What moved out is only the expensive, reproducible output.

Existing history is **not** rewritten. The 23.7 MB stays where it is: SHAs are
cited in notes and pull requests, and rewriting them to reclaim ~20 MB is a
bad trade. The rule changes going forward only.

A **content-addressed index (DVC or equivalent) was deliberately deferred.**
Its advantage over a naming convention is dedup and per-branch push, neither
of which one maintainer and twenty result files needs yet. Recording sha256
from the first commit is what keeps that door open: adopting DVC later becomes
a script over the existing pointers rather than a re-upload. Revisit when a
second person publishes concurrently, or when the same artifact is visibly
duplicated across many SHA prefixes.

## Known weakness

A pointer records a **commit** SHA, and a commit SHA does not survive history
rewriting — a rebase or a squash-merge orphans it. This bit us immediately:
the first published O₂ manifests recorded a commit that a routine rebase then
rewrote away. Until pointers also record the **tree** hash — which is
content-addressed and survives re-parenting — a branch whose artifacts are
already published must be merged with a merge commit or a fast-forward
rebase-merge, never squashed.
