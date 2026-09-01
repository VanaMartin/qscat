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

Objects are addressed **by content**. The sha256 goes into the filename:

```
https://data.qscat.org/<experiment>/cross_section.830cffb8a044.csv
```

The folder is a human label, not an identity — it keeps a URL readable enough
to paste into a paper, while the hash decides what the URL means. The
committed `artifacts.json` carries the full digest per file and the fetcher
derives the address from it, so there is no second copy of the address to
disagree with the first.

Two properties follow, neither of which is a policy anyone has to enforce:

1. **A URL means one thing, by construction.** Different content hashes
   differently and lands elsewhere, so a link in a note cannot quietly come to
   mean something else, and the one-year `immutable` cache header is simply
   true. Re-publishing is idempotent rather than forbidden.
2. **A reproducible re-run costs nothing.** Identical bytes produce an
   identical key, so publishing the same result twice stores it once.

`index.json` is the deliberate exception: it is addressed by name, because it
is what you look up *before* you know any digest, and it carries a short cache
lifetime because it may legitimately be replaced.

Retention is **keep everything**, for now. Content-addressed blobs are shared
between pointers, so deleting by age is unsafe once two runs reference the
same bytes; the correct mechanism is reachability — sweep for blobs no pointer
on `main` names — which is what git does for the same reason. At roughly a
megabyte per sweep with duplicates collapsing, that is not yet worth building.

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

## Why not address by commit

The first version of this did, and it failed twice over.

**A branch commit is not a stable address.** Rebasing rewrites it, and the
manifests were orphaned three separate times — once by a rebase onto main,
once by the rebase-merge that landed the branch, and once because the
corrective commit never reached main at all. Each fix restored a SHA that the
next rewrite invalidated.

**It stored the same bytes repeatedly.** Publishing one experiment from a
branch and then from main left **751 KB of 3009 KB — 25% — as exact
duplicates**, because the commit had changed and the content had not. Since
these sweeps reproduce bit-identically (measured three times, across a host
restart, three images and three commits), that is the common case, not the
exception.

Content addressing removes both. What it gives up is the URL no longer naming
the code that produced the file — a real loss, but a small one against a
signal that was demonstrably unreliable. `git_sha` stays in the manifest as a
**record**: it answers "what produced these bytes", while `git blame` on the
pointer answers "when did this repository start citing them". Different
questions, and only the second is well served by history.

One consequence, recorded because it reverses an earlier decision here: an
undeterminable `git_sha` **warns** rather than failing the run. The hard
failure was justified while the address depended on the SHA; once it does not,
a local figure should not fail over forty bytes of metadata. Publishing is
where provenance matters, and the publisher still refuses a manifest that
cannot name its commit.
