# Coupled precision and synthetic dense-cost public evidence

The [result](../../NUMEROV_COUPLED_PRECISION.md) records89/89 passing checks,
32 typed propagations and original exits0/0. Original interior false/exits1/1
and CO coverage are retained. The finite synthetic dense probe identifies a
large arithmetic cost; it does not authorize a CO sweep.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/numerov-coupled
```

The [artifact index](artifacts.json) pins original results/resources, frozen
source at`4a9baa1`, actual prior arithmetic report, mpmath wheel, independent
saved verifier and replay/profiling contract. Separate support at`d4662d3` adds
one verifier without replacing frozen numerical modules.

Unpack`source.tar.gz` to`source/`, verification support into that source tree,
`results.tar.gz` to`collected/`, and the wheel to`dependencies/`. Retain
`prior-report.json` at the fetch root, verify all hashes, and from`source/` run:

```sh
PYTHONPATH=../dependencies:. python -m projects.ukrmol_co.numerov_coupled_replay .. --output ../verified-decisions.json
```

The saved verifier independently recomputes propagation arithmetic, reference
and hypothesis decisions; stored preflight diagnostics are explicitly audited
rather than claimed to reconstruct missing algebra output matrices. The numerical
worker reexecutes all32 propagations,9 discrete oracle pairs,6 continuous
integrations and actual preflight algebra/tails. Only its two existing298-channel
rows are instrumented with cProfile; no extra propagation or input change.
The worker expects read-only`/public` and empty writable`/outputs`, single-threaded
BLAS and a host with at least64 longdouble significand bits. The replay owner
records its exact pinned image and resources.

Complete replay passes exits0/0: all32 decimal K matrices,89 decisions,
references and preflights reproduce exactly. The14 final payloads include
`public-replay.tar.gz` (10 files, including`dense-binary64.prof` and
`dense-longdouble.prof`), owner, summary and sampled execution/resources.
The extended dense profile attributes61.92% to three Taylor changes and26.11%
to ten mixed solves. Setup is substantial; an eight-interval average is not a
CO per-interval cost measurement. Local public tests pass80 with two
wider-format skips. See [verification](../../verification-numerov-coupled-public-20261019.json).
