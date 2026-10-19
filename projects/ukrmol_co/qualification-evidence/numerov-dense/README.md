# Dense scaling public evidence

The [result](../../NUMEROV_DENSE_SCALING.md) passes166/166 checks / original
exits0/0:24 dense298-channel synthetic propagations, two potentials, three
interval counts, two arithmetic formats and two repetitions. No CO solve is
added; original interior false/exits1/1 and coverage132/594 remain.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/numerov-dense
```

The [artifact index](artifacts.json) pins original source at`1d8cfac`, results,
resource samples, prior coupled report and mpmath wheel. Separate support at
`ff2d258` adds one independent verifier without replacing frozen numerical
modules. The separate replay contract names the exact numerical worker.

Unpack`source.tar.gz` into`source/`, verification support into that source tree,
`results.tar.gz` into`collected/`, and the wheel into`dependencies/`. Keep
`prior-report.json` at the fetch root, verify the hashes and from`source/` run:

```sh
PYTHONPATH=../dependencies:. python -m projects.ukrmol_co.numerov_dense_replay .. --output ../verified-decisions.json
```

The saved verifier recomputes accuracy/invariant decisions, independent discrete
references, solve-role order, returned dtypes, disjoint timer partitions and
affine fits. Original logs omit per-call array shapes; dense298 dimensions follow
the fixed numerical source, not independent shape measurements.

The numerical worker requires read-only`/public`, empty writable`/outputs`,
single-threaded BLAS and a host with at least64 longdouble significand bits.
It reruns all24 propagations and six oracle pairs, requiring exact decimal and
binary64 K components and166 decisions. Timing fits/hypotheses are fresh
observations, separately retained rather than required to reproduce elapsed time.
See [verification](../../verification-numerov-dense-public-20261019.json).

Complete replay passes exits0/0: every decimal and binary64 K component, all166
decisions and six oracle pairs reproduce exactly. Timing hypotheses also repeat
as observations. All14 final public payloads fresh-verify; these include the
8-file`public-replay.tar.gz`, replay summary, owner and sampled resources.
Local public tests pass85 with two wider-format skips, exercised remotely.
