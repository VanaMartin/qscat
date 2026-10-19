# Fixed-grid extended-arithmetic public evidence

The [result](../../NUMEROV_ARITHMETIC_CONTROL.md) records35 control and10
qualification checks passing, original exits0/0. The original interior analytic
prerequisite remains false/exits1/1. The precision-gain hypothesis passes all six
coordinates; this qualifies only the declared free constant-potential problem.

Fetch and verify content-addressed payloads:

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/numerov-arithmetic
```

The [artifact index](artifacts.json) pins source/results archives, actual prior
report, dependency wheel, independent saved-decision verifier, numerical worker
and separately frozen replay contract. Sources are frozen at`44ed195`, separate
verification support at`a372e8b`; support adds no replacement numerical module.

For saved decisions, unpack`source.tar.gz` to`source/`,
`verification-support.tar.gz` into that same source tree,`results.tar.gz` to
`collected/`, and the mpmath wheel to`dependencies/`. Keep`prior-report.json` at
the fetch root, verify all manifest hashes, then run from`source/`:

```sh
PYTHONPATH=../dependencies:. python -m projects.ukrmol_co.numerov_arithmetic_replay .. --output ../verified-decisions.json
```

The public numerical worker expects`/public` as its read-only root and an empty
writable`/outputs`. It runs six tiny short-grid oracle tests and14 full campaign
reexecutions. A host with at least64 longdouble significand bits is required;
unsupported precision stops unset before the campaign. Container provenance,
limits, original exits and sampled resources accompany the replay.

Complete public numerical replay passes exits0/0: all12 typed decimal matrices,
two original baselines and45 decisions reproduce exactly. The14 final payloads
include`public-replay.tar.gz`, original owner, resource/execution record and
summary. Six additional short-grid preflight tests pass in both dtypes; remote
longdouble has64 significand bits. Local public tests pass74 with one
host-dependent wider-precision test skipped. See the
[verification record](../../verification-numerov-arithmetic-public-20261019.json).

No new CO/native/inner solve or unique energy is added by this evidence.
