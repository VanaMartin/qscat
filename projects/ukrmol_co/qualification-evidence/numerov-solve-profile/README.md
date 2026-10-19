# Mixed-solve profile and operand evidence

The [result](../../NUMEROV_SOLVE_PROFILE.md) passes76/76 checks/exits0/0:
eight exact reexecutions of previously qualified dense32-interval rows.
Extended recurring residual evaluation accounts for94.6% of solve time.
Six captured A/B/X triples support later kernel differential testing.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/numerov-solve-profile
```

The [index](artifacts.json) pins source at`bbbc3d8`, original results, six NPZ
operand captures, source/engine/resource checks, actual prior report and mpmath
wheel. Separate support at`f396ed3` adds one independent saved verifier without
replacing numerical modules.

Unpack`source.tar.gz` into`source/`, support into that source tree,
`results.tar.gz` into`collected/` and the wheel into`dependencies/`. Keep
`prior-report.json` at the fetch root, verify all hashes and from`source/` run:

```sh
PYTHONPATH=../dependencies:. python -m projects.ukrmol_co.numerov_solve_profile_replay .. --output ../verified-decisions.json
```

Saved verification reconstructs all76 numerical/instrumentation decisions,
392 solve-call layouts, source-statement partitions, summaries and hypotheses;
it audits capture hashes/metadata. Full numerical replay requires a wider-format
host, single-threaded BLAS, read-only`/public` and empty writable`/outputs`.
It reruns all eight propagations, compares all K/diagnostics/76 decisions
exactly and checks loaded A/B/X numeric identity/dtype and32-eps backward errors
for six captures. Timings/hypotheses remain fresh observations, not exact replay
gates. Capture padding bytes are not portable numeric identity.

See [public verification](../../verification-numerov-solve-profile-public-20261019.json).

Complete replay passes exits0/0: all eight K/diagnostic records and76 decisions
reproduce exactly. Six captured A/B/X triples reproduce loaded values/dtypes
and pass recomputed backward-error gates. All14 final public payloads verify;
these include the14-file replay collection, summary, owner and resource record.
Local public tests pass89 with two wider-format skips, exercised remotely.
