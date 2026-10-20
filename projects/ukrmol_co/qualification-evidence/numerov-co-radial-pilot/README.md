# Long-range saved-CO radial pilot evidence

[Design and result](../../NUMEROV_CO_RADIAL_PILOT.md): **131/131 producer checks**
over five full propagations /141,200 intervals /282,400 mixed solves, original
command/container/owner exits **0/0/0**. Independent saved-byte reconstruction
reproduces all131 decisions,12 full captured residual matrices and36 paired
80/110-digit residual entries, verifier exits **0/0**. All26 remote tests pass.
Finest200/400 native64 K discrepancies1.56617e-7 /1.81307e-7 pass their separate
1e-6 gates; S, phase, radial and invariant checks also pass. This reuses one
previously corner-qualified stretched B1 identity, adding no unique-energy coverage.

The full400-bohr propagation completes and its scientific result reproduces
exactly, but the original wrapper fails **command/container/owner1/1/1** on a
tuple/list coordinate-schema mismatch. A separate saved-only diagnosis reproduces
that bug and verifies the complete fresh record and all three full captured
A/B/X sets, **exits0/0/owner0, zero new propagations**. All **26 published payloads**
fresh-download-verify locally. [Final verification](../../verification-numerov-co-radial-schema-20261020.json)
and [dated result/status](../../progress-numerov-co-radial-replay-result-20261020T1704.json)
retain both outcomes. The earlier15:03UTC snapshots remain historical launch-state
records. The replay owner authenticated the initial10-payload publication.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/numerov-co-radial-pilot
```

Extract `source.tar.gz` and `verification-support.tar.gz` into `source/`,
`results.tar.gz` into `collected/`, `raw.tar.gz` into `raw/`, `native-wheel.whl`
into `extension/`, and the mpmath wheel into `dependencies/`. Keep
`collection-sha256.json` and `replay-contract.json` in the extraction root.
Authenticate each artifact against `artifacts.json` first. Producer source stays
at `37c503f`, separate verification support at `42bb02f`; the support archive
contains three additional files, leaving all107 original sources unchanged.
Verify the107 source hashes, three support hashes,24 raw hashes, two prior files,
wheel and its two components, and the producer's original42-file collection.
`manifest.json` pins support hashes; `source-manifest.json` pins original sources.

From `source/`, on qualified little-endian x86-64 Linux with the pinned runtime:

```sh
PYTHONPATH=../extension:../dependencies:. python -m projects.ukrmol_co.numerov_co_radial_replay --root .. --contract ../replay-contract.json --output ../saved-verification-fresh.json
```

The verifier checks exact identity/energy/header/EOF and independently expresses
the decisions using right-solve Cayley and eigenvalue-arctangent phase. It checks
all12 captured triples against actual sampled A/B/X hashes, reevaluates their full
native/NumPy residual matrices, and independently sums36 paired80/110-digit
entries. Reported96 exact sampled residual evaluations and48 sampled solve gates
are preserved observer records; only the12 captured final residuals can be
independently reconstructed from saved operands. Every propagation solve retains
its original32-extended-epsilon acceptance gate.

## Full400-bohr replay — scientific result exact, original wrapper fails

The [separate contract](../../numerov-co-radial-pilot-public-replay-contract.json)
freezes exactly one complete18→400-bohr native reexecution,61,120 intervals /
122,240 mixed solves. Its original command adds `--replay-output ../reexecuted`.
The complete fresh scientific record, canonical extended K hash, all sampled
operand hashes and every entry of all three captured A/B/X triples must be exact.
Only positive finite auxiliary rcond observations have the prior32-binary64-epsilon
relative budget. Timings and NPZ storage-byte hashes are excluded from numerical
identity; each saved payload's own digest is authenticated independently.

It finishes **16:59:02UTC**,7164.645088s /1.99018h wall /.17221GiB,35,760 resource
samples. The complete fresh row passes22/22 checks. Original scientific K,
extended K hash, complete diagnostics, sampled operand hashes, audit aggregates
and paired residual references all match exactly. The failure is in Python
comparison: the unchanged oracle helper emits nine tuple coordinate pairs,
while archived JSON decodes those same coordinates as lists. Original wrapper
source stays `42bb02f`. Its `verification.json` is the earlier saved-audit
checkpoint, with a true saved verdict and null numerical verdict, not a final
wrapper success. Original command/container/owner exits remain1/1/1.

`failed-public-replay.tar.gz` contains all30 original failure files, pinned by
`failure-collection-sha256.json`, including the full fresh row, three captures,
logs, exits and sampled resources. `saved-tests-and-launch.tar.gz` retains the
earlier20-file snapshot and original downloaded index. `public-owner.py` records
the unchanged detached commands, resource limits and stop rule.

## Saved-output schema diagnosis — no repeated propagation

The [separate diagnostic contract](../../numerov-co-radial-schema-diagnostic-contract.json)
and corrected support freeze at `ca9a5d0` before execution. To reproduce, extract
the failed replay archive into the extraction root, then overlay
`schema-support.tar.gz` into `source/`. Verify the three new support hashes in
`schema-support-manifest.json`; the107 producer sources remain unchanged. The
overlay corrects coordinate-container normalization and adds a saved-output mode.

From `source/`, with the same qualified runtime and dependencies:

```sh
PYTHONPATH=../extension:../dependencies:. python -m projects.ukrmol_co.numerov_co_radial_replay --root .. --contract ../replay-contract.json --output ../schema-diagnosis-fresh.json --existing-replay-output ../reexecuted
```

This reexecutes **no propagation**. It reconstructs131 original and22 fresh row
decisions, calls the original oracle helper on the three saved fresh captures to
reproduce the old container mismatch, normalizes coordinate containers only, and
checks every entry of the three original/fresh full A/B/X sets exactly. Auxiliary
rcond drift6.77226e-16 is within the unchanged7.10543e-15 relative budget; there is
no scientific tolerance change. Original failed exits remain separately checked.

The saved diagnostic takes **12.287171s /.07103GiB**,60 samples, exits0/0/owner0;
five focused regression tests pass locally and remotely. Its16-file collection
is `schema-diagnosis.tar.gz`, pinned by `schema-collection-sha256.json`; support,
contract and `schema-diagnosis.json` are separate public payloads. The original
110 effective sources and111 diagnostic sources,242 read-only engine artifacts
and28 raw/prior/component pins verify before/after. There is no second full-row
replay after the fix. No scientific container remains at17:04UTC.

The original producer ran4.31638h /.27837GiB, with before/after107-source,
242-engine-artifact and28 raw/prior/component checks. Timing comparisons include
unequal audit work and host contention; they are not isolated speed benchmarks.
The maximum recorded all-solve backward error3.24329e-18 is close to its
3.46945e-18 gate; neither it nor minimum sampled rcond9.45956e-4 is a global
condition/forward-error bound. Other CO identities, native298-channel800-bohr
execution, full interior coverage and the electronic/physical model remain open.
Original interior false/exits1/1 and historical132/594 coverage are retained.
