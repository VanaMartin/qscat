# Long-range saved-CO radial pilot evidence

[Design and result](../../NUMEROV_CO_RADIAL_PILOT.md): **131/131 producer checks**
over five full propagations /141,200 intervals /282,400 mixed solves, original
command/container/owner exits **0/0/0**. Independent saved-byte reconstruction
reproduces all131 decisions,12 full captured residual matrices and36 paired
80/110-digit residual entries, verifier exits **0/0**. All26 remote tests pass.
Finest200/400 native64 K discrepancies1.56617e-7 /1.81307e-7 pass their separate
1e-6 gates; S, phase, radial and invariant checks also pass. This reuses one
previously corner-qualified stretched B1 identity, adding no unique-energy coverage.

All **16 published payloads** fresh-download-verify locally; the replay owner
separately authenticates the original10-payload publication before execution.
[Verification summary](../../verification-numerov-co-radial-pilot-public-20261020.json)
and [dated result/status](../../progress-numerov-co-radial-pilot-result-20261020.json)
retain this distinction and leave the final numerical replay verdict/exits null.

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

## Full400-bohr numerical replay — running

The [separate contract](../../numerov-co-radial-pilot-public-replay-contract.json)
freezes exactly one complete18→400-bohr native reexecution,61,120 intervals /
122,240 mixed solves. Add `--replay-output ../reexecuted-fresh` to the command.
The complete fresh scientific record, canonical extended K hash, all sampled
operand hashes and every entry of all three captured A/B/X triples must be exact.
Only positive finite auxiliary rcond observations have the prior32-binary64-epsilon
relative budget. Timings and NPZ storage-byte hashes are excluded from numerical
identity; each saved payload's own digest is authenticated independently.

**At15:03UTC the replay remains running**, after saved verification passes in
8.03s /.113GiB. Its runtime starts14:59:37UTC; expected7800s /upper10800s /
hard cap14400s, cores4–11 /one actual single-threaded worker /4GiB. Final
command/container/owner exits and numerical verdict remain pending.
`saved-tests-and-launch.tar.gz` preserves the completed tests/audit, owner and
bootstrap, original downloaded10-payload index, runtime launch and dated progress
snapshot. `saved-tests-launch-sha256.json` pins those20 snapshot files. Its running
snapshot is not a final runtime/owner collection. `public-owner.py` records the
exact detached commands, hashes, resource sampling and stop rule.

The original producer ran4.31638h /.27837GiB, with before/after107-source,
242-engine-artifact and28 raw/prior/component checks. Timing comparisons include
unequal audit work and host contention; they are not isolated speed benchmarks.
The maximum recorded all-solve backward error3.24329e-18 is close to its
3.46945e-18 gate; neither it nor minimum sampled rcond9.45956e-4 is a global
condition/forward-error bound. Other CO identities, native298-channel800-bohr
execution, full interior coverage and the electronic/physical model remain open.
Original interior false/exits1/1 and historical132/594 coverage are retained.
