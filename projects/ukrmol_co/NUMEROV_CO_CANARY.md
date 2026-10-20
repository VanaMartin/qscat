# Saved-CO recurrence canary

## Decision and finite design

The [native integration](NUMEROV_NATIVE_INTEGRATION.md) qualifies synthetic
accuracy and about 2.26× marginal interval improvement. Actual 298-channel CO
matrices can have different conditioning, correction counts and throughput.
The [prospective contract](numerov-co-canary-contract.json) tests those differences
using existing saved inputs and the unchanged propagator/native wheel.

All six identities are the previously selected interior high-|K| peaks:
compressed B1/B2 index87, equilibrium B1/B2 index53, stretched B1/B2 index18.
Their exact archived Rydberg bytes remain authoritative. Each gets reference
and audited-native N8 prefixes. Three cases (compressed B1, equilibrium B2,
stretched B1) get N32/128 timing pairs, two opposite-order repetitions. This is
**36 prefix executions / 3,996 mixed solves / 1,980 recurrence updates**. The
nominal endpoint-interval sum is 2,016; no matching solve is executed.

The wrapper lets the first Taylor call run unchanged, then stops at the second
call only after the expected `2N-1` solves. This captures exactly the original
initialization and N-1 recurring updates. Final recurrence radii are
18.0109375 / 18.0484375 / 18.1984375 bohr for N8/32/128. There is no tail-based
short-radius K interpretation. Both hooks restore on success or failure.

Accuracy pairs must reproduce all canonical A/B/X hashes, binary64 LU one-norm
reciprocal-condition estimates, final-delta hash and stats exactly. Each of the
90 native accuracy solves independently audits every residual/decision against
NumPy. Twelve actual accepted A/B/X triples (initialization and last-defect per
case) are retained with three residual entries independently checked at paired
80/110 digits. This qualifies residual arithmetic and differential behavior,
not a full forward-error bound. Reciprocal-condition estimates diagnose the
binary64 preconditioner; they do not bound the extended solve or later radii.

Timing rows omit residual observers, condition estimation, per-solve hashing and
captures. The final hash is computed outside the timer. Marginal N32→128 and
direct N128 recurring rates are observations. Near-boundary rates cannot release
a full 18→400-bohr or 60-solve calculation: long-range conditioning remains open.

## Resources and stopping rule

Run on Sadaharu in a fresh persistent `/home` directory: reserved physical cores
0–3, one worker, single-threaded BLAS, 4 GiB. Forecast **600 s**, upper estimate
**1,200 s**, hard cap **1,800 s**. Prior synthetic rates imply about 375 s of
timing-prefix work plus accuracy, capture and independent-reference overhead.
Record existing unrelated host workloads; do not modify them. Compare numerical
records independently of timing fluctuations.

Preflight pins the source archive, original passing integration evidence,
selection, all 24 raw inputs, whole wheel/two components and 242 installed engine
artifacts. Tests run before scientific work. A false accuracy pair/oracle stops
before further work; an exception remains unset with original exits retained.
There is no retry, tolerance change or automatic follow-on. Original interior
false/exits1/1 and historical 132/594 coverage remain. No UKRmol invocation or
new inner solve is required.

## Reproduction

On the qualified Linux x86-64 scalar format, with the already published native
wheel and mpmath on `PYTHONPATH`:

```sh
python -m projects.ukrmol_co.numerov_co_canary \
  --contract projects/ukrmol_co/numerov-co-canary-contract.json \
  --raw RAW --output NEW_OUTPUT_DIRECTORY
```

RAW is the checksum-pinned six-case corner raw directory. The output must be new;
the report checkpoints after each row and retains complete accepted operand
triples under `operands/`. Numerical source modules remain frozen at `2022311`.

## Independent verification

The separate `numerov_co_canary_replay` reader reconstructs the complete decision
map, exact energy/order/role coverage and timing arithmetic. It checks captured
triples against the actual solve hashes, recomputes full native/NumPy residual
matrices and independently repeats the paired high-precision residual references.
Saved residual counters remain observations until a fresh prefix is executed.

```sh
python -m projects.ukrmol_co.numerov_co_canary_replay \
  --report OUTPUT/report.json --operands OUTPUT/operands \
  --output verified.json
```

Adding `--raw RAW` authorizes exactly six audited native N8 prefixes (90 solves),
with canonical output/operand hashes, all stats and residual aggregates required
to reproduce exactly. It does not repeat timing rows. Reexecution retains its
own original exits and resources in a fresh directory; no numerical producer
module is replaced by verification support. Public evidence pins both sources.
