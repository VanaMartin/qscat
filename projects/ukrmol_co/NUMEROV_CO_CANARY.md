# Saved-CO recurrence canary

## Producer result — 20 October, 00:05 UTC

Frozen producer **`8520bbbad8879e450066b07aceb48e7b9e42b51d`** passes all
**303/303 checks / 36 prefixes / 3,996 mixed solves**. Original tests/runtime
exits are **0/0**; **17 remote tests pass**. Every reference/native final-delta
hash and diagnostic matches; all N8 A/B/X and condition-estimate records match.
The 90 audited native solves produce **180 full residual evaluations**, all
numerically identical to NumPy, with **zero accept/reject disagreements**.

At most **one correction** is needed. Maximum independently accepted backward
error is **5.11240321e-20**, below the unchanged 32-epsilon threshold. Minimum
sampled binary64 reciprocal one-norm condition estimate is **.0042606510**
(stretched B1 initialization). These early systems show no severe preconditioner
conditioning; that observation does not bound later radii or forward errors.
The twelve captured accepted operand triples provide **36 paired 80/110-digit
residual references**. Their maximum scaled arithmetic error is
**1.43844435e-19**, below 32 extended epsilons. Paired reference values agree.

| Case | Reference marginal s/interval | Native marginal s/interval | Marginal gain | Full N128-prefix gain |
|---|---:|---:|---:|---:|
| Compressed B1 | .279799469 | .124993730 | **2.23851×** | **2.14843×** |
| Equilibrium B2 | .278763455 | .124918698 | **2.23156×** | **2.14826×** |
| Stretched B1 | .277354928 | .123918143 | **2.23821×** | **2.15079×** |

Direct N128 native recurring rates are **.12296122 / .12258468 / .12190433 s
per interval**. Marginal rates exceed the prior synthetic varying-potential
native rate by only **3.42–4.32%**. These are **near-boundary prefix** rates;
the N128 prefix stops at **18.1984375 bohr**, without tail matching. Its roughly
2.15× total gain is not a complete exterior propagation gain. No long-range
runtime forecast or numerical-accuracy claim follows from this local sample.

Producer runtime finishes **00:05:41 UTC**, **468.927574 s wall / 468.456538
sampled CPU s / 192,290,816-byte kernel peak (.1791 GiB)**, 2,340 samples at .2 s.
This is below the 600/1200/1800-s expected/upper/hard budgets. It uses reserved
physical cores0–3, one worker, single-threaded BLAS and 4 GiB. Unrelated ComfyUI
work remains active and untouched; timings include that host context.
All **99 source hashes / 242 installed artifacts / 24 raw inputs / two wheel
components** verify before/after. The 25-file collection preserves twelve full
operand captures, checkpoints, original owner/commands/exits and resources in
`prepared/numerov-co-canary-20261020/`.

Separate verification support freezes at **`646257e`**. It reconstructs saved
decisions and replays exactly six audited native N8 prefixes from public bytes;
the [replay contract](numerov-co-canary-public-replay-contract.json) pins the
producer, full artifacts and 90/180/300-s resource/stopping rules.

### Preserved strict-replay failure and targeted diagnosis

The original public replay finishes **00:07:46 UTC / exits1/1 / owner1** after
14.439877 s / 13.834249 sampled CPU s / 100,003,840-byte peak. Its saved-decision
and captured-reference audit completes, but the first fresh native prefix fails
the original all-fields-exact comparison. No completed verification JSON is
produced. This remains a failed replay, not a passing result reclassified later.
Its seven-file archive and original logs/exits are retained publicly.

A separately frozen [one-prefix diagnostic](numerov-co-canary-replay-diagnostic-contract.json)
at `858b370` completes at **00:09:27 UTC / exits0/0**, 4.976144 s /
105,230,336-byte peak. It saves the fresh row, Q and first/last operand triples.
Only **eight auxiliary binary64 reciprocal-condition estimates** differ, by at
most **8.88e-16 absolute**. Every other numerical field, final delta, all actual
A/B/X hashes, solver diagnostics and residual decisions match. Both full captured
triples have **zero differing entries**. This isolates the discrepancy to the
auxiliary estimates; an elementary-operation cause has not been established.

The [separate successor supplement](numerov-co-canary-condition-replay-supplement.json)
at `5da71f3` keeps every scientific record exact and compares only reciprocal
condition observations with a **32-binary64-epsilon relative rounding budget**
(7.10543e-15), retaining positive/finite estimates. Original303 saved decisions,
all operand hashes and 32-extended-epsilon / three-correction solve gates remain.
The original producer and strict verifier are unchanged. This is a separate
explicitly scoped successor, with the original failure preserved.

### Completed public successor — 20 October, 00:13 UTC

The successor completes **00:12:25 UTC / exits0/0 / owner0** and independently
reconstructs all **303 saved decisions**, twelve full captured residuals and
36 paired high-precision reference entries. All **six audited native N8 prefixes /
90 solves** reproduce scientific records exactly, including actual A/B/X hashes,
final deltas, complete solver stats and residual aggregates. Maximum auxiliary
rcond relative drift is **4.55336e-16**, below **7.10543e-15**; all estimates
remain positive/finite. Timing rows are not repeated.

Runtime is **35.919531 s wall / 35.478702 sampled CPU s / 137,981,952-byte peak**,
178 samples, below 90/180/300-s budgets. All **103 effective sources / 242 installed
artifacts / 24 raw inputs / two wheel components** verify before/after. Originals
are in `prepared/numerov-co-canary-condition-successor-20261020/`.

All **19 final public payloads** and the final index fresh-download and verify,
including the failed original replay, diagnostic and exact-science successor.
The [00:13:25 UTC capture](progress-numerov-co-canary-result-20261020.json)
reverifies all **40 collected producer/failed-replay/successor hashes**, the
99/101/103-source trees and 242 installed artifacts. No scientific container
remains; unrelated ComfyUI stays active, with 57.92 GiB free. Public evidence:
[`qualification-evidence/numerov-co-canary/`](qualification-evidence/numerov-co-canary/README.md).

**Next:** cost a bounded long-range radial pilot that resolves actual CO
conditioning/accuracy and throughput beyond this near-boundary sample. Do not
extrapolate these prefixes into an absolute native K bound or automatically
release the full 60-solve interior scope. Original interior false/exits1/1,
14 corner-qualified identities, historical132/594 coverage and all model,
extraction, continuity and neutral-domain gates retain their prior scopes.

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
to reproduce exactly; the original strict command additionally requires auxiliary
condition estimates exact and preserves its observed failure. The supplementary
successor command is:

```sh
python -m projects.ukrmol_co.numerov_co_canary_condition_replay \
  --report OUTPUT/report.json --operands OUTPUT/operands --raw RAW \
  --output successor-verification.json
```

Extract the separately pinned `condition-support.tar.gz` after original source
and verification support. Only auxiliary rcond observations use the stated
relative rounding budget. Neither reader repeats timing rows. Reexecution retains its
own original exits and resources in a fresh directory; no numerical producer
module is replaced by verification support. Public evidence pins both sources.
