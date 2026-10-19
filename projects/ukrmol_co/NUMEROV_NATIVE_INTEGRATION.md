# Native residual inside coupled Numerov

## Producer result — 19 October, 15:53 UTC

Producer source freezes/pushes at **`9ed0123b62345d49439999e83fe5aa412816f463`**.
All **538/538 checks** pass over 50 propagations; original test/runtime exits are
0/0. The qualified-host tests exercise both native-format skips from the Mac:
**29 tests pass**. Independent saved verification and public numerical replay
are separately frozen at **`09a4e83`** and pass as recorded below.

The 867,632 audited native solves produce **1,735,263 residual evaluations**.
Every full residual matches NumPy numerically, with **zero gate disagreements**.
The maximum independently accepted backward error is **1.95195573e-18**, below
32 extended epsilons; at most **one correction** is used. Every native/reference
pair has exact decimal K, full numerical diagnostics and invariants. All 20
dense accuracy/timing rows reproduce archived K and diagnostics exactly.

| Model | Reference marginal s/interval | Native marginal s/interval | Marginal speedup | Full N128 speedup |
|---|---:|---:|---:|---:|
| Free | .134776679 | .059406739 | **2.26871×** | **1.87680×** |
| Varying | .270785560 | .119819247 | **2.25995×** | **1.84298×** |

Both predeclared >=1.5× marginal-speedup hypotheses pass. The full propagation
gain is smaller because Taylor/tail setup and other matrix algebra are retained.
These are unaudited timing rows; independent residual checks are confined to
accuracy rows. Conditional synthetic recurring-only estimates for the old
60-solve scope are about **10.09–20.34 single-worker days**, replacing the
23–46-day unoptimized synthetic estimates within that scope only. Actual CO
conditioning, refinement counts and throughput still require a costed pilot.

Short-model discrete arithmetic error is at most **6.66460e-19**. Exact
inverse-square S errors are **5.83160e-7 / 3.64743e-8 / 2.28025e-9**, showing
approximately fourth-order reduction. The 200/400/800-bohr stress arithmetic
errors are **3.10828e-11 / 4.99488e-10 / 6.62407e-10**; the largest continuous
K difference is **4.11497e-8**. The four-channel mixed path may differ from the
old <=3-channel arithmetic (largest archived stress difference 1.37272e-9),
while native and retained NumPy residual paths within that mixed algorithm
agree exactly. All original budgets remain.

Runtime completes at **15:53:22 UTC**, **415.725390 s wall / 415.197557 sampled
CPU s / 147,685,376-byte kernel peak**, with 2,075 samples. This is below the
600/1200/1800-s prospective budgets. Tests cost 1.272381 s and 117,735,424-byte
peak. All **94 frozen sources / 242 installed artifacts / three prior reports /
two wheel components** verify before/after. The owner's legacy
`captured_operand_hashes_verified_*` fields count those five prior/component
hashes; they do not identify new A/B/X captures. The 16-file collection preserves
original commands, exits, sampled resources and checkpointed report in
`prepared/numerov-native-integration-20261019/`.

Original interior false/exits1/1 and 132/594 historical identity coverage remain;
this adds no CO identity, UKRmol invocation or inner solve.

## Independent public result — 19 October, 16:05 UTC

The independent saved verifier reconstructs all **538 decisions** in the pinned
runtime, exits 0/0 in .721554 s. A fresh download of the published original
source, result archive, support, three prior reports and wheel then repeats the
complete 50-propagation scope. **Every numerical row, decimal K, full diagnostic,
residual-audit aggregate, reference error and numerical decision reproduces
exactly**. New timings are observations: marginal speedups **2.27237× / 2.25897×**
again pass both hypotheses, while complete N128 gains are 1.87808× / 1.84304×.

Replay completes at **16:04:32 UTC / exits 0/0**, **416.168333 s wall /
415.640450 sampled CPU s / 147,382,272-byte kernel peak**, 2,077 samples, below
450/900/1800-s budgets. All 96 effective sources, 242 installed artifacts and
five prior/component hashes verify before/after. Its fresh persistent directory
is `prepared/numerov-native-integration-public-replay-20261019/`.

All **14 final public payloads** and the final index fresh-download and verify,
including explicit whole-wheel/initializer/binary checks. Public evidence:
[`qualification-evidence/numerov-native-integration/`](qualification-evidence/numerov-native-integration/README.md).
The [16:05:32 UTC capture](progress-numerov-native-integration-result-20261019.json)
reverifies all 31 producer/replay collected hashes, both frozen source trees and
242 installed artifacts. No scientific container remains; 59.92 GiB is free.

The next decision-changing experiment is a separately costed saved-CO-matrix
conditioning/throughput canary. The synthetic result does not establish actual
CO refinement counts, native K accuracy or full-scope runtime, and does not
release the original withheld interior successor.

## Decision and scope

The native extended residual primitive passes its independent numerical and
performance controls in [NUMEROV_NATIVE_RESIDUAL.md](NUMEROV_NATIVE_RESIDUAL.md).
Integration must establish that substituting it inside the bounded mixed solve
preserves accuracy, accept/reject decisions and useful interval throughput.

The [finite contract](numerov-native-integration-contract.json) retains the
original `numerov_typed_coupled` propagator, NumPy residual and binary64 LU
preconditioner. `numerov_native_coupled` substitutes only residual evaluation;
it retains three correction attempts and the 32-epsilon backward-error gate.
It is mixed precision, not extended LAPACK. The single-worker binding hook
restores the original solver on normal or exceptional exit, including around
the existing dense timing wrappers.

## Smallest decision-changing comparison

The prior short/physical/stress controls have at most three channels and use
the small direct solver. Append uncoupled **closed** spectators to make four
channels: original coefficient blocks are unchanged, spectator Q is .3 bohr,
angular momentum and potential coefficients are zero, and each spectator
threshold is incident energy plus 2+i Rydberg. No open observable is added.
This forces the integrated mixed path while retaining independent original
discrete/continuous/analytic references. Tests compare the embedding's open K
against the original small model and check actual mixed solve counts.

The fixed scope has **50 propagations / 1,730,352 intervals / 1,737,200 mixed
solves**:

| Stage | Propagations | Purpose |
|---|---:|---|
| Short | 18 | Three coupled models, three steps, retained discrete and continuous references |
| Physical | 6 | Exact inverse-square S, invariants and step convergence |
| Dense accuracy | 4 | Free/varying 298-channel N8, archived K/diagnostics and independent oracles |
| Stress | 6 | All original 200/400/800-bohr finest-step endpoints, actual archived binary64 Q |
| Timing | 16 | Free/varying N32/128, two opposite-order reference/native repetitions |

Every native accuracy residual is independently recomputed as `b-a@x` with
NumPy, including rejected refinement evaluations. The auditor records exact
numeric residual counts, accept/reject disagreements, accepted backward errors,
iterations, dtype and accepted solve counts. It covers **867,632 native solves**.
Dense accuracy auditing is separate from unaudited performance rows, so the
independent NumPy work does not inflate native timing. Exact paired K and full
diagnostics are required, as are archived dense outputs. Timing hypotheses
are observations, separate from numerical qualification.

The exact models, grids and references come from three checksum-pinned prior
reports. Reusing them avoids repeating independent high-precision calculations
merely to change the residual implementation. Whole-wheel and both nested
component hashes pin the already qualified Linux binary; no rebuild is needed.

## Cost and stopping rule

Use Sadaharu's reserved physical cores 0–3, one worker, single-threaded BLAS and
4 GiB memory. Prior dense timings imply under about 180 seconds for the timing
portion under the 1.5x native speedup hypothesis. Four-channel mixed-solve and
per-evaluation audit overhead at the long stress endpoints is unmeasured;
the complete prospective budgets are **600 s expected / 1,200 s upper /
1,800 s hard cap**. Preserve sampled CPU/memory, original command/container
exits and fresh persistent `/home` directories.

Failed preflight stops before scientific runtime. A false numerical gate
completes its fixed stage, then stops before later stages. Unsupported arithmetic
or an exception is unset. There are no automatic retries, extra grids or relaxed
gates. This synthetic qualification does not itself release the withheld CO
interior successor; that requires a separately costed pilot.

## Reproduction

With the pinned wheel extracted on the qualified host, the pinned mpmath wheel
on `PYTHONPATH`, and the prior reports named `arithmetic.json`, `coupled.json`
and `dense.json` in `PRIORS`:

```sh
python -m projects.ukrmol_co.numerov_native_integration \
  --contract projects/ukrmol_co/numerov-native-integration-contract.json \
  --priors PRIORS --output report.json
```

The output is checkpointed after every completed row and stage. Frozen source,
wheel/component/prior verification and resource ownership belong to the launch
record; the scientific report preserves numerical rows and gate decisions.

## Independent reconstruction and numerical replay

`numerov_native_integration_replay` reconstructs decisions from decimal K,
checksum-pinned original references and saved residual observations. It checks
the declared row order, embedding coordinates, solve counts, continuous/discrete
errors, physical invariants, stage summaries, timer partitions and precision
trends. Residual counters are saved observations until numerical reexecution;
the saved audit alone does not reproduce individual matrix products.
Exact derived K-to-S/invariant arithmetic is checked in the pinned Linux runtime;
different host BLAS arithmetic can change last bits of those metrics even when
the saved decimal-to-binary64 K conversion is exact.

```sh
python -m projects.ukrmol_co.numerov_native_integration_replay \
  --report report.json --priors PRIORS --output saved-decisions.json
```

Adding `--reexecute` with a fresh output directory repeats the frozen scope and
remeasures every native accuracy residual. Decimal K, full numerical diagnostics,
residual-audit aggregates, reference errors and decisions must match. Wall times
and timing hypotheses are fresh observations; they do not need to equal the
producer. A public replay uses downloaded source/result/dependency/wheel bytes,
verifies both nested extension components and retains its own original exits
and sampled resource record.
