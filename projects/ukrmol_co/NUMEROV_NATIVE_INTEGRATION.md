# Native residual inside coupled Numerov

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
