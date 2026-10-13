# CO sweep readiness — 13 October, 08:03 UTC

**Matched-engine reproduction passes in both sectors; the separately declared
4096-root scientific comparison is running on Sadaharu.**

## Completed prerequisite

The [2048-root control](progress-engine-baseline-20261013T0726.json) finishes at
07:10:36 UTC with all native/outer stages, independent reader and owner exits 0.
Both sectors reproduce all exported energies, all 99 phase points and identical-
window fitted positions/full widths for background orders 1–4 exactly.
The maximum native residual is **1.3276447836219887e-10 Hartree**; Pi differs by
**1.000000082740371e-9 rad**. Native cost is **11.565930 hours / 21.964687 GiB**
kernel peak, excluding target construction and the separate independent reader.
This qualifies the selected audited engine/direct-MPI/telemetry runtime for
reproduction. Earlier crash causes and an engine correction remain unresolved.

## Active scientific comparison

[Owner 1246007](progress-audited-4096-20261013T0803.json), service
`co-cas12-audited-engine-4096-20261013.service`, starts native execution at
**07:59:37 UTC** after fresh runtime-library/telemetry and stricter baseline-
to-itself reading gates pass. At **08:03:21 UTC**, all eight B1 ranks pass graph/
retained-context checks, complete Class1 target 1 and begin target 2. No pointer-
failure, signal, bounds or native-error marker appears. B1 exports/residual
telemetry, B2, outer stages and the scientific pair verdict remain pending.
Sampling is **799.65% CPU / 10.73 GiB**, not final cost.

The [new contract](cas12-audited-engine-4096-contract.json) pins source commit
`c6e6f3b5f48d61851be706c82e23b186dde59edb`, exact installed engine, numerical
libraries/image and telemetry. It changes only requested scattering roots
2048→4096 in the scientific calculation. The completed audited 2048-root run
is the reference. All **43 source / 500 input-evidence hashes** and **242
installed artifacts** reverify; B1's actual four data inputs match that
reference byte-for-byte, and decks differ only in `nstat`. Outer decks match.
The launch observation pins **562 immutable remote records**, including the
completed prerequisite, freeze, fresh gates, launch commands and rank logs.

The native driver/dependencies are inherited byte-for-byte. Shared owner/reader
changes route the declared root comparison, audit actual inputs and preserve
false/unset results before exiting nonzero. **40 focused tests pass**, including
rejection of altered integrals/decks, failed/unset prerequisites and failed/unset
reader verdicts. The [execution note](AUDITED_ENGINE_4096.md) describes the route.

Observer **1246546**, service
`co-cas12-audited-engine-4096-observer-20261013.service`, is observation-only,
with its first immutable capture at 07:59:41 UTC and hourly captures thereafter.
Startup verifies a readable owner record and matching live service PID.
Persistent source/run/observer outputs are under
`prepared/cas12-audited-engine-4096-source-20261013/`. `STOP` in that fresh
source directory ends the active step. Allocation remains 8 physical cores
4–11, 8 MPI ranks, one library thread/rank, 80-GiB cgroup, 8-GiB shared backing
and 16-GiB native memory input. Stage/native/owner bounds remain 24/50/72 hours.

## Next decision

Finish B1/B2 and the declared independent comparison, or preserve the first
failure. The reader checks residuals, ordered finite exports, solver coverage/
selection, shared-prefix energies, unchanged inputs/maps/interfaces, complete
grids and Pi. Scientific budgets remain **0.05 eV position**, **5% full width
with a 0.001-eV floor**, and **0.05 rad phase moduloπ**, over the same window
and all background orders 1–4. False/unset verdicts stop without retry or
relaxed budgets. Completion schedules no additional native job.

Pair stability does not bound the asymptotic omitted spectrum. Further resolutions
or independent omission tests require separate declarations. Model/basis/ensemble/
explicit-A2 controls, three-anchor comparisons, continuity/extraction, neutral
correlation/domain and frozen restartable-pilot gates still precede a production
sweep and held-out potential fitting.
