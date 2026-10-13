# CO sweep readiness — 13 October, 00:13 UTC

**The matched2048-root audited-engine control has reached post-Class4 pair
blocks on all eight B1 ranks.** Each rank records Class1,Class3 and Class4
completion, followed by three Class5,three Class6 and24Class7 pair-block exits.
The latest flushed marker is **Class5 BEGIN, target pair7–8**. Its exit,
Class567/Class2+8 completion and L2 update are not yet recorded.

Owner1228004 remains active at **00:13:56UTC**, after **4.6257hours** of native
B1 since19:36:24UTC on12October. Sampling is **800.04% CPU /13.13GiB current
memory** under the eight-core/80-GiB allocation. These are ongoing samples,
not final cost. No pointer-failure, signal or Fortran bounds-error marker
appears. No CI/boundary exports or solver residual telemetry exist yet;
B2 and the outer stages have not started.

The [separate control](cas12-audited-engine-baseline-contract.json) continues
to use the exact completed address-audit engine and unchanged scientific
driver/input bytes. Runtime-library/telemetry loading and original-baseline
self-reading remain passed. The current run's engine-reproduction verdict
remains unset until both sectors, outer stages and independent comparison
complete. Same-root budgets remain **0.1meV position/full width /1e-6rad phase**,
with original residual/export/interface/grid/Pi checks.

Post-Class4 progress passes the earlier debug Class1 fault stage and the release
run's last confirmed Class4 milestone. It does not identify either failure's
cause or qualify an engine correction. Engine reproduction, root-pair
sensitivity and omitted-spectrum qualification remain separate.

The [00:13 capture](progress-engine-baseline-20261013T0013.json) reverifies
**41source/343input-evidence hashes**, **242installed engine artifacts** and
B1's five native inputs against the original baseline. It pins67immutable
remote records plus its frozen observation, including all eight rank logs,
the passed pre-launch gates and five hourly observer captures.

Owner service `co-cas12-audited-engine-baseline-20261012.service` and
observer1228557/service `co-cas12-audited-engine-baseline-observer-20261012.service`
are active. Observer captures run19:36–23:37UTC on12October. Source/run/observer
outputs remain under`prepared/cas12-audited-engine-baseline-source-20261012/`.
Stage/whole-native/owner limits remain24/50/72hours; `STOP` in that source
directory ends the active step.

**Next:** finish the remaining B1 pair blocks, diagonalization and outer stages,
then B2 and the independent same-root reader, or preserve the first failure.
A passing engine control precedes a separately declared scientific4096-root
comparison. No4096job is automatic. Model/channel/identity/extraction/neutral
and frozen restartable-production-pilot gates still govern the sweep and fit.
