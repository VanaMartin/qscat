# CO sweep readiness — 13 October, 14:04 UTC

**Scientific 4096-root B1 finishes matrix construction and the L2 update on all
eight ranks.** Owner 1246007 remains active on Sadaharu at **14:04:49 UTC**,
after **6.0868 hours** of native execution. Rank0 prints SLEPc/KrylovSchur setup
for dimension 86352, 4096 requested roots, tolerance 1e-12 and 1000 iterations.
Buffered tails do not identify the subsequent live solve operation. No pointer-
failure, signal, bounds or native-error marker appears on any rank.

B1 eigenpair/boundary exports and residual telemetry do not yet exist. B2,
outer stages, Pi and the independent scientific pair verdict remain pending.
Current sampling is **798.68% CPU / 14.46 GiB**, within the unchanged 8-core /
80-GiB allocation; these are ongoing samples, not final cost. The generic
SLEPc symmetry warning is preserved and does not replace the pending Pi check.

The [immutable milestone capture](progress-audited-4096-20261013T1404.json)
reverifies **43 source / 500 input-evidence hashes**, **242 installed artifacts**
and B1's actual native/outer inputs against the qualified audited2048baseline.
It pins **568 immutable remote records**. Observation-only owner 1246546 remains
active with seven hourly captures, 07:59–13:59 UTC.

Matched2048engine reproduction remains passed, with **11.565930 h /
21.964687 GiB** completed native cost. The
[scientific declaration](cas12-audited-engine-4096-contract.json) retains the
same engine, direct-MPI driver, numerical libraries and residual telemetry;
only requested roots change. Persistent outputs and `STOP` are under
`prepared/cas12-audited-engine-4096-source-20261013/`. Stage/native/owner limits
remain 24/50/72 hours. Other host workloads remain separate.

**Next:** finish B1 solve/export/outer stages, then B2 and the declared independent
raw/selection/shared-energy/interface/grid/Pi/fixed-window comparison. Scientific
budgets remain 0.05 eV position, 5% full width with a 0.001-eV floor and 0.05 rad
phase moduloπ, at every background order 1–4. A failure or unset verdict ends
without retry or relaxed limits. Completion schedules no additional native job.

Pair stability and absolute omitted-spectrum qualification remain separate.
Earlier crash causes/engine correction and model/channel/three-anchor/continuity/
extraction/neutral/restartable-production-pilot gates remain unresolved. These
still precede the sweep and held-out potential fitting.
