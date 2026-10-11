# CO sweep readiness — 11 October, 22:50 UTC

**The CAS12 4096-root comparator is blocked by a native B1 SCATCI failure.**
Owner 1163485 ends at **22:48:11 UTC**: PETSc reports **SEGV signal 11 on MPI
rank 4**, followed by MPI_ABORT/native exit **59**. Container and controller
exits are **1**, with no timeout. The last completed logged matrix class is 4.
CI/boundary exports and solver telemetry are absent; B2 and the independent
comparison never start. No scientific container remains active. The hourly
observer is still waiting for its next scheduled completion capture.

The [failure capture](progress-failure-20261011T2250.json) pins nine immutable
owner/input/provenance/log/resource records. All **43 source / 97 input-evidence
hashes** reverify. Measured failed-attempt cost is **2.6112 hours / 15.6082 GiB**
inclusive cgroup peak, below the 80-GiB allocation. This is a failed B1-attempt
cost, not a completed two-sector cost. The root cause remains unknown; this
SEGV is distinct from the earlier 2048-root shared-memory SIGBUS.

The completed [2048-root cost/export pilot](SWEEP_READINESS_20261011T1952.md)
remains passed. The declared [4096-root comparison](sparse-scattering-cas12-4096-comparator-contract.json)
has produced no phase/fit result, so root-pair and omitted-spectrum verdicts
remain unset. The model, target, numerical budgets and original failures retain
their existing status.

**Next:** diagnose the preserved native failure, then separately declare any
bounded diagnostic successor needed to resume the unchanged scientific comparison.
No automatic retry or limit/tolerance change follows from this failure. Electronic,
channel, identity, extraction, neutral-domain and restartable-production-pilot
gates continue to govern release of the full sweep.
