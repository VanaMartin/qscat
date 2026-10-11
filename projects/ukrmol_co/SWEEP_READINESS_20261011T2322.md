# CO sweep readiness — 11 October, 23:22 UTC

**The separately declared B1 native diagnostic is running on Sadaharu.**
Owner **1168755**, service `co-cas12-native-debug-20261011.service`, launches
SCATCI at **23:22:00 UTC** after its diagnostic build and all twelve original
serial/MPI engine reference checks pass. Normal-exit, first-fault SIGSEGV/source
stack and Fortran bounds-error/source-backtrace debugger smokes also pass.
Eight per-rank GDB logs show native execution. The **23:22:43 UTC** capture
samples **801% CPU / 11.33 GiB current memory**; these are ongoing samples,
not final diagnostic costs.

The [diagnostic contract](cas12-native-debug-contract.json) and
[implementation/triage note](NATIVE_FAILURE_DIAGNOSTIC.md) explain why this
run rebuilds the Hamiltonian from saved inputs. All **47 source / 119 input-evidence
hashes** reverify, as do **242 installed diagnostic artifacts** before native
launch. The exact failed B1 deck, target CI/properties, integrals and CONGEN
bytes are reused. The build adds bounds/backtrace/debug information and six
flushed Class 5/6/7 entry/exit rank/target-pair markers over the original pinned
native sources and numerical libraries.

Native allocation remains eight physical cores **4–11**, **80-GiB memory**,
**8-GiB shared backing**, one library thread per rank and **16-GiB `memp`**.
The native limit is 24 hours, container 25 hours and whole owner 30 hours.
Hourly observer **1169951** records progress. Build/check/native failure,
changed frozen inputs, `STOP` or timeout ends the owner with preserved exits.
No second sector or outer-region calculation is scheduled by this diagnostic.

The original [4096-root scientific comparator failure](progress-failure-20261011T2250.json)
and failed service remain preserved. Its observer **1163637** finishes normally
at **23:11:55 UTC** after four captures, including its final owner-completion
reading. The successful 2048-root cost/export/Pi pilot remains passed.

**Next:** inspect the diagnostic's first fault or completion evidence before
declaring a correction or a matched-engine baseline control. A diagnostic
completion is not scientific root-count qualification. Pair sensitivity and
omitted-spectrum verdicts remain unset. Electronic/channel convergence,
state/pole identity, width extraction, neutral correlation/domain and a frozen
restartable production pilot continue to gate the sweep and potential fitting.

The [dated launch capture](progress-diagnostic-start-20261011T2322.json) pins
19 immutable remote freeze/build/reference/debugger/previous-observer records
and an immutable copy of its live observation. All earlier plotted numerical
evidence remains unchanged.
