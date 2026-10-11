# CO sweep readiness — 11 October, 20:12 UTC

**The separately declared CAS12 2048→4096-root comparator is running on Sadaharu.**
B1 is constructing the 86,352-dimensional distributed Hamiltonian; no eigenpairs
have been exported yet. The completed [2048-root cost/export pilot](SWEEP_READINESS_20261011T1952.md)
remains passed. The [new contract](sparse-scattering-cas12-4096-comparator-contract.json)
and [launch capture](progress-observation-20261011T2012.json) fix the experiment.

## Declared comparison

Only the selected scattering root count changes. The native driver/dependencies,
target, CAS12-specific map, retained integrals/CONGEN, model, outer decks, solver
tolerances and energy grid retain the successful pilot's bytes/settings. Both
Pi sectors run on eight physical cores, CPUs 4–11, with an 80-GiB memory limit,
8-GiB shared memory and one library thread per rank. CPUs 0–3 remain available
for small controls. Launch gates pass at 107.62 GiB available RAM and 90.21 GiB
free disk. All 43 source and 97 input/evidence frozen hashes verify.

The independent reader compares shared-prefix energies, native channel/threshold/
multipole metadata, full 99-point phases and identical 1.8–3.3 requested-eV-window
fits with background orders 1–4. Budgets remain **0.05 eV position, 5% full width
with a 0.001-eV floor, and 0.05 rad phase modulo π**. Fits use the pinned native
requested-eV→Hartree conversion and identical initial seeds. Empty/MAXFIT/failed/
bound fits leave the pair verdict unset. Analytic phase-factor tests pass, and the
new reader reconstructs the completed baseline in the pinned image before launch.

Owner **1163485** has 24-hour native-stage, 50-hour native-container, one-hour
reader and 72-hour owner bounds. Observer **1163637** captures hourly for up to
72 hours or until the owner ends. Failure, timeout, changed frozen inputs or
`STOP` terminates the owner and preserves original exits. Successful native
completion triggers only the declared reader.

## Decision after completion

Review the pair's phase/fit sensitivity, expanded emitted energy coverage and
measured runtime/memory. A passing pair establishes stability between these
two tested root counts; an asymptotic omitted-spectrum error bound or stronger
convergence claim requires separately declared further evidence. Electronic-model,
explicit-A2-channel, pole/state identity, extraction, neutral correlation/domain
and restartable production-pilot gates still govern release of the full sweep.
