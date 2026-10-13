# CAS12 scientific root comparison on the reproduced audited engine

The [completed matched2048control](progress-engine-baseline-20261013T0726.json)
passes both sectors' independent engine/runtime reproduction. Native execution
takes **11.565930hours**, with **21.964687GiB** kernel peak. All native/outer
stages, the reader and owner exit0. Exported energies, all99phase points and
fixed-window fits at background orders1–4 match the original baseline exactly.
Maximum native action residual is1.3276447836219887e-10Hartree and the Pi
difference is1.000000082740371e-9rad. The observation pins170remote records,
including final raw exports, resource measurements and completed owner evidence.

This satisfies the prerequisite for the separately declared
[2048-to-4096scientific control](cas12-audited-engine-4096-contract.json).
The comparison uses the completed audited2048run as its reference, so both
sides share the exact242installed engine artifacts, numerical libraries,
direct-MPI scientific driver and residual telemetry. Only requested scattering
roots change. The earlier failed4096release run and successful GDB-only4096B1
diagnostic remain preserved as distinct evidence.

## Execution and reading

Run the shared `projects.ukrmol_co.engine_baseline_owner` from the fresh frozen
source directory on Sadaharu, with `--campaign /home/kooza/ukrmol/co-sa-20261006`.
The contract's `control_kind=root-comparison` selects the scientific route.
The source provenance pins the committed owner/reader/declaration and inherited
byte-identical `sparse_pilot` and dependencies from the successful baseline.

The owner performs completed-baseline/engine/freeze/resource gates, a fresh
runtime numerical-library/telemetry gate, and baseline-to-itself independent
reading using the completed baseline's original **0.1meV/1e-6rad** contract.
It then invokes one4096root B1/B2 calculation, followed by the independent
`engine_baseline root-compare` reader. The reader verifies actual native
data-file hashes, decks differing only in`nstat`, unchanged outer inputs/map/
interfaces, shared-prefix energies, selection, raw exports/residuals, complete
grids and Pi symmetry. It measures full-grid phase changes and identical-window
fits at every declared background order.

Scientific budgets retain **0.05eV position**, **5% full width with0.001eV
floor**, and **0.05rad phase moduloπ**. Empty/MAXFIT/failed/bound fits leave the
pair verdict unset. False/unset results are written before the reader exits
nonzero, ending the owner. The allocation retains8physical cores4–11,
8MPI ranks, one library thread/rank,80GiB cgroup,8GiB shared backing and16GiB
native memory input. Stage/native/owner limits remain24/50/72hours, with
persistent `/home` outputs and observation-only hourly captures. `STOP` in
the fresh source directory ends the active step. No retries or further native
calculations are scheduled by completion.

## Interpretation

A passing comparison measures **pair stability**, without bounding the
asymptotic omitted pole sum. Further resolutions or independent omission
tests need separate declarations. The earlier crash causes and any engine
correction remain unresolved. Electronic-model/explicit-A2/three-anchor,
continuity/extraction, neutral correlation/domain and frozen restartable-pilot
gates still precede a production sweep and held-out potential fitting.
