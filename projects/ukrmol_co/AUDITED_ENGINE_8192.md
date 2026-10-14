# CAS12 third-resolution control: 2048 / 4096 / 8192 roots

The completed4096control independently passes both-sector pair stability against
2048, with phase change1.67e-5rad and maximum same-window fitted position/full-
width changes0.00224/0.04195meV. Its measured native cost is14.924481h/21.736774GiB.
That motivates the [separately declared8192control](cas12-audited-engine-8192-contract.json):
test persistence through a second doubling under the same numerical/electronic
inputs, rather than treating one small successive difference as convergence.

## Decision and scope

Run one fresh B1/B2 native calculation on Sadaharu, using8physical CPUs4–11,
8MPI ranks, one library thread/rank,80GiB cgroup,8GiB shared backing and16GiB
native memory input. Stage/native/owner limits remain24/50/72hours. The completed
2048/4096costs justify this bounded test; larger Krylov workspaces/iteration
counts may change cost, so the limits are not a runtime forecast.

The science driver/dependencies, exact242installed audited artifacts, numerical
libraries/image and residual telemetry are inherited byte-for-byte. Only `nstat`
changes. The frozen source records the committed shared-owner/reader changes.
Fresh engine/freeze/resource/prerequisite checks, runtime symbol-loading and
strict2048baseline-to-itself reading precede native execution.

## Three-resolution reading

`engine_baseline sequence-compare` independently reads all three root counts and
compares **2048→4096**, **4096→8192** and **2048→8192**. The endpoint comparison
is mandatory: two adjacent budget passes can coexist with an endpoint failure.
Every comparison verifies actual native data-file hashes, decks differing only
in`nstat`, original stages, raw poles/coefficients/residuals, native selection,
target mapping, channel/threshold/multipole interfaces, shared-prefix energies,
complete grids, Pi and identical-window fits at background orders1–4.

Budgets remain0.05eV position,5% full width with0.001eV floor and0.05rad phase
moduloπ. Width budgets use each pair's own baseline fit. Empty/MAXFIT/failed/bound
fits stay unset; the sequence verdict passes only if all three pairs pass.
False/unset results are preserved before the reader/owner exits nonzero.
Completion schedules no additional calculation. `STOP` in the fresh source
ends the active step; an hourly observer records immutable status captures.

## Interpretation and independent omission gate

Small differences at three resolutions establish empirical sequence stability,
not a rigorous error bound on the uncorrected omitted pole sum. No decay-rate or
monotonicity assumption is imposed. Absolute omission remains separately unset.
The current saved native data contain no assembled-Hamiltonian checkpoint; a
full-spectrum oracle or independently validated remainder/resolvent test needs
separate design and validation. Repeated root doubling alone cannot supply that
bound, and this declaration schedules only8192.

Engine reproduction, earlier crash causes/correction, electronic-model/explicit-
A2/channel/three-anchor/continuity/extraction/neutral and frozen restartable-
pilot gates remain distinct. Neither this control nor its pass releases a sweep
or qualifies a resonance pole/potential fit.
