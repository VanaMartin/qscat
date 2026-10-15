# Resolving the A2 phase rejection

The [completed channel control](A2_CHANNEL_SCATTERING.md) rejects phase
insensitivity at matched2048roots: **0.0679062rad versus0.05rad**, near1eV, while
position/full-width changes pass. The
[saved-output diagnosis](diagnosis-a2-phase-rejection-20261015.json) explains why:
at background order3, changing the fitted center contributes **0.06585rad**, then
the width contributes **0.00202rad** and the background **0.000067rad** at1eV.
The resulting fitted change0.06794rad is close to the observed0.06791rad.
This is an order-dependent fit decomposition, not unique physical attribution
or a resonance-pole identification. All background orders preserve the basic
center-shift explanation. Existing position and phase budgets measure different
sensitivities; neither should be changed to make them agree.

## The remaining model-specific question

The economical five-A2baseline passes against its saved dense observables, but
the nine-A2model has no same-model resolution comparison. Its highest retained
pole at2048 is **−110.310269Hartree**, versus **−109.951825Hartree** for five-A2:
adding channels lowers the same-count cutoff by **0.35844Hartree**. Equal root
counts therefore do not imply equal spectral coverage. The measured phase
exceedance margin is **0.0179062rad**.

Declare one [nine-A2only2048→4096control](a2-nine-root-control-contract.json).
It answers whether the rejection persists when the expanded model's spectrum
is refined. It can change our working resolution/attribution for this model,
which satisfies the [root-resolution stopping rule](ROOT_RESOLUTION_DECISION_20261014.md#prospective-stopping-rule).
The accepted CAS12decision remains in force.

Only `nstat` changes in the actual completed nine-A2decks. Reuse the same
target vectors/properties, integrals, scattering CONGEN expansion,48-state map,
outer decks and exact completed audited engine. Keep eight physical cores/ranks,
one library thread/rank, native6-GiBmemp,48-GiBcgroup and8-GiBshared backing.
The original failed scientific comparison and its exits remain immutable.

## Decision and finite cost

First pass the usual original-exit, input-byte, map/interface, prefix-energy,
export, residual, Pi and fixed-window/full-grid checks for nine2048versus4096.
Use unchanged **0.05rad / 0.05eV / 5percent full width with0.001-eVfloor**
scientific budgets and0.4–1.6requested-eV fits at background orders1–4.

For the additional attribution decision require **≤0.005rad** root-refinement
phase change in both sectors, comfortably below the0.0179062rad original
exceedance margin, and require the remeasured five2048versusnine4096phase
comparison still to reject0.05rad in both sectors. A true
`resolution_decision_verdict` means this named rejection is stable under the
declared refinement. It does not turn channel insensitivity into a pass or
bound the absolute nine-A2omitted tail. False/unset outcomes stop for a new
decision. Empty/MAXFIT/bound/failed fits stay unset.

Nine2048cost **1.849952hours / 6.306644GiB**. Earlier equilibriumCAS11doubling
cost ratio3.335/2.454≈1.36 on another engine/scope is a conditional cost anchor,
not a matched forecast. Budget **2–4hours for both sectors**, with hard caps
**8hours per stage / 12hours native / 14hours owner**. No repeat five-model,
target optimization or dense build is needed. No third resolution or later
channel/geometry follows automatically. Reserve four physical cores for small
jobs and retain all original limits/exits/telemetry.

## Launch — 15 October, 06:50 UTC

Owner **1334196** starts at06:49:53UTC on Sadaharu from frozen source
`dde183672d5930f2fa8f0d976462424a25bc94c9`. Fresh pinned-runtime prerequisites
reconstruct the original phase-only rejection while retaining its ownerexit1;
the new numerical-library/telemetry and reference gates pass. The
[06:50capture](progress-a2-root-control-20261015T0650.json) records B1matrix
assembly active, dimension27818 and continuum1682. Live **801.66percentCPU /
2.649GiB** describes this initial stage, not final cost. B2 and the independent
resolution decision remain pending.

Frozen **61source / 181input hashes** and **242installed engine artifacts**
reverify; the capture pins254immutable remote records. Observation-only observer
**1334564** launches after readable-owner/live-service PID agreement and captures
every30minutes within the same14-hour bound. **94focused tests pass**, including
actual298-channel unformatted exports, exact count-only deck transformation,
false/unset stopping and preserving the original rejection when the additional
resolution decision passes. Services are
`co-a2-nine-root-control-20261015.service` and
`co-a2-nine-root-control-observer-20261015.service`; persistent source/output root
is `/home/kooza/ukrmol/co-sa-20261006/prepared/a2-nine-root-control-source-20261015/`.
