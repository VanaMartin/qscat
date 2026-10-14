# Matched retained-A2 scattering control

The [native target/property preflight](A2_CHANNEL_TARGET_PREFLIGHT.md) passes on
the stretched CAS11 R=2.5-bohr/DZ frame. The measured coverage gap motivates one
**five-versus-nine A2 channel comparison**, with the original 40-component orbital
average, two frozen A1 cores, ten active electrons and all other target sectors
fixed at five roots per spin. The [declaration](a2-channel-scattering-contract.json)
fixes the experiment and stopping rules before launch.

## What changes in the scattering Hamiltonian?

Both models use the **same checked 48-state target eigenvector/property files**,
original integrals and original uncontracted scattering CONGEN expansions.
Only contraction counts/set/root sequences and the matching outer-interface
target map change. Five-A2 retention requires explicit per-sector selection:
its 40 states are not the lowest40 energy-ordered IDs of the 48-state inventory.
The expected dimensions are **27546 / 27818**, continuum coefficients **1410 /
1682** and boundary angular channels **250 / 298**. The increase is eight A2
components times34 continuum orbitals / six angular channels, with unchanged
L2 functions. Actual native telemetry and exports must confirm these sizes.

Both models request **2048 scattering roots**, using the exact completed audited
engine, eight MPI ranks on physical CPUs4–11 and one library thread/rank. Native
memp remains6GiB; cgroup48GiB/shared8GiB. Radius18bohr, maxl4, deletion1e-6,
propagation100bohr and the uncorrected selected-pole sum remain fixed.

## Baseline gate before paying for expanded scattering

The stretched geometry has saved **complete27546-pole dense boundary data**.
Replay only its outer region on the common99-point0.1–5.0requested-eV grid;
no dense Hamiltonian solve is needed. Then run the five-A2 selected model in
both Pi sectors and independently compare it to this dense replay. A false or
unset baseline verdict stops before nine-A2 scattering. This qualifies the new
geometry's economical baseline rather than assuming equilibrium root evidence
transfers unchanged.

If it passes, run the nine-A2 model once in both sectors and compare raw phases
modulo pi across all99energies, plus identical-window position/full-width fits
at every background order1–4 on **0.4–1.6requested-eV**. This fixed window is
already used in the stretched extraction controls. Keep the original budgets
**0.05rad / 0.05eV / 5percent full width with0.001-eV floor**. Energies passed to
the fitter are requested_eV×0.0735/2Hartree. Empty,MAXFIT,failed or bound fits
stay unset; raw residual1e-7Hartree/Pi0.002rad gates precede comparisons.

Changing channels changes poles, dimensions and interface tables. The comparison
checks each model's own declared target map and input bytes rather than demanding
cross-model equality of Hamiltonian spectra or interface digests. Shared
target/integral/CONGEN bytes and matched outer decks must agree exactly.

## Marginal value, cost and stopping

Existing equilibriumCAS11selected2048cost **2.454h / 7.7GiB**; original stretched
dense B1/B2 solves cost **1.296 / 1.211h** on four ranks. These differ in scope,
engine and allocation, so they are cost anchors rather than matched forecasts.
The expanded contracted dimension rises0.99percent and continuum/boundary sizes
about19percent; the same uncontracted CSFs dominate assembly. Each selected model
is bounded to18hours, each native stage to8hours and the owner to40hours, with
measured profiles for dense replay/five/nine separately. This is a bounded
several-hour-scale experiment under conditional cost expectations; the declared
caps are not forecasts. Two fresh selected models provide the matched comparison;
the cheap saved-dense replay can stop a weak baseline before expanded cost.

No retry, numerical-budget change or further root/channel ladder follows a
failure. A passing pair establishes observable insensitivity for this one
fixed-orbital channel control. A rejection identifies a channel or baseline
resolution discrepancy worth deciding on. Full channel convergence, global
identity, CAS12 qualification and the production sweep remain separate gates.

## Launch — 14 October, 13:58 UTC

Owner **1304346** starts on Sadaharu at13:58UTC from
`cdcf721a9472334845291c951462fd3bae9d0f6f`. The pinned-runtime dry prerequisite
checks the completed same-frame target/property import and both complete dense
boundary datasets: **27546poles × 250channels** per sector. The fresh numerical
library/telemetry gate and channel-input gate pass. The
[14:00capture](progress-a2-channel-scattering-20261014T1400.json) records dense B1
outer stages exit0, R_SOLVE129.31seconds, with B2outer replay still active.
No selected-model or channel-sensitivity verdict exists at this cutoff.

Frozen **58source / 58input hashes** and **242installed engine artifacts**
reverify. Observation-only observer **1304833** starts after readable-owner/live
service PID agreement, captures the nested dense/five/nine model directories
every30minutes and has the same40-hour bound. **82focused tests pass**, including
noncontiguous explicit selection, wrong-frame/CSF/deck rejection, the298-channel
export shape and false/unset baseline stopping before expanded-run creation.

Persistent source/output root:
`/home/kooza/ukrmol/co-sa-20261006/prepared/a2-channel-scattering-source-20261014/`.
Services: `co-a2-channel-scattering-20261014.service` and
`co-a2-channel-scattering-observer-20261014.service`.

## Dense replay complete; selected baseline active — 14 October, 14:03 UTC

The [14:03capture](progress-a2-channel-scattering-20261014T1403.json) records
the complete dense outer replay, all six stage exits0, **4.322minutes / 0.447GiB
kernel peak**, and printed B1/B2 phase difference0. The selected five-A2model is
now active on eight ranks; the nine-A2directory is absent. Dense-baseline and
channel-sensitivity verdicts remain unset. This cutoff pins138immutable remote
records and reverifies the same58source/58input hashes and242installed artifacts.
Live **797.72percent CPU / 2.494GiB** describes the active baseline, not final cost.

The [14:06independent local reading](verification-a2-dense-outer-20261014.json)
verifies six checksum-fetched raw files, both complete99-point phase grids,
original dense outer exits and **Pi phase difference0rad modulo pi**. Native fits
are nonempty and have no MAXFIT diagnostic in either sector. This independently
qualifies the dense outer replay's raw output; it does not yet compare the
selected baseline or channel models.

## Five-A2 baseline passes; expanded control active — 14 October, 16:34 UTC

The [16:30capture](progress-a2-channel-scattering-20261014T1630.json) records the
completed five-A2selected2048run and passing independent dense-baseline gate in
both sectors. All native/reader exits0; measured native cost is **1.7070hours /
5.894GiB kernel peak**. The
[local30-file reconstruction](verification-a2-five-baseline-20261014.json)
confirms2048poles,2048×1410continuum vectors and2048×250boundary amplitudes per
sector; residual **≤1.26308e-10Hartree**, Pi **1e-9rad** and exact boundary/CI
energy agreement. Maximum selected→dense phase difference is **1.02e-5rad**;
position/full-width fit differences over all background orders are at most
**4.57937e-7 / 1.65765e-6eV**, passing the unchanged scientific budgets.
Local reconstruction uses the public readers/fitter on fetched raw files,
independently of producer summaries.

The declared nine-A2model starts only after that pass. B1matrix assembly is
complete and its2048-root eigensolve is active at16:34UTC, with native contracted
dimension **27818** as declared. Live16:30CPU/memory **802percent / 3.43GiB** is
partial telemetry. B1exports/outer stages and B2 remain pending; no owner error
or timeout is observed. **Channel-sensitivity verdict remains unset.**
