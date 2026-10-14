# Fixed-orbital A2 target/property preflight

The next decision-changing control concerns explicit channels, following the
[accepted root-resolution decision](ROOT_RESOLUTION_DECISION_20261014.md).
The [A2 audit](A2_CHANNEL_COVERAGE.md) finds that relevant CAS11 states leave
the five-root retained sector; root6 alone misses substantial root7/root9 weight.

## One bounded prerequisite experiment

The [declaration](a2-channel-target-preflight-contract.json) uses the already
qualified **stretched CAS11, R=2.5 bohr, cc-pVDZ** orbital frame. This is the
model on which the gap was measured. Two frozen A1 cores, ten active electrons,
CAS11 [5,3,3,0] and the **40-component orbital average remain unchanged**.
The native target representation becomes **[5,5,5,9] per spin**, or **48 retained
components**. This is channel coverage, not a 48-component orbital average.

Native five-root target vectors cannot supply the added A2 states. Existing
eight-root native imports are distinct orbital runs and omit root9. Reuse the
original integrals, target CONGEN expansions, Molden export and checkpoint;
compute eight small selected-target sectors into a fresh native CI file, then
one 48-state DENPROP property interface. There is no orbital optimization or
scattering Hamiltonian build in this preflight.

Already independently reconstructed twelve-root PySCF vectors provide the target
energies and action/spin/convergence prerequisites on the **exact same checkpoint**.
Fresh native per-root telemetry must pass 1e-9-Hartree action residuals; energies
must match within 1e-7 Hartree and the ground dipole within 1e-5 a.u.
The nonuniform channel map must include both A2 root9 states, and emitted
dipole/quadrupole records must reference and cover the actual 48 target IDs.
Native interface coverage is not an independent proof of every added transition
moment or of channel convergence.

## Cost and stopping rule

The original stretched target solves total about **8 minutes**, while DENPROP40
takes **39 minutes**. These are measured on the original engine/allocation;
the audited O1 engine and expanded property table can change cost. Use eight
physical cores/MPI ranks, one library thread/rank, 24-GiB cgroup and 8-GiB shared
backing. Target sectors are limited to 30 minutes each, DENPROP to two hours,
the whole native run to four hours and the owner to five hours. Every failure
preserves its original evidence and stops without retry or budget relaxation.

## Decision after completion

Review the native energy/residual and property coverage, exact target map and
actual cost. A passing preflight enables a separately declared **fixed-orbital
five-versus-nine A2 scattering comparison**, using each model's own target IDs
and multipoles. It schedules no scattering calculation itself. Other geometries
and CAS12 require their own coverage checks; neither global continuity nor
resonance-pole identity follows from this target import.

## Launch — 14 October, 11:13 UTC

The [dated capture](progress-a2-target-preflight-20261014T1113.json) records
owner 1282177 and observer 1282512, with a passing fresh runtime gate and
same-frame independent prerequisites. The original 40-state native map covers
656 dipole/quadrupole records and agrees with the fixed-orbital reference to
4.81677e-10 Hartree. The first two target stages finish with exit 0; A2 root9
exports and the complete 48-state interface remain pending at this cutoff.

Implementation is pinned to `aac698f310783b0f760ad197546e0a2b1f478b59`;
53 source / 50 input hashes and 242 installed engine artifacts reverify.
All 120 immutable launch-evidence hashes and 66 focused regressions pass.
The [launch queue](SWEEP_READINESS_20261014T1113.md) records the initial decisions.

## Completion and independent reading — 14 October, 13:19 UTC capture

The [completed capture](progress-a2-target-preflight-20261014T1319.json) records
native/owner exit 0 at **12:03:51 UTC**, all nine stages successful, and both
services inactive. Cost is **51.686 minutes / 5.466 GiB kernel peak**: eight target
solves total **4.291 minutes**, DENPROP48 takes **47.390 minutes**. The frozen
53 source / 50 input hashes and 242 installed engine artifacts reverify.

The [24-file local raw reconstruction](verification-a2-target-preflight-20261014.json)
independently reads all eight CIDATA sets and **48 complete vectors**, with no
undeclared trailing sets. Native energies agree with the independent fixed-orbital
reference to **5.25802e-13 Hartree**; property-table rounding contributes up to
**4.81663e-10 Hartree** relative to binary energies. Maximum native action residual
is **1.29215e-10 Hartree**, Pi splitting **3.69482e-13 Hartree** and vector Gram
error **7.54952e-15**. All original numerical budgets pass.

The model's own property map resolves **921 finite, unique dipole/quadrupole
records**, with valid state/irrep references covering all 48 states. Singlet and
triplet A2 root9 are native energy-ordered IDs **48 / 45**. Selecting the lowest
40 energy-ordered states would omit both; explicit nonuniform sector selection
is required. Ground-dipole difference is **1.88292e-11 a.u.**

A supplementary reader initially rejected an incorrect assumption that both
multipole irrep fields were always zero. Its original exit 1 and source hash are
preserved separately; the resolved raw format contains the actual C2v irrep
codes 0–3, which the completed reader validates against each state label. No
native retry or numerical-budget change occurred.

**Target/property input qualification passes.** The next decision is a separately
budgeted matched five-versus-nine A2 scattering comparison using this model's
own maps and multipoles. Every added transition moment is not independently
validated by interface coverage alone, and channel convergence/global continuity
remain open. No scientific calculation is active at the capture cutoff.
