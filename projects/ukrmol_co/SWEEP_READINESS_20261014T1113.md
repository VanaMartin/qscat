# CO next gate — 14 October, 11:13 UTC

**The fixed-orbital A2 target/property preflight is active on Sadaharu.**
It addresses the measured retained-channel gap, following the user's acceptance
of the completed CAS12 2048→4096 comparison as sufficient resolution evidence
for the next controls. The 8192 experiment remains administratively cancelled.

## Why this experiment?

The independently reconstructed CAS11 tracking audit finds relevant A2 states
outside the five-root sector. Root6 alone misses substantial root7/root9 weight.
An explicit-channel model must have those native eigenvectors, their properties
and its own target map before a scattering comparison can be meaningful.

The [declaration](a2-channel-target-preflight-contract.json) uses the passing
**stretched CAS11, R=2.5 bohr, cc-pVDZ** frame, on which the A2 coverage issue was
measured. It retains the 40-component orbital average, two frozen A1 core
orbitals, ten active electrons and CAS11 [5,3,3,0]. Retained native target counts
become **[5,5,5,9] per spin**, or **48 components**. This changes explicit target
coverage, not the orbital average. CAS12 and other geometries need their own maps
and coverage evidence; a CAS11 target result cannot establish CAS12 convergence.

## Launch evidence

The [11:13 UTC capture](progress-a2-target-preflight-20261014T1113.json) records
owner **1282177** and observation-only observer **1282512**. Same-frame twelve-root
independent action/spin/convergence prerequisites pass, and the original native
40-state map resolves all 656 multipole records. Maximum original native/reference
energy difference is **4.81677e-10 Hartree**. The fresh pinned-library/telemetry
runtime gate passes.

The first native singlet/triplet A1 target stages exit 0 in **17.40 / 34.81 s**;
per-root residual and selection gates pass before the next sector starts.
Singlet B1 is active at the cutoff. Native A2 root9 export and the 48-state
DENPROP interface remain pending. Frozen **53 source / 50 input hashes** and
**242 installed engine artifacts** reverify; the capture pins **120 immutable
remote records**. Sampled **798.72% CPU / 4.69 GiB** is live telemetry, not final cost.

Use eight physical CPUs 4–11, eight MPI ranks, one library thread/rank, 24-GiB
cgroup and 8-GiB shared backing. Target sectors are limited to 30 minutes each,
DENPROP to two hours, the complete native preflight to four hours and owner to
five hours. Original target solves total about eight minutes and DENPROP40 about
39 minutes; the expanded model/audited engine can change cost. The preflight
reuses saved integrals, CONGEN expansions and orbitals and builds no scattering
Hamiltonian. Failures stop with original evidence and no retry or gate relaxation.

## Next decision and ordered gates

1. Independently read the completed target energies, native action telemetry,
   exact nonuniform 48-state map, property coverage, ground dipole and frozen
   input hashes. Record actual native/DENPROP costs.
2. If that input qualification passes, separately declare the fixed-orbital
   **five-versus-nine A2 scattering comparison**, with each model's own channel
   IDs, thresholds and multipoles. Use economical root guidance and the existing
   phase/position/full-width budgets. Completion of this preflight launches no
   scattering calculation.
3. Complete matched three-anchor basis/active-space/ensemble/channel controls,
   continuity and background/window-sensitive extraction qualification. Existing
   compressed/stretched CAS12 and supported aug-TZ imports are numerically checked;
   they are not electronic-model convergence results.
4. Qualify neutral correlation/domain uncertainty and finish remaining raw public
   reconstruction, then freeze and pass a restartable production pilot before
   the sweep and held-out potential fitting.

Further scattering root expansion requires a specific decision-changing
discrepancy. Absolute omitted-spectrum error, electronic-model convergence,
global continuity and resonance-pole identity remain distinct claims. The
six-panel overview retains its latest completed-decision cutoff, 10:34 UTC;
new substantive native/property evidence will trigger the next figure release.
