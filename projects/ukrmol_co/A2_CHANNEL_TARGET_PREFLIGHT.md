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
