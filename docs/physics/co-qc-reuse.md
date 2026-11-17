# Reusing accepted CO orbitals for native target preparation

The target orbital objective and native retained target space are different
objects. For example, native retention of nine A2 states need not reoptimize a
forty-component equal-weight orbital ensemble. Geometry-specific integrals and
native eigenvectors are still required; QC exports do not replace them.

`projects.ukrmol_co.run --target-only --reuse-qc-manifest PATH
--reuse-qc-sha256 SHA256` authenticates an accepted state-averaged CAS bundle
before preparing the native run. It excludes `--qc-only` and an initial-checkpoint
request. The manifest contains `source_directory` and a `files` mapping from each
of `config.json`, `target.json`, `target-diagnostics.json`, `resources.json`,
`co.molden`, and `co.casscf.chk` to its SHA256. Relative directories resolve beside
the manifest; the owner must pin that manifest and source capsule externally.

Replay requires the same geometry, basis, two-core active space and spin/irrep
orbital objective, successful original exit, convergence/fresh-CI/spin/energy
audits, atom/origin and core/active inventory agreement. It checks the Molden
digest and container format markers. These checks authenticate the accepted
record; they do not parse the complete HDF5/Molden wavefunction. Before a new
native release, independently open both files and compare molecule, basis, MO
coefficients and orthogonality in the pinned runtime, without an SCF/CI call.

The upstream quantum-chemistry slot runs the replay branch of `target.main`,
not `build_target`. Original target/diagnostic/checkpoint and Molden bytes are
copied without overwriting existing exports. `qc-reuse.json` distinguishes replay
from the original QC timing/source and records zero new QC calls. Failure has
no recomputation fallback. The ordinary optimization path remains the default.

Native root/residual checks, independent energies, properties, channel ordering,
boundary qualification and multi-geometry identity remain separate gates.
Replay alone adds no resonance geometry and authorizes no scattering calculation.

## Native audit output is per invocation

The retained `scattering_telemetry.c` interposer deliberately opens
`UKRMOL_EPS_AUDIT` with exclusive-create mode (`"ax"`). Give **each target-sector
solve a distinct, absent absolute filename**. Reusing an aggregate path aborts
the next rank-zero audit after `EPSSolve` returns, before native vector export;
a solve-site backtrace alone can therefore
misidentify an instrumentation failure as an eigensolver defect. Preserve the
original failure, inspect the actual interposer and test its file-opening behavior
before changing numerical settings or rerunning completed sectors. Never remove
an old audit file to make an unchanged attempt appear successful.

## Recover the failed stage, not successful peers

Distinguish the QC process exit from downstream consumer/controller exits.
`authenticated_qc` requires the original six-file bundle's
`resources.json["exit_code"] == 0`; a later log-reader failure does not change
that QC exit. Conversely, converged orbitals inside a timed-out QC process do
not meet this API. Separately qualified split completion needs its own explicit
contract and provenance, never an edited resource record or a bypass of this gate.

Reported orbital convergence with a false CI flag is weaker still. A fixed-orbital
completion cannot inherit that flag as a success: explicitly select the numerical
CI driver, retain the failed diagnostics and exit, and audit converged eigenvectors
against both physical and spin-penalized Hamiltonian actions. Keep the same
Hamiltonian, orbital frame, spin/irrep sectors, root counts and ensemble weights.
Root energies resembling the failed printout do not prove stationarity.

At the unchanged MO frame, reconstruct the equal-weight one- and two-particle
active densities from the qualified vectors and **recompute the orbital gradient**.
Use exactly the optimizer's independent, symmetry-allowed orbital-rotation mask
and gradient norm convention; do not include redundant active rotations or omit
allowed core/active/virtual rotations. Check the gradient assembly against an
independent integral/RDM contraction and a generated finite-difference energy
control. In PySCF's one-step convention the packed generalized-Fock difference
`g` is half the energy derivative for its antisymmetric exponential rotation:
compare a directional energy difference with `2*g.dot(direction)`, not `g`.
The optimizer's stationarity limit still applies to `norm(g)`. Preserve densities,
integrals, mask and gradient for saved-only replay.
If the original stationarity, residual, spin, root-coverage or occupation gates
fail, stop; neither a changed variable representation nor exported metadata repairs
the orbitals. A successful new bundle identifies its own completion source and
exit separately from the original failed optimizer process. Later independent
root-coverage probes and native property checks remain separate prerequisites.

Read the native command's declared output unit rather than assuming `fort.921`
for every symmetry. Preserve accepted sectors and real guards byte-identically;
a consumer-only error does not justify another native solve. A search that leaves
its frozen box supplies no pole, not evidence of absent resonance: diagnose
localization and qualify safeguards before a separately released remainder,
without silently reseeding, widening the box or relaxing numerical budgets.

QC, target-only and expanded-continuum layouts are distinct stages. Reuse
geometry-independent signed CSFs only after actual inventory/deck compatibility;
integrals, CIDATA, properties and native vectors remain geometry/frame-specific.
Verify destination phases and target-integral equivalence before boundary reuse.
