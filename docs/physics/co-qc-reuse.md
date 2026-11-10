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
