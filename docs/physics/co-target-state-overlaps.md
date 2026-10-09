# Frozen-core CO target-state overlaps

The coarse three-anchor orbital screen cannot identify many-electron states.
This toy-stage diagnostic calculates overlaps between determinant-CI target
wavefunctions at fixed saved orbitals, with both frozen core orbitals included.
All lengths/energies retain atomic units; overlaps are dimensionless.

## Determinant metric and frozen-core reduction

Each geometry's spherical molecular orbitals are orthonormal in its own AO
metric. For a physical cross-geometry comparison form
`S = C_left.conj().T @ S_AO(left,right) @ C_right`, restricted to core plus
active orbitals. The overlap between two Slater determinants is the determinant
of their occupied-orbital cross block. Alpha and beta determinant overlaps
multiply, with bra CI coefficients conjugated.

Both cores are doubly occupied in every determinant. Partition the cross metric
as `S = [[Scc, Sca], [Sac, Saa]]`. If `Scc` is nonsingular, block determinant
factorization gives

`det(S_occupied) = det(Scc) det((Saa - Sac Scc^-1 Sca)_active_occupied)`.

Therefore the total overlap is `det(Scc)^2` times the alpha/beta active-CI
contraction with the **Schur-complement metric**, not with `Saa` alone. This
includes core/active cross terms and physical overlap lost when atom-centred
orbitals move. Singular/ill-conditioned core blocks reject this implementation's
declared regime; they are not regularized or replaced by a normalized score.

The pure-Python API `projects.ukrmol_co.ci_overlap.frozen_core_overlaps` accepts
left/right arrays shaped `(roots, alpha_strings, beta_strings)`, the core-plus-
active cross metric, core count and active alpha/beta electron counts. Strings
are in ascending integer-bitstring order, matching PySCF determinant arrays.
Minor determinants define each spin's metric; two-sided matrix multiplication
contracts all root pairs. Memory scales with the spin determinant metrics and
CI vectors, rather than a full Kronecker determinant-space matrix.

## Validation and interpretation

Validation precedes the real-anchor screen: compare against a brute-force sum
of **full occupied Slater determinants** on complex, nonorthogonal toy frames,
with zero/one/two cores and unequal alpha/beta populations. Require `1e-12`
absolute/relative agreement, reverse-pair conjugate reciprocity, same-geometry
orthonormal-root identity and root-phase covariance. A core-only Slater control
checks the squared core factor. These are algebraic checks, not orbital or
electronic-model convergence.

The finite CO experiment re-solves eight roots per spin/irrep at each of the
three unchanged CAS11/cc-pVDZ checkpoints, keeping all original energy/spin/
physical and penalized residual gates. It compares physical overlaps and the
separately labelled AO-following Lowdin transport. The latter identifies
coefficient frames and is not a physical wavefunction overlap. Root-order
agreement, finite-manifold singular values, finite-manifold projection loss and
one-to-one maximum-overlap assignments are diagnostics; all candidate roots
and ambiguous assignments remain recorded. Coarse geometry spacing and eight
roots cannot certify a full continuous path or establish resonance-pole identity.
Projection deficits combine geometry/orbital-frame loss and excluded target
states; they are not isolated omitted-state populations.

The [completed eight-root companion](../../projects/ukrmol_co/qualification-evidence/target-state-overlap/README.md)
records state-order changes, retained-root assignments outside the first five,
and nearly singular finite manifolds. The [twelve-root closure contract](../../projects/ukrmol_co/cas11-target-state-overlap-roots12-contract.json)
keeps the same checkpoints and enlarges only the diagnostic root/space request.
First-eight energies are compared within 1e-7 Hartree. For each original vector,
sum squared overlaps with **all twelve new roots within 1e-7 Hartree of its
energy** and require a weight ≥0.99999999. This compares the matching eigenspace
without imposing an arbitrary individual vector in a near-degenerate manifold
at the root8 cut. Full-root CI/residual/spin/orthogonality gates remain mandatory.
Root assignment triggers and final continuity/pole verdicts retain their distinct
diagnostic and qualification roles.

This method stays at the validated Python toy stage. Saved CI vectors and
metrics permit independent reconstruction using direct full occupied minors;
no Rust optimization or QSCAT promotion is released.
