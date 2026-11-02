# Stretched DZ/TZ target correspondence from saved native wavefunctions

## Decision and bounded experiment

The user authorizes continued experiments after the completed finite-pole review.
The first saved-data review (2November2026, started17:00UTC) selects **one physical
target-correspondence screen**, not another larger-basis/scattering calculation.
The103.370meV/16.54% local pole shifts already prevent a production-insensitivity
claim. The useful remaining question is whether the two retained64-state target
manifolds represent predominantly corresponding physical subspaces, or whether
ranked-root exchanges/truncation obstruct that interpretation.

Both native exports retain64 CSF coefficient vectors, all eight CSF definitions,
their orbital checkpoints/Molden files and exact energies. Independent TZ
determinant-CI vectors for the actual imported checkpoint are publicly archived
in the existing feasibility proof; DZ determinant arrays are not saved in the
inspected coverage directories. Do not repeat CI to fill that cache: re-expand
the saved native vectors algebraically. Existing CSF decoders are CAS11-specific
(7/3/3), so they cannot be applied unchanged to the CAS12(8/3/3) target.

No QC, CI eigensolve, Hamiltonian/spin action, native engine, property export,
scattering map or protected holdout read is permitted. Checkpoints/CSFs/CIDATA
remain immutable. A failure stops without a conversion fallback or new solve.

## Mathematics and representation

Use the ordinary Hermitian electronic inner product, Hartree/bohr internally.
Native spatial orbitals are ordered by C2v symmetry, alpha then beta within
each spatial orbital; two doubly occupied A1 cores and12 active orbitals give
14 electrons. Each packed CSF term substitutes electron slots in the saved
reference. Multiply its coefficient by the fermionic parity that sorts those
slots into alpha-spatial order followed by beta-spatial order. Remove only the
four always-occupied core spin orbitals. Active determinants use ascending
integer-bitstring order, matching the existing `ci_overlap` contract.

For native eigenvector V and sparse determinant/CSF expansion D, C=D V.
Require source phases to match CIDATA and independently recover each CSF's
saved first-term canonical phase. Do not multiply an extra arbitrary phase or
normalize converted roots. Gram/norm preservation and agreement with all eight
saved TZ determinant-CI sectors (up to individual overall root phase) validate
the real conversion before interpreting any cross-basis overlap. In particular,
matching energies alone cannot validate a coefficient-address or sign mapping.

Read the actual checkpoint orbital coefficients and target's explicit
(irrep,ordinal) inventory. Check self AO-metric orthonormality, nuclear geometry,
two core/12 active occupation, and orbital symmetries. Reorder each frame into
the same symmetry-grouped core-plus-active indexing. Form the rectangular AO
cross integral S_AB and M=C_A^dagger S_AB C_B; check reverse reciprocity.
Different AO dimensions are expected and are not padded/projectively normalized.

Evaluate all8x8 physical root overlaps per spin/irrep using the existing
frozen-core Schur-minor formula, and independently compare full occupied
Slater minors, reverse overlaps and own-frame Gram matrices. Core/active cross
terms and the doubly occupied core determinant factor must remain included.
Save the complete matrices, singular values, directional capture weights and
maximum-weight assignments; do not renormalize lost physical overlap.

## Predeclared gates and interpretation

Numerical gates: complete record framing/EOF/counts, finite coefficients,
source/CI phase identity; CSF norm1e-12, root/self-MO Gram1e-9; independent full-
minor/Schur and forward/reverse agreement1e-10; AO reciprocity1e-12. Saved TZ
CI energies must match the same native roots within1e-7Ha, physical residual/
spin evidence must remain qualified, and corresponding vector squared overlap
must be>=1−1e-8. No new residual action or eigensolve is substituted for that
authenticated saved evidence. All controls precede cross-basis interpretation.

The existing finite **dominance screen** is extended from retained-five to the
actually retained-eight: min singular-value squared>=0.5 in every sector.
This asks whether every retained direction has a majority projection into the
opposite retained manifold; it is not an accuracy bound or a theorem of unique
root identity. Report all individual squared overlaps/assignments without an
additional invented root-identity threshold. No ninth-root buffer is available;
low capture cannot distinguish rotation into omitted roots from changes of
electronic description. Neither a pass nor rejection proves N+1 resonance-state
identity or absolute basis/correlation convergence.

On a numerical pass, report the full physical screen even if some sectors reject.
A rejection blocks a clean same-manifold interpretation and redirects retained-
model selection; do not increase roots automatically. A pass permits a later
explicit balanced electronic control to be designed, not a production release.
No B2 rescue, geometry, curve, fit, sweep or114-holdout analysis follows.

## Resource and effort envelope

Reuse one set of byte-identical DZ/TZ CSFs where hashes prove equality, both
native coefficient sets and the eight saved TZ reference sectors. About0.66GB
of CSF records,72MB native coefficients and278MB saved determinant references
are inputs, not repeated calculations. Four physical cores/8GiB/no network/
one BLAS thread; expected120–900s (not measured at CAS12), hard1800s science/
2100s owner,30s read-only observations. At most60minutes narrow implementation/
freeze effort, runtime known-truth controls and exact immutable source/contract/
owner push before any molecular re-expansion or overlap calculation.

This remains project-stage validation; no Rust optimization or qscat promotion.
The original40degree/capsule/adapter failures, qualified finite poles,2048 baseline,
nine working reference, parked new50 and equilibrium exclusion all remain intact.
