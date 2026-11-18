# Full interior comparison with one continuum electron

**Synthetic-stage prototype.** The factorization below has a bounded Python
oracle, not a qualified native decoder, cross-geometry orbital integrator or
molecular state-assignment implementation.
It addresses the continuum information omitted by the
[bound-sector enclosure](co-interior-state-comparison.md), without constructing
the entire cross-geometry N+1 CI Gram matrix.

## Four blocks of one observable

Write the retained interior state as `Psi = Q + P`, where Q contains only bound
orbitals and

    P = sum_(I,j) A_Ij a†(u_j) Phi_I.

`Phi_I` is an N-electron target component and `u_j` a continuum spin orbital.
Spin coupling and native phases belong in A and the target components. The
full overlap is the coherent sum `QQ + QP + PQ + PP`; sector overlaps cannot be
normalized separately, squared separately or compared as ordinal profiles.

All orbitals must share a declared physical domain and coordinate frame. Denote
bound orbitals by b and continuum orbitals by c. Four one-electron cross metrics
are needed: `S_bb`, `S_bc`, `S_cb`, `S_cc`. Reversal swaps the geometry labels and
takes the conjugate transpose; the mixed blocks are not generally zero across R.

Define coefficient-removal states by exterior-algebra deletion:

    D_mu |i_1 ... i_N> = sum_k delta_(mu,i_k) (-1)^(k-1)
                                  |i_1 ... omit(i_k) ... i_N>.

In a self-orthonormal frame D is canonical annihilation. Otherwise this is a
coefficient-removal operation, not an assertion that the physical annihilator
of a nonorthogonal basis orbital has Kronecker anticommutators.

For target overlaps `O_IJ=<Phi_I^a|Phi_J^b>` and removal-state overlaps
`R_IJ,mu,nu=<D_mu Phi_I^a|D_nu Phi_J^b>`, the exact product-state block is

    <a†(u_j) Phi_I^a | a†(v_k) Phi_J^b>
      = S_cc[j,k] O_IJ
        - sum_(mu,nu) S_bc[mu,k] R_IJ,mu,nu S_cb[j,nu].

For determinants this follows by expanding
`det([[B,b],[c,d]]) = d det(B) - c adj(B) b`. The removal overlaps are its
signed cofactors, so no inverse of B or of a cross-core block is required.
In orthonormal self frames it is also the canonical anticommutator identity.
The exchange term is essential even when each own-geometry continuum is
orthogonal to its target orbitals.

The mixed block uses target/anion removal amplitudes:

    <a†(u_j) Phi_I^a | Q_b>
      = sum_nu S_cb[j,nu] <Phi_I^a | D_nu Q_b>.

The reverse mixed block is obtained by swapping geometries and conjugating, not
by assuming real coefficients or equal phase conventions. QQ is the ordinary
N+1 bound-determinant overlap, including frozen-core/active cross terms.

## A tractable representation, not a measured cost forecast

O, R and the mixed removal amplitudes involve saved target states and at most
one orbital removal. They can be evaluated in batches using determinant minors
and CI contractions. Removing an occupied core orbital must be included: at
different geometries a core-hole state need not be orthogonal to the other
geometry's frozen-core space. Preserve explicit core-hole sectors rather than
expanding every target into a huge unrestricted full-orbital FCI tensor.

A model with singlet and triplet targets coupled to a doublet needs the relevant
triplet spin projections and spin-changing removal blocks. A saved highest-spin
projection can generate its spin-lowered partner algebraically; that does not
authorize a new target solve or establish native Clebsch–Gordan signs. Compare
against directly decoded coupled CSFs before using the product representation.

Contract the coefficient matrices A into these blocks; do not materialize every
pair of expanded continuum configurations or a full inner-basis Gram matrix.
This removes a major representation explosion, but target-removal contractions
still need an actual bounded benchmark. No laptop or four-core runtime is
established by the algebra alone.

## Domain and native admission gates

Saved native integral files can retain AO and MO basis data as well as integrals.
GBTO's `read_basis_and_integrals` and `get_orbital_coefficient_matrix` provide a
prospective read-only export path. The existing project `native_target_io.f90`
demonstrates saved integral reading, but does not export all continuum orbitals.
Orbital evaluation must preserve the stored normalization; do not normalize the
same continuum basis again while reading it.

Native Gaussian overlap conventions can use all-space target/mixed integrals
and finite-sphere continuum–continuum tails. Neither a generic all-space AO
cross integral nor an arbitrary common sphere is automatically that metric.
Choose a justified common physical domain, quantify its difference from native
conventions and compute self and cross blocks consistently. If self frames are
not orthonormal there, evaluate the complete state's self norm with the same
four-block formula. Do not use `c†c` as an unverified physical normalization or
whiten a failed metric into agreement.

Required before molecular use: independent small determinant/minor tests of all
four many-electron blocks; complex phases and spin coupling; singular cross
metrics; core removal; native target/continuum contraction and phase agreement;
self Gram positivity and reverse consistency; three-level orbital-metric
convergence; numerical uncertainty and implementation/work caps. Qualify the
saved-basis export and scalar block kernels before any full-vector recovery.

The missing Q coefficients remain distinct from the saved P coefficients. A
native eigensolve to recover them is a separately authorized calculation, not a
read-only export. Even a qualified full interior overlap would be an endpoint
character comparison, not proof of unique continuation throughout an unsampled
interval, a continuous crossing location or electronic accuracy.

## Prototype and read-only export interfaces

`four_block_overlap.py` exposes `InteriorExpansion`, signed
`coefficient_removal` and `four_block_overlap`. Sparse determinant dictionaries
and target-by-continuum amplitudes define the two states. It returns complex
QQ/PQ/QP/PP/total amplitudes without normalization or phase projection. The
direct oracle is deliberately capped at five target electrons, ten bound spin
orbitals, six targets/continuum orbitals and one million total determinant pairs;
it is not a molecular FCI implementation. Tests compare all blocks with explicit
antisymmetrized exterior coordinates for orthonormal/nonorthogonal frames,
vacuum through three-electron targets, singular metrics, spin coupling, core
removal, phase/order changes and a complete-state self norm.

`native_orbital_io.f90` is a prospective saved-basis reader only. It copies MO
coefficients, stored native overlaps, AO shell/function addresses, geometry and
CGTO exponents/contractions/normalizations; it does not normalize or evaluate
orbitals/integrals. Its generated `--self-test` checks the double/int64 writer
layout without reading a molecule. `orbital_export.read_orbital_export` validates
the complete export and returns unchanged coefficients/metadata. Stored overlaps
are explicitly not a newly computed common-domain self or cross metric. A
matching compiler/module/library runtime and actual saved read still require
their own bounded execution gate; unsupported BTO input is rejected.
