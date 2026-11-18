# Finite-interior state comparison with bounded determinant truncation

This experimental method compares electronic **interior character**, not unique
analytic continuation between geometries. A resonance has no ordinary full-space
probability norm. Restrict both states to a declared common electronic domain and
use the Hermitian product there; do not substitute the outgoing c-product.
All coordinates and the sphere centre/axes must share one physical frame.

## Sparse determinant representation

Write an approximation as `p = sum_I d_I |I>`, using normalized antisymmetrized
Slater determinants of spin orbitals. Each occupation list is strictly increasing
in a declared common spin-orbital labelling convention. Include frozen occupied
core spin orbitals explicitly; neither core/active nor bound/continuum blocks may
be deleted from the metric. The orbital cross metric is

    S_ij = integral_domain phi_a,i(x)* phi_b,j(x) dx,
    <I_a|J_b> = det S[I,J],
    <p_a|p_b> = sum_IJ d_a,I* det S[I,J] d_b,J.

This permits different orbital counts at the two geometries and does not require
an invertible core or active block. Spin-orthogonal orbitals have zero mixed-spin
entries. Native CSF phases, target contractions and antisymmetrization must be
resolved before this API is given determinant amplitudes. Reordering occupied
orbitals without its determinant parity changes the state, not just its storage.

`sparse_state_overlap` in `projects/ukrmol_co/interior_state_overlap.py` evaluates
occupied minors in bounded tiles. Pair-count and tile-size guards apply before
allocating minors. This is a small sparse Python oracle, not an unrestricted
full-CI algorithm, native decoder, quadrature builder or automatic approximation.
It does not normalize either input or repair a nonphysical metric.

## Conditional residual enclosure

Let `a,b` be normalized physical interior states and `p_a,p_b` their orthogonal
projections onto retained determinants in their respective **self-orthonormal**
orbital frames. These need not be the same subspace. With captured squared norms
`A,B`, define `n_a=sqrt(A)`, `r_a=sqrt(1-A)` and likewise for b. If
`s=<p_a|p_b>`, Cauchy–Schwarz and three telescopings give

    |<a|b> - s| <= min(r_a + n_a*r_b,
                       r_b + n_b*r_a,
                       r_a*n_b + n_a*r_b + r_a*r_b).

For example, the first follows from
`<a|b>-<p_a|p_b> = <a-p_a|b> + <p_a|b-p_b>`.
Add an explicitly supplied absolute numerical-error allowance to this radius;
then square `max(0,|s|-radius)` and `min(1,|s|+radius)` to enclose the squared
overlap. The upper cap is a consequence of normalized physical states, not a
repair of data that violates the separate partial-overlap Cauchy bound.

`projected_overlap_interval` implements this conditional algebra. It does **not**
derive a numerical-error enclosure from successive quadrature differences or
infer captured norms from nonorthogonal coefficients. If self-orbital metrics
are not identities, `sum |d_I|^2` is not the captured norm: a different projection
and residual analysis is needed. Apparent normalization or lost norm must not be
renormalized away. Large omitted fractions correctly leave an inconclusive
interval; surface localization does not make those fractions small.

If the entire bound-only subspaces are known, not just the two retained
vectors, the enclosure can tighten further. Let `P_a,P_b` project onto those
complete subspaces, `x=||(I-P_a)p_b||` and `y=||(I-P_b)p_a||`. Orthogonality
of each omitted component to its own bound subspace gives the additional bound

    |<a|b>-s| <= r_a*x + r_b*y + r_a*r_b.

Use the minimum of the bounds, with the same explicit numerical allowance.
`interior_state_overlap.bound_sector_projection` evaluates this complete-space
projection with full occupied minors, independently of the core Schur reduction
and without an identity-sized collection of CI roots.
The complementary norms must be justified upper bounds, not fitted values.
Setting them to zero merely because two geometries are close is invalid.

## Validation and molecular admission

Synthetic tests compare tiled minors against explicit permutation determinants
and an independent exterior-power coordinate expansion. They cover complex
phases, rectangular frames, reverse metrics, self Grams, frozen-core/active cross
terms, bound/continuum exchange, orbital and determinant permutations, pair caps,
nonfinite/malformed inputs and three truncation sizes with known exact overlaps.
Null and ambiguous comparisons remain unassigned. Units are dimensionless.

Before molecular use, freeze the accepted coefficient combination and physical
domain, qualify native orbital/CSF mappings and self norms, and converge or bound
the cross metric independently. Fix an overlap decision and competitor/subspace
rule before reading its scores. Recovering full native vectors alone meets none
of those metric gates. A passing interior-character test does not locate a
crossing, prove uniqueness in an unsampled interval, establish electronic
accuracy or extend a fitted potential. There is no Rust optimization or QSCAT
promotion without separate profiling and independent physics review.

## Bound-orbital sector without continuum reconstruction

The uncontracted L² tail of a scattering CSF file contains actual N+1-electron
configurations. `native_l2_determinants.l2_determinant_transform` decodes only that
tail into alpha/beta determinant amplitudes, preserving substitution ordering and
the full occupied-permutation sign. The tail boundary and its ordering in the
contracted Hamiltonian must be verified against the pinned native source; a file
header dimension alone is insufficient. Frozen core occupation, spin projection,
bound-only occupations, every packed term and each CSF norm are checked. The
continuum contraction and its phases are not guessed.

In a CAS11 model with two doubly occupied frozen orbitals, a doublet N+1 state
has active occupations `(6,5)`. Its bound-sector tensor has shape `(462,462)`.
The existing `ci_overlap.frozen_core_overlaps` can therefore evaluate its exact
cross-orbital overlap, with core/active cross terms, without pairwise expansion
of the continuum. A synthetic fifteen-electron comparison against full
spin-orbital minors tests that use of the existing API.

One possible observable is an enclosure over **all norm-preserving completions**
of this bound-sector embedding. It fixes the physical bound-orbital overlap and
allows any orthonormal complementary continuum sectors. The residual inequality
above is valid for every such completion. This is stronger than a normalized
bound-only similarity: large omitted weights can prevent any assignment. It is
not a computation of the actual cross-geometry continuum metric. Applying it
to native states requires checking that their bound-sector metric agrees with
the declared embedding and that the native orthonormal-complement convention is
consistent. GBTO tail approximations and finite-sphere versus all-space target
conventions must not be silently exchanged. If these conditions are unresolved,
the output remains conditional model-character evidence, not a physical overlap.
