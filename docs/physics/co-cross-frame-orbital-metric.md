# Loader-aware cross-frame orbital metric on one physical ball

`native_cross_metric` extends the
[one-frame loader-aware prescription](co-loader-aware-self-metric.md) to two
separately tagged frames. Positions/radii are bohr. The cross Gram is a one-body
ingredient, not a many-electron state overlap or a continuation assignment.

## Rectangular assembly and fixed scales

`all_space_cross_metric(left, right)` assembles every directed shell pair using
the retained analytic Cartesian moments or independent positive Gaussian-product
Hermite integration. It does not impose reversal by mirroring the upper triangle.
The same per-basis shell/AO/primitive/angular caps and a rectangular pair-work
budget apply before allocating the matrix. Compute the reverse independently
when checking orientation, complex phases and rounding error.

For each frame, `generator_scales` uses its own fixed native generation sphere
and read-reconstructed shell amplitudes. Keep these diagonal AO-row factors
tagged with that export; do not replace either by a new normalization on each
cross-comparison radius. Unchanged MO coefficient arrays are contracted only
after applying those factors exactly once.

`cross_ball_witness` replaces origin-centered continuum–continuum shell entries
by their exact equal-angular-degree radial ball integrals, with explicit zeros
for angularly orthogonal components. Target-containing entries retain all-space
integrals and a conservative, contraction-aware exterior uncertainty. Let
tau_left/right be the scaled square-root tail bounds from the self prescription.
Then E_ij=tau_left_i tau_right_j for entries with at least one target AO, and
E_ij=0 for explicitly integrated continuum pairs. For complex MO coefficients,

    S_ab = C_a† D_a H_ab D_b C_b,
    |S_ab - exact_uniform_ball_cross| <= |C_a|^T E_ab |C_b|.

No tail is declared zero merely because it appears small; no Gram whitening,
eigenvalue clipping or fitted phase correction repairs a failed metric.

## Domain, reversal and block admission

Both exports must establish the same origin, molecular-axis orientation, charge
inventory and physical comparison sphere. A finite all-space Gaussian extension
beyond a native interior sphere is not a physical outgoing resonance exterior.
Changing the comparison domain is a separate physical decision, not an automatic
radius convergence ladder.

Validate rectangular moment/Hermite agreement, three polynomial-exact Hermite
levels, raw reverse assembly, identical-basis reduction to the self prescription,
and small direct finite-ball controls with propagated complex-coefficient bounds.
At actual saved-input admission, independently qualify the newly acquired self
metric and combine it with the already qualified other frame's self matrix.
The joint Gram [[G_a,S_ab],[S_ab†,G_b]] is positive semidefinite for exact common-
domain integrals. Its numerical spectral allowance must include matrix-dimension
amplification of entrywise rounding bounds plus the tail-bound row sums; an
observed negative eigenvalue must not be silently clipped.

Use explicit native continuum flags to select bound/continuum coordinates and
obtain bb, bc, cb, cc while preserving the original symmetry/orbital order.
Mixed blocks generally survive across frames. Their magnitudes, singular values
or orbital-subspace overlaps do not establish an N+1 state overlap. A coherent
native P/Q adapter and complete physical self norms remain separate requirements.
