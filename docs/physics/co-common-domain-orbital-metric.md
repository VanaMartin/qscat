# A common-domain orbital metric before state comparison

This is a Cartesian-Gaussian toy construction and a native admission contract,
not a qualified GBTO shell adapter or a molecular overlap implementation. It
supplies a candidate one-electron ingredient for the
[coherent four-block observable](co-full-interior-metric.md). Units are bohr.

## Define the observable before integrating

Choose one body-fixed coordinate frame, directed molecular axis, origin O and
physical ball D = {x: |x-O| < a}. Both orbital sets and every bb/bc/cb/cc block
must use that same D. A declared nucleus-aligned frame can compare geometries;
independently translating one orbital set or its integration ball changes the
observable. Moving both sets and D by the same rigid motion does not. Native
symmetry labels, m-order, real-harmonic signs and origin must be explicitly mapped.

The one-electron matrix is

    S_ab[mu,nu;D] = integral_D chi_a,mu(x)* chi_b,nu(x) d^3x.
    S_ab^MO = C_a† S_ab^AO C_b.

Spin adds the spin-coordinate inner product, not a geometry-dependent orbital
label match. Compute self matrices with precisely the same domain/conventions.
Reverse comparison is S_ba = S_ab†. Self Grams are positive semidefinite; they
need not be identity on a restricted domain. Do not whiten, truncate eigenvalues,
individually rescale sectors or renormalize orbitals to force native agreement.

For M-electron Slater expansions this matrix induces the physical product-domain
metric on D^M: each determinant pair integrates to det(S_occ). The four-block
factorization still holds without self-orthonormality; evaluate the complete
state's self norm with it. Only then form the normalized coherent interior
amplitude or its squared magnitude, with uncertainty propagated through both
self norms. Squaring four blocks separately loses interference.

An outgoing resonance is not generally an all-space normalizable state. An
interior overlap on D^15 is a declared finite-domain character diagnostic, not an
all-space probability or a proof of unique continuation through unsampled R.
Changing a changes that diagnostic. Its stability must be tested and reported;
an all-space Gaussian integral is not the missing physical exterior solution.

## Positive quadrature and a small oracle

The toy function is explicit, unnormalized and Cartesian:

    chi(x) = product_d (x_d-A_d)^k_d
             sum_p c_p exp(-alpha_p |x-A|^2).

All amplitudes, including any desired primitive/contraction normalizations, must
already be in c_p; complex coefficients are allowed. No amplitude is inferred
from the native export. This differs deliberately from a native spherical CGTO
whose solid harmonic, primitive norm and contraction norm need separate decoding.

With r in [0,a], t=cos(theta) in [-1,1] and phi in [0,2pi), use Gauss-Legendre
rules for r and t, and a periodic trapezoid rule for phi. Weights include r^2,
dr, dt and dphi exactly once. Every weight is positive. If V_a is the pointwise
basis matrix, the discretized metric is V_a† W V_b, so reverse consistency and
self positivity follow even before convergence. These algebraic properties are
necessary, not convergence evidence. Finite-order quadrature need not be
rotation invariant, and localized displaced functions require angular refinement.

`projects.ukrmol_co.gaussian_ball_metric` exposes `CartesianGaussian` and
`gaussian_ball_metric(left, right, radius, center=..., orders=...)`. It returns
an unnormalized complex spatial cross Gram. Bounded point/basis/primitive/work
counts and chunked values make it a toy oracle, not a molecular performance
claim. It does not load saved orbitals, build spin functions, recover vectors,
compute a Hamiltonian or choose an integration order from the desired answer.

## Independent controls and error separation

For s primitives, Gaussian multiplication gives p=alpha+beta and
P=(alpha A+beta B)/p. The ball integral is

    exp[-alpha beta |A-B|^2/p] (pi/p)^(3/2)
        F_ncx2(2 p a^2; df=3, noncentrality=2 p |P-O|^2).

This follows by viewing the normalized product as a three-dimensional Gaussian
with covariance I/(2p); it is an independent distribution-CDF oracle, not the
quadrature above. Complex contracted s overlaps follow by summing c_a* c_b.
Centered Cartesian products have analytic angular moments and incomplete-gamma
radial integrals; odd powers vanish. An independent all-space Cartesian moment
expansion checks the large-domain limit. Tests need displaced and off-axis
centers, s/p/d polynomials, cancellation, complex phase, deliberately wrong
normalization, and separate three-level quadrature and domain families. A
single small refinement difference or a PSD matrix cannot certify accuracy.

Keep these errors separate:

1. **Quadrature:** refine radial and angular orders independently; compare
   analytic controls and three levels, with a stated roundoff floor. Actual
   native exponents/centers require their own resolution test later.
2. **Representation/native conventions:** scalar analytic agreement does not
   identify real spherical harmonic order/phases or the native use of stored
   primitive and shell normalization. A Gram can miss sign/permutation mistakes.
3. **Domain:** compare identical quadrature prescriptions at several physical
   radii after each is converged. Small orbital tails alone do not bound a
   normalized many-electron comparison with large/cancelling coefficients.
4. **State/model:** saved target/continuum products, core-hole/spin conventions,
   missing Q coefficients and electronic-model errors are distinct unresolved
   requirements. No orbital-metric test resolves them.

For genuinely L2 functions on a larger domain, the omitted cross integral obeys
|<f,g>_outside| <= ||f||_outside ||g||_outside. For normalized M-electron states,
the probability of any electron outside D is bounded by the integral of the
physical one-body density outside D (union bound); this needs state coefficients
or a valid coefficient-aware upper bound, not one typical orbital's tail. Then
the omitted state cross amplitude has the analogous Cauchy-Schwarz bound. These
statements cannot be applied to a nonnormalizable outgoing state on all space.

An orbital entry tolerance is not a state-overlap tolerance. For one occupied
M-by-M matrix with row norms <=B and perturbation row norms <=eta, determinant
multilinearity and Hadamard give absolute error <=(B+eta)^M-B^M. Sum with
absolute CI coefficient products, or use a sharper justified contraction bound,
and propagate self-norm errors as well. A nearly singular occupied cross block
does not justify an inverse or an unbounded relative determinant tolerance.

## Native admission: what one future self-metric test would decide

The exported schema provides centers, angular degree, primitive exponents,
contractions, primitive norms, shell norm, addressing and unchanged MO
coefficients. The [source-resolved native contract](co-native-orbital-conventions.md)
specifies the point-value factors and signs, but also identifies a critical
read-versus-generator normalization distinction: the convenience reader
reconstructs all-space CGTO norms. Its output must not be assumed to preserve
raw stored sphere factors. The source's hybrid domain applies to **AO shell
pairs**, not directly to MO sector labels. A uniform physical D metric need not
reproduce the stored native identity, even after the AO scale is restored.

A separately released saved-basis self-metric sentinel is justified only after:

- pinned, read-only convention evidence maps the primitive/shell normalization,
  real solid harmonic scale/order/sign, AO/MO addressing and coordinate frame;
- the saved-reader reconstruction is distinguished from the coefficient
  generator's AO scale, with any continuum restoration tagged and independently
  validated rather than silently reapplied to the old export;
- independent pointwise analytic s/p/d controls test absolute amplitude and
  signs, not only an identity Gram (which is invariant under orbital signs);
- native reference quantities are classified by integration domain; expected
  finite-domain differences are predicted, not fitted or whitened away;
- a frozen three-level quadrature family, radius-family diagnostic, numerical
  error budget, work/memory cap and first-failure owner are provided for the
  **one existing export only**. No second export or other-R data are necessary
  to test this convention/domain gate.

That sentinel could reject an adapter or qualify a self-metric ingredient. It
cannot alone admit a cross-R state comparison. Native spin/phase reconstruction,
four-block molecular contractions and their measured costs remain prerequisites
before any separately authorized full-vector recovery. No qscat promotion or
Rust optimization is implied by synthetic success.
