# Loader-aware one-body metric implementation

This is a project-local Gaussian self-metric ingredient, not a cross-R
many-electron state comparator. Follow the
[native point-value/normalization contract](co-native-orbital-conventions.md)
and [common physical-domain definition](co-common-domain-orbital-metric.md).
All positions/radii are bohr. No Hamiltonian, eigensolver, native program or
electronic optimization is called.

## Interfaces and immutable representation

`native_orbital_metric.SphericalGaussian` contains one center, angular degree,
exponents, explicit complex primitive amplitudes and an AO continuum flag.
Amplitudes already include shell/primitive/contraction factors. The conversion
`shells_from_export` copies post-read amplitudes without interpreting them as
raw file norms, and rejects unsupported nonsequential AO/shell addressing.

`generator_scales(shells, a0)` returns the source-derived diagonal AO-row factors
restoring origin-centered continuum self norms on the fixed generation sphere.
The original export and C stay unchanged; record the factors separately and
apply them to the AO metric or C rows exactly once. This operation is neither
Gram whitening nor physical-domain-dependent orbital renormalization.
Nonpositive or numerically tiny sphere norms are a conditioning rejection, not
permission to clip them. Keep the same factors when varying the physical radius.

`solid_harmonics` evaluates the source-resolved real solid harmonics through
SciPy's normalized complex spherical harmonics with explicit Condon–Shortley
phase removal. Independent Cartesian coefficients follow from the explicit
Legendre polynomial. Absolute s/p/d order/sign controls are required; self Gram
agreement alone cannot test phases.

## Two independent all-space algorithms

For two primitives at A,B, Gaussian multiplication gives p=alpha+beta,
P=(alpha A+beta B)/p and prefactor exp(-alpha beta |A-B|²/p).

`shell_overlap(..., method="moments")` expands the real solid harmonics into
Cartesian polynomials. Its one-dimensional moments integrate powers of x-A
and x-B against exp(-p(x-P)²) using integration-by-parts recurrences, starting at
sqrt(pi/p). The three Cartesian directions separate, then primitive contractions
and angular coefficients are applied. Same-center angular orthogonality and
analytic radial gamma functions supply an exact shortcut.

`method="hermite"` instead evaluates the spherical-harmonic point values on a
positive tensor Gauss-Hermite rule centered at P, scaled by 1/sqrt(p). This
integrates every polynomial of degree l_a+l_b exactly when each order n obeys
2n-1 >= l_a+l_b. The implementation uses the minimum such n plus zero, one or
two extra orders. Three-level agreement here tests polynomial exactness and
roundoff, not an empirically decreasing discretization error. This Gaussian
product grid resolves narrow displaced core functions without a global real-space
mesh attempting to hit every Gaussian peak.

`all_space_metric` assembles a bounded AO self Gram. Raw reverse shell-pair
comparisons are tested separately from the Hermitian assembly; no eigenvalue
clipping, renormalization or fitted symmetrization is used to pass a gate.
Caps are 128 shells, 512 AOs, 256 primitives total, l<=10 and a precomputed
point-function work budget. These are admission limits, not a molecular cost
forecast; costs must be measured for each actual basis.

## Ball, hybrid domain and omitted-tail witness

Origin-centered equal-l continuum pairs have diagonal-in-m ball integrals with
lower incomplete gamma radial factors. `radial_overlap` also offers independent
positive Gauss-Legendre radial quadrature. In this helper a finite radius is
measured from the shared shell center; native continuum callers enforce that
center is the coordinate origin.

`uniform_ball_witness(shells, G_all, a, scales)` returns the restored native-style
hybrid AO matrix H(a), plus an entrywise upper bound E on its difference from
the uniformly integrated physical ball matrix. Continuum-continuum entries use
the exact ball integral; target-containing entries retain all-space values and
have a **nonzero bounded uncertainty**, not assumed zero tails.

For shell-center distance |A|, points outside the common ball imply distance
from the shell center greater than d=max(0,a-|A|). Each primitive's exterior
norm on that larger radial region is known from the upper incomplete gamma.
The contraction-aware triangle bound is

    T_s(a) <= [sum_p |amplitude_sp| sqrt(J_sp(d))]^2,
    J_sp(d) = 2pi/(2l+1) Gamma(l+3/2) (2alpha_sp)^(-l-3/2)
              Q(l+3/2, 2alpha_sp d²).

Q is the regularized upper incomplete gamma. This is a bound for actual Gaussian
basis functions, not an outgoing resonance's all-space probability. It is
conservative even for cancelling contractions. AO entries with at least one
target have E_ij=sqrt(T_i T_j), with generator scales included; CC entries have
E_ij=0 because their finite-ball values were computed explicitly. Floating-point
evaluation error must still be budgeted; this is not interval arithmetic.

For complex MO coefficients, the propagated tail bound is |C|^T E |C|. Thus

    G_MO(a) ≈ C† H(a) C,
    |G_MO(a) - exact_uniform_ball_MO(a)| <= |C|^T E C_abs.

Here C_abs=|C| in the second expression. Do not replace this by a typical orbital
tail or compare different radii after changing the generator normalization.
At a=a0 compare the hybrid prediction with the stored native overlap column;
separately inspect the uniform-domain uncertainty and physical PSD/reversal.

## Admission boundary

Use generated absolute points, analytic polynomial moments, independent Hermite
rules, incomplete-gamma/radial quadrature, complex covariance, direct finite-ball
controls and malformed/work-cap rejections before saved-basis use. A saved
sentinel additionally needs exact export/source pins, all-space read-shell norm
checks, fixed-a0 restoration, three Hermite/radial levels, a radius family,
propagated target-tail bounds, actual time/RSS and original exits. Reject source,
normalization, conditioning or cost failures without widening a ladder or hiding
them by whitening.

Even a passing one-R self metric says nothing by itself about native coupled-CSF
spin/phases, target/removal contractions, cross-R state identity, missing interior
coefficients, electronic accuracy or unsampled branch continuation. It is not
authorization for native recovery or promotion to qscat.
