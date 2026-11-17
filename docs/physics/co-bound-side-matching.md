# Finite CO bound-side matching: project-stage Python prototype

## Purpose and scope

Test for a normalizable fixed-nuclei N+1 state below the native neutral ground
threshold using the same finite inner spectrum and closed-channel exterior as
the local resonance calculation. The project-level Python prototype
`projects.ukrmol_co.bound_matching` has analytic/synthetic regression controls
and separately frozen finite-CO applications in the
[experiment catalogue](../../projects/ukrmol_co/experiments.md). The API alone
is not a generally qualified molecular solver or promoted QSCAT capability. It does not
establish electronic convergence, unique multi-R identity or attachment dynamics.
The [local representation](co-results-and-rejected-routes.md#local-representations-accepted-scope)
and its real-R domain remain unchanged.

## Units, boundary condition and finite eigenproblem

Let epsilon be incident Hartree energy, E0 the geometry's native neutral ground
total, D the diagonal exported inner eigenvalues in Hartree, and W the real
pole-major/channel-minor stored amplitudes. At zero Bloch parameter:

    T = E0 + epsilon,
    Q(T) = W.T @ inverse(D - T I) @ W.

Q has units of bohr. Stored W already includes the native 1/sqrt(2) factor:
do not add a half, radius division or complex conjugation. See the
[native convention](co-native-boundary-conventions.md) and
[outgoing target definition](co-effective-resonance-target.md#reference-hamiltonian-energy-and-sheet).

For epsilon below every physical detachment threshold, use decaying physical-sheet
solutions in **all channels**, not the resonance sheet continued through zero.
With native incident thresholds Delta_Ry and radial multipoles C2/C3:

    U''(r) = [diag(Delta_Ry - 2 epsilon) + C2/r^2 + C3/r^3] U(r),
    kappa_j = sqrt(Delta_Ry,j - 2 epsilon) > 0,
    Y_dec(a; epsilon) = U'(a) @ inverse(U(a)).

The leading large-r derivative is -diag(kappa); the coupled inverse-power tail
and its error must still be controlled. The factor two converts incident Hartree
to the native Rydberg radial equation. Y_dec is in inverse bohr and all finite
closed-channel couplings remain present. No electron-polarization term absent
from the original model may be silently appended.

Away from inverse-coordinate singularities, matching is

    [I - Q(T) Y_dec(epsilon)] v = 0.

For root verification, a complementary finite inner-coordinate problem is

    F(epsilon) c = 0,
    F = D - (E0 + epsilon) I - W Y_dec(epsilon) W.T.

For real symmetric Y_dec, F is real symmetric. The determinant identity
`det(F)=det(D-T I)*det(I-Q Y_dec)` holds where D-T I is invertible, but raw
determinant signs are not acceptance tests. F avoids a spurious "root" caused
only by Q diverging at an inner eigenvalue; it does not remove poles of Y_dec.
Check residuals and a simple isolated eigenvalue, not determinant magnitude alone.

If U(a) is singular, neither inversion nor a sign change across that chart pole
qualifies a bound state. An inverse-free boundary residual can use

    (D - T I) c - W U'(a) b = 0,
    W.T c - U(a) b = 0.

Column scaling of U and its derivative must remain consistent. A singular chart,
unresolved nullspace or missing stable representation stops the bounded attempt;
it is not permission for an open-ended matching-framework project.

## Implemented API and boundary

### Finite-radius D/N support

`projects.ukrmol_co.bound_box.finite_box_dtn` takes the same incident Hartree
energy, native Rydberg thresholds/C2/C3 and inner radius as `decaying_dtn`, with
an explicit `endpoint_bohr=R` and `boundary="dirichlet"` or `"neumann"`.
The radial equation is

    U'' = V_epsilon(r) U,
    V_epsilon = diag(thresholds - 2 epsilon) + C2/r^2 + C3/r^3.

For Neumann data, integrate `Y'=V_epsilon-Y Y` inward from `Y(R)=0`.
For Dirichlet data, integrate the inverse chart `Z=U(U')^-1` from `Z(R)=0`
using `Z'=I-Z V_epsilon Z`, then solve `Z(a) Y(a)=I`. Matrix multiplication
order matters for coupled channels. These are prescribed separate problems,
not automatic chart alternatives. Invalid inputs, asymmetric maps, singular or
ill-conditioned inverse charts, divergent integration and work caps stop the
calculation without switching boundaries or projecting away transpose defects.
Inputs are copied and checked for changes; this does not preempt arbitrary
concurrent mutations, so external ownership remains required.

These maps are tagged `finite_box=True`, `decaying=False` and with their boundary.
`find_bound_state(..., finite_box_boundary="dirichlet" | "neumann")` requires
matching tags and returns **finite-box-root-only**. Its default path still rejects
explicit nondecaying/finite-box maps. Neither mode supplies a total uncertainty.
Legacy untagged callbacks remain accepted by the default API; callers, not the
driver, are responsible for not misrepresenting a finite box as a physical tail.

For exact D/N bracketing, use the **same fixed finite inner model** in all three
problems. Its quadratic form has inner term `c.T (D-E0) c`, exterior term
`1/2 integral[|u'|^2 + u.T V_0 u] dr`, total norm `|c|^2+integral |u|^2 dr`
and trace constraint `u(a)=sqrt(2) W.T c`. Variation gives precisely
`F=D-(E0+epsilon)I-W Y W.T`, without another radius or sqrt(2) factor.
Dirichlet truncation at R restricts the form domain and gives ordered upper
bounds. Neumann splitting at R enlarges it and gives a direct sum of the finite
box and a detached Neumann half-line. **Only if the latter is nonnegative**
do its eigenvalues not contaminate the ordered negative spectrum: then

    epsilon_N,k(R) <= epsilon_infinite,k <= epsilon_D,k(R).

Matching the same ordered index also requires regular DtN charts and a common
self-adjoint operator/state ordering, not merely nearby numerical roots. A
nonnegative exterior throughout `[a,infinity)` is a sufficient (stronger than
necessary) condition. A negative Neumann box level alone proves no binding;
a short Dirichlet box can miss a shallow physical bound level altogether.

`nonnegative_exterior_witness(thresholds, C2, C3, a)` supplies a conservative
sufficient test, with nonnegative thresholds T. In exact arithmetic,

    r^2 V_0(r) >= a^2 T + C2 + C3/r,       r >= a.

The right-hand side is a convex combination of `B0=a^2 T+C2` and
`B1=B0+C3/a`. Nonnegative Gershgorin row lower bounds at both matrices suffice.
The implementation requires exact input symmetry and applies a floating safety
allowance `64*(n+1)*eps*row_sum(a^2 T+|C2|+|C3|/a)`, retaining uncancelled
input scales. Free exact-zero rows are allowed.
This allowance is **not interval arithmetic or certified roundoff**, and failure
is inconclusive, including for positive matrices that are not diagonally dominant.
The witness covers only the supplied finite-channel native multipole operator;
it adds no polarization tail or omitted-channel/electronic error bound.

An explicitly selected alternative, `positive_exterior_witness`, tests the same
two endpoint matrices without requiring diagonal dominance. It does not replace
or automatically retry the Gershgorin API. For each symmetric endpoint B, compute
a full real eigendecomposition `(Q, lambda)` and reconstruct

    R = B - Q diag(lambda) Q.T,
    delta = ||Q.T Q - I||_infinity,
    rho = max(||R||_infinity, ||R||_1).

If `lambda_min>0` and `delta<1`, the square Q is nonsingular and
`lambda_min(Q diag(lambda) Q.T) >= lambda_min*(1-delta)`. The spectral norm of R
is at most `sqrt(||R||_1 ||R||_infinity) <= rho`. Thus, in exact arithmetic,

    lambda_min(B) >= lambda_min*(1-delta) - rho.

Require this lower bound to exceed a floating safety allowance based on
uncancelled input and reconstruction row scales, and report all terms. Zero or
unresolved margins reject rather than inheriting the free-zero exception of the
Gershgorin API. This is an algebraic positive-form test with numerical padding,
**not interval-certified roundoff**. Small individual eigenpair residuals alone
do not certify the full minimum eigenvalue; the full reconstruction and
orthogonality defects matter. Convexity then gives positive B(1/r) and the same
nonnegative exterior form needed for detached-Neumann D/N ordering. Negative or
inconclusive endpoints stop this criterion; no model term or tolerance is changed.
Analytic rotated spectra, non-diagonally-dominant positive matrices, corrupted
decompositions and explicit full-channel synthetic controls test this alternative.
Its availability does not release a molecular calculation.

Analytic tanh/coth and modified-Bessel controls, independent second-order coupled
collocation, three-level tolerance/radius refinements and a shallow finite-inner
well test exercise these APIs. D/N ordering is an exact quadratic-form statement,
not an ODE error enclosure: molecular acceptance still needs separate controlled
numerical padding/refinement, a resolved negative interval, localization and
state witnesses. No molecular matching prescription or execution follows from
these generated tests. This remains project-stage Python, not promoted QSCAT code.

`projects.ukrmol_co.local_pole_inputs.load_closed_boundary_model(source)` loads
only the same neutral 48-target/298-channel/2048-pole/a=18-bohr native inputs,
without opening a K/S grid. It retains six zero ground thresholds, strictly
positive finite excited thresholds and the native target/threshold energy identity.
Unlike the unchanged positive-grid `load_boundary_model`, it imposes no 0.3675-Ry
excited-threshold floor. This is input support for **strictly negative incident
Hartree energies**, not qualification of a positive grid or molecular bound state.
Both loaders share the D/W/multipole decoder and preserve native normalization.
The decaying solver still enforces its own energy domain and numerical gates.

A positive search must not reinterpret that negative-only policy.
`load_positive_boundary_model(source, positive_interval_hartree)` separately
qualifies an explicit finite ordered interval with `0 < lower < upper` and every
excited native Rydberg threshold strictly above `2*upper`. Threshold touching
fails. The neutral 48/298/2048/a18 view, six zero ground thresholds and D/W
normalization are unchanged; this is not a variable-open-channel exporter or
qualification of the historical wider positive grid. Real and complex consumers
must enforce the declared real-energy projection for every evaluation.
The explicit `local-cas11-six-open-window-v1` pole-runner view requires its search
box within that interval and exactly its own two native input files. It retains
the mandatory measured two-point real-axis guard, canary, root, contour and
residue gates; input-policy success alone neither releases a calculation nor
qualifies a resonance. Historical model views and the default positive-grid
floor remain unchanged.

That threshold window classifies open/closed channels; it does **not** guarantee
complex-contour exposure. On the continued outgoing/decaying sheet, every chosen
angle must also satisfy `Im(exp(i theta) k_i) > 0` for every channel, since
`exp(i k_i exp(i theta) t)` must decay along the rotated ray. Check the joint
sheet/contour/search domain, not just seeds and real-axis controls: a rectangular
complex box can include points outside an exposed wedge, and a root-search trial
can reach them. The existing outgoing map rejects such a point before integration.
That rejection qualifies neither a root nor absence of a resonance. Revising the
domain, angle or seed after a frozen first-failure stop requires a separate
experiment decision; it is not grounds to remove the exposure guard.

Exposure alone also does not control **inverse-Jost row contamination**. The
[exact dual-transport control](inverse-jost-dual-transport-control.md) distinguishes
outgoing-subspace decay from backward attenuation of an unwanted closed dual in
an open row. With native radial energy `E=x-i eta`, write the ground-open momentum
`k_o=a-i b` and a closed momentum `k_c=-c+i d`, where all four quantities are
nonnegative. Since `a*b=c*d=eta/2`,

    Im(exp(i theta) (k_c-k_o)) = (d+b) (cos(theta) - a/d sin(theta)).

A positive value damps that asymptotic relative error during backward transport;
a negative one amplifies it. This is a sensitivity test, not an error bound on
the coupled tail. For a whole rectangle `0<x0<=x<=x1`, `0<=eta<=eta1` below the
lowest closed threshold `T`, sufficient bounds are
`a<=sqrt((sqrt(x1*x1+eta1*eta1)+x1)/2)` and `d>=sqrt(T-x1)`.
Thus `sqrt(T-x1)*cos(theta)-a_max*sin(theta)>0` is a conservative joint-domain
condition, alongside outgoing exposure. Apply it separately to the actual
real-axis controls and every complex-amplitude evaluation. A contour suitable
for a low-energy pole need not be suitable for higher-energy real controls.
An exponential change of amplitude variables cannot remove finite-tail
contamination already present in the initial data. Keep the full flux/Wronskian
gates, including closed columns; never project or whiten them into agreement.

    decaying_dtn(incident_hartree, thresholds_rydberg, c2, c3, radius_bohr,
                 *, endpoint_bohr, tail_order=8, rtol=1e-11, caps=...)
    find_bound_state(poles_hartree, amplitudes, ground_hartree,
                     negative_interval_hartree, decaying_map, *, mode_index=0,
                     energy_tolerance_hartree=1e-12,
                     residual_tolerance_hartree=1e-8,
                     isolation_gap_hartree=1e-6, caps=...)

`decaying_dtn` integrates the real-axis Riccati equation inward with DOP853 and
returns `(Y, diagnostics)`, including raw transpose defect, tail-order/endpoint
settings and RHS cost. It shares the inverse-power recurrence with the outgoing
code, not that code's sheet/angle interface. It rejects nonnegative/complex
incident energies, invalid dimensions/symmetry, unresolved tail/DtN charts and
time/RHS exhaustion. No projection repairs a failed symmetry check.

The default `tail_representation="fundamental"` retains that original recurrence
and behavior. An explicitly selected `"log-derivative"` alternative instead
expands the physical-sheet Riccati solution directly:

    Y(r) = -K + sum_(q=1)^N A_q/r^q,    K=diag(kappa), A_0=-K,
    (kappa_i+kappa_j) (A_q)_ij = [-(q-1) A_(q-1) - V_q
                                      + sum_(p=1)^(q-1) A_p A_(q-p)]_ij,
    V_2=C2, V_3=C3, otherwise V_q=0.

This follows by equating powers in `Y'=K^2+C2/r^2+C3/r^3-Y^2`; in particular
`A_1=0` and `A_2,ij=-C2_ij/(kappa_i+kappa_j)`. Coefficients are symmetric in
exact arithmetic because each product is paired with its transpose. The code
does not symmetrize coefficients, initial data or the propagated result: raw
symmetry checks still apply. It additionally requires a negative-definite initial
Y to resolve the decaying branch. Neither construction repairs a previously
failed map. A new choice needs its own frozen endpoint/order/tolerance family.

Unlike the free p-wave fundamental solution, whose inverse-r series terminates,
its logarithmic derivative has an infinite expansion. Low kappa*R or excessively
high orders may therefore be poor initial data. Move to prospectively justified
farther endpoints and test convergence; the last term is not a tail enclosure.
Analytic free/p-wave and independent second-order coupled BVP comparisons, with
three orders/endpoints, qualify the alternative only in their tested regimes.

`find_bound_state` brackets one prospectively selected **ordered eigenvalue**
of F using Brent's method. It checks the raw matrix residual, local isolation,
unchanged inner inputs and identical-energy map replay. It returns the energy,
inner coefficients, numerical diagnostics and `qualification="finite-model-root-only"`.
Neither a successful bracket nor this local index proves global root count or
state assignment. Callback maps must enforce their own work caps; cooperative
before/after-call time checks cannot preempt an arbitrary blocked callback.

`BoundCaps` defaults to 30 seconds, 20000 RHS calls, 96 map evaluations, at most
32 channels and 256 inner poles. These are synthetic-stage limits, not a
qualified 298-channel/2048-root execution envelope. There are no molecular
readers, implicit downloads, native calls or recomputation fallbacks.
DtN-chart failures stop; the inverse-free representation described above is not
implemented. `total_uncertainty_hartree` remains **unset**: search tolerance,
residual and the last tail term are not a rigorous total error enclosure.

`outgoing_map.sheet_momenta` already selects decaying roots for closed channels
inside a negative reference interval, but the full all-closed bound workflow is
not thereby validated. The current map also requires a positive contour angle;
do not pass angle zero or remove its threshold guard to pretend it is a qualified
real-axis bound solver. Test any reused rotated construction against an independent
real-axis decaying boundary-value solution before trusting it.

The resonance `contour_witness` requires positive width and nonzero open-channel
residue. A bound solution has no open detachment channels, so it needs a separate
bound-state acceptance contract, not a bypass that reports the resonance tests
as passed. Real-axis K/S guards cannot supply a below-threshold scattering oracle.

## Synthetic controls and remaining molecular gates

`test_bound_matching.py` exercises the controls below using derived toy inputs,
not measured fixtures. Three tolerance levels test exact decaying p waves down
to incident energy -5e-9 Hartree; three tail orders/endpoints test a two-channel
map against a separate real second-order boundary-value solve. The latter shares
the finite differential operator but not the Riccati propagation or tail
recurrence, and uses a farther leading-Robin endpoint. A shallow radial square
well is checked against its analytic secular equation at three inner-spectrum
sizes, including the expected omitted-tail trend. Rotated two-mode closed-form
roots and opposite-order coupled boundary matching check the F construction.
Subcritical/short-box and singularity controls retain rejected outcomes.

These controls do not establish a CO binding sign, molecular localization,
full-size cost, omitted-spectrum error, or the adequacy of a molecular exterior.
Before a molecular release, retain these controls and additionally qualify the
actual input/state/domain and uncertainty contract:

- Free closed-channel DtN against exact -kappa, and nonzero-l decaying solutions
  against analytic logarithmic derivatives, including shallow negative energies.
- An analytically solvable finite well with a genuine bound state and an
  unbound/virtual-state negative control. Verify the physical sheet and binding
  sign; do not confuse a finite-box continuum level with binding.
- Coupled-channel comparison against an independently assembled second-order
  boundary-value/eigenproblem, plus convergence in tail order, endpoint and
  integration tolerance. Agreement on the same tail approximation is not an
  independent tail-error bound.
- Q/F or inverse-free matching agreement away from singularities; synthetic
  inner-pole and DtN-chart singularities must not pass as states. Include lost
  isolation, unresolved sign, timeout and changed-input rejection controls.
- Resolved energy error within the separate 0.1-meV same-model budget and a
  binding-sign margin exceeding justified numerical uncertainty. A width floor
  or failed outgoing search cannot set Gamma to zero and declare binding.
- Physical decay/localization and exterior-domain sensitivity of the finite
  solution. Missing full inner CI vectors limit wavefunction interpretation;
  D/W coefficient similarity across different geometries is not an N+1 overlap.

All molecular inputs, state-assignment rules, effort/compute caps and the tested
owner require a separately frozen release. A successful local bound candidate
and an earlier resonance do not alone prove unique continuation between them,
absolute attachment energy or a stable nuclear anion.

## Norm from the physical-sheet energy derivative

For a regular all-closed map, normalize its exterior fundamental solution by
`U(a,E)=I`. Differentiating `U''=(V0-2E)U` in incident Hartree energy and
integrating the Lagrange identity with fixed boundary value gives

    dY(a,E)/dE_Ha = 2 integral_a^infinity U(r,E).T U(r,E) dr.

The boundary term at infinity vanishes for the decaying solution. The boundary
term at a is the derivative of Y, since `dU(a,E)/dE=0`. Thus with `v=W.T c` and
native `u(a)=sqrt(2)*v`, the exterior probability norm is
`v.T (dY/dE_Ha) v`; the total norm is `c.T c` plus this quantity. Equivalently,
`-c.T F'(E)c` equals that total. Neither an extra factor two nor a radius factor
belongs in the final expression. A negative exterior norm rejects the diagnostic.

Finite differences of independently converged Y maps can estimate this norm;
compare at least three energy steps and retain their differences and map-error
sensitivity. They are not interval derivatives or rigorous norm enclosures.
Free channels give exactly `dY/dE=diag(1/kappa)`; nonzero-l analytic radial
solutions provide a quadrature comparison. The inner probability fraction
`(c.T c)/(c.T c+norm_exterior)` is a useful finite-model localization witness,
not an electronic-state assignment or an N+1 overlap across geometries. If it is
unresolved, no localization claim follows from boundary amplitudes alone.

## Comparing transverse bound counterparts

A same-geometry B1/B2 comparison first requires the same target energies, neutral
zero, geometry, retained space and exterior model. Compare qualified energies
with their numerical uncertainty intervals and compare localization estimates
under the same normalization. Reuse the accepted component's saved evidence;
do not repeat its exterior solve merely to compare the missing counterpart.

Raw inner coefficients or ordinal channel profiles across B1/B2 are not invariant
state overlaps: target components, real-harmonic signs and channel ordering can
differ. A limited trace-character check aggregates `|W.T c|^2` by `(l,abs(m))`,
both for all channels and separately for ground-target channels, then normalizes
each profile. This is invariant under channel permutations and phases, but loses
target-resolution information. Agreement is a necessary compatible-counterpart
witness, not a proof of unique electronic identity or continuation across R.
Never reinterpret those surface weights as full inner/exterior probabilities.
