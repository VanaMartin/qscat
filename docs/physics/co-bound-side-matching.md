# Finite CO bound-side matching: design, not a qualified solver

## Purpose and scope

Test for a normalizable fixed-nuclei N+1 state below the native neutral ground
threshold using the same finite inner spectrum and closed-channel exterior as
the local resonance calculation. This is a proposed project-level Python
capability, not implemented/validated here or promoted to QSCAT. It does not
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

## Candidate API and implementation boundary

A small project adapter could take incident Hartree energy, native thresholds/
multipoles, matching radius, a prescribed exterior domain and numerical caps,
and return Y_dec plus tail/endpoint/residual diagnostics. A separate bounded
root driver would consume D/W/E0 and a fixed **negative-energy** interval,
returning an energy, null/residual witnesses and a numerical uncertainty estimate.
Candidate interfaces (not implemented or frozen public APIs):

    decaying_dtn(incident_hartree, thresholds_rydberg, c2, c3,
                 radius_bohr, exterior_settings, caps) -> (y, diagnostics)
    find_bound_state(poles_hartree, amplitudes, ground_hartree,
                     negative_interval_hartree, decaying_map, caps) -> result

It must reject nonnegative/boundary-threshold inputs and nondecaying solutions.
No input should trigger QC, native exports or implicit downloads.

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

## Validation required before a molecular release

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
