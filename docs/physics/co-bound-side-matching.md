# Finite CO bound-side matching: synthetic Python prototype

## Purpose and scope

Test for a normalizable fixed-nuclei N+1 state below the native neutral ground
threshold using the same finite inner spectrum and closed-channel exterior as
the local resonance calculation. The project-level Python prototype
`projects.ukrmol_co.bound_matching` is tested on analytic/synthetic cases only;
it is not a qualified molecular solver or promoted QSCAT capability. It does not
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

`projects.ukrmol_co.local_pole_inputs.load_closed_boundary_model(source)` loads
only the same neutral 48-target/298-channel/2048-pole/a=18-bohr native inputs,
without opening a K/S grid. It retains six zero ground thresholds, strictly
positive finite excited thresholds and the native target/threshold energy identity.
Unlike the unchanged positive-grid `load_boundary_model`, it imposes no 0.3675-Ry
excited-threshold floor. This is input support for **strictly negative incident
Hartree energies**, not qualification of a positive grid or molecular bound state.
Both loaders share the D/W/multipole decoder and preserve native normalization.
The decaying solver still enforces its own energy domain and numerical gates.

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
