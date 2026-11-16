# CO scalar reduction: an outgoing-pole target, not a full scattering model

## Decision and claim boundary

Use **one continuously assigned outgoing resonance pole per Pi sector** as the
scalar factory's prospective electronic target. This is a physically defined
reduction: it preserves a coupled Hamiltonian's resonance energy and total decay
width, not its channel amplitudes. It replaces neither the native exterior nor
the full physical S matrix. The bare full-S route remains
[rejected](../../projects/ukrmol_co/DIRECT_SCATTERING_SCREEN.md).

Compatible finite CAS11/nine/2048 B1/B2 pairs now support a prospectively tested
scoped local representation; the [current handoff](../../projects/ukrmol_co/CONTINUATION.md)
links its immutable data and validation membership. Separate stretched CAS12
results remain distinct. Common N+1 identity and global/factory-model qualification
are not supplied by this local representation.
The [two-goal plan](../../projects/ukrmol_co/RESEARCH_PLAN.md)
recommends local finite-model collection, without releasing a job or fit.
A pole-only fit plus the reduced validation below would
**not complete the campaign's full held-out scattering goal**. In particular,
pole equivalence is not equivalence of vibrational-excitation cross sections.

## Reference Hamiltonian, energy and sheet

The saved feasibility reference is CAS11/nine-A2 at R = 1.9, 2.1323, 2.5 bohr,
separate B1/B2 sectors, radius a = 18 bohr, 298 boundary channels and the accepted
2048 exported inner eigenpairs. Six ground-target channels are open on the saved
grid, with l = [1, 2, 3, 3, 4, 4]. The others must not simply be dropped: their
virtual coupling contributes to the pole. These files do not define a CAS12
curve, a channel-converged model or a complete inner spectrum.

At each R define the **incident Hartree pole**

    z(R) = E_total,pole(R) - E_ground,native(R)
         = E_res(R) - i Gamma(R)/2,       Gamma > 0.

Gamma is the full width; it is not the half width, an eigenphase-fit parameter,
an inner-region R-matrix eigenvalue or a generic real-axis peak width. Preserve
both native total energy and incident energy in any later evidence record.
Requested energies in the saved decks map to incident Hartree by
`E_inc = 0.0735 * E_requested_eV / 2`. For physical eV reporting and tolerances
use `qscat.units`, not that historical requested-grid conversion.

For channel thresholds Delta_j in incident Hartree, continue
`k_j = sqrt(2*(z - Delta_j))` from the physical real-energy interval: outgoing
open waves have positive real k there, continued onto the resonance sheet with
negative imaginary k; closed waves continue the decaying branch k = +i kappa.
All six degenerate ground-target channels cross the cut together. Do not choose
square-root signs afresh at every point or call a lower-half-plane principal
square root a complete sheet specification. Work in a local domain not crossing
a threshold or another singularity; a threshold collision ends this reduction's
scope rather than triggering a single-channel Wigner-law extrapolation.

With native Q in bohr and outgoing/decaying logarithmic derivative Y_out in
inverse bohr, the matching condition is

    det[I - Q(E_total)*(Y_out(z) - beta I)] = 0,   beta = 0,
    Q = W.T @ diag(1/(E_n - E_total)) @ W.

There is no extra half, radius division or conjugation of W. The
[native convention derivation](co-native-boundary-conventions.md) supplies the
normalization and energy reference. This determinant is a definition, not a
qualified pole finder. Later evidence must exclude matching-coordinate
singularities/cancellations and demonstrate a nonzero S-matrix residue at a
simple outgoing solution. A zero of a numerically scaled determinant alone is
insufficient. Null-space residuals and an independent pole construction are
needed, with conditioning and numerical uncertainty reported.

## What survives the reduction

Retain **in the native reference problem** the supplied channel thresholds,
closed-channel coupling and coupled inverse-square/inverse-cube exterior CF
matrices. They are already Rydberg-scaled in the radial equation; do not double
them again. These exports do not establish a complete electron polarization
`-alpha/(2 r^4)` exterior. Numerical pole qualification of the finite 2048-term
problem would still not bound omitted-spectrum or electronic-model error.

The scalar potential absorbs their net shift of E_res and Gamma. Do **not** add
the native polar exterior to an already dressed pole fit as an unaccounted
second correction. Nor is the scalar Gaussian's free exterior a claim to
reproduce CO's physical polar threshold behavior. `ell=1` is the proposed
effective resonance coordinate, not an assertion that all native open channels
are p waves. Its adequacy and potential nonuniqueness remain testable assumptions.

Deliberately omitted from the scalar target: the residue matrix and its channel
branching, background phases, angular mixing, interference, other resonances,
and the energy dependence of off-shell couplings. A fixed-R pole curve does not
determine nonadiabatic/nuclear derivative couplings. None may be reconstructed
from an arbitrary energy-by-energy eigenchannel or declared negligible by fiat.

B1/B2 are independent checks of the Pi pair, not extra geometry nodes. Require
agreement at same-root numerical budgets before treating them as one physical
curve; do not average away disagreement or double-count their poles in a
single-sector determinant. Assignment across R must use continuous pole/residue
subspaces in an identified common channel frame, in addition to pole proximity.
If branches exchange, overlap is ambiguous, a pole is not simple or there is no
isolated continuation domain, report unset/reject rather than reseeding a smooth
width curve. Phase-fit candidates can locate a search domain, not certify it.

## Connection to the existing factory

### Saved local evaluator (not a factory model)

`projects.ukrmol_co.local_representation.load_representation(capsule)` authenticates
the immutable local-delivery capsule by its pinned full SHA256 before reading
selected JSON members in memory; it never extracts or executes archived source.
There is no implicit download or quantum/scattering call. The returned
`LocalRepresentation` accepts finite **real** scalar or array radii in the closed
interval **1.9–2.5 bohr** only (no clipping, extrapolation or complex/ECS extension).
Its methods preserve input shape and return energies in Hartree:

- `v0(R)`: the saved degree-6 5Z Chebyshev series, evaluated at
  `x=(R-2.2)/0.3`; zero is its own molecular total at **R=2.1323 bohr**,
  not the minimum, native target total or separated fragments.
- `er(R)` and `full_gamma(R)`: the saved not-a-knot piecewise cubics in
  `Er` and `log(full Gamma)`, with width exponentiated after evaluation.
  Stored powers descend from cubic to constant in `R-left_knot`; coefficients
  are used literally, not reconstructed by another interpolation/fit.
- `incident_pole(R) = er(R) - 0.5j*full_gamma(R)`.
- `vres_hybrid(R) = v0(R) + er(R)`: an explicitly mixed-method view;
  it is not an absolute correlated anion total.

The loader requires the separate post-validation acceptance sidecar. It does not
rewrite or reinterpret the original frozen candidate's historical false acceptance
flag. Eight development nodes and two prospective validation nodes retain distinct
roles; QZ sensitivity and separate MR endpoints are not substituted into this API.
The CLI renders saved neutral/resonance curves and validation residuals outside
Git, with separate development/validation markers and explicit hybrid labelling.

Synthetic controls check known polynomial/log-width functions, stored-power
orientation, C2 continuity and domain rejection. Opt-in authenticated saved-data
regression must reproduce all 26 primary 5Z predictions and both frozen resonance
predictions within **1e-12 Hartree/eV**, and retain the original **5 meV** position
and **max(1% predicted full width, 1 meV)** sentinel gates. These test arithmetic
and the frozen validation design, not a uniform error bound, electronic accuracy,
unique N+1 continuation, global potential, full S matrix or nuclear dynamics.
No hot path or promotion is claimed; this remains a Python project evaluator.

Let V0_ref(R) be the correlated neutral curve on a declared common zero,
distinct from the native absolute target energy. A local dataset may use a
fixed molecular reference; a global dissociation-constrained factory model
requires a qualified fragment zero. The eventual adapter is

    ResonanceTarget.v_ion(R) = V0_ref(R) + Re z(R),
    ResonanceTarget.gamma(R) = -2 Im z(R).

`fit._pole_target_fn` then targets `v_ion - model.v0 - i gamma/2`. A neutral-fit
error therefore shifts its incident target; report that error rather than
silently mixing the native absolute ground energy with the correlated curve.
No numerical values or `ResonanceTarget` object are supplied by this note.

Existing constraints are real blockers, not reasons to invent targets:

- `fit_resonance` requires at least **four successfully tracked geometry nodes**.
  Three anchors, B1/B2 duplication or spline-generated pseudo-data do not suffice.
- Width classification uses `Tolerances.gamma_floor`; a below-floor positive-
  energy slice is not a resolved resonance. Future settings must resolve the
  declared CO width budget, not inherit defaults that mask it. Do not represent
  an unresolved width as zero or weaken the scientific gate to match a default.
- `ea`, `R_inf` and `tail` require justified dissociation data. The factory's
  `polarisation_tail` is an **ion--atom tail in nuclear R**, not the native
  electron exterior in r. A local 1.9–2.5 bohr dataset cannot justify those
  asymptotic inputs. No CO asymptote is guessed to make the API run. Store local
  data first; a local-fit interface must explicitly defer global constraints.
- `Target.eigenphase` is reserved metadata, not an implemented T2 scattering
  loss. `CouplingTarget` is not determined by Gamma(R) alone.

## Falsifiable validation without a fitted-pole tautology

The primary reduced validation is **prediction of an excluded geometry's native
outgoing pole**, not agreement with a phase curve manufactured from a training
pole. Prospectively reserve the equilibrium R = 2.1323 bohr *pole target* from
factory fitting, parameter selection, smoothing and branch-seed selection.
Historical equilibrium real-axis results are already known: this is not a claim
of globally blind data. Qualification-method controls should use analytic cases
and the compressed reference, not tune against the excluded equilibrium pole.

The original minimum design required four genuine training nodes bracketing
that exclusion and one validation node. The delivered local representation now
has eight development nodes and two different prospective validation nodes;
this minimum is no longer a missing-geometry queue. The equilibrium outgoing
pole stays reserved and needs its own later release. A neutral input at
equilibrium may be supplied openly: this is conditional resonance validation,
not prediction of the neutral curve.

Before seeing the validation pole, freeze the scalar parameters, numerical
settings and assignment rule. Compare its predicted incident pole with an
independently qualified native one in each sector, with these **pointwise** gates:

- |Delta E_res| <= 0.05 physical eV;
- |Delta Gamma| <= max(0.05 Gamma_native, 0.001 physical eV).

Reference and scalar numerical errors must be resolved separately. Same-root
numerical controls retain 0.1 meV in position/full width and 1e-6 rad phase;
include their uncertainty bounds in the comparison before declaring a pass;
unresolved slack is inconclusive. Numerical stability cannot establish
electronic accuracy. Wrong-branch,
wrong-width and overlapping-pole controls must reject; if the required
uncertainty/isolation cannot be established, the result is inconclusive.

A secondary reduced scattering observable can be stated exactly. For one simple
isolated pole, define on the real axis

    b(E;z) = (E - conjugate(z))/(E - z),
    d_res(E) = |Arg[b(E;z_scalar)*conjugate(b(E;z_native))]|/2.

Require d_res <= 0.05 rad over a predeclared, sufficiently resolved validation
window within the pole's isolation domain. This is the resonant determinant
phase contribution, not the full eigenphase sum or maximum full-S eigenphase
distance. For a simple uncancelled determinant pole,
`det S_native(E) = g_native(E) b(E;z_native)` locally; g contains the remaining
background. The reduction predicts b, **not g**. Test the observable at the
excluded geometry only after freezing the model. At training geometries it is
merely a restatement of fitted pole agreement. Even at the excluded geometry it
is derived from the independent pole, not an independent real-axis S experiment.
Sparse energy samples alone cannot certify the whole window.

**Do not compute g from held-out S and then count g*b_scalar as a prediction.**
Nor may an arbitrary smooth background refit turn a wrong-width factor into a
pass. A genuine total-scattering test needs a separately justified, frozen
channel/background embedding. That is an additional missing capability, not a
hidden loss function in this scalar map. Until supplied and tested, the original
held-out physical-scattering deliverable remains open even if reduced gates pass.
All 114 energy-route holdout S matrices from the direct-S screen remain
unanalyzed; none is consumed or reclassified by this design.

A fixed right-background replacement is insufficient even when it preserves
unitarity: if the background and pole factors do not commute, their product
need not be symmetric and reciprocity can fail. A separately specified
channel-aware construction would need its own reciprocity and predictive tests;
no CO embedding or local-potential equivalence is qualified here.

## Capability boundary and stopping decision

At the original capability checkpoint, `native_boundary.surface_matrix`
supported complex denominators, but had no
missing-spectrum correction. `exterior_numerov.numerov_k` is **real-axis only**:
float casts, real threshold decisions and real/imaginary column extraction make
complex-energy reuse invalid. `combined_exterior` provides complex outgoing maps
for its two-channel toy, not the 298-channel CO reference. Saved real-axis
reconstruction, finite controls and the 2048-to-4096 phase difference do not
qualify complex continuation or bound complex pole error.

The missing-node question for the declared local representation is now complete;
see the [canonical accepted scope](co-results-and-rejected-routes.md#local-representations-accepted-scope).
Common N+1 identity, electronic/global accuracy and a channel/background embedding
for full-scattering predictions remain separate unresolved claims. A qualified
local pole interpolation does not supply them.

**Current disposition:** reuse the saved evaluator and stop new compute. No
automatic complex-propagation development, full-S campaign, new50, repeated
QC/properties, dense inner solve or root/channel ladder. A different observable
needs a separate scope decision and bounded release. The historical follow-through
below records the original capability question, not a current missing-node queue.

**23 October follow-through:** the separately frozen
[one-anchor canary](co-outgoing-map-canary.md) now demonstrates finite298-channel
outgoing-map feasibility on compressed/nine/B1. It does not find a CO pole or
qualify complex continuation throughout a search region. The saved references,
2048-root omission and full-scattering limitations above are unchanged. A bounded
pole-reference proposal can now use measured map costs instead of assuming that
a general solver campaign is necessary.

**Later23 October result:** the separately frozen
[compressed pole experiment](co-compressed-pole-reference.md) now witnesses one
outgoing B1 pole of the saved finite2048-root model with a visible channel
residue. It is not a qualified multi-geometry CO width curve or absolute model
accuracy; paired-sector consistency, electronic-model choice, continuity and
the factory/held-out scattering requirements above remain open. No factory
target is populated from that single datum.

**Later paired outcome:** the compressed finite-reference
[B1/B2 open-Pi support/counterpart test](../../projects/ukrmol_co/PAIRED_POLE_CHARACTER.md)
passes with independently reconstructed saved witnesses. This is not merely
degeneracy, but also not exact Lambda conservation, single-l sufficiency,
electronic accuracy, a CO width curve or full-S equivalence. Matched electronic-
model evidence and all factory/held-out validation requirements remain open;
the equilibrium pole-target exclusion and114 route holdouts stay protected.
