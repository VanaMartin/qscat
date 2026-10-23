# CO scalar reduction: an outgoing-pole target, not a full scattering model

## Decision and claim boundary

Use **one continuously assigned outgoing resonance pole per Pi sector** as the
scalar factory's prospective electronic target. This is a physically defined
reduction: it preserves a coupled Hamiltonian's resonance energy and total decay
width, not its channel amplitudes. It replaces neither the native exterior nor
the full physical S matrix. The bare full-S route remains
[rejected](../../projects/ukrmol_co/DIRECT_SCATTERING_SCREEN.md).

This note completes the target-definition stage only. No CO pole has been
qualified, no factory target is populated, and no numerical successor or geometry
sweep is authorized. A pole-only fit plus the reduced validation below would
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

Let V0_ref(R) be the correlated neutral curve on a declared common dissociation
zero, distinct from the native absolute target energy. The eventual adapter is

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
  asymptotic inputs. No CO asymptote is guessed to make the API run.
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

The minimum later geometry design needs four genuine training nodes bracketing
that exclusion and one validation node. Keeping compressed/stretched anchors
therefore requires at least two additional training geometries. Their identities,
electronic prescription, cost and execution require a separate pilot decision;
this note launches none. A neutral input at equilibrium may be supplied openly:
this is conditional resonance validation, not prediction of the neutral curve.

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

## Capability boundary and stopping decision

`native_boundary.surface_matrix` supports complex denominators, but has no
missing-spectrum correction. `exterior_numerov.numerov_k` is **real-axis only**:
float casts, real threshold decisions and real/imaginary column extraction make
complex-energy reuse invalid. `combined_exterior` provides complex outgoing maps
for its two-channel toy, not the 298-channel CO reference. Saved real-axis
reconstruction, finite controls and the 2048-to-4096 phase difference do not
qualify complex continuation or bound complex pole error.

The immediate missing input is a **qualified, assigned native outgoing pole**,
including its sheet, residue witness and finite-model numerical uncertainty.
It cannot be obtained honestly just by rerunning the rejected phase extractor.
The separate production gap is a channel/background embedding for full-scattering
predictions; pole targeting deliberately does not supply it.

**Disposition:** retain this map as the scalar integration specification; stop
before a CO fit or automatic complex-propagation development. The bounded design
cycle is complete with zero scientific compute. Next decision is whether to
authorize a strictly one-anchor pole-reference feasibility contract, or redirect
to a coupled effective model. The former would need explicit implementation and
compute caps, known-truth controls, qualified prerequisites and a stop on failure
before extension to other anchors. No runtime estimate for a nonexistent qualified
298-channel solver is claimed here. New50, QC/properties repetition, dense inner
solves and automatic root/channel ladders remain out of scope.

**23 October follow-through:** the separately frozen
[one-anchor canary](co-outgoing-map-canary.md) now demonstrates finite298-channel
outgoing-map feasibility on compressed/nine/B1. It does not find a CO pole or
qualify complex continuation throughout a search region. The saved references,
2048-root omission and full-scattering limitations above are unchanged. A bounded
pole-reference proposal can now use measured map costs instead of assuming that
a general solver campaign is necessary.
