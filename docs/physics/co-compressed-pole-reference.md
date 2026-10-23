# A finite compressed-anchor outgoing-pole experiment

## Decision, scope and lifecycle

The [outgoing-map canary](../../projects/ukrmol_co/OUTGOING_MAP_CANARY.md) supplies
measured 298-channel map costs, not a pole. Test whether **one outgoing pole of
the saved compressed/nine-A2/B1, 2048-root reference** can be located and verified
without a general solver campaign. Preserve all channels and the native exterior.
No new native call, QC/property work, inner Hamiltonian, root/channel ladder,
geometry, potential fit or other anchor is included. All114 energy-route holdout
S matrices and the equilibrium pole exclusion remain protected.

This is design → project Python → analytic/independent algebra controls → one
costed application, not Rust optimization or qscat promotion. Implementation and
diagnosis are limited to this working cycle. A failed search, witness, dependency
or runtime gate stops the experiment; no automatic repair, reseeding or successor.

## Energy, sheet and finite search region

Use incident **Hartree** z relative to the saved native ground target. Native
surface denominators receive `E_ground + z`; the exterior receives `2*z` in
Rydberg. Preserve total energy, incident z, `E_res=Re z` and full width
`Gamma=-2 Im z`. Convert reporting/budgets with `qscat.units`, not the historical
requested-eV grid conversion. Keep the accepted2048 exported eigenpairs; no
omitted-complex-spectrum bound is inferred from real-axis root controls.

The prospective search box is

    0.10 <= Re z <= 0.18 Ha,       -0.070 <= Im z <= -0.015 Ha.

Seeds are `0.140-0.040i` and `0.142-0.040i Ha`, chosen from the already declared
non-pole canary neighborhood, not a qualified CO width interval. This is a finite
search allocation, **not** a claim that a physical pole must lie inside it.
Leaving the box, a singular secant step or exhausted evaluations is inconclusive;
do not enlarge the box or add seeds. Ground-target channels continue outgoing
from positive real momentum; all other channels continue the closed decaying
sheet from real reference `0.28 Ryd`. No threshold is crossed. Require contour
exposure for every channel on every evaluation, including contour quadrature.

## First construction: a bordered boundary zero

With native spectral matrix Q in bohr and outgoing Y in inverse bohr, form

    M(z) = I - Q(E_ground+z) Y(2*z).

At the starting point choose the last left/right singular vectors and freeze
them as bordering vectors b,c within that search. Solve

    [ M(z)   b ] [ x     ] = [ 0 ],
    [ c.T    0 ] [ alpha ]   [ 1 ].

Alpha is a locally holomorphic scalar zero function; no determinant underflow,
matrix symmetrization or conjugated spectral product is used. Apply complex
secant iteration with a fixed maximum step0.025Ha and the declared box. A point
can terminate only with |alpha|<=1e-10, smallest singular value of M<=1e-8 and
second-smallest>=1e-4. These test one simple null direction; they are not by
themselves evidence of a physical S pole.

Base setting: endpoint80bohr, angle30degrees, tail order8, ODE rtol1e-10.
Two refinements restart only from the previous candidate:120/30/order8/rtol1e-11
and120/40/order10/rtol1e-12. These are the complete three settings, not an adaptive
ladder. Each refinement uses the previous candidate and that candidate plus
2e-5Ha as its secant seed pair. Root locations must agree within **0.1meV position and full
width**, retaining same-root budgets. A candidate must also lie at least0.001Ha
from the search-box edges before a contour witness is allowed.
This joint configuration family measures finite stability; it does not isolate
an endpoint or ODE convergence rate or supply a rigorous absolute error enclosure.

## Second construction: an interior resolvent contour

An independent pole construction must not merely minimize the same scalar zero
function twice. Let D contain the2048 native total-Hartree inner eigenvalues,
P(z)=diag(1/(D-E_ground-z)), and W the stored pole-major surface amplitudes.
The nonlinear effective interior operator is

    H_eff(z) = D - W Y(z) W.T.

This follows from the same boundary spectral expansion, with its absorbed half
factor. No dense2048-dimensional Hamiltonian is constructed or diagonalized.
For eight fixed real orthonormal probe vectors X (seed314159, generated before
CO evaluation), compute the projected physical resolvent by a Schur identity:

    B = W.T P X,
    R_X(z) = -[ X.T P X + B.T Y M^-1 B ].

Here `R_X = X.T [E_ground+z - H_eff(z)]^-1 X`; X is only an observable projection,
not a truncation of W, Q or the298 channels. Evaluate32 equally spaced points
on a radius0.001Ha circle about the final candidate, at the final map setting.
The16-node rule reuses its even indices; no second set of16 propagations.

    C0 = mean[(z_j-z_center) R_X(z_j)],
    C1 = mean[z_j (z_j-z_center) R_X(z_j)],
    z_contour = vdot(C0,C1)/vdot(C0,C0).

For one enclosed simple pole, C0 is its projected residue and C1=z_pole C0.
The conjugated Frobenius scalar product above is only a moment estimator, not
an ECS inner product. Require nonzero residue norm>1e-8, rank-one ratio
sigma2/sigma1<=1e-3, relative16→32 moment change<=1e-5, and contour pole agreement
with the bordered root within0.1meV in position/full width. Agreement of the two
quadratures must also meet those same-root energy budgets. The circle stays in
the lower half plane and contains no real inner eigenvalue. This helps separate
a true resolvent pole from a boundary-coordinate singularity or cancellation.
Also require the discrete phase winding of det M to be+1 on both16/32 rules,
with every adjacent principal phase jump<=pi/2. This is a finite resolution/
isolation check, not a rigorous global pole count. Do not add quadrature points
or shrink the circle if this gate fails.

## A genuine open-channel residue, not a hidden interior state

A pole in an interior projection alone does not establish a scattering pole.
Continue the **outgoing Jost inverse rows A** analytically on the same declared
sheet. This needs an explicit option on `outgoing_map`; it does not allow the
real-axis `physical_s` formula or its conjugations at complex energy. The ODE
and tail recurrence are unchanged; only the existing A equations are exposed
for complex z, with continuous open normalization `sqrt(k) exp(-i(kr-l*pi/2))`.

The real-axis relation is `S=S_D+2i A M^-1 Q A.T`. Where Y and A are regular,
its resonant term has an analytic continuation

    T_out(z) = A(z) M(z)^-1 Q(z) A(z).T,
    Res S = 2i mean[(z_j-z_center) T_out(z_j)].

S_D is not reconstructed from complex conjugates. All pointwise maps and A must
be finite/regular throughout the witness; no exterior coordinate pole may be
accepted as an S pole. Test the6x6 contour residue for rank one and nonzero
open signal, norm(Res S)/Gamma>1e-6, with relative16→32 change<=1e-5.

Independently check its normalization at the center using right/left null
vectors v and ell.T of M. Approximate M' with fixed centered Hartree steps
1e-5 and5e-6. For each,

    R_boundary = v (ell.T Q)/(ell.T M' v),
    R_S,derivative = 2i A R_boundary A.T.

Require finite nonzero derivative denominator, simple singular direction,
boundary residue transpose symmetry and agreement of both derivative residues
with the32-point contour residue to relative1e-3. This is a nonzero flux-channel
witness, not a total-S reconstruction or a branching-ratio production fit.

**Independence limit:** bordered root, interior-resolvent moments and channel
residue are different constructions, but share native W/D and the same outgoing
ODE. This is not an independent298-channel exterior engine or a proof of absolute
CO accuracy. Success means a numerically witnessed pole of the **declared finite
reference**. Omitted spectrum, electronic model, geometry continuity, neutral
correlation, scalar fitting and full-scattering embedding remain separate gates.

## Controls, cost and stopping

Before CO use, test the Schur resolvent against a direct small dense inverse;
the bordered root/contour/residue against the exact positive-weight one-level
boundary model `Q=w^2/(epsilon-z)`, `Y=i sqrt(2z)`; and a decoupled interior-state
counterexample with zero open residue. For the one-level model,

    z_pole = epsilon - w^4 - i w^2 sqrt(2 epsilon-w^4).

Use known-truth wrong-width and off-root counterexamples, complex Jost free-wave
normalization and invalid-sheet/box/evaluation-cap checks. Standard-library test
harnesses use NumPy/SciPy only, with the repository unit helpers included in the
source capsule; no incidental mpmath/pytest installation or image rebuild.

Reuse authenticated canary output and all123 original inputs. Freeze/push code,
tests and exact contract before any CO root search. Base search has at most16
map evaluations; each of the two refinements at most12. At most five center/
derivative evaluations plus32 contour points: **77 map evaluations maximum**.
Use eight physical cores0–7, one BLAS thread per worker, at most eight independent
contour workers,8GiB. Root searches are sequential. Prior3.53–8.27s/map supports
a conditional **600s target /1200s upper /1800s hard cap**,90s per map. Conditioning
near a pole is unmeasured; these are allocations, not guaranteed completion times.

Any failure preserves all partial matrices, points, resources and original exits.
On a contour-worker failure, cancel queued work and collect already in-flight
points before stopping; no new batch or replacement point is submitted.
Accuracy/witness rejection exits1; no candidate, rank/conditioning/ODE failure,
timeout, input failure or evaluation cap exits2 with qualification unset. All
declared gates passing exits0, qualifying only the finite-reference pole scope
above. No retry, expanded region, new model, second sector/anchor, fit or sweep
is released by this contract. Review immediately on any outcome.

## Completed finite application

Freeze `9b4f588` uses47 maps and finishes23 October14:45:31UTC, exits0/0/0.
The finite compressed/nine/B1 reference pole is3.7604eV with2.0846eV full width.
Boundary, interior contour and visible open-channel residue gates all pass.
The tighter settings accept the original candidate without new iterates;
identical reported locations are not a zero uncertainty claim. Container
180.990698s /cgroup peak1.94329GiB. Public saved reconstruction verifies every
witness without another propagation. See the
[result and claim boundary](../../projects/ukrmol_co/COMPRESSED_POLE_REFERENCE.md):
this is one pole of the finite declared model, not independent exterior or
absolute CO accuracy, a width curve or a potential. No automatic next job.
