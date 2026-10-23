# A bounded complex outgoing-map canary, before any CO pole search

## Production decision and scope

The [effective target](co-effective-resonance-target.md) needs a coupled outgoing
pole, but existing CO propagation is real-axis-only. Test whether a small Python
outgoing-map implementation can reproduce **three development S matrices at the
compressed anchor** and give contour-stable maps at one fixed complex energy.
This is a feasibility gate, not a resonance search or a qualification of widths.
One working-cycle implementation/diagnosis cap; fail or timeout stops this route
for reassessment, not an automatic algorithm/precision/refinement campaign.

Lifecycle: design here, project-only Python toy, analytic and independent BVP
controls, one saved-data canary. No optimization or qscat promotion. No new native
calls, inner solves, geometries, QC, properties or potential fits. Equilibrium
pole exclusion and all 114 route holdout S matrices remain protected.

## General-channel mathematics (derived here)

Use the saved Rydberg radial equation `F'' = V(r) F`, with
`V = diag(thresholds - energy) + C2/r^2 + C3/r^3`. C2 includes centrifugal terms;
native CF coefficients are not doubled. Radius and derivatives are in bohr.
At `r=a+(t-a) exp(i theta)`, integrate **backwards** from a finite tail:

    Y = F' F^-1,              dY/dt = exp(i theta) (V - Y Y),
    A = P_open.T F^-1,        dA/dt = -exp(i theta) A Y.

P_open selects the physical open-channel indices, not eigenchannels. Y is complex
symmetric, not Hermitian. Its real-axis imaginary part is the outgoing flux map.
For complex energy, open/closed identities and square-root sheets are continued
from a declared real reference energy, not selected by comparing complex numbers.
Require `Im(exp(i theta) k_j)>0` for every outgoing/decaying channel; thresholds
and invalid contours reject. Riccati coordinate poles may still cause failure.

Initialize from `F = C(r) D(r)`, `C = sum C_p/r^p`, C_0=I, with the explicit
inverse-power recurrence already derived in `exterior_numerov.tail_coefficients`,
but with the explicitly continued column roots. Its original principal-root
selector is not valid for closed channels below the complex cut. Keep that
real-axis implementation unchanged; retain its real threshold grouping algebra.
`Y=(C' + C diag(i k)) C^-1`. D cancels from Y; closed-channel exponential growth
need not be represented. For a real-energy S check also initialize A from the
open rows of `D^-1 C^-1`, with native outgoing normalization
`D_jj=exp(i(k_j r-l_j pi/2))/sqrt(k_j)` for open columns. No real/imaginary split
of a complex-energy solution is used. A is not needed for complex map probes.

For real energy and real symmetric zero-Bloch Q, recover physical open S without
propagating unstable closed Jost columns. The Wronskians imply

    Im Y = A.conj().T A,
    A.conj().T S_D = A.T,
    S = S_D + 2i A (I - Q Y)^-1 Q A.T.

S_D is the exterior Dirichlet-wall S, solved by a rank-revealing least-squares
system, not by normal equations. The transpose in the final term is bilinear.
For one free s wave, `A=sqrt(k) exp(-ika)`, Q=tan(ka)/k, and this gives S=1;
for Q=0 it gives the physical wall factor exp(-2ika). These fix the sign and
normalization independently. All real-axis S are measured unsymmetrized, with
unitarity, reciprocity and flux defects reported, not repaired.

## Independent controls and frozen finite application

Tests must cover free open/closed channels, known radial Hankel limits, invalid
sheets/contours, and real coupled S against the existing independent Numerov
boundary propagation. A coupled unequal-threshold C2/C3 complex case is compared
to a separately assembled second-order Chebyshev BVP, rather than another
Riccati solve. These validate toy regimes only; they do not certify the CO map.

Use unchanged compressed/nine/B1 data: all 298 channels/2048 roots, radius 18,
saved energy indices **60/70/80**, none satisfying i%5==4. Authenticate the original
123-input manifest and prior qualified report before/after; their original false/
unset outcomes stay unchanged. Only those three S values are numerically used.

For each real index, use endpoint 80 then 120 bohr, angle 30 degrees, tail order 8,
ODE rtol 1e-10 /atol 1e-12. Require maximum entrywise S error <=1e-6 against the
saved K-to-S reference, unitarity/reciprocity <=1e-6 and flux defect <=1e-7, as
well as finite full-rank matching. Stop at the first failed row; no complex
calculation is released by partial real-axis agreement.

Then evaluate **one declared non-pole complex probe**, incident energy
`0.14 - 0.04i Hartree`, hence `0.28 - 0.08i Rydberg`, with ground channels
outgoing and all others on the continued decaying sheet. It is a probe, not a
measured resonance or an admissible CO width interval. Four settings use endpoints
80/120 and angles 30/40, order 8, rtol 1e-10. One fifth setting uses endpoint 120,
angle 40, order 10 and rtol 1e-11. Require
`a max(abs(Y_i-Y_final)) <=1e-7` and `a max(abs(Y_i-Y_i.T)) <=1e-7`.
These are finite invariance checks, not a three-rung convergence proof or a
complex-pole error bound. Save every Y and real S, RHS counts and original failures.

At most **11 propagations**, one worker, four reserved physical cores, single-
threaded BLAS, 4 GiB. Runtime allocation is 300s target /900s upper /1800s total
hard cap, with **90s per propagation** and a fixed RHS-evaluation cap. These are
prospective allocations, not measured forecasts; timeout is inconclusive and
preserves partial records. Source/tests/contract freeze and push precede CO use.

**All gates pass:** only a finite outgoing-map feasibility result, enabling a
separately bounded pole-reference proposal. **Accuracy gate fails:** reject this
canary implementation/scope, not all possible outgoing maps. **Timeout, rank or
integration failure:** inconclusive feasibility; no retry or automatic successor.
No pole finding, CO width, omitted-spectrum bound, electronic-model choice or
production sweep follows from any outcome.

## Completed application

The [recorded result](../../projects/ukrmol_co/OUTGOING_MAP_CANARY.md) preserves
the initial test-dependency stop2/2/2 at `c00c09a`. A goal-reviewed dependency-only
successor `5a919f0` keeps the propagator, runner and scientific contract byte-exact
and passes all11 propagations /40 row gates in55.850916s. Worst real-S error
1.51042e-7 and radius-scaled complex map change2.06671e-11 pass the original gates.
Fresh saved public reconstruction verifies both owners and all40 decisions,
with no new propagation. This supplies finite feasibility only, not a pole,
width or all-anchor qualification; any root search still needs its own contract.
