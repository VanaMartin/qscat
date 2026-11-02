# Coupled canonical-Jost boundary control

## Source-only method decision

The exact diagonal [dual transport control](inverse-jost-dual-transport-control.md)
establishes sensitivity to an injected closed-dual error, not a CO tail error
bound. A smaller angle or larger endpoint is therefore not yet a justified CO
repair. Reject flux factorization/projection as a substitute for independent
amplitude normalization: it would make the physical-flux gate tautological.
Retain canonical outgoing-wave data and their derivative as the proposed boundary
prescription, and independently resolve finite-tail error before any CO change.

For a real energy, damping requires Im(k_j exp(i theta))>0. Attenuation of a
closed-dual error relative to open row o additionally depends on
Im((k_c-k_o)exp(i theta)). These are distinct conditions. The magnitude and
endpoint dependence of the omitted tail also matter; the diagonal crossover
alone cannot set an actual coupled error budget. No angle is selected for CO.

The minimum next evidence is one exact, genuinely coupled, real-symmetric
two-channel exterior. It tests the same Riccati/dual equations using both exact
canonical endpoint data and free-tail truncation, not a projection or fitted
normalization. This is a short-range control, not an inverse-power CO oracle.

## Self-contained exact construction

Use native radial Rydberg units: H=-d²/dz²+V, lengths bohr. Fix
T=diag(0,2), E=1, factorization energy e=-1, K=sqrt(T-e)=diag(1,sqrt(3)),
g=0.4, v=(1,1)^T. Products below use transpose, not Hermitian conjugation.

    b(z)=exp(-K z)v
    d(z)=1/g + sum_j (1-exp(-2 K_j z))/(2 K_j)
    B=b b^T/d; W=K+B
    W'=-K B-B K-B²
    V=T-2W'

Then W'+W²=T-e. For real z>=1, W and V are real symmetric and V has nonzero
off-diagonal coupling; V does not commute with T. For Re(z)>=1,
Re(d)>1/g, so this denominator has no contour pole. Define S=diag(i,-1),
N=(S-K)^-1 and the exact canonical wave matrix

    F=(S-W) exp(S z) N
    F'=(-W'+(S-W)S) exp(S z) N.

Direct differentiation of W'+W²=T-e verifies F''=(V-E)F. At infinity F tends
to exp(S z), so these are canonical open-outgoing/closed-decaying columns.
At a finite point Y=F'F^-1 and A=(1,0)F^-1 (sqrt(k_open)=1). Thus Y'=V-E-Y²
and A'=-AY are exactly the executed transport equations. The real-boundary
flux and bilinear Wronskian identities follow from these canonical solutions;
neither is enforced in the numerical result.

An independent oracle can avoid the rank-one formula entirely. Set C=K+g vv^T,
U=cosh(Kz)+K^-1 sinh(Kz)C, W=U'U^-1 and W'=T-e-W², then construct F,F'.
The constant seed Wronskian U^T U'-U'^T U=0 establishes symmetry. This independent
matrix route must agree before and after the numerical run.

## Frozen finite comparison

R=1; angles20/40/50degrees; lengths4/8/12bohr; rtol1e-8/1e-10/1e-12,
atol=0.01rtol; exact and free-tail endpoint data:54 integrations. Integrate
Y and A inward on z=R+t exp(i theta), t=L to0. Exact start uses F,F'; free start
uses Y=S, A=exp(-i z)(1,0), deliberately omitting the nonzero exponential tail.
All final outputs compare with the exact coupled real-boundary map.

Retain all three ODE resolutions and all three endpoint lengths. Finest exact
starts must satisfy maximum absolute Y/A error1e-7, flux1e-7 and Wronskian1e-6.
Finest free starts at length12 must satisfy the same gates. At each angle the
finest free-start maximum Y/A error at length4 must exceed1e-8 and fall by at
least100fold at length12, distinguishing a resolved tail error from the ODE floor.
Do not require strict pointwise monotonicity at roundoff. Any failure stops the
allocation without retries, adjusted parameters or CO fallback.

Freeze/push design, contract, code, synthetic controls and exact owner before any
ODE execution. Implementation60minutes;4 physical cores/2GiB/one library thread,
no network, seven source files only/no CO mount. Conditional1–10s;5s/10000 RHS
per integration,300s science/360s owner caps. Preserve original exits/resources.
Independent saved-only reconstruction and fresh public evidence precede a new
production decision. No native call, CO map, TZ guard, pole, B2, fit or holdout read.

## Interpretation limits

A pass establishes canonical normalization and finite-endpoint error separation
for this noncommuting exponential-tail model only. It cannot validate the CO
inverse-power recurrence, certify the coupled CO tail remainder or release a new
contour. Even above the diagonal rate crossover, a sufficiently rapidly vanishing
endpoint defect may converge: error injection and actual truncation are different.
After this one control, decide whether a separately justified CO method trial is
authorized or a different extraction route is needed; do not accumulate toys
without that decision. Original guard exits2 and both native-export receipts stand.

## Observed result and method decision —2 November2026

[Independent public reading](../../projects/ukrmol_co/qualification-evidence/coupled-jost/README.md)
authenticates all54 integrations from pushed `16f76fb`, all original exits0.
Exact-start maximum map errors are1.40548e-8,1.40328e-10,1.39615e-12 at the three
rtols. Free-start errors fall across the three lengths from2.23e-5 to6.22e-12
at20degrees,8.75e-5 to4.19e-10 at40degrees, and2.35e-4 to8.04e-9 at50degrees.
Every preregistered finest/longest map/flux/Wronskian and tail-separation gate
passes. The independent fundamental-seed reader performs no ODE calls.
Science0.576349s/owner1.236477s/0.052368GiB recorded cgroup high-water.

Canonical endpoint data are a sound normalization prescription for this exact
model. The coupled error still converges above the diagonal45degree crossover,
because the true omitted endpoint tail vanishes sufficiently rapidly. Therefore
neither a diagonal damping factor nor an angle choice substitutes for actual
coupled finite-tail convergence. This explicitly limits the earlier CO mechanism
diagnosis rather than upgrading it to causal proof.

The next useful calculation would now consume actual CO data. An
[unreleased three-map proposal](../../projects/ukrmol_co/stretched-contour-study-proposal-20261102.json)
would test one fixed20degree contour and three endpoints at the failed own-DZ
development energy, retaining original tolerance/order/flux/S budgets. It is
an explicitly new numerical prescription, not an original40degree guard retry;
it requires a user decision before any implementation/execution release.
Until then the original CO stop remains binding. No more generic toys are
selected: authorize a bounded actual-CO diagnostic or park this extraction route.

### Explicit separate authorization

On2November the user selects **Authorize study** for the fixed three-map CO
diagnostic. Its [new execution contract](../../projects/ukrmol_co/stretched-contour-study-contract-20261102.json)
does not alter the historical unreleased proposal or original40degree failure.
At the one own-DZ development energy80, evaluate20degrees/endpoints60/90/120
with unchanged order10/rtol1e-12/original real gates, plus90/120 full Y/A
stability1e-7. No empirical angle/tolerance ladder. Short-endpoint diagnostics
are retained as convergence evidence, not gate successes; input/integration/cap
failure aborts. Source/tests/owner freeze before maps and independent actual
CO evidence after execution remain necessary. No pole or peer release follows.

### Actual-CO outcome within that separate authorization

The [independent public three-map study](../../projects/ukrmol_co/qualification-evidence/stretched-contour-study/README.md)
qualifies the new20degree prescription at own-DZ development80 only. Endpoints
60/90/120 give flux9.85840e-8/4.70332e-9/1.01230e-9, all below the original1e-7
gate; Wronskian residuals1.51706e-9/7.95953e-13/4.74598e-15 pass1e-6. Full90/120
map changes Y2.32024e-13/A6.69041e-9 pass1e-7. Native S/phase and other invariants
also pass without altered tolerances or amplitude projection. All original exits0,
profile41.861348s/0.375618GiB kernel cgroup peak,33 executed Git sources and
210 resource samples independently authenticate. The first reader's mistaken
six-versus-five control count is preserved separately; no scientific retry.

This is actual evidence for a useful lower-angle prescription at the failed
development identity, not a rigorous CO tail bound or unique attribution of the
old40degree failure. It does not qualify DZ60/TZ60/TZ80 at20degrees or complex
continuation. Original40degree exits2 and both native-export receipts stand.
Complete the remaining own-model real prerequisites only under a separately
approved/frozen finite allocation; no automatic peer or pole use follows.

### Separate real-axis completion authorization

The user subsequently authorizes the [fixed five-map completion](../../projects/ukrmol_co/stretched-real-completion-contract-20261102.json):
reuse passed DZ80; DZ60/TZ60 at120 and TZ80 at60/90/120, unchanged20degree/
order10/rtol1e-12/real gates plus TZ80 endpoint stability. Only original own input/
method pins and physical budgets are inherited; historical40degree settings and
four-map routing are not inherited. Standalone/final physical gate failures,
input errors, integration errors or cap exhaustion stop without fallback. Short
TZ80 diagnostics remain part of its predeclared convergence dataset. No original
failure requalification or complex/pole/model release follows from this permission.
