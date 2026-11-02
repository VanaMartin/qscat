# Exact two-channel inverse-Jost transport control

## Decision served and boundary

The [saved CO diagnosis](co-stretched-real-guard-diagnosis.md) identifies a weak
separation of open and closed inverse-Jost rates, but does not prove why the CO
tail is contaminated. The smallest next control tests that transport mechanism
against an exact solution, without any CO/native file access. This is a numerical
method toy under `projects/`, not a CO calculation, correction or promoted solver.

## Analytic oracle

Use the native radial convention E=k², with one open ground channel E=1Ry and
one closed threshold T=2Ry, so k=κ=1bohr⁻¹. Set radius zero solely for this
constant-coefficient mathematical control, not the molecular boundary radius.
Let z=t exp(iθ), t decreasing fromL to0. The exact logarithmic derivative is
Y=diag(ik,−κ). Its single open dual row satisfies

    dAo/dt = −ik exp(iθ) Ao,
    dAc/dt = κ exp(iθ) Ac.

Prescribe Ao(L)=exp(−ikL exp(iθ)) and inject Ac(L)=ε Ao(L), where ε=0 or1e-6.
At zero the independent analytic result is

    Ao(0)=1,
    Ac(0)=ε exp[−(κ+ik)L exp(iθ)].

The closed/open contamination magnitude divided byε is exp[−L(κcosθ−ksinθ)].
Its crossover is exactly45degrees. Below it errors attenuate; at it their
relative magnitude is unchanged; above it they amplify. This mathematical
condition is distinct from requiring all outgoing/closed waves to decay on the
rotated contour; both decay for every tested angle20–50degrees.

For Q=0 the saved physical-S algebra reduces to the least-squares wall phase.
With ε=0 its exact S=1 and flux defect0. With injectedε≠0, S may stay within1e-6
while the full flux identity fails1e-7, reproducing the logical distinction that
blocked CO. Do not enforce or project the identity; retain the intentional error.

## Frozen test and bounds

[Contract](../../projects/ukrmol_co/dual-transport-toy-contract-20261102.json).
Exactly72 DOP853 integrations: angles20/40/45/50degrees, lengths10/20/40bohr,
ε=0/1e-6, and rtol1e-8/1e-10/1e-12 with atol=rtol*0.01. The analytic oracle
contains no ODE integrator. Record the full two-component complex output, analytic
output, scaled error, contamination ratio, S and flux defects and RHS counts.
Require every finest scaled error<=1e-6 using denominator1e-12+|exact|; require
zero-injection closed output<=1e-14. Record all three resolutions rather than
asserting strict monotonicity across roundoff. Require at least one finest
contaminated case whose S passes1e-6 but whose flux fails1e-7.

Freeze/push code, exact parameters, synthetic controls and owner before execution.
Four physical cores0–3/2GiB/one library thread, no network, source read-only and
output-only writable mount; **no CO campaign/input mount**. Per integration5s/
10000 RHS and combined300s cap, conditional1–10s. Implementation bounded30minutes.
No dependency install, image build, native quantum executable or paid resource.
Preserve any error with no retry, new tolerance, angle, length or parameter ladder.

## Interpretation

A pass validates this exact transport mechanism and known negative control only.
It does not validate the CO inverse-power tail, coupled Y integration or complex
pole residue, nor choose a replacement CO angle. A failure stops method work for
review. Any pass returns to a separate numerical-method decision; original CO
exits2, all budgets, qualified native exports and114 holdouts remain untouched.
No new CO map, TZ guard, pole, B2, fit or sweep follows automatically.

## Observed exact-control result —2 November2026

[Independent public evidence](../../projects/ukrmol_co/qualification-evidence/dual-transport-toy/README.md)
authenticates all72 original integrations from pushed `708d3eb`. Maximum scaled
errors at the three declared resolutions are3.51990e-8,3.55123e-10,3.56998e-12;
the approximately hundredfold reduction per hundredfold rtol refinement is
resolved well above roundoff. The finest1e-6 gate and zero-injection control pass.
Seven finest contaminated cases pass S while failing flux. The largest injected
output,50degrees/40bohr, has flux1.38417e-4 but S error7.43918e-9. The independent
reader uses real scalar trig/exponential formulas rather than the integrator or
its complex-exponential oracle. All original exits0;0.350555s science/1.236827s
owner/0.070160GiB recorded cgroup high-water. No CO inputs mounted or new CO maps.

This resolves the toy-level sign/transport mechanism, not actual coupled-CO
error attribution. The original CO stop and all gates remain unchanged. Next
review a separately bounded source-only amplitude/tail method proposal; do not
turn the toy's50degree example into a molecular contour recommendation.
