# Fixed-grid attribution of the interior analytic veto

## Question and finite scope

The [interior prerequisite](../../projects/ukrmol_co/NUMEROV_INTERIOR_CONTROL.md)
rejects scalar roundoff and rotated open/closed K gates at the compressed peak
energy. Finer spacing worsens one K comparison. Before returning to CO, separate
discretisation from floating-point drift on those **unchanged grids** and test
whether equivalent channel representations alter the drift.

This is a free-exterior analytic attribution, not CO qualification. Energies in
the native radial equation are Rydberg coefficients, lengths bohr, and Q maps
u' to u at18 bohr. All old verdicts and exits remain. No inner/native solve.

## Matrix discrete oracle

For diagonal constant W=diag(threshold−energy), define B=I−h²W/12 and
A=2I+h²B⁻¹W. The exact-arithmetic discrete equation is

`u[j+1] = A u[j] − u[j−1]`.

The existing sixth-order Taylor maps are diagonal P(h),D(h), with
`u[1]=P(h)u[0]+D(h)u'[0]`. A final value/derivative pair must satisfy
`u[N−1]=P(−h)u[N]+D(−h)u'[N]`. Thus this oracle retains the same finite-step
initialisation, recurrence and matching, not just a continuous free solution.

For each open channel write Aii=2cos(theta); for closed Aii=2cosh(gamma).
Backward transfer of an endpoint basis column gives

`v[0]=(S[N] v[N−1]−S[N−1]v[N])/S[1]`,

`v[1]=(S[N−1]v[N−1]−S[N−2]v[N])/S[1]`,

where S is sin(j theta) or sinh(j gamma). Scale a closed column by
exp(−N gamma); this is an arbitrary nonzero column normalization and keeps the
initial-radius column finite without eliminating closed-channel effects.
Then `v'[0]=D(h)⁻¹(v[1]−P(h)v[0])`. Enforce `v[0]−Qv'[0]=0` on open f,g
and closed-decaying columns and solve for **unsymmetrized K**.

Independent checks compare this backward closed-form construction with a direct
high-precision forward recurrence on short mixed-threshold grids and with the
retained scalar discrete oracle. On long grids compare80/110-digit values.
These checks do not use the binary64 ratio recurrence as their reference.

## Separating errors and representations

For every twelve archived compressed rows, compute signed matrices
`K_float−K_discrete` and `K_discrete−K_continuous`. Their sum is the continuous
error. Both discrete and continuous oracles use the **actual binary64 Q**, not
the nominal tangent320. This separates floating-point drift from discretisation
for the free matrix problem. It does not identify one elementary operation or
bound a native CO error.

For a free decoupled closed channel, its physical decay derivative is−kappa u.
Analytic elimination at the boundary yields
`Qeff=Qoo−Qoc kappa/(1+Qcc kappa) Qco`.
Run the retained increment propagator on this two-open boundary in its original
rotated frame, and on its orthogonal eigenchannel representation. Compare each
against its own high-precision continuous/discrete oracle; transport diagonal K
back using the recorded eigenvectors. Binary64 elimination/rotation changes the
boundary slightly: require continuous equivalence<=1e-8 in K/S/phase before using
the comparison. Finite-step elimination is not presumed identical across models.

## Frozen calculation and stopping rules

Reuse the twelve original compressed rows. Reexecute scalar finest800 and full
rotated open/closed finest400 once to require **exact stored K reproduction**.
New reduced-boundary propagations are the two frames at four coordinates:
h=.003125/400, h=.0015625/200/400/800, order8. Total **ten new solves**; no new
spacing or endpoint. Pin the prior report hash and all source bytes before launch.

Diagnostic acceptance requires complete declared coverage, exact baseline replay,
all discrete-oracle80/110-digit and scalar consistency drifts<1e-20, continuous
representation equivalence<=1e-8, finite matrices and invariants<=1e-6. Oracle
failure/unset stops propagation. Preserve complete outputs and original exits.
No retry, tolerance change or CO release follows automatically.

Record two hypotheses independently of diagnostic acceptance: (1) floating-point
drift at the rejected coupled200/400 points exceeds discrete truncation by at
least5×; (2) diagonal-frame roundoff is<=one tenth of the reduced rotated-frame
drift and<=1e-7 at **all three** finest endpoints. A false hypothesis is a result,
not grounds to change a gate or repeat the calculation.

Four reserved physical cores0–3, one worker/one BLAS thread,4 GiB. Based on the
prior24-row72.95-s measurement, forecast120 s, upper240 s, cap600 s. Sample
resources at0.2 s, verify pinned sources and242 installed artifacts before/after,
and retain dated persistent `/home` evidence. Publish and independently verify
the saved decisions; public numerical replay scope is frozen separately.

## Recorded finite result

Freeze`aa917aa` completes all ten new solves and twelve reused decompositions,
**46/46 diagnostic checks / exits0/0**. At rejected200/400 points, matrix
floating-point drift exceeds discrete truncation by **120.46×/57.31×**, so the
>=5× dominance hypothesis passes. The all-finest-endpoint tenfold diagonal-frame
improvement hypothesis is **false**:800-bohr improvement is4.37×, although all
diagonal drifts are below1e-7. All representation equivalence checks pass with
maximum4.01e-12. Original interior false/exits1/1 remains.

Measured61.050941 s /60.487653 sampled CPU s /47,837,184-byte kernel high-water,
303 samples at0.2 s;51 sources/242 installed artifacts and prior report verify
before/after. See the [result and finite interpretation](../../projects/ukrmol_co/NUMEROV_PRECISION_ATTRIBUTION.md).
This measurement motivates separately costed increased-arithmetic-precision
attribution on fixed grids, not further spatial refinement or immediate CO use.
