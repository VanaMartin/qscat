# Dtype-preserving coupled exterior: finite qualification and cost control

## Decision

The [free-exterior arithmetic control](../../projects/ukrmol_co/NUMEROV_ARITHMETIC_CONTROL.md)
qualifies64-significand-bit arithmetic on the previously rejected grids. CO
requires noncommuting inverse-power matrices and298 boundary channels. This
control asks whether that arithmetic extends to the coupled equation, and
measures its dense-channel cost before any new CO submission.

The real-axis equation is `u''=W(r)u`, with
`W=diag(threshold_Ryd−E_Ryd)+C2/r²+C3/r³`, lengths in bohr. Rydberg coefficients
are retained at the native-data boundary (twice the Hartree radial energy).
The zero-Bloch boundary is `u=Q u'`; no inverse of Q or K symmetrization is used.

## Algorithm and arithmetic

`typed_coupled_k` has the retained increment propagator's physical signature,
plus a dtype argument. It uses its ordered identity-defect algebra and sixth-order
Taylor endpoints; radial coordinates, potential derivatives, increments, complex
outgoing/decaying tail series, phases and matching remain in the declared dtype.
Use `acos(-1)` in that dtype, not a binary64 pi literal. Construct B differences
from potential differences, rather than subtracting two near-identity B matrices.
Constant free W is cached, including its solved Numerov increment.

The existing one-to-three-channel partial-pivot solver is retained unchanged.
Larger real systems use binary64 LU as a preconditioner, with extended residuals
and at most three binary64 correction solves. Require scaled infinity-norm
backward residual `||B−AX||∞/(||A||∞ ||X||∞+||B||∞) <=32 eps(dtype)` on every
large-system solve. Unsupported wider precision, singular LU, nonfinite values,
or exhausted refinement stops unset. This is explicitly mixed-precision linear
algebra, not an assertion that LAPACK itself runs in longdouble. Small solves
are checked by independent discrete oracles. Diagnostics retain maximum measured
large-solve residual and refinement count. No rejected check is waived by S.

Tail grouping retains the existing absolute pair-gap rule. Tests here use exact
degeneracies or well-separated thresholds; approximate-group error is not
qualified. The complex dtype must retain the real dtype's significand. Closed
columns are scaled at the matching radius, retaining their forced open components.

## Independent controls

1. Small-system pivoting/singularity tests and larger rank-one systems whose
   inverse is known by Sherman–Morrison; test12 and298 channels with two RHS.
2. Free centrifugal outgoing/decaying columns against spherical Bessel functions;
   degenerate inverse-square columns against the retained independent Hankel
   oracle, including analytic derivatives. General unequal-threshold C3 tail
   coefficients are compared against an independently expressed mpmath recurrence.
3. A short-grid mpmath direct **fundamental-value** Numerov recurrence at70/100
   digits. It propagates values, not ratios or identity defects. Initial values
   and final derivatives come from independently assembled sixth-order Taylor
   maps; match using the high-precision finite-order tail. This isolates discrete
   arithmetic from spatial and tail truncation on three noncommuting models.
4. Continuous short-interval DOP853 fundamental-value integration at two fixed
   tolerances, independently matched with mpmath tail derivatives. Reference
   differences must be below1e-9 before judging convergence; this is a bounded
   reference stability control, not an absolute CO error bound.

## Frozen finite calculation

Three short models use radius2, endpoint4, energy1 Rydberg, angular zeros,
tail order8, h=.1/.05/.025. The two-channel models use thresholds[0,.3] and
[0,1.8], Q=[[.3,.12],[.12,−.2]], C2=[[2,−.4],[−.4,6]],
C3=[[.1,.3],[.3,−.2]]. The three-channel model uses thresholds[0,.2,1.8],
Q=[[.3,.12,.04],[.12,−.2,−.03],[.04,−.03,.25]],
C2=[[2,−.4,.15],[−.4,6,.2],[.15,.2,1]],
C3=[[.1,.3,−.08],[.3,−.2,.12],[−.08,.12,.05]].
Run both dtypes:18 typed solves,9 discrete oracle pairs,6 continuous solves.
Require discrete-oracle stability<1e-30; arithmetic error<=1e-10 binary64 and
<=1e-13 extended. Finest continuous K/S/phase<=1e-6 and final raw K error
reduction>=8 unless both successive errors are below1e-10. The short-radius
finite-order tail is not exactly flux conserving: development against the
continuous reference gives open/open nonunitarity about.0067 at4 bohr. Record
that finite-tail defect and require its numerical invariants to agree with the
reference within1e-9; do not assert physical unitarity for this truncated problem.

Separately run the exact degenerate inverse-square hard-wall benchmark with
dipole1.2, radius1, endpoint20, energy1, angular[1,2], order12, and
h=.05/.025/.0125 in both dtypes:6 propagations. The independent noninteger-Hankel
S oracle and physical unitarity/reciprocity retain absolute1e-6 finest gates;
final S-error reduction>=8 unless both errors are below1e-10. This isolates
physical qualification from short-radius finite-tail truncation.

Reexecute the full three-channel compressed free stress at h=.0015625 and
endpoints200/400/800 in both dtypes:6 solves, using the actual archived Q/threshold
and energy bits from the arithmetic report. Require extended drift<=1e-8,
continuous K/S/phase<=1e-6 and invariants<=1e-6 at every endpoint. Compare extended
K against the saved free implementation at1e-8, without requiring identical
roundoff. Record all-three >=100× precision gain separately as a hypothesis.

Finally run a synthetic298-channel free problem in both dtypes, radius18,
endpoint18.125, h=.015625 (eight intervals), order0, energy.7, six open292 closed
channels. Q has a nontrivial3×3 open/open/closed block, all other entries diagonal;
the full298 matrices are still propagated densely. The independent matrix discrete
oracle evaluates only the active block plus four uncoupled scalar open columns.
Require K error<=1e-9 binary64/1e-12 extended, invariants<=1e-6 and every declared
large-solve residual<=32 eps. This qualifies the exercised dense solve path,
not general298-channel conditioning. Report extended/binary64 cost and a linear
per-interval extrapolation for the withheld60-solve CO proposal, labelled a
rough synthetic forecast rather than a measured CO cost. No cost threshold
changes the numerical verdict or automatically releases a successor.

Total32 typed propagations,9 high-precision discrete-oracle pairs and6 continuous
reference integrations, plus fixed algebra/tail preflights. No CO/native/inner
solve. One worker, reserved physical cores0–3/4 GiB, single-threaded BLAS;
forecast300 s, upper600 s, cap1200 s,0.2-s cgroup samples. Freeze source and
actual prior-report hash before submission. Any prerequisite control failure
stops later stages; retain failed/unset reports and original exits. Accuracy
failure completes its fixed stage, then stops before the next stage. No retry,
new precision/grid rung, relaxed gate or automatic CO invocation.

## Recorded result

The frozen control passes89/89 checks, original exits0/0, in45.782150 s. All32
typed propagations complete. Short-grid extended arithmetic error is<=6.67e-19;
exact inverse-square finest S error is2.281e-9, with fourth-order convergence.
All long-grid extended drift/continuous gates pass; the all-three>=100× gain
hypothesis is true. Synthetic298-channel K error is1.48e-18, but its extended
cost4.658 s is79× binary64, motivating a separate timing-only profile during
complete public replay. See the [result and limitations](../../projects/ukrmol_co/NUMEROV_COUPLED_PRECISION.md).
