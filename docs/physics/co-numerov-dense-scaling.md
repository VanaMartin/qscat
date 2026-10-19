# Dense extended-arithmetic setup versus interval scaling

## Decision and retained numerical reference

The [coupled precision control](../../projects/ukrmol_co/NUMEROV_COUPLED_PRECISION.md)
passes89 checks, but its eight-interval298-channel extended row takes4.66 s.
Profiling attributes61.92% to three Taylor setup/endpoint evaluations and26.11%
to ten mixed solves. The naive99-day extrapolation includes setup in every
interval and is not a CO throughput forecast. Measure recurring work before
choosing an optimization or costing the withheld interior calculation.

This benchmark changes no numerical module. Retain the exact coupled propagator
from`4a9baa1`, native Rydberg/bohr coefficients, actual binary64 input values,
zero-Bloch Q boundary, sixth-order endpoints and unsymmetrized K. Timers wrap
the existing solve, Taylor and tail function bindings and restore them after
each call; they change neither inputs nor returned arrays. Every matrix solve
still enforces the existing32-eps backward-error gate.

## Finite scope

Use the previous synthetic298-channel boundary: six open thresholds0,292 closed
thresholds2, energy.7, radius18, h=.015625, angular zeros, Q diagonal.3 except
the saved active[0,1,6] block and the four other open diagonals. Two potentials:

- Free: C2=C3=0, order0, identical to the previous dense probe.
- Varying: put C2=[[2,−.4,.15],[−.4,6,.2],[.15,.2,1]] and
  C3=[[.1,.3,−.08],[.3,−.2,.12],[−.08,.12,.05]] on that same3×3 block,
  zeros elsewhere; tail order8. These matrices do not commute. The full298
  matrices are propagated densely; no block reduction is used by the producer.

For each potential run N=8/32/128 intervals, endpoints18.125/18.5/20, both
binary64 and64-significand-bit longdouble, **two repetitions**:24 propagations.
Run repetition0 in increasing N and repetition1 in decreasing N, reversing
backend order on repetition1. This provides a bounded repeat spread rather than
an adaptive timing ladder. These are synthetic performance coordinates, not a
new CO discretisation family. Every original row is retained independently.

## Controls and independent references

The independent oracle reduces only the *reference* to the coupled3×3 block and
four free scalar open channels. For varying W use the tested direct fundamental-
value Numerov recurrence at80/110 digits with the same finite-order tail; for
free W use the independent discrete Jost oracle. Require oracle stability<1e-30
before propagation at each coordinate. Binary64/extended K arithmetic errors
must be<=1e-9/1e-12; S unitarity/reciprocity<=1e-6 in this tested tail regime.
Retain extended decimal K before binary64 presentation/error assessment.

Require both repetitions' K decimals to reproduce exactly; require both free
N=8 rows to reproduce the old stored dtype-specific K decimals exactly.
Unsupported precision stops unset. A failed control completes its two dtype
rows at that coordinate, then stops later coordinates; no retry or rung change.

## Timing attribution and prospective hypotheses

Free W has N+2 large solves: one initialization, one cached potential solve,
N−1 defect updates and one matching solve. Varying W has2N: initialization,
N−1 potential solves and N−1 defect updates, plus matching. Timer wrappers assign
roles by this frozen call order and require actual diagnostic counts to match.
All rows must have exactly three Taylor calls and one tail call. Preserve each
call's elapsed wall time, shape, role and returned dtype without copying arrays.

Report direct recurring-solve time separately from initialization/cache/matching,
Taylor and tail time. The sum of timed disjoint functions cannot exceed total
kernel wall time by more than timer-resolution slack1e-6 s. Untimed remainder
contains potential construction, endpoint products and other overhead; it is
not silently called setup or propagation.

For each potential/backend fit `T(N)=alpha+beta*N` to the two-repeat mean wall
times at all three N. Also report the32→128 marginal rate and repetition spread.
The predeclared **timing hypotheses**, separate from numerical qualification:

1. Every affine fit has positive slope and predicts all three mean times within
   20% relative error.
2. Free extended128-interval recurring solve time exceeds Taylor+tail time.
3. The varying/free extended marginal-rate ratio is between1 and3.

False hypotheses are retained; no extra repeats or interval counts follow.
Extrapolations use both fitted and direct measured recurring rates and explicitly
name synthetic conditioning/potential/tail scope. They are not native CO cost
or absolute-error bounds. Selection of a later optimization or CO pilot needs
a separate contract; this benchmark itself submits neither.

## Cost and evidence

One worker, reserved physical cores0–3/4 GiB, single-threaded BLAS,0.2-s cgroup
samples. Forecast180 s, upper300 s, cap600 s. Scope contains1344 intervals across
24 dense propagations; each varying interval performs two dense solves. Preserve
source/actual-prior-report hashes, all242 installed artifact checks, original
exits, dated persistent `/home` results and separately frozen public replay scope.
No CO/native/inner invocation or unique energy is added.
