# Coupled extended-arithmetic qualification and dense-channel cost

**Pass:89/89 checks, original exits0/0,32 typed propagations.** The formulation
retains higher precision through noncommuting inverse-power potentials,
derivatives, tail recursion and matching in the tested regimes. Dense-channel
cost now limits a reasonable CO successor; the synthetic probe is79× slower
in extended arithmetic. It does not authorize a new CO sweep.

## Frozen question and implementation

The [design](../../docs/physics/co-numerov-coupled-precision.md) and
[contract](numerov-coupled-precision-contract.json) are frozen/pushed at
**`4a9baa17606bbbc11403eeba1a53fe0c4b2f4f35`** before submission. Inputs retain
the native Rydberg/bohr radial equation and length-valued zero-Bloch Q map.
The [experimental implementation](numerov_typed_coupled.py) preserves ordered
identity-defect algebra, sixth-order Taylor endpoints, unsymmetrized K and
scaled closed columns. Complex tails retain53/64 significand bits in the two
formats; approximate threshold grouping is outside the qualification.

One-to-three-channel systems reuse the unchanged partial-pivot solver. Larger
systems use binary64 LU preconditioning with extended residuals and up to three
corrections; every solve must pass a32-eps scaled backward-error gate. This is
mixed precision, not extended-precision LAPACK. Tests qualify only the exercised
conditioning, not every298-channel surface matrix.

## Separate arithmetic, truncation and physical gates

The short noncommuting two/three-channel models use2→4 bohr, h=.1/.05/.025 and
tail order8. An independently expressed70/100-digit **fundamental-value Numerov
recurrence**, not a ratio recurrence, solves their discrete problem. Continuous
DOP853 fundamental integrations at two fixed tolerances have reference differences
<=2.47e-13. Short-radius finite-order tail nonunitarity is real in the reference
(about.0067 for unequal open thresholds), so it is recorded and compared to that
reference rather than incorrectly asserted zero. Physical unitarity is tested
separately against exact inverse-square scattering.

| Control | Binary64 | Extended | Gate / interpretation |
|---|---:|---:|---|
| Worst short-grid discrete arithmetic K error | 8.73342e-16 | 6.66460e-19 | 1e-10 /1e-13 |
| Worst short finest continuous K error | 5.03487e-9 | 5.03487e-9 | 1e-6, separate S/phase also pass |
| Exact inverse-square finest S error | 2.28039e-9 | 2.28025e-9 | 1e-6, physical invariants pass |
| Synthetic298 discrete K error | 3.25774e-15 | 1.47948e-18 | 1e-9 /1e-12 |

All short families show the required final>=8× raw error reduction; their measured
reduction is approximately16×. Exact inverse-square hard-wall scattering also
converges approximately16× on h=.05/.025/.0125, with order12 matching at20 bohr.
Rank-one12/298-channel algebra and independent tail/Hankel preflights pass.

The long-grid stress uses the exact archived compressed Q/threshold/energy bytes,
zero exterior coefficients, h=.0015625 and endpoints200/400/800:

| Endpoint / bohr | Binary64 K drift | Extended K drift | Extended continuous K error | Precision gain |
|---|---:|---:|---:|---:|
| 200 | 4.06816e-6 | 3.87651e-10 | 9.03346e-9 | 10494.37× |
| 400 | 3.91768e-6 | 4.08505e-10 | 1.93231e-8 | 9590.28× |
| 800 | 3.49962e-6 | 7.10315e-10 | 3.97771e-8 | 4926.86× |

Every extended drift passes1e-8 and continuous K/S/phase/invariants pass1e-6.
The all-three>=100× precision-gain hypothesis is **true**. Difference from the
previous extended free implementation is at most1.24719e-9, within1e-8. The
generic binary64 implementation has larger stress drift than the older method;
its stress accuracy is not qualified by the extended pass. It is a diagnostic
comparison, not a replacement for the retained binary64 implementation.

## Synthetic dense cost is the next bottleneck

The298-channel synthetic problem propagates all matrices densely for eight
intervals, including six open and292 closed channels; its coupled boundary has
an independently solvable3×3 block. This exercises large-system solves and
matching, but has a free diagonal exterior. Extended dense backward error is
<=5.73343e-20, one correction suffices, and unsymmetrized K transpose error is
6.78e-20. Rank-one298-channel forward error is6.82e-19.

Dense binary64 costs **0.0589693 s** versus extended **4.6581683 s**, ratio
**78.993×**. Dividing the latter by eight gives.58227 s per interval *including
startup, tail and matching*. A naive60-solve400-bohr/finest-step extrapolation
is **8,541,217 s /98.86 single-worker days**. This is an intentionally rough
synthetic estimate: setup cost, interval count, nonzero potentials, conditioning
and concurrency differ from CO. It is not a measured CO throughput or accuracy
bound. It motivates profiling before scheduling a costly successor.

The separately frozen complete public replay profiles only its two existing
dense rows, with zero extra propagations. It preserves numerical source/input
bytes and exact decimal comparisons. The next decision should use that measured
profile to target arithmetic cost, retaining this Python formulation as the
differential reference and requalifying any optimized kernel before CO.

## Resources and preservation

Finished **2026-10-19T13:16:58.070332UTC**: **45.782150 s /45.227669 sampled CPU s /
140,115,968-byte kernel high-water**,227 samples at0.2 s. One worker, four
reserved physical cores0–3/4 GiB, single-threaded BLAS; below300-s forecast/
600-s upper/1200-s cap. All63 frozen sources,242 installed artifacts and prior
report verify before/after. Original owner, outputs, exits and resources remain
under`/home/kooza/ukrmol/co-sa-20261006/prepared/numerov-coupled-precision-20261019/`.

The [saved verifier](numerov_coupled_replay.py), frozen at`d4662d3`, independently
recomputes propagation-row errors and decisions. It audits stored preflight
diagnostics; complete numerical replay reexecutes the actual preflight algebra
and tails. [Saved verification](verification-numerov-coupled-saved-20261019.json)
passes all89 decisions. [Public evidence](qualification-evidence/numerov-coupled/README.md)
retains source/results and the separately frozen profiled-replay contract.

## Complete public replay and measured profile

The [public replay contract](numerov-coupled-public-replay-contract.json), frozen
at`d4662d3`, reexecutes all32 typed propagations,9 discrete oracle pairs,6
continuous integrations and all algebra/tail preflights from fresh public bytes.
Every extended decimal and binary64 K component, all89 checks, stage verdicts,
hypothesis values, references and preflight diagnostics reproduce **exactly**.
Original and verification exits remain0/0. Timing-only cProfile wraps just the
two existing synthetic dense rows; zero additional propagations.

Extended dense-row profile,4.656866 s total:

| Function | Calls | Cumulative time | Fraction of dense row |
|---|---:|---:|---:|
| Sixth-order Taylor changes | 3 | 2.883755 s | 61.92% |
| Mixed-precision solves | 10 | 1.215995 s | 26.11% |

These functions are disjoint call groups; together they account for88.04%.
Taylor work is setup/endpoint work, not per-interval propagation. The naive99-day
extrapolation therefore cannot be treated as a measured interval throughput or
CO completion forecast. The profile identifies expensive dense extended products
and residual work, but does not by itself isolate a single elementary operation.
**Next:** a separately costed setup-versus-interval scaling benchmark, then target
the measured hot path while retaining this Python implementation as the numerical
reference. Requalify any optimized formulation before the withheld interior CO
contract; no automated precision/spacing ladder follows this pass.

Replay finishes **2026-10-19T13:21:16.353254UTC**: **45.829095 s /45.263771 sampled
CPU s /112,754,688-byte kernel high-water**,227 samples at0.2 s, four reserved
cores0–3/4 GiB, within90-s forecast/180-s upper/600-s cap. All63 sources, one
separate verifier,242 installed artifacts and prior report verify before/after.
The10-file replay collection includes both binary profiles and remains under
`prepared/numerov-coupled-public-replay-20261019/`. All **14 final payloads**
fresh-verify; local public tests pass **80 with two wider-precision skips**.
The remote experiment/replay exercise wider real and complex arithmetic directly.

The [13:23:03UTC capture](progress-numerov-coupled-result-20261019.json) reverifies
20 collected hashes,63 sources,242 artifacts and prior report. No scientific
containers run;60.75 GiB persistent disk remains. The
[public verification record](verification-numerov-coupled-public-20261019.json)
contains exact replay and full top40 function timings; binary profiles are in
the public replay archive.

Original interior prerequisite false/exits1/1, all older rejections and historical
coverage132/594/462 uncomputed remain. Candidate-qualified CO identities remain14;
no new CO/native/inner solve, unique energy, full sweep or pole interpretation.

The subsequent [dense scaling benchmark](NUMEROV_DENSE_SCALING.md) passes166/166
checks on24 unchanged-kernel propagations. Extended marginal rates are
.135263/.271376 s per interval for free/varying W; recurring solves dominate
at128 intervals. This resolves setup-versus-interval cost in the synthetic
regime and prioritizes mixed-solve internals for the next profiling/optimization
contract. It does not measure actual CO throughput or release the interior scope.
