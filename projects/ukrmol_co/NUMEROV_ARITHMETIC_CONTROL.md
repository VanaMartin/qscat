# Fixed-grid extended-arithmetic control

**Pass:35/35 control checks,10/10 qualification checks, original exits0/0.**
Increasing arithmetic precision on unchanged free-exterior grids suppresses the
measured drift. The original interior prerequisite remains **false/exits1/1**;
the new solver is limited to constant diagonal exteriors and at most three
channels. A coupled inverse-power CO formulation still needs separate validation.

## Question, freeze and scope

The [precision attribution](NUMEROV_PRECISION_ATTRIBUTION.md) found analytic
matrix roundoff57–120 times larger than truncation at the rejected200/400-bohr
points. Its all-endpoint tenfold basis-improvement hypothesis was false.
The [design](../../docs/physics/co-numerov-arithmetic-control.md) and
[finite contract](numerov-arithmetic-contract.json) were frozen/pushed at
**`44ed1952858cd471df64721a837353cac45f8b0a`** before execution.

The control preserves actual archived binary64 Q/threshold/energy inputs and the
Numerov difference problem, sixth-order Taylor endpoints and free-tail matching.
The new typed implementation uses explicit partial-pivot small solves, direct
Taylor increments and column-equilibrated matching in both arithmetic formats.
Binary64 and longdouble differ only in dtype; measured significands are53/64
bits. This compares precision within that implementation; it does not isolate
one operation in the older SciPy implementation.

Six fixed coordinates, each in both dtypes, plus two exact original SciPy
baseline reexecutions give **14 analytic solves /3957760 radial intervals**.
No new spatial rung, endpoint or CO identity. The independent backward discrete
Jost oracle and continuous oracle pass80/110-digit stability gates. Extended K
components are retained as decimal strings before binary64 presentation or error
assessment. S and phase comparisons use binary64 presentation; their budgets are
1e-6, far larger than its presentation precision.

## Result on unchanged grids

| Model | h / bohr | Endpoint / bohr | Extended K drift | Continuous K error | Binary64/extended drift |
|---|---:|---:|---:|---:|---:|
| Scalar | .003125 | 400 | 5.16014e-11 | 4.74610e-7 | 983.19× |
| Scalar | .0015625 | 800 | 1.48802e-10 | 6.34103e-8 | 812.39× |
| Rotated open/open/closed | .003125 | 400 | 3.71634e-10 | 3.04153e-7 | 718.67× |
| Rotated open/open/closed | .0015625 | 200 | 4.51905e-10 | 9.87301e-9 | 1348.90× |
| Rotated open/open/closed | .0015625 | 400 | 5.52006e-10 | 2.02835e-8 | 453.27× |
| Rotated open/open/closed | .0015625 | 800 | 5.36875e-10 | 4.10242e-8 | 2458.25× |

All six extended roundoff errors pass **1e-8**; all four declared finest rows
pass separate continuous K/S/phase **1e-6** budgets. Worst finest continuous K
error is6.34103e-8, S error1.23848e-12 and phase error1.24736e-12 rad. All
invariant controls pass. Both older SciPy baselines reproduce every K component
exactly. The separate prospective hypothesis of **>=100× improvement at all six
coordinates is true**. The typed binary64 result need not equal the older SciPy
result; their arithmetic expressions and matching differ.

## Cost and evidence

Finished **2026-10-19T12:51:11.077309UTC**: **64.803874 s /64.275304 sampled CPU s /
46,731,264-byte kernel high-water**,322 samples at0.2 s. One worker, reserved
four physical cores0–3/4 GiB, single-threaded BLAS. Below360-s forecast/600-s
upper/1200-s cap. All57 frozen sources,242 installed engine artifacts and pinned
prior report verify before/after. Original owner, outputs, exits and resources
remain under`/home/kooza/ukrmol/co-sa-20261006/prepared/numerov-arithmetic-control-20261019/`.

The separate [saved-decision verifier](numerov_arithmetic_replay.py), frozen at
**`a372e8b`**, reconstructs all45 decisions and six precision-gain ratios from
pinned bytes without propagation; [verification](verification-numerov-arithmetic-saved-20261019.json)
passes. It shares the tested discrete/continuous oracles and uses separate
scattering/error expressions. Public sources and results are listed in the
[evidence pointer](qualification-evidence/numerov-arithmetic/README.md).

## Complete public numerical replay

The [replay contract](numerov-arithmetic-public-replay-contract.json), frozen at
`a372e8b`, runs all14 campaign solves from fresh public bytes. Every typed
**extended decimal K component** and every original-baseline K component
reproduces exactly. All45 decisions, oracle controls and hypothesis values are
identical; replay exits0/0. Six tiny short-grid preflight solves in both dtypes
also pass: worst extended discrete-oracle error6.50e-18, plus pivot-demanding
and singular small-solver controls. These are additional preflight tests, not
new campaign coordinates. NumPy calls the stored longdouble dtype`float128` on
this host, but its measured significand is64 bits, not IEEE binary128's113.

Finished **2026-10-19T12:58:08.153239UTC**: **64.704617 s /64.082103 sampled CPU s /
44,285,952-byte kernel high-water**,321 samples at0.2 s, reserved cores0–3/4 GiB.
Below90-s forecast/180-s upper/600-s cap. All57 frozen sources, one separate
verifier,242 installed artifacts and prior report verify. The8-file replay
collection remains under`prepared/numerov-arithmetic-public-replay-20261019/`.
All **14 final public payloads** fresh-verify; local public tests pass **74 with
one wider-precision test skipped** on the Mac. The remote precision preflight
and short-grid controls exercise the wider format. The pinned tools image lacks
pytest (an exploratory capability probe returned ModuleNotFoundError); its
standalone public worker runs the stated remote controls directly.

The [12:59:29UTC capture](progress-numerov-arithmetic-result-20261019.json)
reverifies18 collected hashes,57 frozen sources,242 installed artifacts and
prior report. No scientific containers run;60.76 GiB persistent disk remains.
[Public verification](verification-numerov-arithmetic-public-20261019.json)
records the exact replay, resources and original exits.

## Consequence and next finite decision

The free constant-potential arithmetic control is qualified in its declared
scope. The next numerical step is to design a dtype-preserving coupled
inverse-square/inverse-cube propagator, including matrix solves, potential
derivatives, tail recursion and matching, and qualify it against independent
analytic/reference controls under a separate costed contract. This small
free-exterior solver cannot propagate the298-channel CO model.

Only after that qualification should a separately frozen interior CO contract
revisit the deterministic ten-identity proposal. Candidate-qualified CO identities
remain14, historical coverage132/594 with462 uncomputed. All older numerical,
channel/model, continuity, extraction and neutral-domain rejections remain;
no full sweep or complex-pole interpretation is released.
