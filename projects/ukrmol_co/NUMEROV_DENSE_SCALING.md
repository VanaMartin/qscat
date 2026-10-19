# Dense setup-versus-interval scaling

The finite [design](../../docs/physics/co-numerov-dense-scaling.md) and
[producer contract](numerov-dense-scaling-contract.json) were frozen/pushed at
**`1d8cfacb2f86254c8035af5d74bba3c6869bd301`** before execution. All **166/166
checks pass / original exits0/0**. The existing coupled numerical modules remain
byte-identical to`4a9baa1`. This result establishes synthetic timing scaling and
accuracy in the exercised regimes; native CO precision/conditioning remain open.

## Scope and accuracy

Exactly24 dense298-channel propagations /1344 radial intervals: free and
noncommuting inverse-square/inverse-cube potentials, N=8/32/128 at fixed
h=.015625 bohr and R0=18 bohr, binary64 and64-significand-bit longdouble,
two repetitions in the predeclared reversed order. Six80/110-digit independent
discrete-oracle pairs use a reduced active block only in the reference; the
producer propagates the full matrices. No new CO/native/inner solve or unique
energy is added.

Worst binary64/extended K arithmetic errors are **7.66141e-15 /2.04508e-18**.
Maximum extended solve backward error **8.72570e-20**, one correction at most,
passes the unchanged32-eps gate. Physical unitarity/reciprocity errors are
<=1.28728e-10 in this tail regime. Every repeated K decimal reproduces exactly;
both free N=8 arithmetic formats reproduce the older coupled dense rows exactly.

The independent [saved verifier](numerov_dense_replay.py), frozen at`ff2d258`,
recomputes all166 decisions, six oracle pairs, four least-squares timing fits,
2,040 solve roles and24 timer partitions. It also audits returned dtype. Original
per-call logs omit shapes despite the prospective design requesting them; the
fixed numerical source specifies dense298 matrices. Per-call shapes are not
independently measured. This evidence limitation is preserved in the separate
[public replay contract](numerov-dense-public-replay-contract.json).

## Measured scaling

Two-repeat mean kernel wall seconds:

| Potential / arithmetic | N=8 | N=32 | N=128 | Fitted setup intercept / s | Marginal32→128 / s per interval |
|---|---:|---:|---:|---:|---:|
| Free / binary64 | .058512 | .110060 | .325703 | .039754 | .00224629 |
| Free / extended | 4.654522 | 7.897791 | 20.883044 | 3.571403 | **.13526305** |
| Varying / binary64 | .135760 | .244300 | .679242 | .099450 | .00453064 |
| Varying / extended | 9.603472 | 16.116818 | 42.168961 | 7.432564 | **.27137648** |

All three prospective timing hypotheses are **true**: all four slopes are
positive and affine residuals<=20% (actual maximum1.539%); free extended N=128
recurring solves exceed Taylor+tail; varying/free extended marginal ratio is
**2.006287**, within[1,3]. Largest two-repeat extended spread is.017413 s.
The affine intercept is a fitted aggregate, not a direct measurement of setup.

At free extended N=128, recurring solves take **17.107545 s**, versus **2.892569 s**
for Taylor+tail. Solve call counts follow N+2 for free and2N for varying W.
Measured direct recurring rates are.13470508/.26906906 s per interval; they agree
closely with marginal rates. Recurring mixed solves are now the priority over
fixed setup. Their LU/residual sub-operations need explicit profiling before an
arithmetic primitive is declared the proven optimization target.

For the **synthetic**60 finest400-bohr solves, affine/marginal estimates are
**22.96–46.08 single-worker days**, direct recurring-only estimates22.87–45.68
days. These replace the old setup-inclusive99-day guess as synthetic cost
evidence. They are not measured CO completion forecasts: thresholds, conditioning,
potentials, longer intervals and concurrency differ. No expensive CO successor
is released by this result.

## Resources and next decision

Producer finishes **2026-10-19T13:48:02.166561UTC**: **207.106276 s /206.610400
sampled CPU s /147,623,936-byte kernel peak**,1,033 samples at0.2 s. One worker,
reserved physical cores0–3/4 GiB, single-threaded BLAS. The180-s forecast was
exceeded by15.1%; the300-s upper and600-s cap were respected. All67 frozen
sources,242 installed artifacts and prior coupled report verify before/after.
Persistent originals:`prepared/numerov-dense-scaling-20261019/` under the campaign
root on Sadaharu.

Next isolate residual-product versus binary64 LU work within representative
mixed solves, then select a compiled kernel under a separately frozen numerical
and resource contract. Retain`typed_coupled_k` as the differential reference,
preserve64-bit-significand arithmetic, backward-error gates and existing analytic
controls, and requalify performance/accuracy before an interior CO pilot.
Original interior **false/exits1/1**, candidate-qualified14 identities and
historical132/594 coverage/462 uncomputed remain.

Public evidence and replay status:
[hash-pinned payloads](qualification-evidence/numerov-dense/README.md),
[saved verification](verification-numerov-dense-saved-20261019.json),
[fresh-public verification](verification-numerov-dense-public-20261019.json).

## Exact public numerical replay

The separate verifier/replay scope was frozen/pushed at **`ff2d258`**. Full
reexecution from fresh public bytes finishes **2026-10-19T13:55:08.961302UTC /
exits0/0**: every decimal and binary64 K component in all24 rows, all166
qualification decisions and all six oracle pairs reproduce **exactly**. Timing
fits are fresh observations; all three hypotheses also repeat, but their identity
is not a numerical replay gate.

Replay costs **207.088214 s /206.616432 sampled CPU s /146,092,032-byte kernel
peak**,1,033 samples at0.2 s, reserved cores0–3/4 GiB, below210-s forecast/300-s
upper/600-s cap. All67 sources, one separate verifier,242 installed artifacts
and prior report verify before/after. The8-file replay collection remains under
`prepared/numerov-dense-public-replay-20261019/`. All **14 final public payloads**
fresh-verify; local public-source tests pass **85 with two wider-format skips**,
exercised by remote producer/replay.

The [13:56:02UTC capture](progress-numerov-dense-result-20261019.json) reverifies
17 collected hashes,67 sources,242 installed artifacts and prior report. No
scientific containers remain;60.74 GiB persistent disk is free. Public evidence:
`https://data.qscat.org/co-numerov-dense-scaling-20261019/`.
