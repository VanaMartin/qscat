# Bounded long-range saved-CO radial pilot

## Decision and fixed scope

The [near-boundary canary](NUMEROV_CO_CANARY.md) qualifies recurrence behavior only
through18.1984375 bohr. The new [finite contract](numerov-co-radial-pilot-contract.json)
tests cumulative mixed-precision accuracy, sampled conditioning and recurring
cost over18→200/400 bohr, including the unchanged Taylor endpoint and tail matching.

Reuse **stretched nine-A2 B1, energy index18**, exact archived Rydberg
`0x1.2d0e560418936p-4`. It has a large K entry≈141, the early canary's smallest
sampled reciprocal condition estimate and an independently retained native64
reference. The older binary64 increment solver already passes the
[corner qualification](NUMEROV_CORNER_CONTROL.md) at this identity. This is a
comparison of the new mixed/extended execution against existing evidence.

Five full propagations are authorized:

| Stage | Backend | Endpoint, bohr | Spacing, bohr | Intervals | Mixed solves |
|---|---|---:|---:|---:|---:|
| Initial | Native residual |200|.025|7,280|14,560|
| Initial | Native residual |200|.0125|14,560|29,120|
| Initial | Native residual |200|.00625|29,120|58,240|
| Initial | NumPy residual |200|.00625|29,120|58,240|
| Gated endpoint | Native residual |400|.00625|61,120|122,240|

Total **141,200 intervals /282,400 mixed solves**, including five matching solves.
Three spacings establish a finite convergence sequence. One full same-grid
NumPy-residual run tests accumulated kernel substitution. A sole longer endpoint
tests additional radial accumulation and matching stability. All four initial
rows finish before the endpoint is admitted; any false/unset gate prevents it.

## Numerical gates

Keep the original coupled recurrence, binary64 LU preconditioning, extended
residual/refinement, maximum three corrections and32-extended-epsilon backward
gate. This is mixed precision; the LAPACK factorization remains binary64.
Use all298 channels,2048 poles, exact saved energy, order8 analytic tail,
18-bohr boundary and1e-10-Rydberg threshold grouping. Preserve unsymmetrized K.

Each finest row separately passes **1e-6 K /S /phase moduloπ** against native64
and the archived same-grid binary64 increment result. Final step and200→400
differences pass those same gates. S step change must shrink by8 or reach1e-8;
each row passes separate1e-6 unitarity/reciprocity gates. Coarser rows diagnose
convergence. Backend fine200 rows require exact extended-K hash, complete solver
stats and sampled A/B/X hashes. Auxiliary rcond observations alone use the prior
32-binary64-epsilon relative rounding budget, retaining the original strict
canary replay failure1/1 as a separate record.

Sample initialization, both recurring systems at25/50/100/150 bohr (also200/300
for the400 row), final ratio and matching:11 samples per200 row and15 per400.
Independently reevaluate every sampled native residual/decision with NumPy;
require full residual equality and zero gate disagreements. All unsampled solves
retain the same residual acceptance gate. Retain native initialization/final
ratio/matching A/B/X:12 full triples and36 paired80/110-digit residual entries.
Require scaled oracle stability≤1e-60 and arithmetic error≤32×2^-63.

Positive finite one-norm rcond estimates describe sampled binary64 preconditioners.
They establish neither global long-range conditioning nor a full forward-error
bound. Saved native agreement remains conditional on its earlier subrange controls.

## Resources, provenance and interpretation

Sadaharu: **physical cores4–11, four single-threaded processes initially, one for
the endpoint,8 GiB**; cores0–3 remain available for small work. Reuse the qualified
wheel and digest-pinned runtime/test images. All24 raw files, native64, old corner
report, canary/corner public verification, source snapshot, wheel/components and
242 installed engine artifacts verify before/after. Originals persist under
`/home/kooza/ukrmol/co-sa-20261006/prepared/numerov-co-radial-pilot-20261020/`.

Early native/reference interval rates.124/.279 s suggest approximately
**6.12 single-worker CPU-hours**:2.26h initial critical path plus2.11h endpoint,
before overhead and host contention. Forecast **5h elapsed**, upper **8h**, hard
cap **12h**. This is a cost hypothesis that the pilot measures. Record sampled
resources, command/container/owner exits and unrelated host work. No retry or
automatic extra rung follows a failure, exception or timeout.

Full elapsed times include instrumentation. Per-role solve sums and25-bohr bins
exclude sampled solves from recurring throughput. Native/reference audits differ,
so full-runtime ratios are not a clean performance benchmark. One identity cannot
forecast all energies or qualify800-bohr execution.

The old interior analytic false/exits1/1,132/594 historical energy coverage and
all physical-model gates retain their scopes. This pilot adds zero unique energy
identities. Numerical verification and replay require a separately costed contract
after a substantive result.

## Execution

On the qualified Linux x86-64 scalar format, with existing native wheel and mpmath:

```sh
python -m projects.ukrmol_co.numerov_co_radial_pilot \
  --contract projects/ukrmol_co/numerov-co-radial-pilot-contract.json \
  --raw RAW --prior PRIOR --output NEW_DIRECTORY
```

PRIOR contains the hash-pinned `corner.json` and `native64.bin`. The output must
be new. Every row retains progress after each1024 solves and at matching; full
rows and comparison decisions checkpoint as each worker/stage completes.
