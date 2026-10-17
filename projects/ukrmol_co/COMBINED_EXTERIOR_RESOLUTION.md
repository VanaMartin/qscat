# Combined exterior resolution control — 17 October 2026

**The separately frozen degree128/256/512 study passes all sixteen checks on
Sadaharu.** All **432 paired-method cases + four analytic controls** complete,
with original scientific command/container exits **0/0**. The
[parent study](COMBINED_EXTERIOR_BENCHMARK.md) remains a rejection with exits1/1.

The [design](../../docs/physics/co-combined-exterior-resolution.md) and
[contract](combined-exterior-resolution-contract.json) were frozen at
**`db2158a` before execution**. The runner and both numerical methods are
byte-identical to the parent. The only numerical-setting change is the spectral
family, shifted from64/128/256 to128/256/512 across the entire matrix. No budget,
model, sheet, energy, endpoint, angle or ODE tolerance changes.

## Evidence and limits

| Check / diagnostic | Measured maximum | Frozen budget |
|---|---:|---:|
| Spectral128→256 change | 7.9099065827e-7 | Coarse→middle diagnostic |
| Spectral256→512 change | **4.2812918312e-8** | **1e-7**, plus tenfold improvement |
| Finest independent-method difference | 4.3498113784e-8 | 1e-7 |
| Finest analytic-limit error | 1.5963680366e-9 | 1e-7 |
| Spectral domain/contour spread | 5.2061942421e-8 | 1e-7 |
| ODE middle→finest change | 2.4087410709e-8 | 1e-7 |
| ODE domain/contour spread | 1.6016126141e-9 | 1e-7 |

The original three failing configurations, ε=30,L=80,θ=20°, now have
degree256→512 changes **1.139e-12 / 1.514e-12 / 1.044e-10** for
d=0/0.4/1.2. The new maximum is at **ε=15.9,d=1.2,L=40,θ=35°**, rather than
at the original high-energy long-domain miss. Degree512 increases the maximum
condition estimate to **4.69163e8**, and residual to **1.49153e-8**. Finest
method/domain discrepancies are larger than the parent's degree256 values,
although within the frozen budgets. This supports resolution of the measured
failure, not an assertion that arbitrary further degree growth improves accuracy
or that the successive difference is an absolute error bound.

Both real-axis transpose-symmetry and flux-positivity checks pass. Threshold-
straddling and lower-half-plane energies retain the declared outgoing sheets;
all rotated momenta decay. Omitted dipole, missing derivative rotation and wrong
incoming-sheet controls retain differences0.180034/2.14709/5.02831.

See [numerical verification](verification-combined-resolution-20261017.json)
and [checksum-pinned public evidence](qualification-evidence/combined-resolution/README.md).

Fresh [public reconstruction](verification-combined-resolution-public-20261017.json)
verifies **nine payloads / ten source files**, all436 cases, all sixteen checks
and scientific exit0. Maximum all-rung/finest host-public Y-component drifts are
**4.18531e-8 / 2.83234e-8**. The replay spectral final change is3.17634e-8 and
domain/contour spread4.01043e-8, both passing. Six focused combined-exterior and
native-boundary analytic tests pass before the run.

## Resources and next decision

Sadaharu costs **33.6106 seconds / 0.1292 GiB sampled kernel high-water**,
166 samples at0.2 seconds; sampled CPU usage reaches33.064 seconds. Sampling
is not a post-exit peak. Physical cores0–3 are reserved, one BLAS thread,4-GiB
limit,300-second timeout. This is below the prospectively estimated90 seconds.
Persistent source/output directory:
`/home/kooza/ukrmol/co-sa-20261006/prepared/combined-resolution-source-20261017/`.
Source archive SHA256:
`0288572953b78abc8ab7b63ba0506b074baa31b19b532af56605d6f41a8090a1`.
The existing pinned image and checksum-verified mpmath wheel are reused.

The next decision is a prospective **independent saved-native open-open K
reconstruction**. It must use the source-derived zero-Bloch/Hartree/Rydberg
conventions and test Gailitis amplitude normalization, all saved open-channel
matrix entries, finite-radius/step/tail sensitivity and closed-channel effects.
The defined two-channel inverse-square toy pass does not qualify the native
298-channel inverse-square/inverse-cube exterior, native K or CO poles. No new
inner solve is justified by missing data; production remains gated.
