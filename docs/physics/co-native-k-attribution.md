# Stretched 1-eV K-discrepancy attribution

## Decision and preserved failure

The [six-case study](../../projects/ukrmol_co/NATIVE_FULL_RECONSTRUCTION.md)
stopped on absolute K vetoes at stretched 1 eV in B1/B2, despite passing S,
phase and invariants. Public reexecution reproduces both rejections exactly.
The large K eigenvalue, about −185, explains coordinate sensitivity but does
not identify the underlying propagation/reference error. The original study
stays false, exits 1/1, with 462 energy identities uncomputed.

The saved 0.025→0.0125-bohr step pair has K changes 1.05931e-5/1.22856e-5.
A diagnostic fourth-order Richardson extrapolation still differs from native32
by 7.46229e-6/5.28680e-6. This is a prediction under an unverified asymptotic
model, not an error bound. It argues against an automatic spacing ladder.

This new fixed experiment isolates native subrange, independent spacing,
endpoint and finite-tail-order sensitivities at **only these two identities**.
The [contract](../../projects/ukrmol_co/native-k-attribution-contract.json)
pins original inputs, saved references, failed bundles and public prerequisite.
The 298-channel / 2048-pole / 48-retained-component model, original orbital
average, two multipoles, thresholds, units and Numerov equations stay fixed.
No new inner calculation, pole interpretation or full-grid continuation occurs.

## Native-reference gate

For each sector create a one-energy NRANGE32 deck, then NRANGE64 only after
the NRANGE32 slice passes. The original native energy record is
0.07349999999999998 Rydberg. Its requested-eV literal is
`0.99999999999999978`, whose double multiplication by 0.0735 reproduces those
bits. Change only NESCAT 99→1, the first EINC literal, and the inserted NRANGE
line; preserve the 0.05 increment, comments, output units (B1 921, B2 922),
100-bohr native endpoint, NLEG10 and Gailitis defaults.

The slice must have exactly one complete record/EOF, identical energy bits,
unchanged interface metadata, printed 32 subranges and K/S/phase differences
≤1e-10 against the saved full-grid native32 record. Any failure stops that
case before NRANGE64. Preserve the original subprocess exit and partial log.

The resulting 16/32/64 native sequence uses the saved16/32 records and new64.
Require each final native32→64 K/S/phase difference ≤1e-6 before independent
controls are authorized. Report both successive changes. Small changes are
conditional resolution evidence, not an absolute native-error bound.

## Independent controls

After both native gates pass, perform exactly five additional solves per sector:

| Step / bohr | Endpoint / bohr | Tail order | Isolated control |
|---:|---:|---:|---|
| 0.00625 | 400 | 8 | Third spatial resolution |
| 0.0125 | 200 | 8 | Shorter propagation and endpoint |
| 0.0125 | 800 | 8 | Longer propagation and endpoint |
| 0.0125 | 400 | 4 | First tail-order rung |
| 0.0125 | 400 | 6 | Second tail-order rung |

Reuse the checksum-pinned 0.025/400/8 and 0.0125/400/8 failed-study rows.
Thus spatial, endpoint and tail-order families each have three rungs, with
one factor changed in each family. The third-step point is compared with
native64. Also report each other control's native64 comparison so endpoint
and tail sensitivities cannot be hidden by choosing the best result afterward.

Separate final-step, endpoint400→800, order6→8 and finest/native64 K/S/phase
vetoes remain **1e-6**. Every newly computed matrix must pass separate unitarity
and unsymmetrized S reciprocity gates of 1e-6. Report S spatial improvement
≥8 subject to the existing 1e-8 floor, both successive changes in every family,
signed large-eigenchannel projections and exact Cayley-identity residuals.
No matrix is symmetrized or fitted to the reference.

An attribution experiment can finish false: its measurements still bound which
tested sensitivities matter. A passing point-resolution verdict qualifies only
these two identities and controls. It does not replace the original rejection
or establish a unique error source. Numerical nonfinite/exception/incomplete
records leave the affected verdict unset.

## Cost and stopping

The saved fine solves cost 114.72–115.73 s each for 30560 intervals. The ten new
solves total 398720 intervals, forecasting about **1510 CPU-seconds**; four
single-energy native invocations conservatively add about 80 CPU-seconds.
Four-worker ideal throughput is 398 s. Allow **600 s expected / 1200 s upper /
1800 s hard cap**, plus a 180-s per-native-stage timeout. These are forecasts
from saved measurements, not a benchmark of this new batch.

Use four single-threaded workers on Sadaharu physical cores0–3, a **4-GiB**
container limit and a fresh persistent `/home` directory. Native case bundles
run two at a time; independent controls use at most four concurrent solves.
Verify topology/workloads, the pinned image, all 242 engine artifacts and source/
input/reference digests before launch and afterward. Store every native output,
unsymmetrized K, S components, diagnostics, verdict and per-solve journal.

On a failed/unset native gate or individual invariant gate, stop new submissions;
already executing fixed bundles finish and remain evidence. The family gates
are evaluated once the predefined controls are present. There are no retries,
further rungs, relaxed budgets or automatic successors.
