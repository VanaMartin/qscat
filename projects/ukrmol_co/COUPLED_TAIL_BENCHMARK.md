# Separate coupled-tail and threshold-sheet controls pass — 17 October 2026

**All twenty frozen checks pass:** 189 dipole-tail S matrices (135 real and
54 complex), 324 unequal-threshold shell S matrices and 27 outgoing pole solves.
All three accuracy levels and matching radii are retained; every root succeeds
without reseeding. Original scientific command/container exits are **0/0**.

The [design](../../docs/physics/co-coupled-tail-benchmark.md),
[contract](coupled-tail-contract.json), Python sources and three analytic tests
were frozen at **`047524e` before the complete matrix**. No native CO data is
read. These are two complementary physical controls, not a combined CO propagator.

## Exact angular tail and matrix propagation

The first model is the l=1,2, |m|=1 projection of a centrifugal plus
`−d cosθ/r²` potential outside a hard wall at a=5 bohr. The two channels have a
common threshold. Dipoles are d=0/0.4/1.2; the angular element is independently
checked as **1/√5**. A real orthogonal eigenbasis gives an exact noninteger-order
Hankel solution. The oracle directly evaluates the wall condition at 60/90
digits without radial integration. The numerical route propagates the full
fundamental solution matrix in the original angular basis and matches to exact
tail functions at B/a=1/2/4. The zero-dipole hard-sphere limit is checked separately.

This validates full matrix propagation, angular/asymptotic phases and matching
against analytic radial evolution. Both routes share the derived tail model;
this does not validate arbitrary unknown tail coefficients or a CO channel cutoff.
At complex q, incoming waves are analytic Hankel2 functions rather than the
complex conjugate of outgoing Hankel1. No complex-energy unitarity test is imposed.

## Unequal thresholds, closed elimination and outgoing poles

The second model has radial s channels, thresholds `T=diag(0,16)` and a real
shell jump matrix `Λ=[[10,c],[c,12]]`, c=0/1/3, at x=r/a=1. Threshold energies
are in ε=2a²E. T and Λ do not commute for nonzero c. Twelve real ε include
**15.999 and 16.001**, immediately below/above the second threshold; exact
threshold is excluded because unit-flux open waves are singular there.

The independent oracle uses `L=diag(k cot k)+Λ` and explicit analytic Schur
elimination of the closed channel. The numerical route instead integrates
regular solution columns, applies the jump and solves a boundary matching system
with incident/outgoing open waves and only a decaying closed wave. Relative
`1/√k` flux factors matter when two unequal-momentum channels open.

For poles below the second threshold, the continuation fixes
`k0=q`, `k1=+i sqrt(16−q²)`. Numerical root solves use only the propagated Jost
determinant from the fixed seed **q=π−0.3−0.2i**. An independent high-precision
analytic Jost determinant supplies reference roots, and c=0 matches the scalar
benchmark's Lambert-W pole. Flipping the incident-channel momentum at an outgoing
pole gives a nonzero oracle determinant (minimum **1.587**), excluding numerator
cancellation. This is the declared two-channel sheet, not a general multichannel
sheet search.

Representative toy poles at a=5 bohr (these are **not CO predictions**):

| Shell coupling c | Position / eV | Full width / eV |
|---|---:|---:|
| 0 | 4.50404409 | 0.41663813 |
| 1 | 4.50024271 | 0.42016035 |
| 3 | 4.46866711 | 0.44991970 |

The c=3 closed coupling changes the width by about **8%** relative to c=0
despite only one asymptotic channel being open near this pole. Closed channels
therefore cannot be discarded on the basis of open-channel count alone.

## Three-rung numerical convergence

| ODE relative tolerance | Maximum tail S entry error | Maximum threshold S entry error | Maximum pole q error |
|---|---:|---:|---:|
| 1e-5 | 3.11520e-5 | 4.34323e-5 | 2.18550e-6 |
| 1e-8 | 3.13007e-8 | 7.65938e-8 | 3.68403e-9 |
| 1e-11 | **3.13693e-11** | **6.89440e-11** | **3.02832e-12** |

The frozen finest S/q budgets are each 1e-7. Errors are against independent
oracles, not successive differences. Finest real-axis unitarity/reciprocity
errors are **2.140e-13 / 2.192e-13**, below 1e-8. The largest matching-system
condition number across all rungs is **90,536**, reported rather than hidden.
Finest position/full-width errors are **9.126e-12 / 5.328e-12 eV**, below
the 0.1-meV reproduction budgets. Finest Jost residual is **7.261e-12** and
matching-radius q spread **1.336e-15**. The 60/90-digit reference difference is
**1.018e-60**; analytic Jost and decoupled Lambert errors are below 6e-90/2e-91.

## Consequential negative and tail-cutoff controls

At c=3, selecting a **growing** closed-channel wave at ε=8 changes S by
**0.161019**. Since the resulting effective real boundary can still yield a
unitary scalar S, unitarity does not select the physical decaying sheet.
Omitting relative flux normalization at ε=16.1 changes S by **0.100306**.
Both exceed the fixed 0.001 negative-control discrepancy threshold.

For d=1.2,q=1.3, propagating the correct dipole tail only to B and then matching
to uncoupled centrifugal waves gives:

| B/a | Error versus uncut polar-tail S | Unitarity error |
|---|---:|---:|
| 8 | **0.0488767** | 2.561e-16 |
| 32 | **0.0121008** | 1.568e-15 |
| 128 | **0.00301035** | 5.735e-15 |

The truncation effect decreases by a factor **16.2**, yet remains visible even
at B=640 bohr for this toy model. These are absolute S-entry errors, not phases
or native CO uncertainty bounds. They demonstrate why a unitary matrix and a
large free-wave matching radius do not by themselves qualify a polar boundary.

## Evidence, cost and remaining gate

[Nine public payloads](qualification-evidence/coupled-tail/README.md) preserve
**eleven frozen source files**, the dependency wheel, original resource/exit
records and all matrices/root statuses/cutoffs. Fresh public replay reproduces
all twenty checks and every case; maximum host/public S-component, q and
position/width drifts are **7.84e-13 / 7.60e-14 / 2.31e-13 eV**. See
[quoted verification](verification-coupled-tail-20261017.json) and
[public replay verification](verification-coupled-tail-public-20261017.json).

Sadaharu cost: **1.766 s wall / 0.050 GiB sampled kernel high-water**, seven
samples at 0.2-second intervals, reserved physical cores 0–3, one BLAS thread,
4-GiB limit. The sample is not a post-exit peak. The existing checksum-pinned
mpmath wheel is reused read-only; its earlier dependency-probe exit 1 is preserved
with explicit original provenance, not asserted as a new failure. Persistent
sources/outputs/logs are in `prepared/coupled-tail-source-20261017/`.
Source archive SHA256:
`96e6abdd41403a82c1a7da816d5a77bad145ee661e9a3360808d35090dffd09b`.
**31 focused tests pass.** No native file is read or native solve launched.

**Next:** combine noncommuting unequal thresholds and polar tails under independent
real/complex-energy controls, and derive native boundary amplitude/Bloch and
Hartree/Rydberg normalization. The two separate passes do not qualify that
combined problem. Full saved-native K-matrix reproduction with radius/step/tail
controls is still required before CO poles. Five-A2 interface failure, extraction,
channel/model and continuity rejections remain preserved; production is gated.
