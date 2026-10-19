# Native extended residual primitive

The [design](../../docs/physics/co-numerov-native-residual.md), kernel and
[finite contract](numerov-native-residual-contract.json) freeze/push at
**`400e7cb78935f3c377c6279bd94960aadbf626df`** before execution. All **69/69
checks pass / original build and runtime exits0/0**. The Rust/PyO3 boundary
and C longdouble arithmetic shim preserve64-significand-bit arithmetic on the
qualified little-endian x86-64 format. The NumPy residual expression is retained
as the differential reference. The coupled propagator remains byte-identical
to`4a9baa1`.

## Finite accuracy and speed

Exactly13 cases: six archived A/B/accepted-X triples (three native extended,
three binary64 NumPy controls), four independent dense rectangular array layouts,
one2^-60 signal, one large-cancellation case and one zero case. Paired80/110-digit
dot references use exact binary inputs via integer ratios:12 preselected entries
per captured triple and every small-case entry, **220 scalar entries /13 pairs**.

All native/NumPy full output matrices agree **exactly in scalar bytes excluding
padding** in the exercised cases, including all298×298 captured residuals. Worst
native independent-oracle scaled error is **2.58889e-20**, below32 extended eps.
Full-matrix differential errors are zero. All captured accepted-X backward gates,
input-preservation, shape/dtype/ownership and finite-output gates pass. The2^-60
signal survives; binary64 coercion would lose it. Fortran-order, transposed and
negative-stride inputs pass, with newly owned outputs.

End-to-end times include Python checks, contiguous packing, byte copies, native
allocation/transposition and output ownership. Exactly one warmup and two measured
repetitions per backend/capture, with backend order reversed in the second:

| Extended captured solve | NumPy / ms | Native / ms | Speedup |
|---|---:|---:|---:|
| Varying potential | 63.8562 | 26.0073 | **2.4553×** |
| Varying defect | 63.5055 | 25.5948 | **2.4812×** |
| Free defect | 63.6205 | 25.6000 | **2.4852×** |

All three prospective>=1.5× timing hypotheses pass. Criterion's separate native
byte-interface benchmarks give95% intervals **261.70–262.60 microseconds** at64
and **25.012–25.106 ms** at298; they omit Python packing/checks and use dense
dyadic fixtures, so they are not end-to-end speedups. Each uses ten samples,
one-second warmup and a two-second measurement target. The298 case's estimated
measurement duration is2.7574 s; all observations remain in the public logs.

## Build and measured resources

The locked release wheel uses Rust1.97.1, GCC12.2, PyO3 0.29.3, `-O3`,
`-fno-fast-math -ffp-contract=off` for the C shim. Rust format, two unit tests,
clippy with warnings denied, release wheel build, both Criterion benchmarks and
**nine Python differential tests** pass. The capability probe builds/imports on
the Mac; its five native-format tests skip there. Local focused tests pass14
with six skips including the existing coupled wider-format check.

Build finishes **2026-10-19T14:34:16.816353UTC**: **31.716034 s /83.366155 sampled
CPU s /2,095,026,176-byte kernel peak**,157 samples at0.2 s. Four Cargo build jobs,
reserved physical cores0–3/4 GiB, below180-s forecast/360-s upper/900-s cap.
The existing layered test image is pinned by full image ID in the contract.

Runtime finishes **14:34:18.987109UTC**: **1.823027 s /1.218827 sampled CPU s /
46,571,520-byte kernel peak**,seven samples, one worker, same cores/memory and
single-threaded BLAS. Below15-s forecast/60-s upper/180-s cap. The original
tools image preserves the previous NumPy2.3.4/64-bit-significand runtime.
All85 frozen sources,242 installed engine artifacts and six operand hashes
verify before/after. Originals remain under
`prepared/numerov-native-residual-20261019/` on Sadaharu. The31-file collection
retains command logs, source/compiler/toolchain provenance, Criterion JSON,
resource samples and original exits. The wheel and extension are separately pinned.

## Interpretation and next qualification

This establishes a useful primitive on the declared operands and cancellation/
layout controls. It does not establish a propagator speedup or a new CO resource
forecast. The next experiment should integrate the primitive behind a separately
qualified mixed-solve path, preserve the Python oracle, bounded corrections and
32-eps gate, and reexecute independent analytic/coupled stress checks plus finite
dense interval comparisons. Only measured integrated precision and cost can
inform a separately frozen interior CO pilot.

Synthetic captured matrices contain a small active block; arbitrary CO
conditioning, cancellation and thresholds remain unqualified. Original interior
**false/exits1/1**,14 candidate-qualified identities and historical132/594
coverage/462 uncomputed remain. No propagation, native UKRmol or inner solve is
added by this primitive study.

[Public evidence](qualification-evidence/numerov-native-residual/README.md) ·
[Saved audit](verification-numerov-native-residual-saved-20261019.json) ·
[Public replay contract](numerov-native-residual-public-replay-contract.json) ·
[Public verification](verification-numerov-native-residual-public-20261019.json).

## Complete public replay and preserved preflight

Saved-metric verifier/full replay freeze at **`1ad062e`**. The first public replay
owner exits **1 before any container or arithmetic**, with`StopIteration` while
looking for a top-level extension binary. Maturin places the binary inside a
package directory. Inspection also establishes that the producer's field
`extension_sha256` hashes the imported package initializer, not the binary.
The whole-wheel pin was correct. Those original fields and failure remain literal.

The [artifact-layout supplement](numerov-native-residual-replay-preflight-supplement.json)
freezes/pushes at **`5b1c1fb`** before a separate successor. It permits exact nested
initializer/binary path checks and adds the true binary SHA256
`d977b8b652a17c4de629138bb22d153bee235cc7a268dc049349e867583be61e`.
It changes no numerical source, wheel, operands or comparison. The initial
preflight's scientific verdict remains **unset**, with zero container submissions.

The successor finishes **14:39:49.293474UTC /exits0/0**: all **13 numerical rows,
13 oracle pairs/220 scalar entries,23 full output numeric hashes and69 decisions**
reproduce exactly from fresh public bytes. Fresh speedups are **2.4517×/2.4964×/
2.4917×** and all three timing hypotheses pass independently of numerical identity.
Replay costs **1.772593 s /1.219100 sampled CPU s /46,776,320-byte kernel peak**,
seven0.2-s samples, reserved cores0–3/4 GiB; below5-s forecast/15-s upper/180-s cap.
All85 producer sources plus one verifier,242 installed artifacts and six operand
hashes verify before/after, with exact initializer and binary pins. Originals
remain at`prepared/numerov-native-residual-public-replay-preflight-successor-20261019/`;
the failed preflight directory is preserved alongside it.

All **17 final public payloads** fresh-verify, including the original failed
owner/log, supplement, successful owner, numerical outputs and resources. The
[14:41:45UTC capture](progress-numerov-native-residual-result-20261019.json)
reverifies41 collected hashes,85 sources and242 artifacts. No scientific container
remains;59.93 GiB persistent disk is free. Public:
`https://data.qscat.org/co-numerov-native-residual-20261019/`.
