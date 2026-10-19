# Fixed-grid Numerov algebra control

**Status: analytic qualification, the ten-solve CO control and public numerical
replay pass.** All24 case checks pass with original
runner/container exits0/0. The original two-point and six-case rejections
remain false/exits1/1. This qualifies only the defined control at two identities.
The [mathematics and prospective gates](../../docs/physics/co-numerov-increment-control.md)
derive an identity-defect recurrence that avoids cancellation in initialisation,
propagation and matching. `exterior_numerov.py` remains the original oracle.

The [analytic verification](verification-numerov-increment-analytic-20261019.json)
passes16 frozen checks and49 focused tests. At the finest scalar grid, original
versus high-precision discrete K drift is1.34450e-4; the increment form gives
9.82188e-10. The separate continuous discretisation error is6.11308e-8.
The rotated two-open/one-closed case and short noncommuting coupled controls
also pass their defined criteria. This validates the tested analytic regime,
not native CO propagation or a complex pole.

Analytic source freezes: `381a201`, dependency-only supplement `e10e387`.
The first attempt exits1/1 before calculation because the tools image lacks
mpmath; its original logs and archive are preserved. The separately pinned
wheel enables the same numerical source in a fresh directory. The passing run
costs10.4293 s /9.826654 sampled CPU-seconds /46,272,512-byte recorded kernel
high-water,50 samples at0.2 s;30 sources and242 engine artifacts verify
before/after. No native data or engine invocation occurs.

The [fixed-grid CO contract](numerov-increment-co-contract.json) permits only
the existing two stretched1-eV identities at the existing three spacings and
three endpoints: ten candidate solves,368160 intervals,500-s forecast/800-s
upper/1200-s cap,four physical cores/4 GiB. All final K/S/phase and invariant
budgets remain1e-6. Stop submissions on false/unset invariants; in-flight
bundles finish. No new energy, inner/native solve or automatic successor.

The [original two-point rejection](NATIVE_K_ATTRIBUTION.md) remains false,
exits1/1, as does the six-case reconstruction.132/594 unique-energy coverage,
462 uncomputed energies and all pole/production gates remain unchanged.

## Completed CO result — 19 October, 07:12 UTC

Source and contract were frozen/pushed at **`a8665f4`** before launch. All ten
candidate controls finish; their invariant stage and all **24/24** case checks
pass. The [separate saved-decision verification](verification-numerov-increment-co-20261019.json)
reconstructs the complete records, every gate and same-grid differences without
calling the candidate evaluator. No native executable or inner solve runs.

| Frozen final comparison | B1 max entrywise K change | B2 max entrywise K change | Budget |
|---|---:|---:|---:|
| Step0.0125→0.00625 at400 bohr | 6.54939498e-7 | 7.59042706e-7 | 1e-6 |
| Endpoint400→800 at step0.0125 | 7.90719326e-7 | 8.14519211e-7 | 1e-6 |
| Finest0.00625/400 versus native64 | 1.48708608e-7 | 1.08836815e-7 | 1e-6 |

All corresponding S/phase and unitarity/reciprocity gates pass. Finest/native64
S differences are2.64492e-10/6.56354e-10, phase differences3.05218e-10/
3.60833e-10 rad moduloπ. Original-form finest/native64 K differences were
7.25612e-5/7.25605e-5; the candidate reduces these measured discrepancies by
about488/667-fold on identical physical inputs, grid, tail and reference.

The two step changes now decrease by factors16.97/14.57 in K and15.99/15.97
in S, consistent with approximately fourth-order behaviour on this finite
sequence. Unlike the previous worsening sequence, raw S changes here satisfy
the eightfold improvement criterion without needing its1e-8 floor.

### Consequential limitation: endpoint800 is not an all-row pass

The frozen native comparison is the finest-step400-bohr row. It does **not**
assert that all ten rows agree with native64 within1e-6. In particular,
step0.0125/endpoint800 gives native K differences **1.59436743e-6 /
1.46472510e-6**, larger than the budget. Coarse-step rows also differ by about
1.2e-5. These are separate diagnostics, not substitutions for the declared
native gate. The endpoint family passes with only about19% margin in B2, and
its change grows with endpoint. Thus a longer endpoint is not automatically
more accurate; no absolute-error or endpoint-independent conclusion follows.

The analytic discrete oracle demonstrates floating-point drift in the original
formulation. The same-grid CO improvement supports cancellation as a major
contributor to the earlier two-point vetoes. This combined initialisation,
recurrence and matching change does not isolate a single offending operation,
prove an absolute native64 error bound, or establish298-channel complex-plane
propagation. The original source/results remain untouched.

## Measured resources and reproducibility

CO runtime **390.672689 s /6.5112 min**, sampled CPU **1334.834105 s**, recorded
kernel memory high-water **252,329,984 bytes /0.235001 GiB**,1950 samples at0.2 s.
Four physical cores0–3/four single-threaded workers/4-GiB cap. All34 archived
sources and242 installed engine artifacts verify before/after; eight raw inputs,
fifteen attribution outputs and eight predecessor references verify in the
runner before/after. The20-file collection retains all ten candidate outputs.

[Public evidence](qualification-evidence/numerov-increment/README.md) includes
the original failed-import and passing analytic attempts, the dependency wheel,
CO source/inputs/references and independent verifier. Fresh public downloads
pass saved-decision verification and all49 focused tests. A separate replay is
restricted to exactly the same ten candidate controls, four physical cores/
4 GiB,500-s forecast/1200-s cap. Replay results and exits remain separate from
the producer and original rejection records.

The public numerical replay finishes at07:20:44.170194 UTC with **exits0/0**.
Every K and stored S component reproduces exactly, with identical case checks
and passing decisions. Runtime390.215354 s /1332.684159 sampled CPU-seconds /
239,140,864-byte recorded kernel high-water,1947 samples at0.2 s. The17-file
replay collection and separately published owner/report/resources preserve
that evidence. All **20 final public payloads** are freshly downloaded and
checksum-verified after supplementation, including the all-configuration
diagnostic that retains the endpoint800/native failure.

The [07:21 UTC result capture](progress-numerov-increment-result-20261019.json)
reverifies55 collected remote hashes,34 sources,242 installed artifacts and
all8/15/8 input/reference pins. **No scientific container is running**;
61.68 GiB persistent space remains. Original failed-import evidence and both
original scientific rejections retain their separate verdicts/exits.

The next decision is a separately frozen cross-anchor/energy qualification of
the increment form, using saved matrices and a meaningful endpoint/precision
check. The present study does not add unique-energy coverage or release a
full-grid run. No automatic successor is submitted by this contract.

## Later finite corner / endpoint qualification

The separately frozen [cross-anchor endpoint study](NUMEROV_CORNER_CONTROL.md)
now passes294 case checks and36 native-reference controls at14 identities,
exits0/0. It preserves this study's coarser h=.0125/800-bohr discrepancies.
At the explicitly contracted h=.00625, all200/400/800-bohr native comparisons
pass; the800-bohr K differences are2.15388e-7/5.59271e-8. This is new finite
qualification at a declared finer spacing, not a rewritten original decision
or an absolute-error/full-grid certificate.
