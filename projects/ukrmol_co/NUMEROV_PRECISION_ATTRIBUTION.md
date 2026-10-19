# Interior analytic precision attribution

**Diagnostic pass,46/46 checks / original diagnostic exits0/0.** The earlier
interior prerequisite remains **false / exits1/1**. This calculation diagnoses
its failure; it does not qualify new CO points or repair the propagator.

## Fixed question and source freeze

Source/contract/design frozen and pushed at **`aa917aa`** before submission.
The [design](../../docs/physics/co-numerov-precision-attribution.md) derives a
high-precision matrix discrete Jost oracle retaining the original sixth-order
Taylor initialisation/matching and Numerov difference equation. Short-grid direct
high-precision mixed-channel recurrence and the independent scalar discrete
oracle pass six tests;17 focused tests pass before launch.

The [finite contract](numerov-precision-attribution-contract.json) reuses the
twelve archived compressed analytic rows, reexecutes two critical baseline rows
exactly, and propagates eight reduced-boundary controls in rotated/eigenchannel
frames. All coordinates were visited previously: no new h or endpoint. Total
**ten new analytic solves /2712320 radial intervals**, one worker, reserved four
physical cores0–3/4 GiB, single-threaded BLAS. No CO/native/inner solve.

## Result: distinguish roundoff from truncation

The high-precision discrete oracle separates signed `Kfloat−Kdiscrete` from
`Kdiscrete−Kcontinuous`, using each actual binary64 boundary Q. All80/110-digit
precision and scalar consistency checks pass. The two critical baseline K
matrices reproduce every component exactly. The full rotated open/closed model
at h=.0015625 gives:

| Endpoint / bohr | Floating-point K drift | Discrete truncation K error | Drift/truncation |
|---|---:|---:|---:|
| 200 | 1.13485700e-6 | 9.42110195e-9 | 120.46 |
| 400 | 1.13086255e-6 | 1.97315201e-8 | 57.31 |
| 800 | 8.86922890e-7 | 4.04873632e-8 | 21.91 |

The predeclared hypothesis that drift dominates truncation by>=5× at the two
rejected200/400 points is **true**. Signed cancellation means continuous error
can be smaller than floating-point drift; a passing total error does not erase
the separate drift comparison. This establishes substantial matrix roundoff in
this analytic regime, not a unique elementary operation or native CO error bound.

## Equivalent representation controls

Physical closed-channel elimination and open eigenchannel rotation preserve
continuous K/S/phase within **4.01e-12**, below the1e-8 equivalence budget.
Each finite-step model has its own discrete oracle; finite-step elimination is
not assumed identical. All eight new controls pass invariants.

Common-frame maximum K drift at h=.0015625:

| Endpoint / bohr | Reduced rotated frame | Reduced diagonal frame |
|---|---:|---:|
| 200 | 5.90927556e-7 | 1.83269706e-8 |
| 400 | 5.27912891e-7 | 3.75160024e-8 |
| 800 | 3.38504404e-7 | 7.73765116e-8 |

The all-endpoint tenfold-improvement hypothesis is **false**: the diagonal800
row is only4.37× better, although every diagonal drift is below1e-7. Retain the
false conjunction as declared. Representation sensitivity supports an arithmetic
contribution; neither an exact basis change nor eliminating the free closed
channel is a qualified general solution for coupled inverse-power CO potentials.

## Resources, preservation and next decision

Finished2026-10-19T12:27:38.384530UTC: **61.050941 s /60.487653 sampled CPU s /
47,837,184-byte recorded kernel high-water**,303 samples at0.2 s. Within120-s
forecast/240-s upper/600-s cap. All51 frozen sources,242 installed artifacts and
the pinned original report verify before/after. The10-file collection retains
owner, original results, exits and sampled resources under
`/home/kooza/ukrmol/co-sa-20261006/prepared/numerov-precision-attribution-20261019/`.

The [separate saved-decision verifier](numerov_precision_replay.py) reconstructs
all oracle/error/check and hypothesis decisions from pinned bytes without
propagation; [verification](verification-numerov-precision-attribution-20261019.json)
passes. It shares the tested oracle and uses separate scattering/error expressions.

The next useful numerical control is increased arithmetic precision on the same
failed analytic grids, with separately costed recurrence/matching attribution.
No spacing ladder or tolerance relaxation is justified. Original interior
false/exits1/1 and all older rejections stay; candidate qualification remains14,
historical coverage132/594/462 uncomputed. No CO successor is submitted.
