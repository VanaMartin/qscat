# Cross-anchor corner and endpoint qualification

**Status: analytic endpoint prerequisite passes; the frozen CO job is submitted
and completion is pending.** Source/contract freeze **`7ea6243`** was pushed
before submission to Sadaharu. The original increment two-point pass, earlier
scientific rejections and their exits remain separate records.

The [design and stopping rules](../../docs/physics/co-numerov-corner-control.md)
answer two finite questions:

1. Does the unchanged increment formulation pass the previously qualified
   energy-extreme settings across all compressed/equilibrium/stretched B1/B2
   sectors?
2. Do all 200/400/800-bohr endpoints at the final spacing agree with the saved
   native reference, including the two stretched 1-eV large-K identities?

## Prerequisite: endpoint truncation versus floating-point drift

The [analytic endpoint study](verification-numerov-endpoint-analytic-20261019.json)
passes **32/32 checks**, original exits0/0. Twelve solves cover both scalar and
rotated open/open/closed boundaries, two spacings and three endpoints. Independent
80/110-digit continuous and scalar discrete Numerov oracles separate error
sources. At800 bohr/h=.00625: continuous scalar K error1.26839e-7, discrete
truncation1.25185e-7, candidate/discrete drift1.65453e-9. The coarser rows are
reported without relabelling them as accurate. The measured analytic cost is
15.0354 s /14.451588 sampled CPU-seconds /46,170,112-byte kernel high-water;
all36 sources and242 installed engine artifacts verify before/after.

## Frozen CO contract

[numerov-corner-contract.json](numerov-corner-contract.json) pins24 raw files,
12 native16/32 references, two native64 references, the original increment report
and two verified prerequisites. No physical input, propagator, tail, decoder or
native engine is changed. Exact binary energies and unsymmetrized K are retained.

- **Fourteen identities:** all twelve geometry/sector low/high corners plus the
  two stretched1-eV identities. No new unique-energy claim follows from them.
- **78 maximum new solves**, six unchanged reused rows. Initial two endpoint
  bundles must pass before the twelve six-solve corner bundles.
- **21 independent case checks:** separate K/S/phase gates for step, endpoint,
  tail order and each of the three finest-endpoint/native comparisons; plus
  unitarity, reciprocity and the declared S-step trend rule. The corners also
  require their separate native16→32 reference-control checks.
- **Eight physical cores4–11**, eight single-threaded workers,8 GiB. Cores0–3
  remain available for smaller jobs. Expected1500 s /upper2400 s /cap3600 s,
  based on2364640 maximum new intervals and the preceding measured timing.
- First false/unset stops new submissions; fixed in-flight bundles finish. No
  retry, relaxed tolerance, new rung or automatic full-grid continuation.

The [saved-decision verifier](numerov_corner_replay.py) uses separate K/S/phase
expressions and decodes pinned native references. It does not call the producer
evaluator or propagator. Numerical replay and public publication follow a
completed finite result, with their scope/cost recorded separately.

The [prospective public replay subset](numerov-corner-public-replay-contract.json)
contains each produced corner's finest400 row and both produced stretched1-eV
finest800 rows, at most14 reexecutions. Selection is independent of scientific
verdict. All other rows receive saved-byte and decision verification; the subset
must not be described as reexecuting all78 solves. The replay has a450-s
forecast/1200-s cap on eight physical cores/8 GiB. Verification support includes
two additional tests of the separate saved expressions and the analytic phase
defect coefficient; numerical producer sources remain pinned to their archive.

## Evidence locations and interpretation

Remote campaign root: `/home/kooza/ukrmol/co-sa-20261006`.

- Analytic source/results/resources: `prepared/numerov-corner-analytic-20261019/`.
- Raw/reference preflight archive: `prepared/numerov-corner-preflight-20261019/`.
- Frozen producer/owner/resources/results: `prepared/numerov-corner-co-20261019/`.

The [previous two-point study](NUMEROV_INCREMENT_CONTROL.md) still passes its
defined24 gates; its step.0125/endpoint800 native discrepancies1.59e-6/1.46e-6
remain preserved. This study explicitly tests every finest-step endpoint rather
than waiving that observation. A pass would qualify only the declared corners
and two endpoint identities, with conditional native reference evidence. Full
594-matrix reconstruction, an absolute native error bound, complex poles,
channel/electronic/continuity/neutral convergence and production remain gated.
