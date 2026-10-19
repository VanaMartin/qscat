# Cross-anchor corner and endpoint qualification

**Status: both finite CO stages and the public numerical subset replay pass,
original exits0/0; every saved decision verifies.** Source/contract
freeze **`7ea6243`** was pushed before submission to Sadaharu. The original
increment two-point pass, earlier scientific rejections and their exits remain
separate records.

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

## Completed finite result — 19 October, 08:50 UTC

Both initial endpoint bundles and all twelve corner bundles complete:
**14 identities, 78 new solves, six unchanged reused rows**. All **294/294 case
checks** and **36/36 native16→32 reference checks** pass, original runner/container
exits **0/0**. No native/inner engine runs. The
[independent saved-decision replay](verification-numerov-corner-saved-20261019.json)
checks all84 rows, exact energies, streamed/disk/report identity, hashes,
reference controls, separate S reconstruction and every case decision.

| Final comparison, worst across the14 identities | K | S | Phase moduloπ, rad |
|---|---:|---:|---:|
| Final spacing change | 7.59043e-7 | 4.73852e-8 | 1.41692e-7 |
| Endpoint400→800 at final spacing | 6.66791e-8 | 3.33170e-9 | 9.99560e-9 |
| Tail order6→8 | 4.34852e-11 | 2.13577e-12 | 2.21253e-12 |
| All final-spacing endpoints versus native | 2.15388e-7 | 6.50401e-9 | 1.94614e-8 |

All budgets remain separately1e-6. The raw S-step reduction is at least12.19
across every tested identity, so the trend criterion passes without its1e-8
floor. This is finite-sequence evidence, not a Richardson error certificate.
Corner native16→32 maxima are1.31247e-10 K /2.62219e-10 S /3.47599e-10 rad;
these are conditional reference controls, not absolute error bounds.

### The prior800-bohr observation is resolved at the declared finer spacing

| Stretched1-eV finest-step/native64 K | B1 | B2 |
|---|---:|---:|
| Endpoint200 | 1.52102e-7 | 6.26673e-8 |
| Endpoint400 | 1.48709e-7 | 1.08837e-7 |
| Endpoint800 | 2.15388e-7 | 5.59271e-8 |

The original h=.0125/800 discrepancies **1.59437e-6/1.46473e-6 remain unchanged**.
At the separately frozen final h=.00625, all three endpoints now pass. The final
400→800 K change shrinks to6.66791e-8/5.29097e-8. This is consistent with the
prospective long-interval discretisation explanation after analytic separation
of truncation and floating-point drift. It does not prove a unique source of
the residual or an absolute native64 error bound.

### Cost and verification

Measured producer **1880.814483 s /31.3469 min**, sampled CPU **8441.309248 s**,
recorded kernel high-water **436,060,160 bytes /0.406113 GiB**,9384 samples at0.2 s.
The wall time exceeds the1500-s expectation but stays below the2400-s upper
estimate and3600-s cap; staged initial two-bundle execution limits parallelism.
Eight physical cores4–11/8 GiB; all41 frozen sources and242 installed engine
artifacts verify before/after, with all24/12/2/1 input/reference hashes.
The38-file producer collection retains every original output and exit.

[Fresh public evidence](qualification-evidence/numerov-corners/README.md) passes
all saved decisions and **57 focused tests**, including separately pinned
verification support (the independent replay module and two additional tests).
The numerical subset is explicitly restricted to14 predeclared rows; the other
70 saved rows receive saved-byte and decision verification, not reexecution.

The public numerical subset finishes at **08:59:13.263552 UTC, exits0/0**.
All14 selected rows reproduce every stored K/S component exactly. Measured
**441.291700 s /2229.002709 sampled CPU-seconds /396,722,176-byte recorded
kernel high-water**,2200 samples at0.2 s, eight physical cores/8 GiB. The21-file
replay collection, owner and separate report/resources are published. All
**15 final public payloads** are freshly downloaded and checksum-verified after
the replay supplement. This remains14 numerical reexecutions, not all78 solves.

The [08:59:58 UTC result capture](progress-numerov-corner-result-20261019.json)
reverifies69 remote collected hashes,41 frozen sources,242 installed engine
artifacts and all24/12/2/1 input/reference pins. **No scientific container is
running**;61.51 GiB persistent space remains. The
[six-panel overview](../../docs/physics/figures/co-progress-20261019-corners.png)
shows the old coarser diagnostic beside the new finest-endpoint passes.

**Qualification scope:** the increment formulation is now qualified at these
14 identities and declared grids. These identities were already in historical
coverage; this adds **zero unique energy identities**. The original132/594
coverage and462 uncomputed identities remain, and untested interior/large-K
energies need a separately costed contract before full-grid release.
