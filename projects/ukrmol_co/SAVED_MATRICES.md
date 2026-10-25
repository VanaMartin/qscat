# Saved matrices support the next physical check — 17 October 2026

**The nine-A2 outputs contain enough information to reconstruct all saved
real-axis S samples at all three anchors.** Six cases pass (both Pi sectors at
each anchor), reproducing **594 six-channel matrices** and native eigenphase
sums to **4.94e-8 rad modulo π**. The full twelve-case gate nevertheless
**rejects**, with the six five-A2 interface readings unset. Original reader/
container exits **1/1** and every case are retained.

This resolves the immediate data-availability question. Saved poles, boundary
amplitudes, thresholds, angular labels and long-range multipole coefficients
are available. No new inner-region solve or export is justified for the next
bounded physical benchmark. Independent outer propagation and native complex
poles remain unqualified.

## What was frozen and measured

The [design](../../docs/physics/co-saved-scattering-matrices.md),
[contract](saved-matrices-contract.json), reader and three analytic/format tests
were frozen at **`343cc1b` before the full reading**. The matrix of twelve cases
is compressed/equilibrium/stretched × five/nine-A2 × B1/B2. Every native K energy
record is parsed; actual EOF/completeness, rather than printed record counters,
is checked. No phase fitting, model selection or complex continuation is done.

For the nine-A2 cases:

| Check | Maximum over six cases | Frozen budget |
|---|---:|---:|
| Native eigenphase sum difference modulo π | **4.939e-8 rad** | 1e-6 rad |
| S†S − I, entrywise | **1.511e-14** | 1e-10 |
| S − Sᵀ, entrywise | **1.508e-14** | 1e-10 |
| Cayley solve versus spectral S | **9.271e-14** | 1e-10 |
| Determinant phase versus summed eigenphases | **6.594e-14 rad** | 1e-10 rad |

The spectral route diagonalizes K, while the Cayley route solves
`(I − iK) S = I + iK`. Both use the same native K, so this is independent
format/algebra verification, **not independent native outer propagation**.
Only one K triangle is stored: reciprocity after unpacking cannot measure the
unsymmetrized native matrix's reciprocity. The largest Cayley-system condition
number is approximately **313**. Known rotated-channel eigenphase tests also
pass; missing energy records are rejected by a separate format test.

Each passing boundary file contains **2048 exported poles**, **298 channels**,
radius **18 bohr**, and two packed **44,551-entry** multipole coefficient
matrices. Both orders have nonzero finite entries. The reader checks thresholds
against the native target energies. The six open channels attach to the ground
target and have `l = 1,2,3,3,4,4`, with the sector's real-harmonic angular labels;
the other 292 channels are closed on this grid. This is not a scalar p wave.

The first excited-target thresholds are approximately **7.49562 / 6.38502 /
5.02524 requested-eV** at compressed/equilibrium/stretched. All 99 saved points
lie at 0.1–5 requested-eV, below those thresholds. Native energy conversion
remains `0.0735 Ryd/requested-eV`; the small physical/requested-eV distinction
is preserved. Boundary pole energies are total Hartree energies. The binary
geometry field is zero for these polyatomic exports; anchor identity is carried
by the original hashed run/deck provenance, not inferred from that field.

## Five-A2 target labels: preserve the failed gate

The five-A2 channel files declare **40 selected target rows**, but channel
target labels retain the original 48-state property IDs. The native writer
densifies the selected target table independently of the channel labels. The
frozen strict reader therefore reports `ValueError: Invalid channel label` in
all six five-A2 cases, leaving their verdicts unset and the aggregate false.

| Five-A2 anchor | Maximum channel target ID | Channels above 40, per sector | Channels needing dense-row remapping, per sector |
|---|---:|---:|---:|
| Compressed | 42 | 6 | 12 |
| Equilibrium | 43 | 12 | 18 |
| Stretched | 43 | 12 | 12 |

A **separate source/map diagnosis** verifies that the actual `IDTARG` decks
exactly match the saved pilot maps. Looking up each original property ID in the
sorted selected-ID list reproduces **every threshold exactly** in all twelve
files. Using an in-range label directly as a dense row instead produces errors
up to **0.015165426 Ryd** compressed and **0.060025268 Ryd** equilibrium; labels
above 40 are invalid dense-row addresses. The reader has not silently remapped
these inputs, and the diagnostic does not replace its failed gate.

The checksum-published external source provides exact anchors:

- `Postprocessing_module.F90::get_channel_info`, lines 748–778: skip unselected
  states but retain original property index `i` in `ICHL` and compute the
  threshold from that state's energy.
- `::write_channel_info`, lines 1070–1083: write selected target rows with dense
  counter `k`, then write unchanged `ICHL` in channel rows.
- `outerio.f::READTC`, lines 1231–1239 and 1295–1305: document/read target pointers,
  target energies and channel thresholds separately.

This identifies a numbering mismatch relevant to independent readers and
index-based consumers. It is **not** evidence of an engine correction, a
quantified error in the historical native phases, or a reason to replace old
channel verdicts. Native-consumer impact remains to be traced before using the
five-A2 interface for a new calculation. Any future mapped decoder needs its own
explicit provenance and tests, while retaining this reading's original exit.

One other diagnostic is misleading: native `rsolve.out` prints zero calculated
K energies although the files contain 99 complete records. The active
`rsolve_compak.f` writes K at lines 741–755, but increments `NRK` inside the
optional `CALCAK` block at lines 761–785. The binary data, not that counter,
establish completeness. No engine patch is made for either observation.

## Physical design requirement

The actual native path propagates coupled multipole potentials from **18 to
100 bohr** and uses Gailitis asymptotic matching. The two exported coefficient
orders correspond to `r⁻²` and `r⁻³` terms in the Rydberg radial operator, in
addition to the centrifugal term. Do not replace them with free spherical waves
or assert a complete physical polarization model from the presence of these
two orders. Closed-channel coupling and the selected target representation
remain consequential. The 2048-pole sum is **uncorrected and truncated**; its
passing real-axis refinement comparison is not an absolute complex-pole bound.

Next freeze a solvable physical scattering benchmark with independent
real-axis scattering and outgoing complex poles, then a coupled polar-tail
benchmark. Validate boundary normalization/Bloch conventions, energy zero,
threshold sheets and outgoing conditions before implementing an independent
outer propagator. Require that propagator to reproduce the saved **full
open-open K matrices**, not just their eigenphase sum, with radius/step/tail
controls. Only then can a separately declared native complex-pole study address
the width ambiguity. A rational interpolation of these real-axis samples would
not by itself provide independent pole evidence.

## Evidence and cost

[Quoted verification](verification-saved-matrices-20261017.json) preserves the
six passing cases, six unset cases, exact failures, original exits, all map
diagnostics and resource provenance. [Ten public payloads](qualification-evidence/saved-matrices/README.md)
reconstruct **123 raw files / six reader-source files / thirteen writer-source
and deck files**. Public replay preserves the false aggregate and exit 1,
matches every case verdict/failure, and independently reruns the published
label-map diagnostic. Maximum host/public matrix-metric drift is **3.88e-14**.

Sadaharu reader cost: **0.664 s wall**, **0.043 GiB sampled kernel high-water**,
reserved physical cores 0–3, one BLAS thread, 4-GiB container limit. The sampler
runs every 0.2 s; its recorded high-water value is not a post-exit peak. Raw
collection/checksumming and publication are separate from this reader timing.
Persistent sources/outputs/logs are in `prepared/saved-matrices-source-20261017/`.
**37 focused tests pass.** No native scattering solve runs. Earlier extraction,
channel and active-space rejections remain intact; production remains gated.

## 25 October scope update: exact mapped-view prerequisite, not a repaired gate

The [compressed matched-pole proposal](../../docs/physics/co-compressed-a2-pole-contrast.md)
now traces original property IDs through the authenticated **zero-Buttle,
fixed-nuclei K path**. Its active matching uses stored channel thresholds; the
relevant target-pointer dereferences are in branches disabled by these inputs.
This is not evidence that all consumers or correction paths tolerate the IDs.
No native source/input is edited and the original false aggregate/six unset
five-case readings above remain unchanged.

The exact next prerequisite is one explicit250-channel compressed-five mapped
view and narrow pole-runner adapter, with selection/threshold/identity negatives
and development full-S gates before complex use. Existing eigenpairs are present;
no QC/property/inner solve or native rebuild is needed merely to supply the view.
Metadata/source inspection is not decoder or pole qualification. The proposal
does not launch a new calculation, consume equilibrium poles/114 route holdout S
values, or establish electronic accuracy or a qualified width curve.
