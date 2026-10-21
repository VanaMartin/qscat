# Immutable-frame continuum and target-CSF preparation

## Completed — 21 October 2026,17:40 UTC

**Full-continuum integrals and all eight target-CSF files are generated on the
qualified new fifty-component orbital frame, without reoptimization.** Target
counts7/3/3/0 and continuum counts51/34/34/22 give **154 total orbitals**, including
141 continuum orbitals. There are no PCOs. All eight CSF dimensions, spin/symmetry
headers, final native write/projection markers and six-record binary EOF checks
pass. The unchanged original forty-component frame remains the existing
scattering baseline.

The stage uses the exact baseline18-bohr /l≤4 /deletion1e-6 continuum input deck
and exact new-frame Molden bytes, not a projected/reoptimized restart. Checkpoint/
Molden core/active physical-overlap singular values retain minima
0.9999999999999499 /0.99999999999998. **No QC/CASCI, target diagonalization,
DENPROP, scattering CONGEN/SCATCI or propagation is performed.** This is a
generation milestone, not target-vector reuse qualification or Hamiltonian accuracy.

## Original wrapper failure and separately frozen completion

Original source **`bc1a0c3`** runs17:31:30–17:31:42UTC. The four-rank integral
call and singlet-A1 CONGEN return **native0/0**. The wrapper then wrongly looks
for DENPROP's `has completed execution` trailer, which CONGEN does not emit.
Its original **wrapper/container/owner1/1/1**, null generation verdict, traceback,
partial files and source freeze are preserved unchanged. Integral runtime10.630s,
first CONGEN0.264s; original whole-container11.530s, CPU42.620s, sampled kernel
peak **2.46218GiB**,56 samples. All13 frozen sources /17 input pins verify.

Saved native output instead contains the final CSF unit/count and
`Wavefn projection (projec()) - completed`; the saved singlet-A1 binary has
15548 CSFs and complete EOF. It is also byte-identical to the old full-continuum
CSF definition. This establishes a program-specific reader defect, not a
numerical or physical rejection. The corrected marker reader has missing,
duplicate, wrong-unit and wrong-count regression controls.

The separate [seven-sector completion contract](immutable-continuum-completion-contract-20261021.json),
freeze **`14c4c61`**, authenticates17 original inputs plus11 saved-failure/input
pins. It copies the saved integral file and first completed CSF into a fresh
directory and invokes **only the seven previously unrun CONGEN sectors**.
It finishes **17:36:12UTC**, original command/container/owner **0/0/0**, in
**3.519s**, CPU2.857s, kernel peak **0.28282GiB**,16 samples. All15 successor
sources /28 input pins verify before/after. Integral generation and the first
CONGEN call are **not repeated**. The failed original is not retroactively
declared a pass.

Across both attempts: **one integral call and eight target CONGEN calls**, all
nine native exits0; total measured container wall **15.049s**, sampled CPU45.478s.
Four reserved physical cores0–3,4GiB. The original60s/180s/300s forecast and
successor10s/30s/120s forecast are not exceeded. No scientific job remains.

## Native layout and independent saved reading

| Sector | Target CSF dimension | Full-layout CSF byte observation |
|---|---:|---|
|Singlet A1|15548|Exact baseline hash|
|Triplet A1|24366|Exact baseline hash|
|Singlet B1|15196|Exact baseline hash|
|Triplet B1|24540|Exact baseline hash|
|Singlet B2|15196|Exact baseline hash|
|Triplet B2|24540|Exact baseline hash|
|Singlet A2|15044|Exact baseline hash|
|Triplet A2|24564|Exact baseline hash|

All eight generated CSF files are byte-identical to the baseline full-continuum
definitions in a **read-only host comparison**. This is an observed basis-definition
agreement, not a comparison of many-electron eigenstates between orbital frames.
The public reader authenticates the new bytes against the recorded baseline
hashes; it does not regenerate or independently download those baseline CSFs.
The new integral file is **12,983,436 bytes**, SHA256
`577873767fd9ef7df58c809a460aea444c8a5fdfcc1307d8149276e4b8db1c89`.
It differs from the old orbital frame's integrals, as expected. New/old identical
CSF definitions do not imply identical electronic Hamiltonians or states.

[Fresh public reading](verification-immutable-continuum-public-20261021.json)
checks all **6 payloads /57 generated-owner files /28 source-file checks across
the two freezes**,17 original inputs and11 successor pins. Eight unchanged copied
files are reconstructed by exact pointers instead of republished. A separate
SciPy Fortran-record reader reconstructs all **48 binary records**, headers and
complete EOF, native inventories, final CSF text markers/counts, all nine native
exits and original failed-wrapper outcome. It also verifies exact native deck and
Molden bytes. This saved-only reading takes1.596s and runs no native calculation.
Integral bytes remain opaque authenticated output, not independently recomputed
or checked element-by-element. Native error counts are zero; all warnings and
floating-point-status notes remain in raw logs, not discarded as irrelevant.

[Public evidence/reproduction](qualification-evidence/immutable-continuum/README.md),
[dated result/resources](progress-immutable-continuum-result-20261021.json),
[prospective design](../../docs/physics/co-immutable-continuum-preparation.md) and
[original contract](immutable-continuum-contract-20261021.json) retain the exact
scope. Native/program-specific readers have34 focused local tests passing,
plus one local PySCF analytic skip exercised in the earlier remote control.

## Next production step

**18:20UTC follow-up:** the [exact signed target-CSF correspondence](NATIVE_CSF_MAP.md)
now passes every term/all8 sectors/all74 native vectors. Only16 singlet-A1 CSFs
flip sign. Native source explains phase-aware reuse with the original two-centre
CIDATA unchanged, so no converted/reheadered file is needed. Actual native
consumer, integral subblock and scattering-channel phase qualification remain.
The original generation result below and physical A1 rejection are unchanged.

The target-only and full-continuum layouts are still not interchangeable by
padding: orbital addresses and signed CSF records differ, even where dimensions
match. These new files allow an exact **signed basis/address correspondence** to
be designed and checked against the saved target-only native vectors, without
another QC optimization or37-minute property calculation. Decode and prove that
map before converting coefficients; otherwise cost a finite native target
reimport explicitly rather than silently assuming reuse.

The physical old40/new50 retained-five screen **still rejects both A1 sectors**.
No retained map is selected by this target-CSF construction. Declare combined
orbital/channel prescription sensitivity versus a justified matched-manifold
selection before scattering, and verify compatible native target/property/channel
maps. No new target eigensolver, properties, sparse/dense scattering Hamiltonian,
phase/width claim or production sweep is automatically released.
