# Native target read-only interface and saved integral comparison

## Completed — 21 October 2026,19:55 UTC

**Native target reuse and the saved target-Hamiltonian comparison pass.** The
original phase-aware CIDATA is reusable unchanged: all74 native vectors and
relative energies match byte-exactly. Saved target-only/full-continuum one-
electron and overlap blocks are exact; the largest ERI difference is
**4.44089e-15Hartree** and the conservative common-basis operator-difference
bound is **6.26708e-13Hartree**, below the frozen1e-7 gate. This removes the
technical reason to repeat QC, target diagonalization or properties.

The ABI completion, freeze **`1a785ee`**, ends19:40:49UTC, command/container/owner
**0/0/0**, **1.066s /.09886GiB**, sampled CPU.588181s. It reuses the completed
native target read and executes only the two missing successful integral reads,
native exits0/0,0.01556s /0.03154s. Compilation0.11395s. All22 source-file /
30 input pins pass before/after, all46 recovered module headers verify, and five
linked engine libraries match the original runtime. No completed native read is
repeated. Original inputs and exact private integral read-copies remain unchanged.

Across the four finite job attempts: measured container wall **3.761s**, sampled
CPU **1.589s**, peak **.12254GiB**. Four physical cores0–3 reserved, one native
rank/thread,4GiB cap. The original compiler-wrapper and link failures plus the
failed ABI read remain separate; the metadata-only header inspection/extraction
is not a native scientific calculation or part of that job timing total.

| Compared target quantity | Entries per file | Maximum difference / Hartree |
|---|---:|---:|
|One-electron h|169|0|
|Ordered ERIs|28,561|4.44089e-15|
|Nuclear repulsion|1|0|

Both169-entry target overlaps have identity defect **6.66134e-16** and exact
interfile agreement (dimensionless). ERI L1 difference **3.13354e-13Hartree**,
h L1 difference0. All required h/ERI symmetries, physical charge/geometry and
CIDATA nuclear-reference checks pass. The finite operator bound compares the
saved numerical Hamiltonians under the common orthonormal-basis identification;
it is **not absolute integral accuracy, physical eigenvector identity across
frames or a scattering phase/width error bound**.

[Fresh public reconstruction](verification-native-target-interface-public-20261021.json)
verifies **7 payloads /3 prerequisite archives /115 generated-owner files /
73 source checks /30 input pins /46 matching double headers**, and reconstructs
all74 native roots exactly. Its independent stdlib reader uses record framing,
typed fields, IEEE sign-bit flips, scalar symmetry checks and `math.fsum` for the
operator bound, importing no producer numerical reader. No native replay occurs.
Six unchanged copied files are referenced rather than republished. Original
failed exits remain1/1/1 for each of three predecessors; completion is0/0/0.

## Preserved partial outcome and failures

The installed original `CIRMAT → READCIP/CIVIO → CPAKPZ` consumer completes all
eight sectors /74 target roots in0.03165s, native exit0. Saved-only reconstruction
checks every coefficient and relative-energy byte against the qualified sign map
exactly, to complete EOF. The source CIDATA stays unchanged; no target
diagonalization, QC, properties or scattering Hamiltonian is run.

Three original failed owners remain separate:

| Freeze | Outcome | Original command/container/owner |
|---|---|---|
|eb092a3|Open MPI points to absent relocated compiler; no compile/read|1/1/1|
|87860b8|Compiler relocation corrected; original utility link dependency missing; no read|1/1/1|
|c1c3026|Link succeeds; target read completes; source integral reader SIGSEGVs, full read unrun|1/1/1|

The last native stage exits are compile0 /target0 /source-integrals−11. Its native
target output independently verifies despite the later failure; the original
producer verdict remains unset, not retroactively changed. Runtime source inputs
remain unchanged and the failed read-copies/logs/executable are retained.

## Diagnosed precision/header mismatch

The inherited installed GBTO modules declare **cfp16 /REAL16 arrays**, while
the source-built `lib.double/libGBTO.a` is **cfp8**. Library calls through the
correct integer/double `CIRMAT` interface succeed; direct GBTO object access with
quad headers does not have a valid ABI. The runtime Docker stage had copied
source-built libraries without their corresponding module headers.
This affects **new clients compiled against inherited headers**, not the earlier
source-built executables, which were built with their own matching modules.
It does not invalidate previously qualified target/scattering outputs or imply
that those executables ran in mixed precision.

Retained source-inspection image
`sha256:6f080f2f678b7048bb67403a6013d5f9a2cfc04305f7080eae3d7b8511816f85`
contains46 matching double headers emitted alongside **the identical GBTO library**,
SHA256 `77c84296264a2a9b3bed16f958acd8ce04db4636cce475ee762f1c9b2d7c686a`,
with the same source archive/compiler provenance. Metadata-only extraction pins
all46 modules; no engine is rebuilt or replaced. The corrected Dockerfile copies
matching headers for future builds, but no new image build/test is claimed here.

The [integral-only ABI completion contract](native-target-interface-abi-contract-20261021.json)
inherits all original scientific gates, five engine libraries,30 input pins and
the60s expected /180s upper /300s cap, four-core/4GiB budget. It reuses completed
target.bin/target.out and **does not repeat the native target read**. It compiles
the I/O harness with matching double headers plus a precision guard, then admits
only the two missing successful integral reads. A new failure stops without
changing numerical/model gates.

## Scope and next production decision

All169 overlap entries,169 h entries and28,561 ordered ERIs per saved integral
file are compared, including the geometry-derived nuclear constant.
The conservative common-orthonormal-basis bound
`|ΔE_nuc| + 2Σ|Δh| + 2Σ|ΔERI|` must be≤1e-7Hartree. This compares the two
saved numerical target Hamiltonians, not absolute integral accuracy, physical
orbital equality to exact arithmetic or a scattering-error bound.

The actual scattering CONGEN target-channel phases and compatible root/property
maps remain unset. Both A1 physical old40/new50 retained subspaces still reject;
declare combined-prescription sensitivity versus a separately justified matched-
manifold selection before scattering. No dense Hamiltonian, root/channel ladder,
production sweep or potential fitting is released by native interface checks.

[Mathematical/native API design](../../docs/physics/co-native-target-interface.md) ·
[Public evidence and reproduction](qualification-evidence/native-target-interface/README.md) ·
[Dated commands/resources/results](progress-native-target-interface-result-20261021.json)
