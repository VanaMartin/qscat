# Immutable CO orbital-frame continuum preparation

## Scope and scientific boundary

Prepare the full-continuum native integral and **target-CSF** layout from the
completed fifty-component CAS(10,11)/cc-pVDZ orbital export at R=2.1323bohr.
This is a construction stage common to a future combined-prescription contrast
and a future matched-manifold design. It neither selects retained target roots
nor runs a scattering calculation. Keep the original forty-component frame as
the existing baseline and the two-A1 physical-subspace rejection unchanged.

The finite stage executes **one four-rank `scatci_integrals` and eight serial
target `congen` calls**. It does not execute SCF/CASSCF/CASCI, target SCATCI,
DENPROP, scattering CONGEN/SCATCI, propagation or extraction. Reuse the immutable
Molden/checkpoint/report, not `run.py --target-initial-checkpoint`, whose adapter
projects orbitals and optimizes them again.

All quantities remain atomic units. The same baseline radius18bohr, Gaussian
continuum through l=4, deletion1e-6, exact target13-orbital selection7/3/3/0 and
native double precision apply. The target-only producer's1e-7 deletion setting
did not construct continuum functions; it is not a continuum convergence result.
Use the authenticated **baseline full-continuum input deck unchanged** with the
new immutable Molden. Do not increase angular cutoff, alter continuum exponents,
retained roots or active/frozen orbital spaces.

## Why the layouts cannot be interchanged by padding

The pinned upstream `main.pl` calls `run_code(qchem)`, integral generation,
target CONGEN/SCATCI and DENPROP before scattering. In `ukrmollib.pm`,
`make_scatci_integrals_input` inserts Gaussian exponents only when the scattering
flag is enabled; simply rerunning a target-only driver will not create them.
The library's `only` execution filter still calls output readers and therefore
does not by itself create a fresh restart. Avoid entering that entire pipeline
for this bounded stage: run the explicit native executables on frozen cards.

The completed baseline's native counts are target7/3/3/0, continuum51/34/34/22,
total58/37/37/22. Target-only counts have zero continuum. Both CONGEN singlet-A1
CSF dimensions are15548, but their binary integer address records differ
(observed initial full-layout entries117/118 versus target-only15/16), and
signed records also differ. Only one record's size differs by48 bytes; **that
does not make the other bytes interchangeable**. Target CIDATA additionally
contains three centres in the continuum layout versus two in the target-only
layout. These are observed binary differences, not a decoded signed-CSF
transformation law. No vector conversion or metadata-only reheader is attempted.

Source trace is against the unchanged campaign `main.pl` SHA256
`c81a67fa5f33b990eb035e279ebb2efca0be29a09930e9c18ffd97e352f8727d`
and instrumented `ukrmollib.pm` SHA256
`f531586e68396c89e18623516fae8f1a5dc1135e60858255fe66a4da109f46e9`.
The relevant qualified source anchors are `run_code`, `run_program_this_time`,
`make_scatci_integrals_input` and `read_scatci_integrals_output`; the runner
adapter is `projects.ukrmol_co.run.enable_state_average` and its
`state_average.pm.read_state_average_output`. No reference-tree code is ported.

## Construction and finite checks

`immutable_continuum.prepare(root, output, contract)` first authenticates every
input, the passing independent new-frame coverage and the preserved physical
subspace rejection. Check checkpoint-to-Molden physical core/active overlaps
against the existing1−1e-8 floor. Copy only the exact Molden bytes; preserve
metadata/checkpoint provenance through input hashes. Native programs cannot
reach an orbital optimizer, and all original inputs mount read-only.

After integral generation, independently parse the unique final native
target/PCO/TGT+PCO/continuum/total inventories. Require eight nonnegative counts,
target7/3/3/0 (remaining zeros), no PCOs, positive continuum in all four C2v
irreps, equal B1/B2 counts and exact total=target+continuum. Native MPI error
counts must all be zero and integral completion/files present. Preserve all
warnings, including the earlier engine's floating-point status flags; do not
reinterpret a zero exit as an accuracy certificate.

Generate each target CONGEN card from its authenticated baseline card, changing
**only `nob` to the actual generated total orbital inventory**. This is native
address accounting, not a changed physical CAS. Keep `nob0=7,3,3,0`, all CSF
groups, spin/symmetry, electron count and workspaces unchanged. Require completion
and six well-framed records to complete EOF in each generated target CSF file.
The pinned64-bit native header must retain14 electrons, its sector multiplicity,
the generated total MO count and the existing dimensions:

| |A1|B1|B2|A2|
|---|---:|---:|---:|---:|
|Singlet|15548|15196|15196|15044|
|Triplet|24366|24540|24540|24564|

Record each binary record's length/hash without claiming to have decoded its
coupling convention. Persist every generated file, actual native commands/exits,
inventory, timing, resource record and before/after integrity checks. A construction
pass means only that these native files were generated with the declared layout
and immutable input frame. It is **not** target-vector reuse qualification,
target-Hamiltonian agreement or a scattering/model pass.

## Cost and next decision

The prior same-basis full-continuum integral stage measured11.020s at four ranks;
its eight target CONGEN calls total about3s. The target-only integral time0.114s
is not an appropriate forecast. Declare **60s expected /180s upper /300s hard
container cap**, four physical cores0–3, four MPI ranks/single-threaded native
libraries,4GiB memory and persistent `/home` outputs. Require8GiB available RAM
and30GiB disk. Individual integral/CONGEN caps are120s/30s. These are forecasts,
not a benchmark or a whole-scattering budget.

Freeze/push source, tests and machine contract before the native stage. Stop on
any failed integrity, inventory, native-error or CSF framing/dimension check;
no automatic retry, larger space, new solver or altered model follows.

After a pass, use these saved full-layout CSFs to design the **signed basis/
address correspondence** for reuse of saved native target eigenvectors, or
justify and cost a finite native target reimport if exact reuse cannot be
established. Do not repeat QC or the37-minute property calculation merely to
change a file layout. Select the scientific contrast and validate compatible
target/property/channel maps before any sparse2048-root scattering Hamiltonian.
The A1 subspace rejection cannot be repaired by this technical construction.

## Saved native completion-marker diagnosis

The first execution (`bc1a0c3`) stops after the first CONGEN because the wrapper
uses DENPROP's `has completed execution` trailer, which this pinned CONGEN does
not emit. Both actual native calls exit0. The original wrapper/container/owner
1/1/1 and partial coverage remain unchanged. Native CONGEN's actual markers are
`Data on CSFs has been written to file (MEGUL)` with the exact unit and CSF count,
and `Wavefn projection (projec()) - completed`; require these **and** complete
six-record binary framing/header checks. This corrects a program-specific reader
defect, not a numerical/physics tolerance. Regression cases reject missing,
duplicated and wrong-unit/dimension outputs.

The separate seven-sector completion contract reuses the saved integral data
and completed singlet-A1 file, after reconstructing their original native exits,
markers, counts and binary EOF. It authenticates17 parent inputs plus11 new pins,
copies saved inputs into a fresh directory and invokes only the seven unrun
CONGEN sectors. Forecast10s expected /30s upper /120s cap on the same reserved
four cores /4GiB. No integral or first-CON GEN calculation is repeated; no
scientific setting is changed. The original failed producer and successor have
separate results and original exits. Any new failure stops without another
automatic call.
