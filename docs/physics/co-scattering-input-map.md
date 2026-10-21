# Saved CO scattering phase and property map

## Scope

Read the existing C2v doublet B1/B2 CONGEN files and cards, the original
74-root CIDATA and the already exported lowest-five-per-sector properties.
No CONGEN, integral, QC, diagonalizer, property or scattering executable runs.
The orbital inventory, target CSF definitions and grouping cards are unchanged
between frames; CONGEN constructs spin/angular algebra, not orbital coefficients.
Existing scattering definitions are candidates for reuse, not new physical
eigenstate identity evidence. Both A1 physical dominance rejections remain.

This project-only Python I/O audit is designed, tested against small fermionic
permutations and malformed records, then compared with native saved metadata.
It is not a method promoted into qscat; neither Rust nor an engine rebuild is
needed.

## Pinned native source convention

UKRmol-in3.3.0, archive SHA256
`1b81e423b7e2a5a530c336cfc7f7ab50d6f2705abf771c71d6119fc065d9cac2`:

- `source/congen/source/congen_driver.f90`: `mdegen` initializes to zero and
  changes only for `symtyp<=1`. C2v uses `symtyp=2`, hence the ordinary `dophz`
  branch applies, not the degenerate-target special case.
- `congen_projection.f90:dophz`: for each target CSF, take the first determinant
  of the first prototype continuum CSF, reconstruct substitutions in electron
  slots, and multiply its sorting parity, coefficient sign, and
  `(-1)**count(occupied > marked_continuum_spin_orbital)`. This last factor places
  the continuum electron after the ordered target electrons. Step through CSFs
  by the **prototype** `notgt`, not the expanded physical continuum count.
- `congen_projection.f90:wrnfto`: six sequential electron-only records. The
  second contains target phases followed by `nctgt,notgt,mcont,gucont`. The
  remaining records hold inventory/reference/counts, offsets, packed
  substitutions and coefficients, as in the target-only decoder.
- `ci-diag/scatci_routines.f:CIRMAT/CPAKPZ`: read source CIDATA phases and apply
  source times destination phase to coefficients. Different destination phases
  are allowed; no converted/reheadered CIDATA is written.

The audit derives the marked continuum orbital from the unique non-target
orbital in each first determinant. Its irrep and first-prototype ordinal must
match `mcont` and the saved cards. All target addresses, 15 distinct occupied
spin orbitals and total spin projection are checked. All recorded target phases
must reproduce exactly. This reconstructs **phase-defining first terms**, not
every N+1 determinant or a full eigenstate correspondence.

## Root and property contract

Target sectors are interleaved singlet/triplet A1,B1,B2,A2. CIDATA stores
13/13/8/8/8/8/8/8 roots; SCATCI selects roots1–5 in each sector by `ntgtf/ntgts`.
The property file is energy-ordered and its IDs are not sector ordinals.
Reconstruct its five-per-sector root map against CIDATA absolute energies
(`relative+nuclear reference`) within1e-7 Hartree, rtol0, preserving all moment
bytes. Require forty unique IDs, sorted finite energies, valid spin/irrep labels,
676 unique finite new-frame dipole/quadrupole records (679 in the saved old-frame
file) and complete state/rank coverage. Different nonzero-record counts are not
themselves evidence of incomplete selection; neither file is silently padded.
Check the old saved `idtarg` against the old property file. Generate the new
`idtarg` from the new property file, never copy the old energy-order permutation.

The prototype continuum multiplicity in CONGEN is two per sector. SCATCI's
expanded continuum counts come from the full inventory minus target inventory,
with continuum irrep `total_irrep XOR target_irrep`; require consistent cards
and `sum(5*expanded_continuum_count)=1410`. An equality of counts does not
prove channel convergence. This audit checks moment selection/indices/energies,
not independent numerical accuracy of transition moments or the downstream
SWINTERF application of every moment.

## Execution and release boundary

Freeze source/tests and hash-pin both original large CONGEN files, all cards,
CIDATA, properties and eight qualified full target-CSF files before the finite
saved-only audit. Four reserved physical cores, one process,4GiB cap,60s expected,
180s upper,300s hard timeout. Hash all inputs before/after; mmap the large files
read-only and verify complete six-record framing/EOF. Keep raw compact witnesses
(each phase-defining determinant and coefficient sign), exact small native
inputs, file hashes, command, exits and measured time/memory in a fresh persistent
directory. Compact witnesses are **derived slices**, not a substitute for the
full original-file hash or a claim to republish1.4GB of raw CSFs. No native call.

A pass supports saved-definition and map preparation only. Actual CIRMAT reads
with the scattering destination phases, downstream property-consumer semantics
and a frozen scientific contrast remain prerequisites for scattering execution.
Any later native check needs its own finite contract and matching double headers.
The lowest-five contrast must be labelled **combined orbital/channel-prescription
sensitivity**; matched-manifold selection requires separate physical evidence.
No dense Hamiltonian, automatic root/channel ladder, sweep or fitting is admitted.

## Preserved parser failure

The initial saved-only audit (freeze `1899892`) reconstructs the B1 first-term
witnesses but stops during card validation: its integer-card reader includes the
digit in the following `nob0` keyword while reading `nob`. Original exits and
partial witness remain in `prepared/scattering-input-map-source-20261021/`.
The property-map stage has already derived25 changed IDs and a maximum
4.91639e-10Hartree discrepancy to CIDATA. These partial observations do not turn
the original failed audit into a pass. A separately frozen parser-only successor
limits the captured assignment to integer tokens/commas/whitespace, with a
`nob` followed by `nob0` regression. Input pins and all physical/numerical gates
remain unchanged; no native calculation is retried.
