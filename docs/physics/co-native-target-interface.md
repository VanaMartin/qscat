# Native target read-only consumer and saved integral-block comparison

## Purpose and fixed scope

The same-frame signed CSF correspondence is already proved for all eight sectors
and74 roots. Resolve two remaining technical questions without repeating QC,
target diagonalization or properties: does the installed native phase-aware
consumer read the original two-centre CIDATA with the full-layout phases, and do
the saved target-only/full-continuum integral files encode the same target
Hamiltonian? This does not select retained roots or assemble any scattering
Hamiltonian. Physical old40/new50 A1 subspace rejection is unchanged.

Lifecycle: source-trace/design, a small native **I/O harness** with a Python
differential checker, analytic failure controls, then a finite measured native
read. The existing Python signed-map proof is the independent coefficient oracle.
No method is promoted, no engine source is ported/modified or recompiled, and no
Rust optimization is warranted. The Fortran driver links the installed engine
static libraries and their installed module interfaces; both enter provenance.

## Native source/API anchors

Source is the checksum-pinned UKRmol-in3.3.0 archive used in the signed-map note.
`ci-diag/scatci_routines.f:CIRMAT` calls `READCIP`, `CIVIO` and `CPAKPZ`; the last
multiplies saved source phases by requested destination phases without changing
root ordering. Call this library entry for all13/13/8/8/8/8/8/8 roots, not a
diagonalizer or the SCATCI executable. Supply phases from each authenticated full
target-CSF file; compare all returned coefficients/relative energies byte-exactly
against saved CIDATA and the already qualified signs. The harness output is an
I/O diagnostic array, **not a converted native CIDATA file**. A final explicit
completion marker and binary EOF reject a native `STOP` that returns exit0.

`gbtolib/source/molecular_basis_mod.f90:read_basis_and_integrals` reads the existing
atomic/molecular basis metadata, symmetric one/two-electron arrays and integral
options. Use `get_number_of_orbitals`, `get_continuum_flags` and
`get_absolute_index` to derive the target embedding, not a guessed global offset.
`integral_index` plus stored array access extracts target matrix elements, as in
the pinned `get_integrals_template.f90` and native `ukrmol_interface.f90:GET_INTEGRAL`.
The latter's scalar one-electron column is T+V, and its two-electron access uses
the symmetric electron-repulsion array. No integral evaluation is requested.

Require target7/3/3/0, no occupied continuum, total13 versus154 and exact expected
full inventories58/37/37/22. Extract all169 overlaps,169 one-electron Hamiltonian
elements and28,561 ordered spatial ERIs per file, including symmetry zeros.
Preserve geometry/charge metadata and derive nuclear repulsion only from it;
the extra continuum centre has zero charge. No ECP is part of this CO model.
All quantities are atomic units: integrals/repulsion in Hartree, geometry in bohr.

## Prospective numerical gates and bound

Both target-overlap matrices must agree with identity and each other within
1e-10 (rtol0). Physical charges/coordinates must agree within1e-12, and the
destination's only extra centre has charge0. Require finite values, h symmetry
and all ERI pair/exchange symmetries within1e-10. Every h/ERI difference must be
≤1e-10Hartree (rtol0). Do not infer a many-electron error from that maximum alone.

In the common orthonormal spin-orbital Fock space, each creation/annihilation
product has operator norm≤1. For spatial h and chemists' ERIs v,

    ||ΔH|| ≤ |ΔE_nuc| + 2 Σ_pq |Δh_pq| + 2 Σ_pqrs |Δv_pqrs| = B.

The factors account for two one-electron spin choices and four two-electron
spin choices with the Hamiltonian's1/2. Sum **all ordered entries**, not just
unique stored integrals. The same bound applies after restriction to the frozen-
core CAS and any orthonormal CSF subspace. Weyl's inequality bounds ordered
target eigenvalue changes by B for the exact finite matrices. Require B≤1e-7
Hartree, the existing native target import tolerance. This bounds the difference
of the **saved numerical integral Hamiltonians under the common-basis
identification**, not absolute integral accuracy, model convergence, eigenvector
stability near degeneracies or a scattering phase/width error. Measured target
overlap defects are also reported; finite overlap tolerance is not a proof of
exact physical orbital equality. Exact Molden provenance and prior frame checks
remain essential. No native target eigenproblem is solved in this check.
Each source CIDATA nuclear-energy reference must also match the source geometry's
repulsion within1e-12Hartree (rtol0), so a changed constant shift cannot hide behind
unchanged relative eigenvalue records.

## Execution and stopping contract

One fresh persistent Sadaharu output directory, original inputs mounted read-only,
source frozen/pushed before execution. Reuse the pinned selected-root engine
image through the existing debug-tools image only for its compiler/binutils; pin
image IDs and confirm library bytes against the original runtime. Native reads
run at one MPI rank/thread on the small-job reservation of four physical cores,
4GiB memory; available RAM/disk floors8GiB/30GiB. Forecast60s expected /180s upper,
300s container hard cap; compilation120s, each of three reader processes30s.
These are bounded forecasts, not measurements. Preserve compiler command/log,
executable and module/library hashes, all native commands/exits, raw binary/text
outputs and sampled resources. Stop on failure with no automatic retry or model
change. This is one compile, one eight-sector target read and two integral-file
reads; no new integral generation, QC/CI, property export or scattering call.
The integral API may open files with read/write access even when only reading:
consume exact private read-copies in the fresh output, checksum them before/after,
and leave original mounted inputs read-only. Do not republish these unchanged
copies; keep exact prior-capsule pointers. Pass the library's declared fixed-length
path character argument, not a shortened temporary string.

Even a complete pass leaves **actual scattering CONGEN target-channel phases,
root/property/channel maps and the scientific contrast** unset. Retain original
CIDATA and40-state properties. Declare combined prescription sensitivity versus
a separately justified matched-manifold selection before any sparse2048-root
scattering calculation; no dense Hamiltonian or automatic root/channel ladder.
