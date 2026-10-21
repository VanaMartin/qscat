# Exact target-CSF address and phase correspondence

## Purpose and boundary

Reuse the already qualified 74-root new-frame CIDATA without reoptimizing CO,
diagonalizing a target Hamiltonian or repeating properties. This is a **same
orbital-frame representation problem**, not the rejected old40/new50 physical
subspace comparison. Units remain Hartree/bohr. No scattering solve or retained
root choice follows automatically.

Lifecycle: design the finite correspondence, implement a readable saved-file
Python checker under `projects/ukrmol_co`, validate analytic determinant controls
and all saved sectors, and stop. No Rust optimization or library promotion is
warranted. This note is the prospective contract for the saved-only check, not
an announcement that it has passed.

## Authenticated native source

UKRmol-in 3.3.0 archive SHA256
`1b81e423b7e2a5a530c336cfc7f7ab50d6f2705abf771c71d6119fc065d9cac2`
is downloaded from the URL pinned in `docker/ukrmol-plus/source-build.sh`.
The running image's `/opt/ukrmolp/source-provenance.json` identifies the same
archive, double precision and eight-byte default integers. Source anchors below
are relative to that archive's `source/`, not to an unpinned online branch:

* `congen/source/congen_driver.f90:getref`: spin orbitals are consecutive within
  each spatial orbital, alpha then beta, and spatial orbitals run by symmetry.
* `congen/source/congen_projection.f90:wrnfto`: six electron-only sequential
  records are header, target phases, orbital/reference/count metadata,
  coefficient/packed-determinant offsets, packed substitutions, determinant
  coefficients. The header's `idiag` is not multiplicity; S and Sz are reals.
* `congen/source/congen_projection.f90:dophz0,iphase`: canonical CSF phase is the
  first determinant's sorting parity times its coefficient sign. Substitutions
  replace reference electron slots before evaluating that parity.
* `mpi-ci-diag/Modules/Hamiltonian/CSF_module.F90:create_csfs` (the determinant-creation
  body): substitutions replace reference slots, then sorting parity multiplies
  the stored determinant coefficient.
* `mpi-ci-diag/Modules/Hamiltonian/SolutionHandler_module.F90:write_header_sol,
  write_eigenvalues,write_eigenvector`: CIDATA stores phases as **integers**, then
  energies and diagonals as reals, and unmodified CSF coefficients per root.
* `mpi-ci-diag/Modules/Hamiltonian/Target_RMat_CI_module.f90:read_ci_mat` calls
  `ci-diag/scatci_routines.f:CIRMAT`, which calls `READCIP`, `CIVIO`, `CPAKPZ`.
  `CPAKPZ` multiplies each coefficient by destination and saved source phases.
  `READCIP` reads its own CIDATA centre count and that many centre records; it
  does not require a target-only file to have the continuum file's third centre.

This source trace explains why a source target file with two physical centres
does not by itself need a third zero-charge centre inserted. It is not a live
test of the eventual scattering driver or property reader. Those consumers and
the scattering CONGEN target-channel phase map still require qualification.

## Mathematics

Let the occupied target spin-orbital address with irrep g, ordinal k and spin
s∈{0,1} be

    a(g,k,s;n) = 2 (sum_{h<g} n_h + k - 1) + s + 1.

The source inventory is (7,3,3), the destination (58,37,37,22); the **used target
space** is (7,3,3,0) in both. Embed only these target addresses by retaining
(g,k,s), rejecting any continuum occupation. This map is strictly increasing
on the used target addresses, so it introduces no permutation sign of its own.

For CSF j, its packed terms specify ordered determinants D_jt by substitutions
into reference electron slots. Require identical CSF and term order, identical
offsets/counts, and every reference/packed address to match under the embedding.
Let p_j be its saved canonical phase and c_jt its stored determinant coefficients.
With q_j = p_j(destination) p_j(source), require for **every term**, not just the
first,

    c_jt(destination) = q_j c_jt(source).

Then Φ_j(destination)=q_j Φ_j(source) after embedding. Both references and
substituted electron-slot order are checked; the determinant sorting parities
are equal because the embedding is increasing. Reconstruct p_j independently
from the first substituted determinant to reject stale/forged phase metadata.
Thus for every root r the destination CSF coefficients are

    v_jr(destination) = q_j v_jr(source).

Q=diag(q) is real orthogonal and Q²=I: norms, Gram matrices, root ordering,
arbitrary overall root phase and energy values are preserved. No normalization,
root relabelling, overlap fitting or changed tolerance is permitted. This is
also the `CPAKPZ` source/destination phase product. Prove this finite exact
signed-diagonal relation; reject any need for a general CSF permutation or
rotation rather than inventing a fallback.

## Finite saved-only verification contract

Authenticate the source archive/provenance, eight target-only and eight completed
full-layout CSFs, and the original `fort.26` hash
`f81547d9776f6b5f23fd8921d4ba0202ddbb912d069d313e23b629a01527ec59`.
All inputs remain unchanged. Use complete record framing/EOF and explicit native
field types. Require finite coefficients, unit CSF norms to1e-12, target-only
occupations, valid substitution references, positive consistent one-based
offsets and no missing/trailing packed or coefficient data. Compare sign-adjusted
coefficient **bytes** exactly: multiplying by ±1 is an exact sign operation, not
a floating arithmetic approximation. Independent analytic controls must cover
fermionic slot-reordering signs, negative CSF orientation, extra empty symmetry,
bad phases, one changed nonleading term, continuum occupation and truncated data.

Read all eight CIDATA blocks to complete EOF. Check native source phases match
their source CSFs, root IDs/order13/13/8/8/8/8/8/8, dimensions, two physical
centres and finite energies/diagonals. In memory only apply Q and Q again;
require exact round-trip coefficient bytes, unchanged energy bytes and Gram
error≤1e-9. Save signs/counts/checksums and original input hashes, **not a new
native vector file**. A second saved-only determinant decoder should verify
the full expansions using canonical occupation keys and sorting parity.

Bound the check to the existing16 CSFs /74 vectors, one CPU,2GiB expected memory,
180s upper estimate (not yet a measured benchmark),300s hard cap. No native
engine, Hamiltonian, eigensolver or property call is admitted. Preserve any
failure in a fresh dated result; do not silently widen the transformation class.

A pass settles the signed target-CSF representation issue only. Reuse original
CIDATA through the native phase-aware consumer rather than padding/reheadering
it. Target integral-block equivalence, full scattering channel/phase and property
compatibility remain explicit prerequisites; the physical A1 rejection remains.
