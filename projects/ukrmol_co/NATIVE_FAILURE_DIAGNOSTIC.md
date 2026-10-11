# CAS12 native matrix-construction diagnostic

The preserved 4096-root comparator fails before eigenpairs are exported. Its
[failure capture](progress-failure-20261011T2250.json) records B1 SCATCI's SEGV11,
MPI rank 4, MPI_ABORT59 and the original container/controller exits. The root
cause is unknown. The [separate diagnostic contract](cas12-native-debug-contract.json)
declares one instrumented B1 invocation using its saved native inputs.

## Why this diagnostic

Read-only inspection of all non-master rank logs places the failure after the
Class 4 completion messages. In the checksum-pinned UKRmol-in 3.3.0 source,
`source/mpi-ci-diag/Modules/Hamiltonian/Contracted_Hamiltonian_module.f90`
constructs Classes 5/6/7 in the following block. `mpi_scatci.f90` builds the
Hamiltonian before entering diagonalization. These facts narrow the investigation;
they do not identify an individual operation or establish root-count causality.

The original and successful B1 `fort.821` files are identical 712-byte files:
one 704-byte Fortran record, SHA256
`68a751002c94bf5f1aeea4053cd5217ee41c0c6eb82138154569b7147bf4981c`.
They contain a header, not an assembled Hamiltonian checkpoint. Saved target
CI/properties, integrals and CONGEN data can be reused; matrix construction
must restart from the beginning.

## Reproducible implementation

- `build_debug_native.sh` checks both retained upstream archives, builds over
  the original numerical libraries with `-O1 -g -fbacktrace -fcheck=bounds`,
  and records installed artifact hashes and debug/source-line sections.
- `native_debug_patch.py` requires the exact original contracted-Hamiltonian
  source hash and adds six flushed Class 5/6/7 entry/exit markers containing
  rank, target-symmetry pair and state counts. Removing these markers restores
  the original source bytes; the patch changes no computational statement.
- `debug_engine_checks.py` requires all twelve original serial/MPI target-HF
  reference checks and debugger smokes for normal exit, first-fault SIGSEGV
  with source stack, and Fortran bounds failure with source backtrace.
- `debug_scatci.py` copies the exact failed `scatci.inp`, `fort.16`, `fort.24`,
  `fort.26` and `fort.721` into a fresh B1 directory. Eight MPI ranks execute
  under individual GDB sessions through `rank_debugger.py`/`rank_gdb.py`.
  An intercepted signal and debugger status are recorded separately from a
  native inferior exit. Zero exit without CI/boundary exports is incomplete.
- `native_debug_owner.py` verifies the freeze and allocation before each
  launch, stops at failure/STOP/deadline and preserves each original exit.
  It schedules only the declared build, checks and one B1 diagnostic.

The native allocation is eight physical cores 4–11, 80-GiB memory, 8-GiB
shared backing, one library thread per rank and the original 16-GiB `memp`.
Native SCATCI is bounded to 24 hours, its container to 25 hours and the whole
owner to 30 hours. A separate hourly observer records progress.

## Interpretation and follow-up

Diagnostic output may identify a source line, bad index, class/pair or rank.
Even successful construction/eigenpair export is diagnostic evidence: it
does not resume the original scientific comparator or establish convergence.
Review any fault before declaring a correction or a matched 2048-root engine
control. A changed engine needs baseline reproduction before a scientific
root-count comparison. CAS12 pair sensitivity, omitted-spectrum, electronic
model, continuity, extraction, neutral-reference and production gates retain
their existing status.
