# CO sweep readiness — 12 October, 07:01 UTC

**The instrumented B1 diagnostic captured a source-level fault and has stopped.**
Owner1168755 finishes at **00:26:31 UTC**, with container/controller exits1.
GDB intercepts **SIGSEGV on rank0** in
`csf_module::compare_excitations_fast`,
`source/mpi-ci-diag/Modules/Hamiltonian/CSF_module.F90:506`:

```fortran
manager => this % csf % manager
```

The call stack passes through `Hamiltonian_module.f90:207` (`slater_rules`),
`Contracted_Hamiltonian_module.f90:494` (`evaluate_class_1`, target-symmetry
index8), its Class1 loop at181, and `mpi_scatci.f90:259`.
This identifies an observed failing pointer-access operation. The underlying
cause and an engine correction are not established by the stack alone.

**This diagnostic fault is earlier than the preserved release-engine failure.**
The original4096-root comparator fails on rank4 after Class4 completion;
the instrumented build stops on rank0 in Class1, before any Class5/6/7 marker.
The two failures may share a cause, but that relationship remains unproven.
Unchanged log tails during the earlier status checks did not locate the active
computation: the backtrace supplies the first direct Class1 location.

The failed diagnostic costs **3867.5975 seconds =1.07433hours**, with
**13.87420GiB inclusive cgroup peak**, under the declared80-GiB allocation.
There is no timeout, Fortran bounds-error report, CI export or boundary export.
GDB and MPI return139 after the intercepted signal; this is a debugger-controlled
status, **not a recorded native inferior exit code**. Peer ranks are terminated
by MPI and provide no independent completed rank-exit records.

All **47 source /119 input-evidence hashes** and **242 installed artifacts**
reverify after failure. The [result capture](progress-diagnostic-result-20261012T0701.json)
pins36immutable remote execution/resource/rank-log/native-source/observer records
plus an immutable copy of the observation. Both original failures and the
successful2048-root pilot remain preserved. Observer1169951 ends normally at
**01:18:52UTC** after three captures. No scientific container remains active.

**Next:** inspect CSF object initialization, pointer lifetime and configuration/
determinant-index provenance around the reported dereference. Any additional
pointer-state diagnostic, matched-engine2048 control or correction needs a
separate declaration. A changed engine must reproduce its baseline before
scientific root-count comparison resumes. CAS12 pair sensitivity and omitted-
spectrum verdicts remain unset. Electronic/channel convergence, state/pole
identity, extraction, neutral correlation/domain and a frozen restartable
production pilot continue to gate the sweep and held-out potential fitting.
