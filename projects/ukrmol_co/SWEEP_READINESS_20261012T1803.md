# CO sweep readiness — 12 October, 18:03 UTC

**The declared CAS12 4096-root B1 address audit completed at15:04:29UTC.**
All eight native inferiors and debuggers exit0 with no intercepted signal.
MPI, container and owner also exit0. All ranks pass the graph/context/address
audits, complete the matrix classes and L2 update, then reach diagonalization
and native export. No pointer-failure or Fortran bounds-error marker appears.

The native measurement is **7.359945hours /20.523655GiB kernel peak** on
eight physical cores4–11, one library thread per rank,80-GiB cgroup and8-GiB
shared backing. The native timer reports17229.2s Hamiltonian construction
and9259.8s diagonalization. SLEPc prints a stopping iteration of159 with
tolerance1e-12/maxiter1000. These are diagnostic B1 costs, not the two-sector
scientific comparator's cost. Original release and debug failures remain
preserved, and their cause or relationship is still unknown.

## Independently read exports

The [completion capture](progress-pointer-audit-20261012T1803.json) includes
checksum-verified independent local binary reading using the current
`sparse_scattering.ci_spectrum`, `boundary_spectrum` and `interface_digest`:

- **4096finite, strictly energy-ordered poles**, from−112.91596570082095
  to−108.96270140373412Hartree.
- **4096×1410continuum coefficients** and **4096×250boundary amplitudes**,
  all finite; CI and boundary pole energies match exactly.
- Native exports are47,683,352bytes (`fort.25`) and8,727,072bytes (`fort.21`).
- `fort.821` remains the original712-byte header-only file; this run did not
  establish a reusable Hamiltonian checkpoint.

Per-root residual telemetry and a full-vector action-residual check are absent.
The saved diagnostic produces no B2 sector, Pi comparison or outer phase/fit
comparison. Thus successful native completion and raw-export reading establish
neither a fix for the earlier faults nor scientific root-count convergence.
The [address-audit contract](cas12-pointer-diagnostic-address-successor-contract.json)
continues to leave engine-correction, pair-stability and omitted-spectrum
verdicts unset. Rank0's normal exit reports IEEE underflow/denormal flags,
with native warning/error message counts both0; the raw notice is preserved.

## Evidence and next gate

All **55source/192input-evidence hashes**, **242installed artifacts** and five
saved native inputs reverify; numerical-library hashes match the prior debug
build. The capture pins60immutable remote records plus its frozen observation,
including debugger exits/logs, exports, resources and the independent reading.
Observer1191504 ends normally after nine captures07:44–15:44UTC. Both services
are inactive and the diagnostic container has exited; unrelated host workloads
remain present.

**Next:** review the successful instrumentation alongside the preserved release
and first-debug failures, then separately declare a matched2048-root control
for any engine/build selected for scientific use. Reproduce the qualified
baseline before resuming the4096-root two-sector comparator with residual,
interface,Pi,full-grid and fixed-window checks. Pair sensitivity still cannot
bound the asymptotic omitted spectrum. Electronic/channel convergence,
state/pole identity, extraction, neutral correlation/domain and a frozen
restartable production pilot continue to gate the sweep and potential fitting.
