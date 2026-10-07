# CAS(10,12) tight-restart SLEPc target import

The equilibrium cc-pVDZ CAS(10,12) target completes on **Sadaharu**, with original
engine/batch and controller exits **zero**. Every one of the **40 native roots**
and the independent DENPROP ground dipole passes the unchanged physical gates.
Wall is **22324.47 seconds / 6.201 hours**, aggregate CPU **9.437 hours**, and
kernel peak **43384504320 bytes / 40.405 GiB** in a 64-GiB container.

The [recipe](../../calibration-target-slepc-cas12.json) retains R=2.1323 bohr,
two inactive core orbitals, ten active electrons, active counts `[6,3,3,0]`,
five equally weighted roots per singlet/triplet C2v sector, forty ensemble
components, five computed target roots per sector and forty configured retained
states. This is a target-only job. It does not measure a scattering Hamiltonian
or a resonance. The four-rank SLEPc Krylov–Schur target path uses `crite=1e-12`,
500 iterations and 16-GiB internal budgets; QC retains the tight orbital/fresh-CI
settings, space 40 and ordinary 200-cycle CI / 150-cycle orbital limits.

## Independent gates

| Comparison | Maximum absolute error |
|---|---:|
| Native roots versus current QC | **4.801e-11 Hartree** |
| Nine-decimal saved roots versus current QC | **5.465e-10 Hartree** |
| Native roots versus all covered seed first-five spectra | **6.815e-9 Hartree** |
| Native roots versus independent RHF-start target | **3.485e-8 Hartree** |
| Current QC versus imported ground dipole | **1.279e-11 a.u.** |
| Imported dipole versus seed | **3.810e-9 a.u.** |
| Imported dipole versus RHF-start target | **5.881e-8 a.u.** |

All are below the **1e-7-Hartree root / 1e-5-a.u. dipole** gates. Seed
reoptimization shifts a root by at most **6.794e-9 Hartree**. Minimum initial/final
core/active overlap singular values are **0.9999999999999997 /
0.9999999999999988**; the separate RHF-start target matches the final active space
to **0.9999999999995508**. QC/fresh-CI solver flags, spin, MO orthogonality and Pi
degeneracy remain mandatory.

The initial checkpoint remains
`9451c7baed762ab674e9d36111cc220202972703ac52346e1ce13dcdc5ac55a0`.
Both exact previously published QC starts and the restart's full **48-probe /
312-eigenpair** coverage record travel with this companion. All five ensemble
roots pass every probe. The **two original extra-root failures** at space 40
remain failures; their spaces 80/160 pass. This import computes five roots per
sector and does not claim an eight-root CAS(10,12) import gate.

The verifier reads the unique newly solved CIDATA set in each native sector
(sets 1–8), checks requested roots, MPI ranks, backend and matrix dimension,
then compares every native root with every covered first-five spectrum. Saved
tables use nine decimal places and native output ten; exact-decimal rounding
errors stay within **5e-10 Hartree**, independently of the physical gate.

## Measured bottleneck

| Stage | Total wall time |
|---|---:|
| QC adapter, upstream slot named `psi4` | 1472.31 s / 24.54 min |
| Integral preparation | 0.095 s |
| Eight CONGEN stages | 12.22 s |
| Eight SCATCI sectors | **2808.79 s / 46.81 min** |
| Serial DENPROP | **18030.97 s / 5.009 h** |

DENPROP is **80.768%** of profiled wall. Selected-root target diagonalization
avoids the rejected all-spectrum dense eigensolver layout but retains a dense
PETSc Hamiltonian: largest dimension **70860**, matrix floor **40169116800
bytes / 37.41 GiB**. The measured complete-job peak is **40.405 GiB**.
These are single-shot timings under recorded concurrent local work, not MPI
scaling or an external-host completion forecast. The original dense failure and
approximately $200 deferred external experiment remain documented in
[`LARGE_HOST.md`](../../LARGE_HOST.md).

## Reconstruct

All **335 payloads** verify by size/digest, with byte-identical repackaging and
public fetch. The executable verifier reconstructs **40 raw current QC states**,
all native roots/dipole, all 48 raw reference spectra/spin/convergence blocks,
both initial/final subspace comparisons and **111339 raw resource samples**.
Original v3 source and console retain the earlier CAS(10,11) retry batch exit
one alongside this successful CAS(10,12) entry; their scientific verdicts remain
distinct. The 46-attempt historical aggregate and its pointers remain frozen.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/cas12-import
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/cas12-import/cas12-import.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run --with pyscf==2.11.0 python \
  "$EVIDENCE_DIR/state-averaged-evidence/verify-cas12-import.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

The [continuation](../../CONTINUATION.md) records completion and release of CPUs
8–11 to the gated stretched CAS(10,11) import/basis queue. CAS(10,12) scattering
still needs its actual contracted dimensions, workspaces, memory/scaling and
electronic-model checks. The numerical target import is now locally established.

A [native scattering preflight](../../scattering-preflight-contract.json) is
now queued on Sadaharu. It reuses qualified target data, reproduces the completed
CAS(10,11) dimension as a control, generates native CAS(10,12) integral/CONGEN
preparation and queries installed-library workspaces at its inferred contracted
dimension. This avoids repeating DENPROP and retains its separate scientific scope.
