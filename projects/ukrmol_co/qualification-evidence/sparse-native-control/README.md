# First native sparse scattering solve and preserved interface failure

The equilibrium CAS(10,10) l=4 B1 control reuses the qualified forty-state dense
target and selects **128 of 8350** scattering eigenpairs. Its native PETSc/SLEPc
telemetry reports `mpisbaij`, Krylov–Schur, largest-magnitude selection, 256
Krylov vectors, 153 converged roots and 215 iterations. All 128 emitted states
are dense roots 1–128. Maximum errors are:

- eigenvalue: **5.969e-13 Hartree**;
- absolute eigenpair residual: **1.340e-10 Hartree**;
- sign-aligned continuum coefficient: **7.877e-8**.

Through the subsequent rejected outer-region replay, wall time is **292.932 s**
and kernel container peak is **1.228 GiB**. This is a scattering-only replay:
reused electronic, integral and configuration generation costs are separate.

The original reduced-coefficient SWINTERF continuation prints its unimplemented
partitioned-R-matrix diagnostic and stops with exit zero. RSOLVE then fails on
the missing boundary data. The driver and owner reject; B2 and larger root
counts under that owner remain unlaunched. **No phase or resonance verdict is
assigned.** The [native-path guide](../../SPARSE_SCATTERING.md) explains the
additional diagonal-return defect in the pinned partitioned reader and the
subsequent MPI-export controls.

## Reconstruct

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/sparse-native-control
mkdir -p /tmp/co-sparse-native-control
tar -xzf projects/ukrmol_co/qualification-evidence/sparse-native-control/sparse-control.tar.gz \
  -C /tmp/co-sparse-native-control
uv run python /tmp/co-sparse-native-control/state-averaged-evidence/verify-sparse-control.py \
  /tmp/co-sparse-native-control/state-averaged-evidence
```

The standalone verifier checks **212 payloads**, **1464 profile samples**, frozen
owner/source hashes, actual solver telemetry, native CIDATA framing and spectrum,
the dense continuum oracle, resource peaks and preserved downstream exits. Its
dense oracle retains all 8350 energies and only the first 128 continuum vectors;
the original native file's digest and extraction code are included. The pinned
inner/outer source archives, native source locators, installed SLEPc header,
compiler-sidecar bytes, downloaded development packages/checksums and original
scheduling/build/MPI-stdin failures are retained. Archives can be repackaged
byte-identically from sorted members and normalized tar/gzip metadata.

The producer's baseline commit is recorded alongside the exact experimental
source snapshots; the experimental files were uncommitted at execution. This
companion preserves the first rejected continuation, independently of later
successful or unsuccessful selected-spectrum experiments.
