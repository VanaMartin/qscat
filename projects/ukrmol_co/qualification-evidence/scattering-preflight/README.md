# CAS(10,12) native scattering dimension and workspace preflight

Fresh native preparation reproduces the completed CAS(10,11) control and
measures the equilibrium CAS(10,12) scattering configuration space, reusing
the qualified forty-root/dipole targets. QC, target diagonalization and
DENPROP are not repeated. Only target integral generation and B1/B2 CONGEN
run; generated SCATCI inputs fix the exact retained inventory at five roots
in each of eight sectors.

| Per Pi sector | CAS(10,11) control | CAS(10,12) |
|---|---:|---:|
| Raw CONGEN configurations | 344124 | 990990 |
| Uncontracted L² configurations | 26136 | 84942 |
| Contracted continuum configurations | 1410 | 1410 |
| Contracted Hamiltonian dimension | **27546** | **86352** |

The two native preparations pass in **46.92 / 172.99 seconds**, with kernel
peaks **2.984 / 12.730 GiB**. The control dimension is also checked against
the original completed native SCATCI output. Missing L² markers and altered
retained-state inventories reject. Target inputs, actual native CIDATA sector
order, exact checkpoints and preserved original source hashes are retained.

One dense CAS(10,12) matrix needs **59653343232 bytes / 55.556 GiB**.
The checksum-pinned installed-library queries give:

| Process grid | Original exit | Valid aggregate array floor | Maximum rank array floor |
|---|---:|---:|---:|
| 2×2 | 5 | Invalid query; no memory estimate | Unset |
| 4×4 | 0 | 238790226304 bytes / 222.391 GiB | 13.956 GiB |
| 4×8 | 0 | 238967087872 bytes / 222.555 GiB | 6.983 GiB |

The valid floors include two local matrices, eigensolver/orthogonalization
workspaces, integer work arrays and eigenvalues. They exclude other live
engine arrays and MPI/library overhead. All rank-local sizes and aggregate
totals reconstruct from raw logs. The 2×2 overflow rejection is retained;
it cannot be used as an estimate. These are preparation/workspace results,
not completed scattering memory or MPI-scaling measurements. The audited dense
full-spectrum CAS(10,12) scattering path is beyond Sadaharu's RAM. A native
[sparse/iterative route](../../SPARSE_SCATTERING.md) is being qualified separately
before considering a large-host allocation.

## Preserved failures and public reconstruction

The original queued preflight exits before native work because its waiting
owner had already advanced. Four fresh preparations preserve subsequent
restart defects: copied `moints` aborts integral generation; enlarged CONGEN
arrays overflow the default workspace; skipped target stages leave retained
sector bookkeeping incomplete; restoring that inventory still exposes the
native `NBMX` workspace floor. The passing successor initializes the inventory
from unique new CIDATA set numbers, preserves skipped DENPROP execution,
and uses `NDIMX=10000000`, `CDIMX=NODIMX=1000000`,
`LNDO=NBMX=100000000`. Native CONGEN's `NBMX` backing workspace is independent
of `LNDO`; the pinned source documentation was inspected before the final retry.

All **814 payloads**, **245 retained source hashes** and **1099 raw resource
samples** verify. The Linux workspace executable and Fortran source travel
with the archive. Repackaging is byte-identical; public fetch checks the bytes.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/scattering-preflight
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/scattering-preflight/scattering-preflight.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run python "$EVIDENCE_DIR/state-averaged-evidence/verify-preflight.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

See the [original preflight contract](../../scattering-preflight-contract.json),
[qualified CAS(10,12) target](../cas12-import/README.md),
[large-host budget](../../LARGE_HOST.md) and [live handoff](../../CONTINUATION.md).
Publication sources and snapshots persist on Sadaharu under
`prepared/publication-cas12-preflight-20261007/`.
