# Fresh-basis competing starts: a distinct aug-TZ orbital solution

Four fixed-geometry CAS(10,11) TZ↔aug-TZ starts complete with **original batch
and supervisor exits one**. Two attempts fail at the projection interface
before QC; two pass QC but disagree with the numerically qualified fresh aug-TZ
lineage. The original failures remain evidence.

| Projection | CI space | Original exit | Runner wall (seconds) | Kernel peak (MiB) |
|---|---:|---:|---:|---:|
| Fresh aug-TZ → TZ | 80 | 1, before QC | 0.403 | 98.17 |
| Fresh aug-TZ → TZ | 160 | 1, before QC | 0.403 | 99.54 |
| Fresh TZ → aug-TZ | 80 | 0 | 1737.00 | 434.28 |
| Fresh TZ → aug-TZ | 160 | 0 | 1705.00 | 491.52 |

Batch wall is **3444.76 seconds / 57.41 minutes**. All entries retain R=2.1323
bohr, two inactive core orbitals, ten active electrons, active irreps `[5,3,3,0]`,
forty equal-weight components and the existing tight tolerances. The ordinary
CI/orbital iteration caps remain **200/150**. Recipe:
[`calibration-sa11-fresh-basis-competing-starts.json`](../../calibration-sa11-fresh-basis-competing-starts.json).

## The two passing aug-TZ targets agree with each other

The space-80/160 pair passes roots, averaged objective, dipole and orbital
subspaces: maximum forty-root difference **5.941e-12 Hartree**, objective
difference **−1.137e-13 Hartree**, z-dipole difference **+1.901e-11 a.u.**,
and minimum active-subspace overlap **0.9999999999999984**. Their ground
energies are −112.89330697453909/−112.89330697453876 Hartree, objective energies
−112.45285158448974/−112.45285158448985 Hartree and z dipoles
0.15596686232972573/0.15596686234872786 a.u.

Against the earlier [fresh aug-TZ repairs](../fresh-augtz-repairs/README.md),
both comparisons **fail** the unchanged restart gate:

| New minus prior fresh solution | Space-80 result |
|---|---:|
| Averaged objective | **−0.538388 eV** |
| Ground energy | **−1.115756 eV** |
| Maximum root-energy difference | **2.126321 eV** |
| Ground z dipole | **−0.06670992 a.u.** |
| Minimum core overlap | 0.99983956 |
| Minimum active overlap | **0.14108460** |

The earlier branch passed QC, fixed-orbital coverage and independent import;
those checks do not establish a unique orbital optimum. The lower averaged
objective motivates qualifying the new branch, not declaring model convergence
or discarding the prior evidence. Its independent coverage and all-root/dipole
import remain separate gates.

## The failed downward projections have an interface repair

PySCF 2.11.0 rejects a 92-column source orbital matrix in a 60-orbital destination:
`RuntimeError: Too many orbitals in mo_init (try passing only the occupied orbitals)`.
The checkpoint adapter now passes the **13 core+active orbitals** when the source
has more columns than the destination. PySCF constructs the destination virtual
complement. Core/active sizes and irreps still must match, and simultaneous
geometry/basis changes retain their unsupported-system rejection.

On the exact failed checkpoint pair, the repaired projection agrees with the
documented truncated-input API and has **1.932e-14** AO-metric orthogonality
error with active irreps `[5,3,3,0]`. Four genuine PySCF regressions check
downprojection, same-basis restart, upward projection and simultaneous-change
rejection. All **80 CO tests pass with PySCF 2.11.0**; these optional tests skip
when that dependency is absent. The control is an interface check, not a passing
CASSCF optimization. Four fresh-name
[downward retries](../../calibration-sa11-tz-projection-coreactive-retries.json)
are queued behind coverage of both new aug-TZ checkpoints.

## Reconstruction

All **277 payloads** verify by size/digest, with byte-identical repackaging and
public download. The verifier reconstructs all **80 raw final QC states**, both
projection exceptions and all three root/objective/dipole comparisons. It also
independently recomputes orbital singular values from the preserved checkpoints.
The projection control retains its patched source, regression log and exact
failed input pair. Original engine/batch sources remain separate frozen payloads.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/competing-starts
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/competing-starts/competing-starts.tar.gz \
  -C "$EVIDENCE_DIR"
ROOT="$EVIDENCE_DIR/state-averaged-evidence"
PYTHONPATH=. uv run --with pyscf==2.11.0 python "$ROOT/verify-competing-starts.py" "$ROOT"
CONTROL="$ROOT/diagnostics/checkpoint-projection-interface-control"
PYTHONPATH="$CONTROL/source" uv run --with pyscf==2.11.0 python \
  "$CONTROL/verify-projection-interface.py" "$ROOT"
```

The [continuation](../../CONTINUATION.md) records the original owner resumption,
new coverage lease and dependent finite downward retry queue on Sadaharu.

The subsequent [new-branch coverage companion](../competing-augtz-coverage/README.md)
now passes all 96 probes/624 eigenpair evaluations, with physical/penalized
residual maxima 6.630e-10/9.977e-10 Hartree, in 10.17 minutes / 0.380 GiB.
Downward retries are running; the independent native import is queued.
