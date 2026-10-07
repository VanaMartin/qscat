# Angular-cutoff controls and competing aug-TZ native import

Both CAS(10,11)/cc-pVDZ angular controls pass the declared position, full-width
and phase gates against the equilibrium **l=4** calculation. The complete
99-point grids and B1/B2 symmetry contributions are retained.

| Control | Position change | Full-width change | Maximum phase change modulo π | Wall | Kernel peak |
|---|---:|---:|---:|---:|---:|
| l=3 versus l=4 | +3.9602 meV | +0.7434% | 0.0166575 rad | 10621.10 s | 24.329 GiB |
| l=5 versus l=4 | +0.07987 meV | +0.2655% | 0.0162174 rad | 14353.38 s | 26.300 GiB |

Independent full-data fits on common windows, with backgrounds 1–4, change
position/width by at most **3.982 meV / 0.7214%** for l=3 and **0.3930 meV /
0.2279%** for l=5. All **36 native background/detection replays** reconstruct
from raw output. Detection thresholds 0.7/1.0/1.3 give identical candidates at
each fixed background; changing the background still changes width by
**12.60% / 12.65% / 11.98%** for l=3/4/5. That exceeds the 5% extraction gate.
Passing angular controls therefore does not qualify a production width.
These complete-grid controls are separate from the earlier fine-grid automatic
fits whose `MAXFIT=100` truncation leaves acceptance verdicts unset.

The competing equilibrium aug-TZ branch also passes its covered **64-root**
native import: maximum seed-reference error **3.8811e-8 Hartree**, native/seed
dipole difference **1.8756e-7 a.u.**, and required forty-root/current-QC error
**4.9431e-10 Hartree**. Cost is **3259.47 s / 5.260 GiB**. Exact-decimal saved
table rounding is checked independently of the 1e-7-Hartree physical gate.
Both reproducible aug-TZ branches remain scientifically distinct; this import
does not establish a unique orbital optimum or geometry/state continuity.

## Public evidence and reconstruction

The companion checks **1014 payloads**, **39 batch-source hashes** and
**219998 raw resource samples**. It retains all three angular runs, their
replays, the competing native import, its seed checkpoint and covered reference
spectra, original completed ownership records, and executable reconstructions.
Archive repackaging is byte-identical and public fetching verifies the bytes.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/continuum-competing
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/continuum-competing/continuum-competing.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run --with pyscf==2.11.0 --with h5py==3.15.1 \
  python "$EVIDENCE_DIR/state-averaged-evidence/verify-continuum-competing.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

See the [competing-branch coverage](../competing-augtz-coverage/README.md),
[downward projections](../downward-projections/README.md),
[angular recipe](../../calibration-sa11-tight-continuum.json) and
[live handoff](../../CONTINUATION.md). Publication sources and snapshots persist
on Sadaharu under `prepared/publication-continuum-competing-20261007/`.
