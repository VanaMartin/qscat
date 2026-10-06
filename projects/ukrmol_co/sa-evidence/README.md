# CO state-averaged continuation evidence

This archive preserves completed **calibration** calculations and failed
attempts from the multi-spin/multi-irrep PySCF continuation. The numerical
engine and target import checks are established; the electronic model is still
unqualified for production potential fitting. The original RHF/SEP and
ground-state-CASSCF archive remains in [`../evidence/`](../evidence/README.md).

The published snapshot contains **37 completed attempts**: 17 validated
pipelines (12 target-only and five scattering), 19 failed runner stages and
one pre-run setup rejection. It includes six paired comparisons and 48 native
RESON replays. Public-client verification checked all **3,831 payload digests**
and exactly reconstructed the aggregate; all 17 successful runs were reanalyzed
from their raw outputs, and repackaging produced identical archive bytes.

Fetch the checksum-verified archive from a cloned QSCAT workspace, then extract
outside the checkout to avoid duplicate historical Python test packages:

```bash
uv run qscat-run fetch projects/ukrmol_co/sa-evidence
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/sa-evidence/state-averaged-evidence.tar.gz \
  -C "$EVIDENCE_DIR"
EVIDENCE_ROOT="$EVIDENCE_DIR/state-averaged-evidence"
uv run python -m projects.ukrmol_co.collect "$EVIDENCE_ROOT" \
  --pairs "$EVIDENCE_ROOT/calibration-sa-pairs.json" \
  --provenance "$EVIDENCE_ROOT/sa-provenance.json" \
  --output "$EVIDENCE_DIR/reconstructed.json"
```

The reconstructed JSON matches the archive's `sa-results.json`. The file index
names every payload's byte count and SHA256. Per-batch source snapshots retain
the exact runner/backend/adapter revisions, including the failed revisions.
The archive also carries current analyzer/collector/packager sources, build/test
provenance, licences, generated inputs, QC/CI/scattering/fitting logs, raw phases
and cross sections, resource/stage records, and binary RHF/CASSCF checkpoints.
Projected starts retain the initial checkpoint copy and its hash.

The small checkpoints make the saved projected starts reusable. Put extracted
`runs/` below a new host root mounted at `/work` when replaying the historical
`/work/runs/...` checkpoint arguments. Give recomputed outputs and batches new
names. Early import probes used a different molecular origin; per-run atom
coordinates and the source snapshots preserve that distinction.

Large integrals, channel files, R-matrix amplitudes and K-matrices remain on the
calculation host. Recompute them from the recipes for new native RESON replays.
Saved replay inputs/outputs are included. Native automatic fit candidates may
overlap or be unphysical; their count is not a physical resonance count.

`projects/ukrmol_co/archive.py` packages only the batches and runs named by a
collector snapshot. It uses sorted regular files, zero uid/gid/mtime, mode 0644,
and gzip with zero mtime and an empty filename. Additional provenance/licence
files are explicit attachments:

```bash
uv run python -m projects.ukrmol_co.archive "$COPIED_ROOT" \
  --snapshot projects/ukrmol_co/sa-results.json \
  --attachments projects/ukrmol_co/sa-provenance.json \
    projects/ukrmol_co/calibration-sa-pairs.json LICENSE \
    projects/ukrmol_co/analyze.py projects/ukrmol_co/collect.py \
    projects/ukrmol_co/archive.py \
  --output /path/outside/checkout/state-averaged-evidence.tar.gz
```

Include the engine test/build logs and upstream licences as additional
attachments for a new publication. The output must be a new file. Fetched
archives and local run/scratch directories are excluded from Docker contexts
and batch source snapshots.
