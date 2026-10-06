# CO state-averaged continuation evidence

This archive preserves completed **calibration** calculations and failed
attempts from the multi-spin/multi-irrep PySCF continuation. The numerical
engine and target import checks are established; the electronic model is still
unqualified for production potential fitting. The original RHF/SEP and
ground-state-CASSCF archive remains in [`../evidence/`](../evidence/README.md).

The published snapshot contains **55 completed attempts**: 29 validated
pipelines (17 target-only and twelve scattering), 23 failed runner stages,
two failed target-import analyses and one pre-run setup rejection. It includes
twelve paired comparisons, 96 native RESON replays and six fixed-orbital CI
coverage probes. The public fetch client verified all **6,616 payload digests**
and exactly reconstructed the aggregate. All 29 successful runs passed raw-output
reanalysis; the six CI probes match independent UKRmol roots, and repackaging
produced identical archive bytes. Earlier snapshots remain available at their
immutable URLs in the manifest.

The CAS(10,12) 6-GiB-per-process retry passed QC/fresh-CI checks and singlet-A1
diagonalization, then failed the triplet-A1 matrix allocation (9.31 GiB per rank).
Its complete failure, resources and source snapshot are included. CAS(10,11)
fresh-CI reoptimization passes every averaged-root import and the dipole check;
its batch also preserves an A1-only control rejected for computed Pi splitting.
Every launched continuation batch is complete and included. Active-space target
properties remain model-sensitive; import agreement does not establish
electronic-model convergence.

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
`diagnostics/cas11-ci-coverage/` retains the six probes, checkpoint, execution
record and exact source snapshot. They recover a missing lowest CI root in
the rejected CAS(10,11) target; they do not validate its orbital optimization.

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
collector snapshot, plus explicitly selected completed diagnostics. It uses
sorted regular files, zero uid/gid/mtime, mode 0644,
and gzip with zero mtime and an empty filename. Additional provenance/licence
files are explicit attachments:

```bash
uv run python -m projects.ukrmol_co.archive "$COPIED_ROOT" \
  --snapshot projects/ukrmol_co/sa-results.json \
  --diagnostics cas11-ci-coverage \
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
