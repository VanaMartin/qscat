# CO calibration evidence

The public artifact pointer beside this file supplies a checksum-verified
`calibration-evidence.tar.gz` archive. Fetch it from a cloned QSCAT workspace:

```bash
uv run qscat-run fetch projects/ukrmol_co/evidence
mkdir -p projects/ukrmol_co/evidence/extracted
tar -xzf projects/ukrmol_co/evidence/calibration-evidence.tar.gz \
  -C projects/ukrmol_co/evidence/extracted
```

The archive contains lightweight **historical calibration evidence**, not a
qualified production scattering dataset: per-run decks, configurations,
templates/instrumented libraries, generated inputs, selected target/SCATCI/
RESON outputs, raw eigenphases/cross sections, stage/resource summaries,
failure logs, refit summaries, batch provenance and available source snapshots.
It also contains the quoted result snapshots, source-build provenance/log,
an archive file-hash index and applicable source licences.

The four early batches predate source snapshotting. Their harness hashes and
generated run inputs are saved, but their original Python revisions are not
all archived. Raw integrals, R-matrix amplitudes and K-matrices remain on the
calculation host and are not included in this lightweight archive. Recompute
them from the committed recipes to perform new RESON replays.

Rebuild the quoted 40-run/21-comparison calibration snapshot using the archive's
raw phase files and summaries:

```bash
uv run python -m projects.ukrmol_co.collect \
  projects/ukrmol_co/evidence/extracted/calibration-evidence \
  --pairs projects/ukrmol_co/calibration-pairs.json \
  --provenance projects/ukrmol_co/campaign-provenance.json \
  --output projects/ukrmol_co/evidence/recollected.json
```

Compare the parsed result with `projects/ukrmol_co/calibration-results.json`;
it should match exactly. This checks reconstruction of the saved evidence,
not electronic-model accuracy. Publication provenance in `manifest.json`
identifies the repository packaging commit; the historical calculations used
the image IDs/source hashes/configurations inside the archive. No claim that
the packaging commit produced the original engine runs is made.
