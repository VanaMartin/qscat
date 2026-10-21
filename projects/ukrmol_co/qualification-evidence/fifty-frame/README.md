# Independently qualified new fifty-component frame

The [result](../../FIFTY_FRAME_COVERAGE.md) passes all48 fixed-orbital CI probes
and all74 saved native target energy/vector/header/EOF checks, original
command/container/owner0/0/0. The original seed-preservation rejection and1/1/1
exits remain separate and unchanged. New-frame qualification does not mean that
the globally selected50-state property file contains the desired five-A2 manifold.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/fifty-frame
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/fifty-import
```

Five new payloads fresh-download-verify; four prerequisite payloads reauthenticate.
`prerequisite-artifacts.json` names the latter's exact content-addressed bytes.
The new `results.tar.gz` preserves108 files /970,503 uncompressed bytes, including
all48 probe logs/results, final/progress/coverage reports, owner/commands/exits,
resource samples and10-file source archive/manifest. All files are pinned in
`collection-sha256.json`; the archive SHA256 is
`7a739e6354ddeb803d969a0a2123a9ab6bfabd2362625a32281791fe8a272d9a`.
The complete original target/seed/binary inputs remain in the prerequisite bundle.

Extract this collection into a fresh directory and its `source.tar.gz` into
`source/`. Extract the prerequisite target collection separately into `prior/`.
Authenticate the10 sources,108 collection files and10 prerequisite input hashes
before interpreting outcomes. The original new-frame computation is frozen at
`55e89a6`, with explicit prospective bounds and controls in the source contract.

For a lightweight saved-native reconstruction, from the extraction root:

```sh
PYTHONPATH=source python - <<'PY'
import json
from pathlib import Path
from projects.ukrmol_co.fifty_import_verification import covered_sectors, native_binary_sectors
report = json.loads(Path('outputs/coverage/result.json').read_text())
reference = covered_sectors(report, report['checkpoint_sha256'])
rows = native_binary_sectors(Path('prior/run/output/CO/geom1/fort.26'), reference)
assert len(report['probes']) == 48 and sum(row['roots'] for row in rows) == 74
assert report['original_seed_preservation_passed'] is False
print('48 saved probe decisions and74 complete native roots pass; old seed rejection retained')
PY
```

[Public saved verification](../../verification-fifty-frame-public-20261021.json)
also authenticates every fetched file and frozen original input, recomputes the
new/seed root changes and ensemble/reference comparisons, and independently
recomputes the native vector Grams. CI residual-action values and exported-frame
overlaps are qualified producer observations; reading JSON is not a fresh CI
solve. A new numerical replay of all48 probes is outside this saved reading.

Measured cost5.144min /.27065GiB; no native engine or orbital optimizer is called.
The [selected40-property successor](../../FIFTY_SELECTED_PROPERTIES.md) is a
separate native property calculation. It does not authorize a new inner solve.
