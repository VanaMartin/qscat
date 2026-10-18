# Resource-corrected two-point K-attribution evidence

The [completed study](../../NATIVE_K_ATTRIBUTION.md) is **false / original exits
1/1**. Four native calls exit0, native32 slices reproduce saved32 exactly and
native32→64 passes. Ten new independent controls and four reused rows produce
**24/30 passing checks**: step, endpoint and native-reference **K** vetoes reject
in each sector. All S/phase, order and invariant gates pass. Smaller steps
worsen the discrepancy; no further refinement is automatic.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/native-k-attribution-resource
```

[artifacts.json](artifacts.json) pins **15 public payloads** and
[manifest.json](manifest.json) records producer freeze
`4a558b35777df51ddedf2e2e2d5cede396708805`. The original producer implementation
is unchanged from`94fccee`; only the separately frozen timing allowances change.
`results.tar.gz` includes39 original files, preserving the cleanup traceback,
partial/final reports and0.2-s resources. Duplicate copied`fort.10/fort.21` are
omitted after checking them against the eight-file`inputs.tar.gz`. The eight-file
`prior.tar.gz` supplies saved native16/32 and the four reused independent rows.
`source.tar.gz` supplies all23 frozen source/prerequisite files.

Extract the archives under a fresh directory as follows: `results.tar.gz` into
`collected/`, `inputs.tar.gz` into`raw-preserved/`, `prior.tar.gz` into
`prior-preserved/`, and`source.tar.gz` into`source/`. Then run:

```bash
PYTHONPATH="$PWD/source" python verify-saved.py "$PWD"
```

The verifier uses the shared binary reader but independently expresses S/phase
and every family veto. It reconstructs four native matrices, fourteen independent
rows and thirty case checks. Verification exit0 means the **false** scientific
decision is reproduced; original scientific exits1/1 remain.

The published`numerical-worker.py` reexecutes **only the ten new independent
controls**, with four physical cores,4 GiB,450-s forecast and1200-s cap. It
reuses the four saved rows and saved native64 references; no native or inner
engine invocation occurs. It requires all reexecuted K and stored S components
to reproduce exactly and the case checks/verdicts to match. The public numerical
report and resources retain their own exits, separate from the original run.

The replay passes **exits0/0**, with all ten K matrices and stored S components
reproducing exactly and all thirty case checks unchanged. Original scientific
false/exits1/1 remain. Measured replay cost392.657 s /1449.435 sampled CPU-seconds /
231,624,704-byte recorded kernel memory high-water. Its17-file archive, separate
report, execution record and checksum-pinned owner are included in the bundle.
All15 final payloads are freshly downloaded and verified after supplementation.

[Fresh public-byte verification](../../verification-native-k-resource-public-20261019.json)
also checks binary/report agreement, all ten disk/report records and all29 tests
from downloaded source. The initial verification helper's grid-vs-matrix shape
error is separately recorded; it did not change producer bytes or scientific
results. The original180-s attempt remains unset/exits1/1; the original six-case
study remains false/exits1/1 with462 uncomputed energies.
