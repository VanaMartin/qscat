# Independent native exterior pilot: original rejection

**The frozen21-case equilibrium-B1 pilot rejects**, preserving original
scientific command/container exits **1/1**. All21 cases complete. Eight of ten
checks pass; finest native S-entry error **1.49310e-6** and phase error
**2.72174e-6rad moduloπ** exceed the unchanged1e-6 budgets at5eV. K-entry,
unitarity,reciprocity and independent step/endpoint/tail-order checks pass.

See [results](../../NATIVE_EXTERIOR_PILOT.md),
[contract](../../native-exterior-pilot-contract.json),
[derivation](../../../../docs/physics/co-native-exterior-reconstruction.md) and
[verification](../../verification-native-exterior-pilot-20261017.json).
Method/source/tests are frozen at **`0a9b0ce`**; the separate native pilot
contract is frozen at **`00c5839` before comparison**. No budget/grid adjustment
is used to replace the false result. The
[method/source controls](../native-exterior-method/README.md) remain separately
passing; no new native engine solve occurs.

Nine public payloads preserve the full report,21-row journal,empty original
log,original exits/container metadata,eleven source/design/contract inputs and
four checksum-pinned native input files with their collection provenance.
[Fresh public numerical reconstruction](../../verification-native-exterior-pilot-public-20261017.json)
verifies all nine payloads,eleven frozen source/design inputs,four native files,
all21 case records and all ten check verdicts. It reproduces **false and exit1**.
Maximum host/public K-entry/S-component drifts are1.61643e-9/3.23265e-9.

## Replay from public bytes

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/native-exterior-pilot
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/native-exterior-pilot"
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/source" "$REPLAY_ROOT/raw"
tar -xzf "$EVIDENCE/source.tar.gz" -C "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/pilot-inputs.tar.gz" -C "$REPLAY_ROOT/raw"
python3 - "$EVIDENCE" "$REPLAY_ROOT" <<'PY'
import hashlib, json, pathlib, sys
evidence, root = map(pathlib.Path, sys.argv[1:])
for name, digest in json.loads((evidence / "source-manifest.json").read_text()).items():
    assert hashlib.sha256((root / "source" / name).read_bytes()).hexdigest() == digest
contract = json.loads((root / "source/projects/ukrmol_co/native-exterior-pilot-contract.json").read_text())
for name, digest in contract["input_sha256"].items():
    assert hashlib.sha256((root / "raw" / name).read_bytes()).hexdigest() == digest
PY
cd "$REPLAY_ROOT/source"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  PYTHONPATH="$REPLAY_ROOT/source" \
  uv run --with numpy --with scipy python -m projects.ukrmol_co.native_exterior_study \
  --contract projects/ukrmol_co/native-exterior-pilot-contract.json \
  --directory "$REPLAY_ROOT/raw" --output "$REPLAY_ROOT/pilot.json"
```

Expected exit **1**, `pilot_verdict=false`, failed checks exactly
`native_s_error` and `native_phase_error_rad_mod_pi`; the complete report and
journal must be retained. Full native K and CO complex-pole verdicts remain
unset. Cost on Sadaharu:735.161s/0.08564GiB sampled kernel high-water,3671
samples at0.2s,physical cores0–3 reserved,one BLAS thread,4-GiB cap.
