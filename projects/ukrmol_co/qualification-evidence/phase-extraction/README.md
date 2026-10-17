# Prospective phase-extraction calibration and saved-anchor evidence

**The 30-case analytic calibration passes its frozen requirements; application
to saved nine-A2 phases rejects compressed/equilibrium and provisionally passes
stretched, in both Pi sectors.** Original calibration/native-reading command
exits are **0/1**, respectively. No native scattering calculation was repeated.

The [design](../../../../docs/physics/co-phase-extraction.md),
[contract](../../phase-extraction-contract.json),
[quoted verification](../../verification-phase-extraction-20261017.json) and
[interpretation](../../PHASE_EXTRACTION.md) explain the selection, untouched-audit
split, width-counterexample veto and limited meaning of acceptance. All 30 cases,
failed fits, alternative widths, folds and starts are downloadable. The 13
payloads contain both source freezes, host/local reports and resource records;
raw scattering archives are reused from their existing pointers.

Calibration source is frozen at `902b76c`; the native reader at `4e089e3` checks
the calibrated source hashes before opening the raw data. Local/host verdicts
agree, with maximum native fitted-parameter drift **7.50e-9 eV** over every
development and inner-fold start. These replay the same algorithm on independent
platforms; the analytic S-factor generator supplies the independent known-pole
oracle. Raw reconstruction is not an independent scattering engine.

## Reproduce outside the checkout

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/phase-extraction
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/compressed-a2-scattering
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/equilibrium-a2-scattering
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/saved-model-controls
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence"
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/source" "$REPLAY_ROOT/raw/compressed" \
  "$REPLAY_ROOT/raw/equilibrium" "$REPLAY_ROOT/raw/saved"
tar -xzf "$EVIDENCE/phase-extraction/reader-source.tar.gz" -C "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/compressed-a2-scattering/raw-inputs.tar.gz" -C "$REPLAY_ROOT/raw/compressed"
tar -xzf "$EVIDENCE/equilibrium-a2-scattering/raw-inputs.tar.gz" -C "$REPLAY_ROOT/raw/equilibrium"
tar -xzf "$EVIDENCE/saved-model-controls/raw-inputs.tar.gz" -C "$REPLAY_ROOT/raw/saved"
cd "$REPLAY_ROOT/source"
python3 - "$EVIDENCE/phase-extraction/reader-source-manifest.json" <<'PY'
import hashlib, json, pathlib, sys
for name, digest in json.loads(pathlib.Path(sys.argv[1]).read_text()).items():
    assert hashlib.sha256(pathlib.Path(name).read_bytes()).hexdigest() == digest, name
PY
PYTHONPATH=. uv run --with numpy --with scipy python -m projects.ukrmol_co.phase_extraction_study \
  --contract projects/ukrmol_co/phase-extraction-contract.json \
  --output "$REPLAY_ROOT/calibration.json"
PYTHONPATH=. uv run --with numpy --with scipy python -m projects.ukrmol_co.phase_extraction_native \
  --root "$REPLAY_ROOT/raw" --calibration "$REPLAY_ROOT/calibration.json" \
  --contract projects/ukrmol_co/phase-extraction-contract.json \
  --output "$REPLAY_ROOT/native.json"
```

Expected final exit is **1** because the three-anchor extraction verdict is
false; `native.json` is written first. Do not suppress that exit in automation.
The calibration command must exit zero before the native reading is launched.
The raw readers check all **270 original manifest entries**, stages, exports,
resources, native fits and numerical prerequisites before new fitting.

The unchanged calibration logic accepts 17/30 cases. The other 13 are retained
rejections/unset outcomes, including every smooth-nonpolynomial case. Passing
the finite calibration requirements therefore does not qualify arbitrary smooth
backgrounds. Audit-only phase mutation and failing-calibration tests enforce the
data-use boundary. The next estimator requires a new physical/analytic design;
this run authorizes no relaxed tolerance, dropped window or automatic successor.
