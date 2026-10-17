# Native radial-subrange attribution: passing bounded control

**All ten frozen checks pass; original scientific command/container exits0/0.**
Native8→16→32 radial-subrange refinement accounts for **99.72%** of the measured
5-eV S mismatch against a fixed finer independent reference. Finest native32
versus independent K/S/phase errors are **2.68388e-9 /4.08632e-9 /1.15175e-8rad
moduloπ**,with unchanged1e-6 comparison budgets.

See [results](../../NATIVE_EXTERIOR_ATTRIBUTION.md),
[verification](../../verification-native-attribution-20261017.json),
[design](../../../../docs/physics/co-native-exterior-attribution.md) and
[contract](../../native-exterior-attribution-contract.json). Source/design/tests
were frozen at **`666c056` before execution**,following fresh public
reconstruction of the [original pilot rejection](../native-exterior-pilot/README.md).
That false verdict and original exits1/1 remain unchanged.

The literal baseline reproduces all99 saved K matrices exactly. Two declared
one-line decks change only NRANGE to16/32; all99 energies,298 channels,2048
uncorrected poles,100-bohr matching radius,NLEG10 and native Gailitis defaults
remain fixed. Two new independent5-eV cases use step0.0125,endpoints100/400
and tail order8. Native16→32 changes at5eV are2.85943e-11 in S and3.84278e-11
in phase; independent step/endpoint S changes are4.71630e-8/2.54094e-9.

Ten public payloads preserve eleven source/design inputs,the full report,
two-case journal,all297 native energy records,all native decks/logs,source
and input digests,original engine/container exits and sampled resources.
Native output archives omit redundant copies of fort.10/21; use the identical
checksum-pinned originals from `inputs.tar.gz` when reexecuting native stages.
That input archive also binds the prior false pilot and its passing public
reconstruction prerequisite. No new inner solve occurs.

[Fresh public reconstruction](../../verification-native-attribution-public-20261017.json)
verifies ten payloads,eleven frozen sources,nine native outputs,four original
inputs and297 energy records. Both independent cases and all ten gate decisions
reproduce,as does the conditional attribution,exit0. Four public tests pass;
maximum host/public K-entry/S-component drift is9.48917e-10/1.52518e-9.
It does not independently reexecute the native engine;
the actual source-pinned three native invocations and original exits are retained
in `execution.json` and each native log. Native-only16→32 refinement over99
energies is not independent full-grid reconstruction.

## Fetch and inspect

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/native-attribution
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/native-attribution"
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/source" "$REPLAY_ROOT/raw" "$REPLAY_ROOT/native"
tar -xzf "$EVIDENCE/source.tar.gz" -C "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/inputs.tar.gz" -C "$REPLAY_ROOT/raw"
tar -xzf "$EVIDENCE/native-results.tar.gz" -C "$REPLAY_ROOT/native"
python3 - "$EVIDENCE" "$REPLAY_ROOT" <<'PY'
import hashlib, json, pathlib, sys
evidence, root = map(pathlib.Path, sys.argv[1:])
for manifest, directory in [(evidence / "source-manifest.json", "source"),
                            (root / "native/native-output-manifest.json", "native"),
                            (root / "raw/input-manifest.json", "raw")]:
    for name, digest in json.loads(manifest.read_text()).items():
        assert hashlib.sha256((root / directory / name).read_bytes()).hexdigest() == digest
PY
cd "$REPLAY_ROOT/source"
PYTHONPATH="$REPLAY_ROOT/source" uv run --with numpy --with scipy --with pytest \
  python -m pytest -q projects/ukrmol_co/test_native_exterior_attribution.py
```

On Sadaharu with the same read-only engine mounts as `execution.json`,the frozen
driver can reexecute the full native/independent study in a **fresh output**:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m projects.ukrmol_co.native_exterior_attribution \
  --directory "$REPLAY_ROOT/raw" --prior "$REPLAY_ROOT/raw/prior-pilot.json" \
  --replay "$REPLAY_ROOT/raw/prior-public-replay.json" \
  --engine /opt/ukrmolp/bin.double/rsolve \
  --contract projects/ukrmol_co/native-exterior-attribution-contract.json \
  --output "$REPLAY_ROOT/new-study"
```

Expected control verdict true,ten checks true,conditional native-subrange
attribution true,exit0. The original pilot verdict remains false; full native K,
CO poles and production stay gated. Measured cost1339.551s/0.56024GiB sampled
kernel high-water,6689 samples at0.2s,physical cores0–3 reserved,one BLAS thread.
