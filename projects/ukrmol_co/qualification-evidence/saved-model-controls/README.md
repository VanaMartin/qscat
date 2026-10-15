# Saved CO model and extraction controls

The [declared study](../../saved-model-controls-contract.json) checks124 original
raw files across equilibrium CAS11/2048,CAS12/2048/4096 and stretched
five/nine-A2/2048. Every original exit,CI/boundary export,residual,grid and Pi
check passes. Both runtimes reproduce the **rejection** of equilibrium model
insensitivity and the four models' extraction-stability gates.

The seven public payloads include raw inputs,the57-file frozen reconstruction
source,and full host/local reports retaining all96 extraction cases,
288 full-data starts and480 held-point fits. This is raw scattering-output
reconstruction; it does not reconstruct the original target QC or full native
Hamiltonians. CAS12 import/production publication remains separate.

## Reproduce from public downloads

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/saved-model-controls
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/saved-model-controls"
# Choose a new directory outside the checkout.
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/raw" "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/raw-inputs.tar.gz" -C "$REPLAY_ROOT/raw"
tar -xzf "$EVIDENCE/reconstruction-source.tar.gz" -C "$REPLAY_ROOT/source"
cd "$REPLAY_ROOT/source"
PYTHONPATH=. uv run --with numpy --with scipy python -m projects.ukrmol_co.saved_model_controls \
  --root "$REPLAY_ROOT/raw" \
  --contract projects/ukrmol_co/saved-model-controls-contract.json \
  --output "$REPLAY_ROOT/reconstructed.json"
```

The raw archive is332166079 bytes,SHA256
`44cbb73ee3d0651c2c295d1df7a6100515e226f5156a2b845378437f6b0446cd`.
The source archive is104412 bytes,SHA256
`e60b3404aef2976c54069e87a087a728bf2cca230d2488b8a6a70afdba19185b`.
`source-provenance.json` identifies the producing commit and every source hash.
The reader verifies the124-file raw manifest before reading numerical records.
Host and macOS share the reader/fitter; this validates reconstruction and
cross-runtime consistency,not an independent scattering engine or absolute error.

See [the conclusions and scope](../../SAVED_MODEL_CONTROLS.md) and the
[quoted verification summary](../../verification-saved-model-controls-20261015.json).
