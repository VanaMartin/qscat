# Endpoint-isolated analytic identifiability audit

**94 known-wrong-width counterexamples remain after explicitly withholding all
window endpoints. All 45 mandatory true-pole smooth-background controls pass.**
The fixed 405-cell matrix and scientific tolerances are unchanged. Corrected
source was frozen at `69f13b6` before repeating the matrix.

The [original producer run](../phase-identifiability/README.md) reported 98
counterexamples, but its endpoint-unused declaration independently fails in
14/45 grid/window combinations. Original source, outputs, exit 0 and producer
qualification remain preserved. They do not qualify the fully unused protocol.
The corrected reading verifies development disjointness against all three audit
sets in all 45 grids, retaining every cell. Counts are compressed-like 27,
equilibrium-like 13, narrow-control 54, resolved/stretched-like 0. These remain
local analytic/unitary counterexamples, not causal finite-range CO models.

See [results and physical limits](../../PHASE_IDENTIFIABILITY.md),
[corrected quoted verification](../../verification-phase-identifiability-endpoint-20261017.json)
and [public verification](../../verification-phase-identifiability-endpoint-public-20261017.json).
Seven payloads contain six frozen source files, all grids, coefficients, cells,
host/local reports and resources. Host/public cell verdicts all agree; maximum
phase-error metric drift is **3.02e-14 rad**. Original corrected command/container
exits are **0/0**, wall time **0.965 s**, sampled kernel high-water **0.060 GiB**.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/phase-identifiability-endpoint
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/phase-identifiability-endpoint"
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/source.tar.gz" -C "$REPLAY_ROOT/source"
cd "$REPLAY_ROOT/source"
python3 - "$EVIDENCE/source-manifest.json" <<'PY'
import hashlib, json, pathlib, sys
for name, digest in json.loads(pathlib.Path(sys.argv[1]).read_text()).items():
    assert hashlib.sha256(pathlib.Path(name).read_bytes()).hexdigest() == digest, name
PY
PYTHONPATH=. uv run --with numpy python -m projects.ukrmol_co.phase_identifiability \
  --contract projects/ukrmol_co/phase-identifiability-contract.json \
  --output "$REPLAY_ROOT/audit.json"
```

Expected exit 0 qualifies representation controls, not width identification or
native use. No native data is read and no native successor is released. Sources,
outputs, logs and original resource records remain under
`prepared/phase-identifiability-endpoint-source-20261017/` on Sadaharu.
