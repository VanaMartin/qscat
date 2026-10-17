# Separate matrix polar-tail and threshold-sheet controls

**All twenty fixed checks pass:** 189 tail matrix samples (135 real, 54 complex),
324 unequal-threshold matrix samples, 27 outgoing pole solves and three finite-
tail cutoff controls. The models are tested separately. See
[design](../../../../docs/physics/co-coupled-tail-benchmark.md),
[contract](../../coupled-tail-contract.json),
[results](../../COUPLED_TAIL_BENCHMARK.md),
[verification](../../verification-coupled-tail-20261017.json) and
[public replay](../../verification-coupled-tail-public-20261017.json).

Source/design/tests/contract were frozen at **`047524e` before the full matrix**.
The exact degenerate-threshold dipole tail has a noninteger-order Hankel oracle;
the unequal-threshold shell uses independent analytic Schur elimination and
high-precision Jost roots. Numerical propagation stays in the original channel
basis. Pole roots start once from the declared seed without oracle input or
retries. The decoupled pole matches the earlier independent Lambert-W oracle.

Nine payloads preserve **eleven source files**, host/local reports, exact
dependency wheel and original scientific exits **0/0**. The same mpmath wheel
is reused read-only from the scalar benchmark. Its metadata preserves that
earlier missing-dependency probe's exit 1 and explicitly identifies its origin;
this study did not have a new dependency failure.

## Replay

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/coupled-tail
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/coupled-tail"
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/source" "$REPLAY_ROOT/dependencies"
tar -xzf "$EVIDENCE/source.tar.gz" -C "$REPLAY_ROOT/source"
python3 - "$EVIDENCE" "$REPLAY_ROOT" <<'PY'
import hashlib, json, pathlib, sys, zipfile
evidence, root = map(pathlib.Path, sys.argv[1:])
for name, digest in json.loads((evidence / "source-manifest.json").read_text()).items():
    assert hashlib.sha256((root / "source" / name).read_bytes()).hexdigest() == digest, name
dep = json.loads((evidence / "dependency.json").read_text())
wheel = evidence / dep["wheel"]
assert hashlib.sha256(wheel.read_bytes()).hexdigest() == dep["sha256"]
with zipfile.ZipFile(wheel) as archive:
    archive.extractall(root / "dependencies")
PY
cd "$REPLAY_ROOT/source"
PYTHONPATH="$REPLAY_ROOT/source:$REPLAY_ROOT/source/libs/qscat:$REPLAY_ROOT/dependencies" \
  uv run --with numpy --with scipy python -m projects.ukrmol_co.coupled_tail_study \
  --contract projects/ukrmol_co/coupled-tail-contract.json \
  --output "$REPLAY_ROOT/benchmark.json"
```

Expected exit **0** and `separate_tail_threshold_benchmark_verdict=true`, with
all twenty checks true. Reports retain every S matrix, linear-system condition
number, all 27 root statuses/seeds, three accuracy rungs and cutoff controls.
Unitarity applies only to real energies; complex tail S values are compared
directly with the analytic oracle. The scalar shell source is included for its
independent decoupled-pole comparison, not as native CO input.

Sadaharu cost: **1.766 s wall / 0.050 GiB sampled kernel high-water**,
seven samples at 0.2-second intervals, reserved physical cores 0–3, one BLAS
thread, 4-GiB limit and 300-second timeout. Sampling is not a post-exit peak.
Persistent sources, wheel, logs and outputs are in
`prepared/coupled-tail-source-20261017/`. The pinned image remains
`sha256:48b9af05d30ce3c289e519433d9fc72a399b8caff6c109ea6c361e8a6b3599b6`.

The combined noncommuting unequal-threshold **and** polar-tail problem, native
amplitude normalization and CO complex poles remain unqualified. Earlier
five-A2 interface failure and all extraction/model rejections remain preserved.
