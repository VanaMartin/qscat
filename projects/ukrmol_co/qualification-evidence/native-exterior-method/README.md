# Native asymptotic normalization, independent Numerov tests and measured cost

The read-only native source fixes **f→sin(kr−lπ/2)/sqrt(k)** and
**g→cos(kr−lπ/2)/sqrt(k)**, with physical columns `f+gK+h c` and all closed
columns decaying. This resolves the sign and flux convention for an independent
real-axis reconstruction; it does not establish native finite-radius accuracy.
Four original source/licence files and twenty literal anchors are preserved.

The independent Python method uses ratio Numerov propagation and a directly
derived coupled inverse-power tail, without native Gailitis recurrences or
continued fractions. Design/source/tests/synthetic contract were frozen at
**`0a9b0ce` before the timing run**. See
[derivation](../../../../docs/physics/co-native-exterior-reconstruction.md) and
[validation](../../verification-native-exterior-method-20261017.json).

**17 focused tests pass**, including eight new solver tests. Two coupled smooth
finite-interval controls improve K errors by factors16.07/16.20 and16.12/16.20
under step halving, consistent with fourth order. A noninteger-Hankel tail
control gives value/derivative errors4.655e-8→1.830e-11→8.539e-14 at orders4/8/12.
Independent 60-digit exact polar hard-wall scattering agrees within3.648e-8 in
S. Other controls test unequal-momentum flux normalization, closed elimination
and invariance under arbitrary closed-column rescaling.

The synthetic **298-channel/100-interval** Sadaharu kernel costs **0.479163s**
(command1.015079s), sampled kernel high-water **82,604,032bytes**, original exits
**0/0**. This is a cost control with no native CO data, not a native accuracy
verdict. Its four0.2-second resource samples are not a post-exit peak. The
separately [frozen21-case pilot](../../native-exterior-pilot-contract.json) costs
an estimated975.5s by linear work extrapolation, with2000s upper estimate and
2400s cap. No native comparison was used to choose its settings.

[Fresh public replay](../../verification-native-exterior-method-public-20261017.json)
verifies **twelve payloads / ten Python-design inputs / four native files /
twenty literal anchors**, reruns all eight solver tests and reproduces the
synthetic command's exit0. The original engine/image are unchanged; the reused
mpmath wheel's earlier missing-dependency exit1 belongs to the scalar preflight.

## Replay

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/native-exterior-method
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/native-exterior-method"
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/source" "$REPLAY_ROOT/native" "$REPLAY_ROOT/dependencies"
tar -xzf "$EVIDENCE/source.tar.gz" -C "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/native-sources.tar.gz" -C "$REPLAY_ROOT/native"
python3 - "$EVIDENCE" "$REPLAY_ROOT" <<'PY'
import hashlib, json, pathlib, sys, zipfile
evidence, root = map(pathlib.Path, sys.argv[1:])
for manifest, directory in [("source-manifest.json", "source"),
                            ("native-source-manifest.json", "native")]:
    for name, digest in json.loads((evidence / manifest).read_text()).items():
        assert hashlib.sha256((root / directory / name).read_bytes()).hexdigest() == digest
dep = json.loads((evidence / "dependency.json").read_text())
wheel = evidence / dep["wheel"]
assert hashlib.sha256(wheel.read_bytes()).hexdigest() == dep["sha256"]
with zipfile.ZipFile(wheel) as archive:
    archive.extractall(root / "dependencies")
PY
cd "$REPLAY_ROOT/source"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  PYTHONPATH="$REPLAY_ROOT/source:$REPLAY_ROOT/dependencies" \
  uv run --with numpy --with scipy --with pytest python -m pytest -q \
  projects/ukrmol_co/test_exterior_numerov.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  PYTHONPATH="$REPLAY_ROOT/source:$REPLAY_ROOT/dependencies" \
  uv run --with numpy --with scipy python -m projects.ukrmol_co.native_exterior_study \
  --contract projects/ukrmol_co/native-exterior-timing-contract.json \
  --output "$REPLAY_ROOT/timing.json"
```

Native full K, CO poles and production remain unqualified. The original
[combined-exterior rejection](../combined-exterior/README.md), its separate
[bounded resolution pass](../combined-resolution/README.md), and the five-A2
strict interface failure retain their original meanings.
