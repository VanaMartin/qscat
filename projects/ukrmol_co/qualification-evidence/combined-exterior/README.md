# Combined polar exterior: original refinement rejection and native source trace

**The frozen combined gate rejects.** All 432 paired-method cases and four
analytic controls complete; **15/16 checks pass**. Spectral middle→finest change
**7.90991e-7** exceeds **1e-7** at ε=30,L=80,θ=20° for all three dipoles.
Original command/container exits **1/1** are preserved. Finest method agreement
**2.654e-9** does not replace the failed refinement check. No successor grid,
threshold relaxation or native solve is run.

See [design](../../../../docs/physics/co-combined-exterior-benchmark.md),
[contract](../../combined-exterior-contract.json),
[results](../../COMBINED_EXTERIOR_BENCHMARK.md),
[verification](../../verification-combined-exterior-20261017.json) and
[public replay](../../verification-combined-exterior-public-20261017.json).

Source/design/contract/tests were frozen at **`a3300b0` before the full matrix**.
The native surface/energy [derivation](../../../../docs/physics/co-native-boundary-conventions.md)
and analytic implementation/tests are independently frozen at **`d621541`**.
The external engine source is read-only, with its licences and manifest retained;
seven Fortran files, build provenance and three COPYING files form eleven inputs.
This source trace establishes the zero-Bloch length-valued surface sum, the
already absorbed inverse-sqrt(2), total-Hartree denominators and already-Rydberg
multipole coefficients. It does not establish full native K reproduction.

## Replay

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/combined-exterior
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/combined-exterior"
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/source" "$REPLAY_ROOT/dependencies" "$REPLAY_ROOT/native"
tar -xzf "$EVIDENCE/source.tar.gz" -C "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/native-sources.tar.gz" -C "$REPLAY_ROOT/native"
python3 - "$EVIDENCE" "$REPLAY_ROOT" <<'PY'
import hashlib, json, pathlib, sys, zipfile
evidence, root = map(pathlib.Path, sys.argv[1:])
for name, digest in json.loads((evidence / "source-manifest.json").read_text()).items():
    assert hashlib.sha256((root / "source" / name).read_bytes()).hexdigest() == digest
trace = json.loads((evidence / "native-trace.json").read_text())
for name, digest in trace["native_source_sha256"].items():
    assert hashlib.sha256((root / "native" / name).read_bytes()).hexdigest() == digest
for anchor in trace["literal_source_anchors"]:
    lines = (root / "native" / anchor["path"]).read_text().splitlines()
    for hit in anchor["matches"]:
        assert lines[hit["line"] - 1].strip() == hit["text"]
dep = json.loads((evidence / "dependency.json").read_text())
wheel = evidence / dep["wheel"]
assert hashlib.sha256(wheel.read_bytes()).hexdigest() == dep["sha256"]
with zipfile.ZipFile(wheel) as archive:
    archive.extractall(root / "dependencies")
PY
cd "$REPLAY_ROOT/source"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  PYTHONPATH="$REPLAY_ROOT/source:$REPLAY_ROOT/dependencies" \
  uv run --with numpy --with scipy python -m projects.ukrmol_co.combined_exterior_study \
  --contract projects/ukrmol_co/combined-exterior-contract.json \
  --output "$REPLAY_ROOT/benchmark.json"
```

Expected scientific **exit 1** after writing the complete report with
`combined_exterior_dtn_verdict=false`, solely `spectral_refinement=false`.
Retain it. Independent replay success means reproducing this rejection.
The exact mpmath wheel is reused read-only; its metadata preserves the earlier
scalar-study dependency probe's exit 1 with explicit historical provenance.

For the separate analytic surface tests, extract `native-python.tar.gz` to
another fresh directory, check its five file hashes against `native-trace.json`,
and run `python -m pytest projects/ukrmol_co/test_native_boundary.py` there.
It tests the free Neumann-box sum's normalization and its **finite omitted-tail
convergence**, not exact equality for a truncated spectrum. No saved CO boundary
file is read by these tests.

Sadaharu cost **11.178 s wall / 0.073 GiB sampled kernel high-water**, 54 samples
at 0.2 s, reserved physical cores 0–3, one BLAS thread and 4-GiB limit. Sampling
is not a post-exit peak. Persistent evidence:
`prepared/combined-exterior-source-20261017/`. Twelve public payloads preserve
the original failure, every complex boundary map, both source freezes and the
unmodified native normalization sources. CO K reconstruction/poles remain unset.
