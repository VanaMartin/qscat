# Bounded combined-exterior spectral-resolution control

This separate degree128/256/512 family passes all sixteen frozen toy-DtN checks
on Sadaharu, with **432 paired-method cases + four analytic controls**, exits
**0/0**, and no tolerance change. The [parent false gate and exits1/1](../combined-exterior/README.md)
remain preserved. Source/design/contract were frozen at **`db2158a` before
execution**; both numerical methods and the runner are unchanged.

See [results](../../COMBINED_EXTERIOR_RESOLUTION.md),
[design](../../../../docs/physics/co-combined-exterior-resolution.md),
[contract](../../combined-exterior-resolution-contract.json) and
[verification](../../verification-combined-resolution-20261017.json).
Nine public payloads include every complex matrix, the original empty log,
container metadata, sampled resources, source manifest and pinned dependency.

Fresh [public replay](../../verification-combined-resolution-public-20261017.json)
verifies all nine payloads, ten source hashes, every case and all sixteen check
verdicts, reproducing exit0. All-rung/finest host-public component drifts are
4.18531e-8/2.83234e-8. The parent public replay still reproduces its original
false verdict and exit1.

## Replay from public bytes

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/combined-resolution
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/combined-resolution"
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/source" "$REPLAY_ROOT/dependencies"
tar -xzf "$EVIDENCE/source.tar.gz" -C "$REPLAY_ROOT/source"
python3 - "$EVIDENCE" "$REPLAY_ROOT" <<'PY'
import hashlib, json, pathlib, sys, zipfile
evidence, root = map(pathlib.Path, sys.argv[1:])
for name, digest in json.loads((evidence / "source-manifest.json").read_text()).items():
    assert hashlib.sha256((root / "source" / name).read_bytes()).hexdigest() == digest
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
  --contract projects/ukrmol_co/combined-exterior-resolution-contract.json \
  --output "$REPLAY_ROOT/benchmark.json"
```

Expected scientific exit **0**, aggregate true and all sixteen checks true.
The parent replay separately expects exit1 and false. No native data is read or
native solve launched; full native K and CO complex-pole verdicts remain unset.
Cost:33.6106 seconds /0.1292 GiB sampled kernel high-water,166 samples at0.2s,
cores0–3 reserved,one BLAS thread,4-GiB limit. The dependency metadata's original
missing-mpmath exit1 belongs to the earlier scalar preflight and is preserved
with that explicit provenance.
