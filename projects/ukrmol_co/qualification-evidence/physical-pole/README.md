# Physical radial-shell real-axis and pole benchmark

**The complete prospectively frozen matrix passes:** 36 real-axis profiles and
81 outgoing shooting solves, with three accuracy levels and three matching
radii. The independent oracle uses exact Jost functions and Lambert W at
60/90 decimal digits. No CO output is read. See the
[design](../../../../docs/physics/co-physical-pole-benchmark.md),
[fixed contract](../../physical-pole-contract.json),
[results](../../PHYSICAL_POLE_BENCHMARK.md),
[verification](../../verification-physical-pole-20261017.json) and
[public replay](../../verification-physical-pole-public-20261017.json).

Source/design/tests/contract were frozen at **`3d94f49` before the full matrix**.
Nine payloads preserve the ten-file source archive, host/local reports,
original exits/resources and the exact mpmath wheel. The initial dependency
preflight exited **1** (`mpmath` absent). A checksum-verified pure-Python wheel
was extracted to an isolated dependency mount; the scientific command/container
then exited **0/0**. The image and native engine were not modified.

## Public replay

Fetch from the checkout, then extract into a fresh directory outside it:

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/physical-pole
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/physical-pole"
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
  uv run --with numpy --with scipy python -m projects.ukrmol_co.physical_pole_study \
  --contract projects/ukrmol_co/physical-pole-contract.json \
  --output "$REPLAY_ROOT/benchmark.json"
```

Expected exit **0**, with `physical_shell_benchmark_verdict=true` and all
nineteen checks true. The report retains all 3,852 real-axis samples, all 81
root statuses/initial seeds, precision checks, convergence sequences and negative
controls. Pole seeds are fixed `nπ−0.3−0.2i`; the shooting solve receives no
Lambert-W answer. No reseeding or relaxed threshold is used.

The pinned execution image is
`sha256:48b9af05d30ce3c289e519433d9fc72a399b8caff6c109ea6c361e8a6b3599b6`,
with reserved physical cores 0–3, one BLAS thread, 4-GiB memory limit and a
300-second timeout. Sadaharu measured **6.722 s wall / 0.051 GiB sampled kernel
high-water**, 32 samples at 0.2-second intervals; the high-water observation is
not a post-exit peak. All persistent sources, dependency bytes, logs and outputs
remain under `prepared/physical-pole-source-20261017/`.

This qualifies the defined scalar-shell boundary convention only. The coupled
polar tail, native Hartree/Rydberg/amplitude normalization and multichannel
threshold sheets remain open. Earlier five-A2 interface rejection, extraction
and electronic-model rejections are preserved; native pole searches and
production remain gated.
