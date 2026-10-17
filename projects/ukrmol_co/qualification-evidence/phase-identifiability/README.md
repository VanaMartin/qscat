# Local analytic width-identifiability evidence

**Historical producer reading: 98 wrong-width phase matches, but its fully
unused-data declaration fails.** Window endpoints coincide with development
points in 14/45 grids. Original command exit 0 and producer qualification remain
preserved; the fully unused protocol is not independently qualified. The
[separately frozen endpoint repair](../phase-identifiability-endpoint/README.md)
retains **94** counterexamples and passes all 45 mandatory true-pole controls.

The original producer reports 45 passing true-pole nonpolynomial controls.
The 405-cell matrix, six source files, exact audit grids and resource records
are published in seven checksum-addressed payloads. Source is frozen at
`c273a71cfe61648637d9fced44f027b470c2886c`, before the full matrix run.

These scalar S factors have analytically known poles and entire nonzero
backgrounds. They are local analytic/unitary counterexamples, **not qualified
causal finite-range CO scattering models**. They show the limit of what those
local checks alone establish. No native data is read, no estimator is selected,
and no pole curve is released. Zero counterexamples in a tested family does not
prove uniqueness. See the [design](../../../../docs/physics/co-phase-identifiability.md),
[contract](../../phase-identifiability-contract.json),
[results](../../PHASE_IDENTIFIABILITY.md),
[quoted verification](../../verification-phase-identifiability-20261017.json) and
[public replay](../../verification-phase-identifiability-public-20261017.json).

## Reproduce

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/phase-identifiability
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/phase-identifiability"
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

Expected original exit: **0**. Its endpoint data-use defect means that this
producer exit does not independently qualify the fully unused protocol or
identify widths. The report retains every degree/window/grid/width cell,
all development/reserved/dense/endpoint errors, coefficients, matrix ranks,
conditioning and maximum-error energies. The true S generator does not call the
candidate predictor. Host/local cell verdicts agree for all 405 cells; largest
phase-error metric drift is **7.72e-15 rad**.

Sadaharu used reserved physical cores 0–3 and one BLAS thread in the pinned tools
image with a 4-GiB memory limit. Original command/container exits **0/0** and
the persistent log are preserved. Wall time is **0.965 s**; the kernel high-water
value sampled during the run is **0.061 GiB**, not a post-exit peak measurement.
The source/output directory is `prepared/phase-identifiability-source-20261017/`.
