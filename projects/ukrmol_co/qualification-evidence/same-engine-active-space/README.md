# Same-engine CAS11→CAS12 supplementary evidence

**The active-space model comparison still rejects phase and position budgets**
using the fresh equilibrium five-A2 CAS11 baseline, the exact same 242-artifact
audited scattering engine, the same image and eight MPI ranks. All four fit
orders in 1.8–3.3 requested-eV reject position insensitivity; widths pass.

| CAS12 roots | Maximum phase change | Center shift | Width change |
|---|---:|---:|---:|
| 2048 | 0.1212798 rad | 63.652–73.447 meV lower | 16.684–26.661 meV smaller |
| 4096 | 0.1212919 rad | 63.653–73.447 meV lower | 16.718–26.677 meV smaller |

The fresh CAS11 phases exactly match the old selected CAS11 phases in both
sectors. This removes the earlier scattering-engine/rank confound; target
diagonalizers/imports, optimized orbital frames and native target-file
representations remain distinct. The result is model sensitivity rather than
a pure active-orbital-count effect or absolute omission bound.

## Original strict preflight and supplementary audit

The [strict preflight](../../verification-same-engine-strict-preflight-20261017.json)
exits 1 with an unset observable verdict because two `rsolve.inp` lines differ
in spacing/comments. It remains failed. The
[supplementary contract](../../same-engine-active-space-contract.json) permits
only the two exact, declared full-line substitutions after checking that their
uncommented tokens are identical. The entire resulting files must then be
byte-identical. All `eigenp`/`reson` bytes remain exact. An injected change in
incident energy is rejected. No numerical budget changes.

The [reader](../../same_engine_active_space.py) reconstructs 197 existing raw
files, six new owner-provenance files, the original numerical prerequisites and
both model comparisons. No native solve is repeated. Host replay costs 2.268
seconds; host/local fitted-parameter drift is 2.717e-10 eV. Shared fitting/export
routines do not provide an independent scattering engine.

## Reproduce

Fetch this bundle and the two existing raw references; raw archives are reused
rather than republished:

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/same-engine-active-space
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/equilibrium-a2-scattering
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/saved-model-controls
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence"
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/source" "$REPLAY_ROOT/provenance" "$REPLAY_ROOT/cas11" "$REPLAY_ROOT/saved"
tar -xzf "$EVIDENCE/same-engine-active-space/reader-source.tar.gz" -C "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/same-engine-active-space/provenance.tar.gz" -C "$REPLAY_ROOT/provenance"
tar -xzf "$EVIDENCE/equilibrium-a2-scattering/raw-inputs.tar.gz" -C "$REPLAY_ROOT/cas11"
tar -xzf "$EVIDENCE/saved-model-controls/raw-inputs.tar.gz" -C "$REPLAY_ROOT/saved"
cd "$REPLAY_ROOT/source"
PYTHONPATH=. uv run --with numpy --with scipy python -m projects.ukrmol_co.same_engine_active_space \
  --cas11 "$REPLAY_ROOT/cas11" --saved "$REPLAY_ROOT/saved" \
  --provenance "$REPLAY_ROOT/provenance" \
  --contract projects/ukrmol_co/same-engine-active-space-contract.json \
  --output "$REPLAY_ROOT/reconstructed.json"
```

Verify the 61 frozen source hashes in `reader-source-manifest.json` before
running. Raw and native-owner manifests are checked by the reader. Source is
frozen at `ff04b2c`; original native scientific runtimes retain their own commits.
Original extraction rejections, unset global gates and the failed strict
preflight remain preserved. Production is gated.
