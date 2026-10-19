# Fixed-grid analytic precision attribution: public evidence

The [precision attribution](../../NUMEROV_PRECISION_ATTRIBUTION.md) passes
46/46 diagnostic checks/exits0/0, retaining the original interior prerequisite's
false/exits1/1. Ten new analytic solves and twelve reused rows identify dominant
floating-point drift; the all-endpoint tenfold basis-improvement hypothesis is
false. No CO/native/inner solve or new unique energy is produced.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/numerov-precision
```

[artifacts.json](artifacts.json) pins public payloads; [manifest.json](manifest.json)
records producer freeze`aa917aa`, verification support and hypothesis decisions.
Extract `source.tar.gz` into`source/` and `results.tar.gz` into`collected/`.
Verify all51 source hashes against `source-manifest.json`. Extract
`verification-support.tar.gz` over`source/` and verify its two module hashes from
the manifest; no frozen module is replaced. Extract the pinned mpmath1.3.0 wheel
into`dependencies/`. Keep `prior-report.json` at the working root.

With NumPy/SciPy available:

```bash
PYTHONPATH="$PWD/dependencies:$PWD/source" python -m projects.ukrmol_co.numerov_precision_replay \
  "$PWD" --output fresh-saved-replay.json
```

The verifier independently expresses scattering/error/check decisions, uses the
tested matrix/scalar oracles, and checks all twelve reused rows, ten new matrices,
signed decompositions and both hypotheses without repropagation.

The initial fresh-public test collection stops with **exit2**, missing the
existing`numerov_corner_replay` test dependency. Its original support archive and
test log are preserved publicly. A dependency-only support supplement adds that
module; numerical producer sources remain unchanged. Fresh downloads then pass
**68 tests** and saved decisions. The separately pinned complete numerical replay
contract covers all ten new solves and twelve reused-row oracle decompositions,
four reserved physical cores0–3/4 GiB/600-s cap.

## Completed numerical replay and preservation

All **ten new analytic matrices** reproduce K components exactly; all twelve
reused-row oracle decompositions,46 check decisions and hypothesis values remain
identical. The twelve archived rows are reused, not repropagated. Verification
exits0/0 preserve both diagnostic acceptance and original interior false/exits1/1.
All **19 final payloads** fresh-verify. Extract `public-replay.tar.gz` separately;
its8-file manifest pins the owner, complete report, replay comparisons, resources
and exits. Standalone owner/execution/replay payloads are byte-identical copies.

The original replay launcher fails before any container/solve with local-variable
shadowing (owner exit1). `public-owner-initial-preflight.py`,
`initial-owner-preflight-failure.log` and `public-owner-supplement.json` preserve
it. A variable-only correction in a separate persistent directory runs the
unchanged pinned worker and inputs. Neither verification supplement changes any
scientific setting or upgrades the original interior rejection.

Measured60.548024 s /60.079134 sampled CPU s /45,404,160-byte kernel high-water,
301 samples at0.2 s, four cores0–3/4 GiB/600-s cap.
[Public verification](../../verification-numerov-precision-public-20261019.json)
records fresh downloads,68 tests and full fixed-scope numerical replay.
