# Fixed-grid Numerov increment control evidence

The [two-identity CO control](../../NUMEROV_INCREMENT_CONTROL.md) passes
**24/24 checks**, producer exits0/0. The original two-point and six-case studies
remain false/exits1/1. The candidate changes floating-point algebra, retaining
the exact-arithmetic Numerov problem, grid, boundary, tail and unsymmetrized K.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/numerov-increment
```

[artifacts.json](artifacts.json) pins **20 public payloads** and [manifest.json](manifest.json)
records source freeze`a8665f4714f5ae2de3f310ca1f42aecb6dd3ccce`.
`initial-import-failure.tar.gz` preserves the missing-mpmath attempt/exits1/1;
`analytic-results.tar.gz` preserves the passing dependency-supplemented run,
16 checks/exits0/0. Each contains its own original source archive and resource
records. The wheel is independently checksum-pinned in the tracked
[dependency declaration](../../numerov-increment-dependencies.json).

For the independent saved-decision replay, extract:

| Archive | Destination |
|---|---|
| `results.tar.gz` | `collected/` |
| `inputs.tar.gz` | `raw-preserved/` |
| `references.tar.gz` | `prior-preserved/` |
| `saved.tar.gz` | `saved-preserved/` |
| `source.tar.gz` | `source/` |

Then, with NumPy/SciPy available:

```bash
PYTHONPATH="$PWD/source" python verify-saved.py "$PWD"
```

The separate verifier decodes saved native64 directly, checks disk/report rows,
input hashes, complete family coverage, S reconstruction and all24 vetoes. It
does not call the candidate evaluator. The verification exit is separate from
producer scientific exits. The bundle also supplies a fixed public numerical
replay worker for exactly ten solves, no native/inner engine call.

**The pass is scoped:** the finest-step400-bohr/native64 comparison passes,
but the step0.0125/endpoint800 rows differ from native64 by1.59437e-6/1.46473e-6.
Those rows are not all-native-reference-qualified. No favourable-row selection,
absolute error bound, full594-matrix reconstruction or complex-pole claim is made.
The original462 uncomputed energies remain uncomputed.

[Public-byte verification](../../verification-numerov-increment-public-20261019.json)
checks34 source files, saved decisions and49 tests from fresh downloads.
The source archive includes the analytic oracle and candidate tests. Extract
the checksum-verified mpmath wheel to a separate dependency directory on
PYTHONPATH to run those analytic tests in an environment lacking mpmath.

The public numerical replay passes **exits0/0**, reproducing every K and
stored S component in all ten controls exactly, with identical24 passing
checks. Its17-file archive and separate owner/report/execution resources are
published. Runtime390.215 s /1332.684 sampled CPU-seconds /239,140,864-byte
recorded kernel high-water, four physical cores/4 GiB. All20 final payloads
are freshly downloaded after supplementation. `configuration-diagnostics.json`
preserves every grid/native comparison, including the endpoint800 discrepancy;
it is explicitly separate from the frozen gate set.
