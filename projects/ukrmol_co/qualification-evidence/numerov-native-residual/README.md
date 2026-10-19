# Native residual primitive evidence

[Qualification](../../NUMEROV_NATIVE_RESIDUAL.md):69/69 checks, build/runtime
exits0/0,13 cases,13 paired80/110-digit references/220 scalar entries. Captured
extended residual matrices agree with NumPy exactly; end-to-end speedup2.46–2.49×.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/numerov-native-residual
```

The [index](artifacts.json) pins original source at`400e7cb`, source manifest,
31-file build/runtime collection, Linux abi3 native wheel, six original operand
captures and mpmath wheel. Verification support at`1ad062e` adds a saved-metric
audit and complete primitive replay without replacing numerical source.

Extract`source.tar.gz` and`verification-support.tar.gz` into`source/`,
`results.tar.gz` into`collected/`, and`operands.tar.gz` into`operands/`.
Extract the mpmath wheel into`dependencies/`. From`source/`:

```sh
PYTHONPATH=../dependencies:. python -m projects.ukrmol_co.numerov_native_residual_replay --report ../collected/result/report.json --output ../verified-decisions.json
```

This independently recalculates metric-derived decisions and timing arithmetic;
input/ownership/signal observations remain saved observations until reexecution.
On the qualified x86-64 host, extract`native-wheel.whl` into`extension/`, add
that directory to`PYTHONPATH` and run:

```sh
PYTHONPATH=../extension:../dependencies:. python -m projects.ukrmol_co.numerov_native_residual_replay --report ../collected/result/report.json --inputs ../operands/result --output ../fresh-replay
```

The output directory must be absent. This reruns all13 cases,13 oracle pairs,
23 full output numeric hashes and69 decisions. Numerical rows must reproduce
exactly; new timings/hypotheses remain observations. The pinned wheel is reused;
a fresh native binary rebuild is outside this numerical replay scope. Build
commands, compiler versions, locks and original Criterion results are retained
for source rebuilds. Python3.11 in the pinned external tools image is explicitly
supported by the published abi3 wheel; ordinary workspace execution uses>=3.12.

The first public owner exits1 in artifact-layout preflight with scientific
verdict unset, before a container is submitted. Its top-level binary lookup
missed maturin's nested package. The producer's`extension_sha256` field also
names the package initializer; that original label is retained. The separately
frozen [supplement](../../numerov-native-residual-replay-preflight-supplement.json)
pins the exact initializer and actual native binary paths/hashes. The successor
replays all numerical rows/hashes/decisions exactly, exits0/0. All17 final payloads
verify; failed owner/log and successful replay evidence are both included.
