# Native K-attribution timeout evidence

The first fixed two-identity attribution attempt is **unset**, runner/container
exits **1/1**. Both native32 slices timed out at180 s and were killed with native
exit−9. No native K file, native64 reference or independent control was produced.
The [result and source-cost diagnosis](../../NATIVE_K_ATTRIBUTION.md) explain why
the native setup forecast was wrong and why no physics decision follows.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/native-k-attribution-timeout
```

[artifacts.json](artifacts.json) pins twelve public payloads. The manifest records
producer commit`94fccee7b79b1f2c125accf19821fd9ffa3657e0`. Archives include original
source and owner, timed-out partial logs/decks, resource samples, eight inputs,
eight saved references, two checksum-pinned engine sources and a separate verifier.
`results.tar.gz` omits only unchanged copied fort.10/fort.21, independently checked
against `inputs.tar.gz`; all other original files and empty logs are retained.

From a fresh evidence directory, extract `results.tar.gz` into`collected/`,
`inputs.tar.gz` into`raw-preserved/`, `prior.tar.gz` into`prior-preserved/` and
`source.tar.gz` into`source/`. Then run with NumPy/SciPy available:

```bash
PYTHONPATH="$PWD/source" python verify-saved.py "$PWD"
```

The saved-decision verifier independently checks deck reversibility, digests,
timeout flags, native signal exits, missing K files, lack of64-subrange folders
and suppression of independent work. Its exit0 means the saved **unset** decision
is reproduced; it does not convert scientific exits1/1 into success. The original
six-case study stays false, exits1/1, with462 energies uncomputed.

[Fresh public verification](../../verification-native-k-timeout-public-20261019.json)
checks every payload, nineteen source files and two engine sources and runs
all29 focused tests from the downloaded source. No numerical job was reexecuted.
