# Saved native scattering matrices and channel interfaces

**The original aggregate reading rejects: six nine-A2 cases pass, six five-A2
cases stay unset on incompatible target labels. Original command/container
exits are 1/1.** This read-only audit uses the saved outputs; no native solve,
outer propagation, fitted width or complex pole search is launched.

The [design](../../../../docs/physics/co-saved-scattering-matrices.md),
[contract](../../saved-matrices-contract.json), reader and three analytic/format
tests were frozen at `343cc1b` before the complete reading. All three anchors,
both Pi sectors and five/nine-A2 were retained. The archive contains **123
checksum-bound input files**: K matrices, channel/boundary data, phase outputs,
decks, owner/model records and external engine source/provenance. The separate
writer archive contains the source and twelve scattering decks used for the
label diagnosis. Ten content-addressed payloads preserve host/local readings,
resources and the original failure. External engine sources retain their own
licence; `engine/COPYING` travels with the raw evidence.

See [results](../../SAVED_MATRICES.md),
[quoted verification](../../verification-saved-matrices-20261017.json) and
[public replay](../../verification-saved-matrices-public-20261017.json).

## Replay the frozen reading

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/saved-matrices
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/saved-matrices"
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/raw" "$REPLAY_ROOT/source" "$REPLAY_ROOT/writer"
tar -xzf "$EVIDENCE/raw-inputs.tar.gz" -C "$REPLAY_ROOT/raw"
tar -xzf "$EVIDENCE/source.tar.gz" -C "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/writer-evidence.tar.gz" -C "$REPLAY_ROOT/writer"
cd "$REPLAY_ROOT/source"
python3 - "$EVIDENCE/source-manifest.json" <<'PY'
import hashlib, json, pathlib, sys
for name, digest in json.loads(pathlib.Path(sys.argv[1]).read_text()).items():
    assert hashlib.sha256(pathlib.Path(name).read_bytes()).hexdigest() == digest, name
PY
PYTHONPATH=. uv run --with numpy --with scipy python -m projects.ukrmol_co.saved_matrices \
  --directory "$REPLAY_ROOT/raw" \
  --contract projects/ukrmol_co/saved-matrices-contract.json \
  --output "$REPLAY_ROOT/matrices.json"
```

Expected **exit 1**, after writing all twelve results. Six nine-A2 cases recover
594 six-channel S matrices with phase agreement ≤4.94e-8 rad modulo π. The
five-A2 rows stay unset with `ValueError: Invalid channel label`; the aggregate
is false. Host/local/public case verdicts and failures agree. This reconstructs
native K outputs, not independent outer propagation or pole identity.

## Recheck the supplementary label diagnosis

This diagnostic uses **sorted original property IDs from the actual frozen
decks/pilots**, not energy-nearest matching. It does not replace the frozen
reader's failed interface gate. From the extracted source directory:

```bash
PYTHONPATH=. uv run --with numpy --with scipy python - "$REPLAY_ROOT" <<'PY'
import json, pathlib, re, sys
import numpy as np
from projects.ukrmol_co.saved_matrices import records
root = pathlib.Path(sys.argv[1])
for anchor in ("compressed", "equilibrium", "stretched"):
    for model in ("five", "nine"):
        mapping = json.loads((root / "raw" / anchor / model / "pilot.json").read_text())["target_map"]
        rank = {old: i for i, old in enumerate(sorted(mapping))}
        for sector in ("B1", "B2"):
            deck = (root / "writer" / anchor / model / sector / "scatci.inp").read_text()
            assert mapping == [int(v) for v in re.search(r"\bidtarg\s*=\s*([\d,]+)", deck, re.I)[1].strip(",").split(",")]
            raw = records(root / "raw" / anchor / model / sector / "fort.10")
            n = np.frombuffer(raw[2], "<i8")[0]
            energies = np.array([np.frombuffer(v[32:], "<f8")[0] for v in raw[4:4+n]])
            labels = np.array([np.frombuffer(v[:32], "<i8") for v in raw[4+n:]])
            thresholds = np.array([np.frombuffer(v[32:], "<f8")[0] for v in raw[4+n:]])
            assert set(labels[:, 1]) == set(mapping)
            expected = 2 * (energies[[rank[v] for v in labels[:, 1]]] - energies[0])
            np.testing.assert_allclose(thresholds, expected, atol=1e-10, rtol=0)
            print(anchor, model, sector, "out-of-table channel labels:", np.count_nonzero(labels[:, 1] > n))
PY
```

The writer `get_channel_info` stores original property IDs in `ICHL`; its
`write_channel_info` writes a densely renumbered selected target table but
unchanged `ICHL`. Explicit old-ID lookup reproduces every recorded threshold
exactly on these files. Reusing labels directly as dense indices is unsafe.
No engine correction or replacement decoding contract is released here.
