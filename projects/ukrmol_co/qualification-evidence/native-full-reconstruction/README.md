# Six-case reconstruction: preserved K-only rejection

**Original verdict false, command/container exits 1/1.** All six native controls
and all 18 independent probes pass. The subsequent remainder stops after 112
passes and two K-only rejections at stretched B1/B2, requested 1.0 eV. Coverage is
132/594 energy identities and 354 independent solves. The other 462 energies
remain uncomputed.

See [results and diagnosis](../../NATIVE_FULL_RECONSTRUCTION.md),
[verification](../../verification-native-full-20261017.json),
[contract](../../native-full-reconstruction-contract.json) and
[design](../../../../docs/physics/co-native-full-reconstruction.md).
The producer is frozen at `fd7c871`, the separate saved-decision verifier at
`64226b9`. Original pilot false/exits 1/1 remain.

Thirteen public payloads preserve the source, verifier, 24 input files, 335 collected
outputs/metadata, original exits, 354 per-solve journal rows and all 1782 native
energy records. `results.tar.gz` contains `result/`, the execution metadata,
owner and sampled resources. Redundant native copies of fort.10/21 are omitted;
the collection verified them against the checksum-pinned originals in
`inputs.tar.gz`. The input archive preserves all six original deck-specific
K units, including B2's fort.922.

The [fresh public verification](../../verification-native-full-public-20261017.json)
checks all 13 payloads and reconstructs every decision and the false aggregate.
The verifier uses a separately expressed Cayley transform and eigenvalue-arctangent
phase sum; it does not rerun native Fortran. Numerical reexecution of the two
rejected independent bundles, four solves, reproduces K and stored S components
exactly and preserves both false decisions. The successful verifier/numerical
replay exits are 0/0; original scientific exits remain 1/1. This reexecutes four
solves, not all 354. Its 190.720-s / 141,156,352-byte recorded kernel high-water
cost is separate from the original stopped study.

`numerical-replay.tar.gz` preserves the four replay rows, per-energy reports,
worker, original fetched pointer and sampled execution metadata. A separate
`test-support.tar.gz` supplies the coupled-tail analytic oracle and predecessor
contract needed by the frozen tests. All 31 focused tests pass. The initial
missing-dependency public test exits (2, then 1) remain in `public-tests.json`;
the original producer source archive and all scientific outputs are unchanged.

## Fetch and reproduce saved decisions

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/native-full-reconstruction
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/native-full-reconstruction"
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/source" "$REPLAY_ROOT/raw" "$REPLAY_ROOT/collected"
tar -xzf "$EVIDENCE/source.tar.gz" -C "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/replay-source.tar.gz" -C "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/test-support.tar.gz" -C "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/inputs.tar.gz" -C "$REPLAY_ROOT/raw"
tar -xzf "$EVIDENCE/results.tar.gz" -C "$REPLAY_ROOT/collected"
cd "$REPLAY_ROOT/source"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m projects.ukrmol_co.native_reconstruction_replay \
  --root "$REPLAY_ROOT/collected" --raw "$REPLAY_ROOT/raw" \
  --output "$REPLAY_ROOT/replayed-decisions.json"
```

Expected verifier exit is 0, with `saved_decision_replay_passed=true` and
`reconstruction_verdict=false`. The successful verification of a rejection is
distinct from the original scientific exits 1/1. `artifacts.json` checks every
download's bytes; the replay also checks collection/input hashes and exact
journal/report parity. Producer and replay source hashes are carried in
`source-manifest.json` and `manifest.json` respectively.

All K/S/phase budgets remain 1e-6. Large tangent coordinates explain the disparate
K/S sensitivity but do not waive the K veto. Full native K reconstruction,
complex poles, electronic/channel convergence and production remain unqualified.
