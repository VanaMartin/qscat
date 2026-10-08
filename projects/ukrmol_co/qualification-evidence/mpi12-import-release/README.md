# Initial twelve-rank target-import release

This is a historical **zero-launched-step declaration**, not a completed
target import or MPI speedup measurement. The [user-directed host policy](../../sadaharu-cpu-allocation-policy.json)
allocates **8 physical cores to large jobs**, **12 to very large jobs**, and
reserves **4 cores (0–3)** for small controls. MPI stages use matching ranks
and one BLAS/OpenMP thread per rank.

The [fresh import contract](../../cas12-native-import-mpi12-successor-contract.json)
administratively supersedes the original waiting package2 owner **1065882**,
after verifying it has launched no steps. Its exact unfinished execution file
and frozen source remain preserved; a successful systemd stop is recorded
without fabricating a scientific controller exit. The active four-rank CAS11
time-budget owner is untouched.

New owner **1072749** waits for that CAS11 owner, then serially imports TZ,
records the original rejected aug-TZ coverage, imports independently supported
aug-TZ, and imports the passing compressed space160 target. Native jobs use
**12 ranks / CPUs4–15**, fresh **`-mpi12-package3`** names, the same **64-GiB
container / 80-GiB available-RAM / 40-GiB disk** gates, and mandatory independent
64-root/dipole/subspace verifiers on CPUs0–3. Only MPI ranks and new run names
change; all scientific arguments, checkpoints, proofs and tolerances remain exact.

Live preflight parses all twelve-rank recipes and intercepts native Docker
launches while exercising real batch snapshot creation. It verifies both
package markers, CPU group4–15, argument parity and **120/123/132 scientific
files**. The companion verifies **847 payloads**, including **832 frozen release
files / 820 unchanged parent file hashes**, exact original-owner supersedure
records, original parent index and the initial current-owner/hourly journal.

The original reporting monitor's startup race is retained: it reads the owner
record before the initial integrity check finishes creating it. A fresh v2
monitor accepts that pending-record state and journals hourly for48hours.
Original release-verifier metadata-parity and nested-index-copy failures are
also retained with their fresh reconstruction repairs. Scientific owner/input
records are unaffected. Current queue progress lives on Sadaharu under
`prepared/mpi12-import-progress-v2-20261009/`; this archive remains the initial
declaration, not a live status feed.

## Public reconstruction

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/mpi12-import-release
REPLAY=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/mpi12-import-release/mpi12-import-release.tar.gz -C "$REPLAY"
EVIDENCE="$REPLAY/state-averaged-evidence"
python "$EVIDENCE/verify-mpi12-release-v2.py" "$EVIDENCE"
python "$EVIDENCE/pack-small-anchor-controls.py" "$EVIDENCE" "$REPLAY/repacked.tar.gz"
```

Public fetch, source/argument/supersedure reconstruction and deterministic
byte-identical repackaging pass. SHA256:
`413e752a05a83ccdfecb400fb6e44c0106ccf0953f44fcf294e9302d7580ff4a`.
Native-import and rank-scaling verdicts remain unset. Every new native import
must pass its unchanged electronic checks before the queue advances.
