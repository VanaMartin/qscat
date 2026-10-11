# Completed CAS11 16K and TZ two-frame evidence

Public download and independent macOS reconstruction pass **352 payloads**:

- CAS11 time-budget successor: both **16384-root** sectors pass raw eigenpair,
  continuum/boundary,full-grid phase and identical fixed-window fit checks.
  Maximum phase error is **1.0e-8 rad modulo π**. Native cost remains
  **19.653 hours / 37.068 GiB**. The separately published
  [2048/4096/8192 controls](../cas11-sparse-refinement/README.md) remain passed;
  **2048 is the smallest qualified CAS11 count**. No further CAS11 root expansion
  is needed for the current discrepancy.
- TZ replay: all **128 widest-space states / 48 probe records** pass independent
  physical/penalized/spin-action,spin,dipole,space and orbital-subspace checks.
  The triplet-A1 root3 frame difference remains **−3.02488e-7 Hartree**,above the
  unchanged **1e-7** seed gate. The original TZ import remains rejected; the
  separately qualified reoptimized target is a distinct contract.

Original24-GiB OOM,six-hour timeout and TZ final-report serialization failure
remain in the archive. This reconstruction certifies those numerical records;
CAS12 scattering,model convergence and state/pole continuity remain separate.

## Reproduce from public downloads

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/cas11-16k-tz-continuation
# Choose a fresh directory outside the checkout; do not collect archive tests.
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT"
tar -xzf projects/ukrmol_co/qualification-evidence/cas11-16k-tz-continuation/continuation-evidence.tar.gz \
  -C "$REPLAY_ROOT"
PYTHONPATH="$REPLAY_ROOT/state-averaged-evidence" \
  uv run --with pyscf==2.11.0 --with h5py==3.15.1 \
  python "$REPLAY_ROOT/state-averaged-evidence/verify-continuation-evidence.py" \
  "$REPLAY_ROOT/state-averaged-evidence" "$REPLAY_ROOT/replay"
```

The archive SHA256 is
`82d8562585c6133eca411a5693d5ac824a9dd3b3baa9ccfe5a30a1c4c885c847`,
**1291953735 bytes**. It includes native raw records,dense oracle,telemetry,
checkpoints/CI vectors,resources,original exits and frozen reconstruction source.
The verifier reads and hashes all indexed payloads before recomputing the physics.

Deterministic repackaging is byte-identical: sort the indexed files plus
`file-index.json` as **pathlib paths**,write only files under
`state-averaged-evidence/` in PAX tar,uid/gid0,empty names,mode0644,mtime0,and gzip
with empty filename/mtime0. Plain string sorting puts `cas11-dense-oracle/` ahead
of `cas11/`; that first local attempt preserved identical payloads but changed
member order. The corrected fresh path-order attempt matches the published bytes.
The initial upload exceeded its120-second client bound; a longer bounded upload
published the same unchanged archive successfully. Both packaging attempts and
the first upload record are retained on the working hosts.
