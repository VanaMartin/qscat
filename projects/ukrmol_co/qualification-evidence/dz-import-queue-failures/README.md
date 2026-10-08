# Stretched CAS12 DZ import and preserved queue failures

The completed stretched CAS12 DZ target import passes **all 64 native roots**,
forty-state dipoles, saved-table rounding and core/active orbital-subspace checks.
Runtime is **21336.38 seconds / 40.371 GiB**. Maximum independent root error is
**5.460e-8 Hartree**; seed-to-reoptimized ensemble-root difference is
**4.262e-8 Hartree**, native dipole difference **8.919e-8 a.u.** Minimum active
subspace singular value is **0.99999999999421**. All **96** independently covered
fixed-orbital probes on both stretched DZ repair checkpoints remain passing.

The companion also preserves two subsequent failures and their dependent stops:

- **CAS11 16384-root, 40-GiB attempt:** B1 SCATCI is interrupted by its declared
  **21600-second timeout**. Total wall time is **21601.48 seconds**, kernel peak
  **31.994 GiB**, controller exit **one**. No completed eigenpairs, boundary
  export, B2 solve or observable verdict exist. The logs enter the Krylov-Schur
  solver; this record establishes neither convergence nor nonconvergence.
- **Original TZ import owner:** batch snapshot creation fails because
  `projects/__init__.py` is missing. No TZ native engine or run is launched.
  The supported-space aug-TZ and compressed owners stop at their predecessor
  guards. Every original exit, traceback and stopped queue remains exact.

The [fresh time-budget contract](../../sparse-scattering-cas11-time-budget-successor-contract.json)
retains the scientific inputs and 40-GiB cap with a 24-hour per-stage limit.
The [package-complete queue contract](../../cas12-native-import-packaging-successor-contract.json)
adds the empty package marker to fresh frozen sources, retains identical recipe
arguments, uses fresh `-package2` run names and preserves TZ → supported-space
aug-TZ → compressed priority. The original aug-TZ coverage rejection remains
explicit. The archive captures their initial declarations and waiting/started
owner states; their final verdicts belong to subsequent evidence.

## Independent public replay

Extract outside the checkout to avoid collecting frozen tests:

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/dz-import-queue-failures
REPLAY=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/dz-import-queue-failures/dz-import-queue-failures.tar.gz -C "$REPLAY"
EVIDENCE="$REPLAY/state-averaged-evidence"
PYTHONPATH="$EVIDENCE/prepared/electronic-successors-source" \
  uv run --with pyscf==2.11.0 --with h5py==3.15.1 python \
  "$EVIDENCE/verify-dz-import-failures.py" "$EVIDENCE"
python "$EVIDENCE/pack-dz-import-failures.py" "$EVIDENCE" "$REPLAY/repacked.tar.gz"
```

Reconstruction verifies **3040 payloads / 225208 resource profiles**, native
final-set and decimal analytic/corruption controls, coverage spin/residual/root/
space gates, source/checkpoint hashes, exact timeout and pre-engine failure
records, and the measured resource peaks. Repackaging is byte-identical; SHA256:
`c920f0836aa93ce195ab76a072b7676729110ec7d90c3d1ece136e83032ecf68`.
Large integrals and the incomplete timeout Hamiltonian remain on Sadaharu,
with their original hashes retained. Full-model convergence, extraction,
threshold identity and CAS12 scattering remain separate gates.
