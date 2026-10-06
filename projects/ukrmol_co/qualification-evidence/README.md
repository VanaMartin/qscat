# CO electronic-qualification evidence

This supplement preserves **eighteen completed attempts in eight batches**:
three tightened DZ CAS(10,11) QC targets, seven passing correlated-neutral
pilots, the initial failed Psi4 Python import and one tightened independent
UKRmol target import, four rejected Davidson controls and two passing small-model
SLEPc controls. It also includes all 48
fixed-orbital CI coverage probes, final checkpoint comparisons and the
CAS(10,12) installed-library workspace audit and guarded utility checks,
including thirteen additional valid 16-/32-rank target-sector queries.

All three orbital starts recover the same fixed-ensemble solution: maximum
root/dipole differences are 7.69e-9 Hartree / 5.13e-8 a.u., with minimum final
active-subspace overlap singular value 0.999999999999924. Forty-six CI probes
converge completely; the eighth singlet B1/B2 root fails at trial-space size
40 and converges at 80/160. Every ensemble root converges in every probe.
The original nonzero diagnostic exit and failed controls remain preserved.

The neutral pilot compares conventional frozen-core CCSD(T) at R=1.9/2.1323/2.5
bohr in aug-TZ/QZ bases. Its dipoles are CCSD lambda-density values. The
independent cc-pVDZ PySCF/Psi4 energy control agrees within 8.3e-11 Hartree;
the sentinel relative curves still change by 76.04/68.47 meV across bases.
These are qualification pilots, not a converged electronic model or neutral
curve. The tightened import passes all 40 averaged roots within 5.33e-10
Hartree and the independent dipole within 7.95e-11 a.u.; it took 97.75 minutes
with an 18.42-GiB kernel memory peak. Live scattering and QC refinement batches
are excluded. The prior eleven-/twelve-attempt supplements remain in
`manifest.json`'s `prior_snapshots` with its immutable URL and digest.

The Davidson controls request 5/8/16/32 roots and all fail required root imports
despite engine convergence. Independent diagonalization of the stored small
singlet-A1 Hamiltonian agrees with QC within 3.69e-13 Hartree. The SLEPc-enabled
source build passes 12/12 upstream checks; five-root and eight-root/tighter
small-model targets pass all 40 required roots within 4.98e-10 Hartree and
dipoles within 5.62e-11 a.u. These controls use about 0.12 GiB in 49.53/46.91
seconds. The larger-space SLEPc differential/import gates are still pending.
The preceding 55-attempt target/scattering snapshot remains available in
[`../sa-evidence/`](../sa-evidence/README.md).

## Fetch and reconstruct

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/qualification-evidence.tar.gz \
  -C "$EVIDENCE_DIR"
EVIDENCE_ROOT="$EVIDENCE_DIR/state-averaged-evidence"
uv run python -m projects.ukrmol_co.collect "$EVIDENCE_ROOT" \
  --provenance "$EVIDENCE_ROOT/qualification-provenance.json" \
  --output "$EVIDENCE_DIR/reconstructed.json"
```

Extract outside the checkout: historical source snapshots contain Python test
modules. The reconstructed JSON matches `qualification-results.json` in the
archive. `neutral-results.json` retains the eight-attempt neutral-only subset.
The file index records every payload's size and SHA256. Publication verification
checked 1786 payloads, 127 batch-source hashes, 17 embedded-image source hashes,
raw reanalysis of all thirteen successful records, all 48 CI spectra/spins and
exact reconstruction of both aggregates. Repackaging is byte-identical.
All four Davidson rejections are reproduced; the stored Hamiltonian is
independently reconstructed and its eigenpair residuals checked.

The archive retains raw logs, configs, resource/stage records, source snapshots,
engine build/test provenance, upstream licences, RHF/CCSD/CASSCF checkpoints,
projected initial checkpoints and all workspace query commands. Valid eight-/
sixteen-rank triplet queries give array floors of 148.92/148.99 GiB; the
undersized four-rank query is explicitly rejected by the shipped utility.
The extended all-sector audit reaches 149.78/149.91 GiB on 16/32 ranks;
these are array floors before engine/library overhead.
Large UKRmol integral/channel/K-matrix artifacts remain on the host.

The archived `verify-qualification-evidence.py` reproduces the checks above:

```bash
PYTHONPATH=. uv run python "$EVIDENCE_ROOT/verify-qualification-evidence.py" \
  "$EVIDENCE_ROOT" \
  --original-archive projects/ukrmol_co/qualification-evidence/qualification-evidence.tar.gz \
  --repacked "$EVIDENCE_DIR/repacked.tar.gz"
```

Exact reanalysis/repackaging uses the analysis source commit recorded in the
manifest. Later analyzer versions deliberately record their own source hash.
To replay the numerical recipes, copy archived `runs/` to a new host root
mounted at `/work` so retained checkpoint arguments resolve, and choose new
output/batch names. See [`../CONTINUATION.md`](../CONTINUATION.md) for the live
jobs and remaining gates, and
[the method contract](../../../docs/physics/co-electronic-qualification.md).
