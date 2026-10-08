# Compressed CAS12 fresh singlet-driver QC: mixed verdict

The two fresh orbital optimizations use the fixed-orbital-qualified alternative
singlet-A1 driver and unchanged forty-state ensemble, basis, geometry and
tolerances. Both orbital optimizers converge, but their CI verdicts differ:

| Trial space | Original runner exit | QC verdict | Wall / peak memory |
|---|---:|---|---|
| 80 | 1 | Singlet-A1 CI flag false; rejected | 1660.89 s / 0.7860 GiB |
| 160 | 0 | All optimizer/fresh-CI/spin/Pi/MO checks pass | 1724.00 s / 0.9610 GiB |

The space160 analyzer also exits zero. The owner retains combined exit one;
the failed space80 run has no exported target or successful analyzer record.
Final gradients are **5.947e-8 / 5.932e-8**, and saved roots agree between the
two checkpoints within **9.813e-11 Hartree**. This numerical agreement cannot
override the failed space80 CI flag. Both objectives are
**−112.23910084461343 Hartree** at printed precision. Space160's fresh roots
agree within **1.279e-13 Hartree** and its ground dipole is **0.1773807 a.u.**

The [finite coverage contract](../../cas12-compressed-singlet-qc-coverage-contract.json)
checks the unchanged new checkpoints independently at five/eight roots and
spaces80/160/240, with physical/spin-penalized residuals, independent singlet
Hamiltonian action and root/orbital continuity. QC self-consistency is not
independent coverage or native-import qualification. Original rejected seeds
and the rejected new space80 optimization remain preserved.

## Public replay

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/compressed-singlet-qc
tar -xzf projects/ukrmol_co/qualification-evidence/compressed-singlet-qc/compressed-singlet-qc.tar.gz \
  -C projects/ukrmol_co/qualification-evidence/compressed-singlet-qc
PYTHONPATH=projects/ukrmol_co/qualification-evidence/compressed-singlet-qc/state-averaged-evidence \
  uv run python projects/ukrmol_co/qualification-evidence/compressed-singlet-qc/state-averaged-evidence/verify-singlet-qc.py \
  projects/ukrmol_co/qualification-evidence/compressed-singlet-qc/state-averaged-evidence
```

The archive independently reconstructs **228 payloads / 16864 profiles**, both
original exits, failed flag/error, successful fresh-CI/spin/Pi/MO/Molden checks,
initial/final checkpoint hashes, immutable source and pinned-image driver controls.
Public fetch, reconstruction and byte-identical repackaging pass.
