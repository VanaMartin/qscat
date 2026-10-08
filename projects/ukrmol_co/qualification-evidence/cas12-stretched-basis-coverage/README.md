# Stretched CAS12 TZ/aug-TZ fixed-orbital coverage

This companion independently reconstructs **354 payloads / 8603 profiles**
from the complete 96-probe scan at the two unchanged, qualified QC checkpoints.
Each spin/irrep uses 5/8 roots, spaces 40/80/160 and 600 cycles, with the existing
energy, physical/spin-penalized residual, spin, orthogonality and checkpoint gates.

**TZ passes all 48 probes. Aug-TZ passes 47 and rejects one:** triplet A2,
eight roots, trial space 40. Root eight's CI flag is false, with physical and
spin-penalized residuals **9.41093e-7 / 1.38254e-6 Hartree** against **1e-9
Hartree**. All forty ensemble roots pass every probe, and the independent
energies remain close; these facts do not override the rejected extra-root flag.

The scan itself completes with exit zero in **1726.94 seconds / 0.7026 GiB**.
The gate-checking verifier exits **one**, and the owner preserves controller exit
**one** before launching the coupled outer-window stage. Original source, logs,
checkpoints, flags/residuals and all profiles travel with the archive.

The [fresh supported-space contract](../../cas12-stretched-augtz-supported-spaces-contract.json)
declares a separate 80/160/240 trial-space family at unchanged orbitals, 600
cycles and unchanged tolerances. It cannot erase this original rejection or
replace the existing native queue's rejected aug-TZ proof automatically.

## Public replay

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/cas12-stretched-basis-coverage
tar -xzf projects/ukrmol_co/qualification-evidence/cas12-stretched-basis-coverage/cas12-stretched-basis-coverage.tar.gz \
  -C projects/ukrmol_co/qualification-evidence/cas12-stretched-basis-coverage
PYTHONPATH=projects/ukrmol_co/qualification-evidence/cas12-stretched-basis-coverage/state-averaged-evidence \
  uv run python projects/ukrmol_co/qualification-evidence/cas12-stretched-basis-coverage/state-averaged-evidence/verify-cas12-basis-coverage-bundle.py \
  projects/ukrmol_co/qualification-evidence/cas12-stretched-basis-coverage/state-averaged-evidence
```

Public fetch, independent reconstruction and byte-identical repackaging pass.
The replay verifies the recorded mixed verdict; successful reconstruction of
the rejection does not qualify aug-TZ for a native import.
