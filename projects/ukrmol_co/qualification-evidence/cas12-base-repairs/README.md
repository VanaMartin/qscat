# Original CAS(10,12) DZ numerical repairs

Supervisor 918331 completes the finite trial-space-80/160 retry pair at
R=1.9/2.5 bohr, using the unchanged CAS(10,12) forty-component ensemble,
two frozen core orbitals, tight tolerances and original checkpoint seeds.

| Geometry / space | QC verdict | Wall / seconds | Peak / GiB |
|---|---|---:|---:|
| 1.9 / 80 | Reject: singlet-A1 CI flag | 1527.87 | 0.7892 |
| 1.9 / 160 | Reject: singlet-A1 CI flag | 1481.43 | 0.9646 |
| 2.5 / 80 | Pass QC; reject original coverage | 1427.43 | 0.7858 |
| 2.5 / 160 | Pass QC; reject original coverage | 1350.82 | 0.9684 |

Both compressed runs converge the orbitals but retain a failed singlet-A1 CI
sector flag. They stop before fresh-root/export qualification. Both stretched
runs pass QC, agree in ground energy within **1.973e-11 Hartree**, and in
ensemble objective within **2.321e-12 eV**.

The original fixed-orbital coverage scan runs **48 probes per stretched
checkpoint**, varying five/eight roots and trial spaces 40/80/160 with
**200 cycles**. All ensemble energies agree within **2.71e-13 Hartree**, but
three fifth roots fail at space 40: singlet A2 and triplet B1/B2. Extra singlet
A1/B1/B2 eighth roots also have failed flags. Energy agreement cannot override
those flags, so neither pair releases a larger-basis successor.

The subsequently launched 600-cycle fixed-orbital physical/penalized-residual
scan is a separately owned control. This companion preserves the original
200-cycle contract and rejection exactly.

## Fetch and independently reconstruct

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/cas12-base-repairs
mkdir -p /tmp/co-cas12-repairs
tar -xzf projects/ukrmol_co/qualification-evidence/cas12-base-repairs/cas12-base-repairs.tar.gz \
  -C /tmp/co-cas12-repairs
uv run python /tmp/co-cas12-repairs/state-averaged-evidence/verify-cas12-repairs.py \
  /tmp/co-cas12-repairs/state-averaged-evidence
```

The verifier checks **437 payloads**, **28844 resource samples**, frozen owner
and executable-source hashes, original seed/checkpoint hashes, original exits,
QC/gradient/spin/orbital/Molden consistency, all 96 fixed-orbital probes and the
two stretched-target differentials. Raw logs, inputs, profiles, diagnostics,
checkpoints and the original failure records travel with the archive.
Prior companion archives remain separate with their original source-snapshot
hashes retained. Sorted normalized tar/gzip repackaging is byte-identical.
