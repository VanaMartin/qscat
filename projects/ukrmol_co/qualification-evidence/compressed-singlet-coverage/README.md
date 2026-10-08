# Compressed CAS12 independent singlet-driver coverage

All **96 fixed-orbital probes pass** on the two new unchanged checkpoints:
five/eight roots, spaces80/160/240, both spins and four irreps, with a 600-cycle
diagnostic cap. Runtime is **1420.60 seconds / 0.8013 GiB**. Independent flags,
physical/spin-penalized/spin-penalty residuals, spins, CI/MO orthogonality,
ensemble-root, space/root-count/Pi and singlet full-Hamiltonian-action checks pass.

Maximum physical/penalized residuals are **8.889e-10 / 9.968e-10 Hartree**;
independent singlet full-action difference/residual are **7.081e-10 / 5.690e-10
Hartree**. Maximum ensemble-root error is **3.837e-13 Hartree**. Space/root-count
differences stay below **2.558e-13 Hartree**. Between-checkpoint roots differ by
at most **9.814e-11 Hartree**; minimum core/active subspace singular values are
**0.9999999999999994 / 0.9999999999999991**.

**Only space160 has passing original QC.** The [fresh optimizer companion](../compressed-singlet-qc/README.md)
retains the failed space80 CI flag. Successful fixed-orbital diagnostics cannot
retroactively qualify that orbital-optimization record. The separately declared
[native import](../../cas12-compressed-singlet-qualified-import-contract.json)
uses only space160, respects the current native queue priority and requires
all64 native roots, dipole and core/active orbital-subspace checks. Scattering
and full model/threshold qualification remain separate.

## Public replay

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/compressed-singlet-coverage
tar -xzf projects/ukrmol_co/qualification-evidence/compressed-singlet-coverage/compressed-singlet-qc-coverage.tar.gz \
  -C projects/ukrmol_co/qualification-evidence/compressed-singlet-coverage
PYTHONPATH=projects/ukrmol_co/qualification-evidence/compressed-singlet-coverage/state-averaged-evidence \
  uv run python projects/ukrmol_co/qualification-evidence/compressed-singlet-coverage/state-averaged-evidence/verify-singlet-qc-coverage.py \
  projects/ukrmol_co/qualification-evidence/compressed-singlet-coverage/state-averaged-evidence
```

The archive reconstructs **438 payloads / 7075 coverage profiles**, every
recorded gate decision, mixed parent QC verdicts, unchanged checkpoint hashes,
analytic/perturbed controls and immutable source. Public fetch, reconstruction
and byte-identical repackaging pass.
