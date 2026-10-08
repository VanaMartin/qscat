# Compressed CAS12 alternative singlet-driver pilot

The separate eight-probe pilot uses PySCF's `direct_spin0_symm` singlet driver
with the same spin penalty and physical Hamiltonian on both unchanged compressed
checkpoints. It requests five/eight roots at trial-space sizes 80/160 with a
600-cycle diagnostic cap. **All eight probes pass** in **252.48 seconds /
0.6062 GiB**. The [original 96-probe diagnosis](../compressed-cas12-diagnostic/README.md)
retains every rejected singlet-A1 result from the original driver.

The driver enforces alpha/beta exchange symmetry, which alone is insufficient
to certify total spin. Separate physical, spin-penalized, spin-penalty action,
spin-square, CI orthogonality and full-Hamiltonian-action checks remain mandatory.
The independent action uses `direct_spin1.contract_2e` on the returned vectors;
maximum action difference is **6.752e-10 Hartree**, and its independently measured
physical residual is **5.514e-10 Hartree**. Maximum physical/penalized residuals
from the driver are **8.708e-10 / 9.988e-10 Hartree**. All CI flags pass.

Common-root count/space differences fall below **2.843e-13 Hartree**, compared
with the original 1.6561-Hartree mismatch. The roots differ from the rejected
optimizer's saved singlet-A1 energies by up to **1.670142 Hartree**. Thus this
pilot qualifies a fixed-orbital numerical driver, not the original orbital
objective. The [fresh QC contract](../../cas12-compressed-singlet-driver-qc-contract.json)
declares two new orbital optimizations, gated by pinned-image analytic and
nontrivial mixed-spin CASSCF controls. Original QC stays rejected, and independent
coverage remains required for any newly optimized checkpoint.

## Public replay

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/compressed-singlet-driver-pilot
tar -xzf projects/ukrmol_co/qualification-evidence/compressed-singlet-driver-pilot/compressed-singlet-driver-pilot.tar.gz \
  -C projects/ukrmol_co/qualification-evidence/compressed-singlet-driver-pilot
PYTHONPATH=projects/ukrmol_co/qualification-evidence/compressed-singlet-driver-pilot/state-averaged-evidence \
  uv run python projects/ukrmol_co/qualification-evidence/compressed-singlet-driver-pilot/state-averaged-evidence/verify-singlet-driver-pilot.py \
  projects/ukrmol_co/qualification-evidence/compressed-singlet-driver-pilot/state-averaged-evidence
```

The archive verifies **49 payloads / 1259 profiles**: complete recorded gate
decisions, analytic/perturbed controls, frozen source, original QC diagnostics
and unchanged checkpoints. Runtime independent Hamiltonian-action checks are
retained in the per-state records; archive replay reconstructs those decisions.
Public fetch, reconstruction and byte-identical repackaging pass.
