# Stretched CAS(10,12) fixed-orbital residual coverage

**Both completed R=2.5-bohr DZ repair checkpoints pass all 96 probes** under the
previously validated **600-cycle** fixed-orbital residual contract. The saved
orbitals, twelve active orbitals, ten active electrons, two frozen core orbitals,
five/eight roots, trial spaces 40/80/160, tight tolerances and spin penalty remain
unchanged. This qualifies lowest-root coverage of these two saved checkpoints.

Maximum physical/penalized residuals are **7.638e-10 / 9.997e-10 Hartree**,
spin-squared error **8.216e-15**, and ensemble-energy difference **2.701e-13
Hartree**. Trial-space and common five/eight-root differentials are at most
**1.706e-13 Hartree**. Controller 1017316 exits zero after **2240.11 seconds /
0.6252 GiB**.

The [original 200-cycle coverage](../cas12-base-repairs/README.md) continues to
reject under its original contract. The fifth singlet-A2/triplet-Pi roots and
extra eighth roots that failed there converge in the longer scan. A separate
initial wrapper exits one before any probe because its extracted imports omit
`hashlib`; that immutable source/log/resource record is preserved within this
companion. Neither accurate energies nor a corrected wrapper erase prior flags.

## Fetch and reconstruct

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/cas12-residual
mkdir -p /tmp/co-cas12-residual
tar -xzf projects/ukrmol_co/qualification-evidence/cas12-residual/cas12-residual.tar.gz \
  -C /tmp/co-cas12-residual
uv run python /tmp/co-cas12-residual/state-averaged-evidence/verify-cas12-residual.py \
  /tmp/co-cas12-residual/state-averaged-evidence
```

The portable verifier checks **627 payloads / 11138 scan-profile samples**, all
96 residual/spin/orthogonality/coverage records, saved checkpoint hashes, the
Hubbard analytic and perturbed-vector controls, passing QC consistency, original
failed 200-cycle probes and the rejected wrapper. Full independent Hamiltonian
action recomputation uses the retained validated PySCF scan in the pinned image;
the portable verifier reconstructs gates from its raw recorded norms and roots.
Owner/producer snapshots, unchanged checkpoint inputs, logs and profiles travel
with the archive. Normalized tar/gzip repackaging is byte-identical.

Native all-root/dipole import and larger-basis target qualification are distinct
successors. This result alone releases no scattering sweep.
