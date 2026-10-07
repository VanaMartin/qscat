# Compressed CAS(10,11) DZ anchor

At **R=1.9 bohr**, the common forty-component state-averaged CAS(10,11)/cc-pVDZ
checkpoint passes all **48 fixed-orbital coverage probes**: five/eight roots,
trial spaces 40/80/160, both spins and all four irreps, with the previously
validated 600-cycle residual contract. Maximum physical/penalized residuals
are **7.034e-10 / 9.716e-10 Hartree**; the maximum ensemble-energy differential
is **1.422e-13 Hartree**. Spin, orthogonality and Hubbard analytic/perturbed-vector
controls pass.

The native selected-root target import agrees with every ensemble root and its
dipole. All **64 roots** also agree with the independent fixed-orbital coverage
oracle within **9.050e-9 Hartree**, including the extra roots. Target import
takes **3035.78 seconds / 5.262 GiB**.

The tight l=4 forty-channel dense scattering anchor completes both Pi sectors,
each dimension **27546**, in **11439.93 seconds / 25.277 GiB**. Its complete
energy grids, nonnegative final-state cross sections and Pi symmetry checks pass;
maximum phase splitting is **1.0e-9 rad modulo π**. The native fitted candidate
has position **3.826384 eV** and full width **2.322242 eV** in both sectors.
These are pipeline/anchor diagnostics, not qualified potential inputs: electronic
active-space/basis/ensemble selection, background sensitivity and continuity
remain open.

## Fetch and reconstruct

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/compressed-anchor
mkdir -p /tmp/co-compressed-anchor
tar -xzf projects/ukrmol_co/qualification-evidence/compressed-anchor/compressed-anchor.tar.gz \
  -C /tmp/co-compressed-anchor
uv run python /tmp/co-compressed-anchor/state-averaged-evidence/verify-compressed-anchor.py \
  /tmp/co-compressed-anchor/state-averaged-evidence
```

The portable verifier checks **1115 payloads / 75913 run-profile samples**, the
frozen source/checkpoint hashes and original completed coverage/import/scattering
step exits. It reconstructs residual/spin/orthogonality gates from the retained
scan records, ensemble/all-64-root agreement from native target data, the
DENPROP dipole, complete phase/cross-section grids, symmetry, native dimensions
and fitted candidates. Raw logs, inputs, checkpoints and resource records travel
with the archive. Residual vectors are produced by the retained validated PySCF
scan; the portable verifier checks their recorded norms, while rerunning that
scan in the pinned image recomputes Hamiltonian actions.

The parent controller's execution record is a historical snapshot after the
compressed steps; its separately owned stretched-scattering successor was still
active at capture. This companion certifies only the completed compressed
steps. Prior companions remain separate with their original snapshot hashes;
normalized tar/gzip repackaging is byte-identical.
