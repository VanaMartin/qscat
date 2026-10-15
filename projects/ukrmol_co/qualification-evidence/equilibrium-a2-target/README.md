# Equilibrium CAS11 nine-A2 target/property evidence

The [input contract](../../a2-equilibrium-target-preflight-contract.json) completes
with all eight target solves,DENPROP48,resource wrapper and owner exits0.
An independent read of **44 checksum-fetched files** reconstructs **48 complete
native target vectors** and **930 unique finite dipole/quadrupole records**.
The original forty-component orbital average and orbital-frame hashes remain
unchanged. This qualifies inputs for the separately declared
[equilibrium channel comparison](../../a2-equilibrium-scattering-contract.json).

The [local report](../../verification-a2-equilibrium-target-20261016.json)
checks CIDATA headers,three centers with charges6/8/0,complete vectors,Gram
matrices,energies against independent twelve-root references and the original
first-five states,Pi degeneracy,telemetry,residuals,property-state irrep codes,
own noncontiguous map and ground dipole. Singlet/triplet A2root9are IDs48/45.
Vector Gram error is at most1.077e-14;reference energy difference5.969e-13Hartree;
native residual1.20921e-10Hartree;dipole difference7.984e-11au.

Native wall cost is **51.817minutes / 5.425GiB** kernel peak. DENPROP takes
47.675minutes. Host and macOS reconstructions pass;their largest Gram-error
metric difference is1.443e-15,with other quoted metrics identical. The shared
reader is a raw reconstruction,not a separate electronic-structure solver.
The multipole inventory does not independently validate every added transition
moment or establish full channel convergence.

## Reproduce

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/equilibrium-a2-target
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/equilibrium-a2-target"
# Use a new external directory, outside the checkout.
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/raw" "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/raw-inputs.tar.gz" -C "$REPLAY_ROOT/raw"
tar -xzf "$EVIDENCE/reader-source.tar.gz" -C "$REPLAY_ROOT/source"
cd "$REPLAY_ROOT/source"
PYTHONPATH=. uv run --with numpy --with scipy python -m projects.ukrmol_co.a2_target_reading \
  --directory "$REPLAY_ROOT/raw" --output "$REPLAY_ROOT/reconstructed.json"
```

`manifest.json` pins all three reader source files and the raw archive digest.
The reader validates every raw manifest hash before reconstruction.
`first-reader.py` preserves the rejected supplementary reader's two-atom center
assumption;CIDATA also contains the zero-charge scattering center. Corrected
source checks all three records. No native outcome or numerical budget changed.
