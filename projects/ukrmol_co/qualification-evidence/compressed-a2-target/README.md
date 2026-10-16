# Compressed CAS11 nine-A2 target/property evidence

The [input declaration](../../a2-compressed-target-preflight-contract.json)
completes at08:51:07UTC,16October,with all eight native target solves,DENPROP,
resources and owner exits0. It preserves the original1.9-bohrCAS11/cc-pVDZ
forty-component orbital average;explicit retention is[5,5,5,9]per spin.

The [44-file independent reading](../../verification-a2-compressed-target-20261016.json)
reconstructs **48complete target vectors/Gram matrices and915unique finite
multipole records**,with reference-energy error8.669e-13Hartree,native residual
1.319e-10Hartree,Gram error1.077e-14,Pi energy split5.684e-13Hartree and dipole
difference1.359e-11au. First-five preservation and own compressed state map pass;
A2root9IDs are48/45. Checkpoint,Molden,integrals/CONGEN remain unchanged.

Measured cost is **51.178minutes /5.477GiB**;DENPROP47.068minutes and target
solves246.065seconds. Host/macOS raw readings pass;quoted Gram metric differs
8.882e-16and other quoted metrics are identical. This qualifies native inputs
for the [separately declared compressed scattering](../../a2-compressed-scattering-contract.json).
It does not establish full channel convergence or independent accuracy of every
added transition moment.

## Reproduce

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/compressed-a2-target
EVIDENCE="$PWD/projects/ukrmol_co/qualification-evidence/compressed-a2-target"
REPLAY_ROOT=/path/to/fresh-external-directory
mkdir -p "$REPLAY_ROOT/raw" "$REPLAY_ROOT/source"
tar -xzf "$EVIDENCE/raw-inputs.tar.gz" -C "$REPLAY_ROOT/raw"
tar -xzf "$EVIDENCE/reader-source.tar.gz" -C "$REPLAY_ROOT/source"
cd "$REPLAY_ROOT/source"
PYTHONPATH=. uv run --with numpy --with scipy python -m projects.ukrmol_co.a2_target_reading \
  --directory "$REPLAY_ROOT/raw" --output "$REPLAY_ROOT/reconstructed.json"
```

The manifest pins all three reader source hashes and raw archive digest.
The reader verifies all44raw hashes before reconstruction. Six published payloads
preserve native exits,raw exports,decks,telemetry,references,reader and reports.
