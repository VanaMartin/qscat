# Near-equilibrium correlated neutral pilot

The completed neutral pilot contains **26 independently calculated geometries
per basis** for aug-cc-pVQZ and aug-cc-pV5Z: the 25-point 1.9–2.5-bohr grid
at 0.025-bohr spacing plus the independent 2.1323-bohr reference. Six existing
anchors are reused byte-exactly and **46 new calculations** finish successfully.

Energies use spherical RHF/frozen-core CCSD(T), two frozen occupied spatial
orbitals and ten correlated electrons. Dipoles use the all-electron CCSD lambda
density, not a CCSD(T) energy derivative. Every point passes RHF convergence,
internal/external stability, CCSD/lambda convergence, finite amplitude diagnostics
and the 14-electron density check. Each basis has its own 2.1323-bohr energy zero.

## Independent held-point and basis checks

Fourteen interpolation knots comprise the 0.05-bohr grid and the additional
reference. Twelve independently calculated interval midpoints are withheld.
The [predeclared contract](../../neutral-near-equilibrium-contract.json) sets
1-meV energy and 0.001-a.u. dipole budgets, without extrapolation.

| Basis / interpolant | Maximum held energy error / meV | Maximum z-dipole error / a.u. |
|---|---:|---:|
| QZ / not-a-knot cubic | 0.11386 | 1.775e-7 |
| QZ / PCHIP | 0.95428 | 7.192e-6 |
| 5Z / not-a-knot cubic | 0.11330 | 1.840e-7 |
| 5Z / PCHIP | 0.91453 | 7.211e-6 |

All held-point budgets pass. The maximum QZ→5Z relative-energy difference is
**17.2433 meV**, within the 20-meV pilot budget. This constrains near-equilibrium
basis and interpolation error; it does not establish the correlation treatment,
the complete neutral potential or a dissociation limit. Four original stretched
references at 3.0/4.0 bohr retain their external-RHF-instability rejection and
stop before CCSD(T).

The 46 new job walls sum to **11557.49 s / 3.210 h**, excluding the reused anchors,
controller gaps and analysis. Maximum kernel container charge is **24 GiB**;
the maximum sampled anonymous memory is **17.138 GiB**. Raw profiles distinguish
container cache/overhead from resident numerical arrays.

## Fetch and verify

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/neutral-near-equilibrium
mkdir -p /tmp/co-neutral-pilot
tar -xzf projects/ukrmol_co/qualification-evidence/neutral-near-equilibrium/neutral-pilot.tar.gz \
  -C /tmp/co-neutral-pilot
uv run python /tmp/co-neutral-pilot/state-averaged-evidence/verify-neutral-pilot.py \
  /tmp/co-neutral-pilot/state-averaged-evidence
```

The standalone verifier checks **959 payloads**, **63867 profile samples**, the
completed owner/source hashes, original anchor digests, raw RHF/CCSD/triples
energies, stability and stage exits, density/convergence flags, geometry/basis
inputs, checkpoint provenance and both interpolation/basis reconstructions.
Linear/cubic analytic interpolation controls verify the held-point machinery.
The original rejected references and their stopped stage sequences are checked.

RHF checkpoints travel with the archive. Larger CCSD checkpoints remain on the
calculation host with byte counts and SHA256s recorded. The original QZ anchor
source is retained separately from the later new-job source; recorded
source checksums resolve each anchor and new calculation. The owner snapshot
includes pointers/hashes for prior public companions, whose
archives remain separate. Deterministic repackaging reproduces the archive bytes.
