# CAS11 many-electron target-state overlap screen

The [finite eight-root contract](../../cas11-target-state-overlap-screen-contract.json)
completes **24 fixed-orbital sector solves / 192 wavefunctions / 48 overlap
matrices** at the saved CAS11/cc-pVDZ anchors **R=1.9/2.1323/2.5 bohr**.
Original native-imported orbitals, model and forty-state ensemble are unchanged.
Eight roots per spin/irrep use space160 and a separately declared600-cycle
diagnostic budget. All CI flags, native first-five roots, physical/penalized
Hamiltonian residuals, spin and orthonormal-root gates pass. This is an electronic
identity diagnostic, with final target continuity and resonance-pole verdicts unset.

The [determinant-overlap method](../../../../docs/physics/co-target-state-overlaps.md)
includes both frozen cores and core/active cross terms. Physical cross-geometry
overlaps and separately labelled AO-following Lowdin transport are retained.
The portable verifier uses **full occupied Slater minors**, independent of the
producer's Schur reduction, to reconstruct every matrix from saved CI vectors.
It additionally rebuilds AO cross integrals and the physical/spin-penalized
Hamiltonian actions from original checkpoints with PySCF2.11.0.

## Measured checks and cost

| Quantity | Maximum / measured value |
|---|---:|
| Original native first-five root error / Hartree | 5.460e-10 |
| Reconstructed physical residual / Hartree | 6.947e-10 |
| Reconstructed penalized residual / Hartree | 9.956e-10 |
| Reconstructed spin-square error | 5.418e-14 |
| Independent full-minor overlap difference | 3.320e-14 |
| Wall / aggregate CPU seconds | 3703.42 / 11265.38 |
| Kernel memory peak / GiB | 0.836895 |

The four-core job completes with exit0 under its2-GiB cap alongside the large
CAS11 solver. Resources include **18455 sampled profiles**. Extra roots6–8 are
independent fixed-orbital CI controls; only the original five roots per sector
are compared with the native forty-state import.

## State-assignment findings

Singlet-A2 roots **1↔2** exchange their strongest matches from compressed to
equilibrium geometry in both metrics: assigned physical squared overlaps are
approximately **0.4325**, while AO-following values are approximately **0.946**.
Therefore energy-order labels alone do not identify a target state across these
anchors.

From equilibrium to stretched geometry, first-five retained roots assign to
extra roots in **singlet.A2 (4→7, 5→8), triplet.A1 (4→7), and triplet.A2
(4→6, 5→8)**, in both metrics. Examples in AO-following transport are
triplet-A1 **4→7 at0.7954** squared overlap and triplet-A2 **4→6 at0.8558**.
These are finite-manifold maximum-overlap assignments, not certified continuous
paths; all competing row scores and assignment triggers are preserved.

| Geometry pair / bohr | Metric | Weak assignments | Ambiguous rows | First-five → extra | Minimum manifold singular value |
|---|---|---:|---:|---:|---:|
| 1.9→2.1323 | Physical | 64 | 2 | 0 | 9.314e-11 |
| 1.9→2.1323 | AO-following | 2 | 2 | 0 | 1.357e-10 |
| 2.1323→2.5 | Physical | 64 | 21 | 5 | 3.995e-11 |
| 2.1323→2.5 | AO-following | 12 | 6 | 5 | 9.720e-11 |
| 1.9→2.5 | Physical | 64 | 64 | 5 | 3.428e-13 |
| 1.9→2.5 | AO-following | 21 | 10 | 5 | 1.352e-12 |

Counts cover eight sectors × eight assignments per pair/metric. Weak means
assigned squared overlap<0.5; ambiguous means the row-best/second-best squared
overlap margin<0.1. Across **384 assignments**, there are **227 weak / 105
ambiguous / 20 first-five-to-extra triggers**. Physical values include loss from
the moving atom-centred core orbitals; AO-following transport identifies
coefficient frames and is not a physical wavefunction overlap. Nearly singular
singlet-B1/B2 eight-root manifolds occur even for the first adjacent pair.
Projection deficits combine geometry/orbital-frame loss and excluded states;
they are not isolated omitted-state populations.

The [fresh twelve-root/space240 closure control](../../cas11-target-state-overlap-roots12-contract.json)
preserves all eight-root findings and separately tests whether extending the
finite root manifold recovers missing matches. Its same-geometry energy-matched
projection gate includes all twelve new roots, allowing a near-degenerate cut
at root8. Enlarged-manifold agreement still cannot establish a full continuous
path, select an electronic model or qualify resonance poles.

## Public reconstruction

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/target-state-overlap
REPLAY=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/target-state-overlap/target-state-overlap.tar.gz -C "$REPLAY"
EVIDENCE="$REPLAY/state-averaged-evidence"
uv run --with pyscf==2.11.0 --with h5py==3.15.1 python \
  "$EVIDENCE/verify-target-state-overlaps.py" "$EVIDENCE"
python "$EVIDENCE/pack-small-anchor-controls.py" "$EVIDENCE" "$REPLAY/repacked.tar.gz"
```

The archive contains **118 payloads**: original checkpoint/native-table inputs,
frozen source/control files, complete CI arrays/energies/spins/residuals, physical
and AO-following metrics, overlaps/assignments, execution records, logs and
resource profiles. Its deterministic packer retains the earlier small-control
archive layout. Fetch/reconstruction and byte-identical repackaging reproduce
SHA256 **`68458d72f565ed527146fc8f47db60f4ef22d94a418e7ad7f46dbd6cf3f766bd`**.
