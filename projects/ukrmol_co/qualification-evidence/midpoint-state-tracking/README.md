# CAS11 two-start midpoint state tracking

The [midpoint contract](../../cas11-midpoint-state-tracking-contract.json)
tests **R=2.31615 bohr** between the equilibrium and stretched CAS11/cc-pVDZ
anchors. Two projected endpoint starts retain the forty-component ensemble,
two frozen core orbitals,ten active electrons and eleven active orbitals.
Both ordinary200-cycle QC calculations pass before the separately declared
600-cycle fixed-orbital diagnostic. All32 probes request twelve roots in
spaces160/240. Existing endpoint vectors are reused.

## Independent and public reconstruction

The [fresh replay contract](../../cas11-midpoint-state-tracking-replay-contract.json)
reconstructs **231 payloads / 576 wavefunctions / 32 probes / 96 overlap matrices /
40963 producer profiles**. It rebuilds AO integrals,physical/spin-penalized/action
residuals,spins,ground dipoles,space-matched projections and full occupied Slater
minors; it also verifies row/global assignments and composed path maps.
All original **11 source / 62 endpoint-input hashes** are retained.

| Quantity | Measured maximum / cost |
|---|---:|
| Reconstructed physical residual / Hartree | 6.9258e-10 |
| Penalized residual / Hartree | 9.9877e-10 |
| Spin-penalty action / Hartree | 6.5193e-10 |
| Saved/native first-five root error / Hartree | 5.4591e-10 |
| Ground-dipole error / a.u. | 2.8356e-11 |
| Host full-minor overlap error | 1.3323e-15 |
| Public-download overlap error on macOS | 2.9976e-15 |
| Producer wall / kernel peak | 8231.19 s / 1.5059 GiB |
| Host capture/pack/replay wall / kernel peak | 164.61 s / 1.8718 GiB |

Public fetch/reconstruction/exact repackaging costs **187.54 seconds** on macOS;
this is not a matched performance comparison with the Linux host replay.

## What the midpoint resolves and leaves open

Two-start ensemble objectives agree within **1.7053e-13 Hartree**,first-five
roots within **1.2259e-8 Hartree**,and ground dipoles within **3.5712e-8 a.u.**
Same-geometry twelve-root manifold minima are about **0.999999999985** in both
metrics. Four higher diagnostic root components nevertheless differ by more
than1e-7 Hartree between independently converged frames,up to **3.8876e-7**.
The mean objective is not a guarantee of every excited-root precision.

| Physical geometry pair | Weak assignments / 96 | Ambiguous rows | Retained → extra |
|---|---:|---:|---:|
| Equilibrium → stretched,direct | 96 | 39 | 5 |
| Equilibrium → midpoint,either start | 21 | 3 | 3 |
| Midpoint → stretched,either start | 38 | 13 | 2 |

These are separate pair diagnostics. Weak means squared overlap<0.5; ambiguous
means best-minus-second squared-overlap margin<0.1. Physical loss includes
moving atom-centred orbitals; separately labelled AO-following transport is
not physical wavefunction overlap.

Two-leg labels differ from direct labels for **16/96 roots**,including one
retained component: **singlet A2 root5 maps directly to stretched root8,but
through the midpoint maps5→6→5**. Its physical squared overlaps are **0.08076**
direct,and **0.41406 / 0.38844** on the two legs. The latter are still weak.
Both midpoint starts and both metrics give the same route disagreement.
The [quarter-interval control](../../cas11-quarter-interval-state-tracking-contract.json)
targets that discrepancy while retaining twelve roots. **Continuous-path state,
resonance-pole and model qualification remain unset.**

## Reproduce from public bytes

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/midpoint-state-tracking
REPLAY=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/midpoint-state-tracking/midpoint-state-tracking.tar.gz -C "$REPLAY"
EVIDENCE="$REPLAY/state-averaged-evidence"
uv run --with pyscf==2.11.0 --with h5py==3.15.1 python \
  "$EVIDENCE/verify-midpoint-state-tracking-20261010.py" "$EVIDENCE"
python "$EVIDENCE/pack-small-anchor-controls.py" "$EVIDENCE" "$REPLAY/repacked.tar.gz"
```

Extract outside the checkout. Complete producer source,ordinary-QC runs,
checkpoint/CI arrays,endpoint inputs,metrics and profiles travel with the
portable verifier. The two starts and original coarse comparisons remain
unchanged. Byte-identical repackaging reproduces SHA256
**`acf7b08581da773b3741ae1faf237afcaf179a5e9e40e295f9146b6971e1e7d2`**.
The [manifest](manifest.json) links the full public independent report.
