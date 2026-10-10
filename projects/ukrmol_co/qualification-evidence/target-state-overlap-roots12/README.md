# CAS11 twelve-root target-state overlap closure

The [twelve-root control](../../cas11-target-state-overlap-roots12-contract.json)
extends the [eight-root screen](../target-state-overlap/README.md) at unchanged
CAS11/cc-pVDZ checkpoints, **R=1.9/2.1323/2.5 bohr**. Its24 fixed-orbital solves
retain twelve roots per sector at space240 with the separately declared600-cycle
diagnostic cap. Five native roots per sector remain the scattering import;
roots6–12 are diagnostic CI controls. No orbital optimization or native import
is repeated by this control.

## Independent reconstruction

The [portable replay contract](../../cas11-target-state-overlap-roots12-replay-contract.json)
reconstructs **171 payloads / 288 wavefunctions / 48 overlap matrices / 21542
profiles**. It checks original-eight energies and energy-matched projection
weights, native first-five roots, physical/spin-penalized actions,spin and root/MO
orthogonality. Full occupied Slater minors independently reproduce the producer's
frozen-core Schur-complement overlaps, including core/active cross terms.

| Check | Maximum error |
|---|---:|
| Native first-five root / Hartree | 5.4591e-10 |
| Physical residual / Hartree | 6.9258e-10 |
| Penalized residual / Hartree | 9.9713e-10 |
| Spin-square | 1.5099e-14 |
| Host full-minor overlap | 1.3323e-15 |
| Public-download full-minor replay on macOS | 2.7756e-15 |

Producer cost is **4322.37 seconds / 1.1766 GiB**, exit0. Public download,
reconstruction and exact repackaging take **109.97 seconds** on the independent
macOS replay. This is not a matched engine-performance comparison.

The larger manifold recovers several nearly missing directions, but still has
**344 weak / 161 ambiguous / 20 retained-to-extra triggers** across the physical
and AO-following comparisons. AO-following transport is a coefficient-frame
diagnostic, not physical overlap. **Continuous-path state and resonance-pole
identity remain unset.** The later [midpoint control](../midpoint-state-tracking/README.md)
tests the larger interval and projected-start dependence.

## Reproduce from public bytes

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/target-state-overlap-roots12
REPLAY=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/target-state-overlap-roots12/target-state-overlap-roots12.tar.gz -C "$REPLAY"
EVIDENCE="$REPLAY/state-averaged-evidence"
uv run --with pyscf==2.11.0 --with h5py==3.15.1 python \
  "$EVIDENCE/verify-target-state-overlaps-roots12.py" "$EVIDENCE"
python "$EVIDENCE/pack-small-anchor-controls.py" "$EVIDENCE" "$REPLAY/repacked.tar.gz"
```

Extract outside the checkout. The archive includes the immutable owner sources,
native/original-eight inputs,all twelve-root vectors,metrics,matrices,assignments,
raw logs and resources. Fetch/reconstruction and deterministic repackaging
reproduce SHA256 **`2a2ce0fa29d9d42b790bc999693f45f777c8423e72acfe80298f4c824d4fb613`**.
The [manifest](manifest.json) links the full public independent report.
