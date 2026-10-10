# CAS11 quarter-interval state tracking

The [quarter control](../../cas11-quarter-interval-state-tracking-contract.json)
tests R=**2.224225 / 2.408075 bohr** between equilibrium,midpoint and stretched
CAS11/cc-pVDZ frames. It targets singlet-A2 retained-root5's direct5→8 versus
midpoint5→6→5 route disagreement. Both new ordinary-QC starts pass before32
separately declared600-cycle,twelve-root probes at spaces160/240. Existing
equilibrium/midpoint/stretched states are reused without an eigensolve.

## Independent reconstruction

[The portable replay](../../cas11-quarter-state-tracking-replay-contract.json)
reconstructs **298 payloads / 672 wavefunctions / 32 probes / 112 matrices /
22916 producer profiles**. Wavefunctions include384 quarter states across two
trial spaces plus288 reused states. AO Hamiltonian actions,spin penalties,
ground dipoles,space-matched projections,orthogonality and Pi checks pass.
Full occupied Slater minors independently reproduce frozen-core Schur overlaps;
all coarse archived comparisons and four/two/direct assignment maps reproduce.

| Quantity | Host maximum / measured cost |
|---|---:|
| Physical residual / Hartree | 6.9258e-10 |
| Penalized residual / Hartree | 9.9852e-10 |
| Spin-penalty action / Hartree | 6.5825e-10 |
| Native/saved root error / Hartree | 5.4591e-10 |
| Ground-dipole error / a.u. | 5.2958e-11 |
| Full occupied-minor overlap error | 1.3323e-15 |
| Producer wall / kernel peak | 4604.51 s / 1.5686 GiB |
| Host capture/pack/replay wall / kernel peak | 208.04 s / 2.3492 GiB |

Host replay owner **1141177** exits0. It verifies the original **11 source /
124 input hashes** plus a fresh **12-source / 293-input** replay freeze. Original
files stay read-only. Public fetch,independent reconstruction and byte-identical
repackaging pass on macOS in **280.15 seconds**; full occupied-minor overlap error
is **2.7756e-15**. This has a different scope/platform from the host replay and
does not establish a speedup. The [manifest](manifest.json) records both report hashes.

## Geometry-path finding

Singlet-A2 root5 follows **5→6→6→5→5**,with physical squared overlaps
**0.81897 / 0.77757 / 0.74913 / 0.82329**. It ends at stretched5 as in the
midpoint route,rather than the direct8 assignment.

| Fine physical geometry link | Weak assignments / 96 |
|---|---:|
| Equilibrium → lower quarter | 2 |
| Lower quarter → midpoint | 1 |
| Midpoint → upper quarter | 6 |
| Upper quarter → stretched | 1 |

Weak means squared overlap<0.5; each pair has its own96 assignments. Fine and
two-leg maps still differ for **9/96 roots**,including **2/40 retained components**,
in both physical and AO-following metrics: singlet-A2 root4 ends at6 rather than9;
triplet-A2 root5 ends at5 rather than9. Fine-versus-direct differences are25/96.
Across seven pairs and both metrics,the replay reproduces202 weak,82 ambiguous
and30 retained-to-extra triggers. AO-following transport is not physical overlap.
**Global state continuity,resonance-pole identity and model verdicts stay unset.**

## Reproduce from public bytes

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/quarter-state-tracking
REPLAY=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/quarter-state-tracking/quarter-state-tracking.tar.gz -C "$REPLAY"
EVIDENCE="$REPLAY/state-averaged-evidence"
uv run --with pyscf==2.11.0 --with h5py==3.15.1 python \
  "$EVIDENCE/verify-quarter-state-tracking-20261010.py" "$EVIDENCE"
python "$EVIDENCE/pack-small-anchor-controls.py" "$EVIDENCE" "$REPLAY/repacked.tar.gz"
```

Extract outside the checkout. SHA256 is
**`17ed7af6894ad1e4b3cafccb853bfa5ba67a9b5ee2ae4aa1d2ca1b29bd7630fc`**.
The archive includes ordinary QC,checkpoint/CI arrays,producer source,parent
evidence,native tables,profiles and independent portable controls.
The [full host report](https://data.qscat.org/ukrmol-co-quarter-state-tracking-2026-10-10/independent-verification.10a55b5992be.json)
records all reconstructed maps and errors.
