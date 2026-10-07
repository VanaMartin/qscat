# Stretched CAS(10,11) TZ/aug-TZ QC

At **R=2.5 bohr**, both cc-pVTZ and aug-cc-pVTZ state-averaged targets pass
tight QC when started from the independently covered/imported stretched-DZ
space-80 checkpoint. Each uses two frozen core orbitals, ten active electrons,
eleven active orbitals and forty equal-weight spin/irrep ensemble components.

| Basis | Wall | Kernel memory peak | Final orbital gradient |
|---|---:|---:|---:|
| cc-pVTZ | 874.04 s / 14.57 min | 0.345 GiB | 7.077e-8 |
| aug-cc-pVTZ | 1219.58 s / 20.33 min | 0.434 GiB | 5.791e-8 |

The sequential four-core batch takes **2094.73 seconds / 34.91 minutes**.
Original engine/batch/controller exits are zero. Orbital, ordinary CI,
fresh-CI, spin, MO orthogonality and Pi-degeneracy checks all pass; numerical
CI space is 80 and the ordinary CI cap remains 200 cycles. No numerical retry
was needed. These are QC/self-consistency results.

## Same-geometry basis differences

The root inventory is identical across the DZ-lineage models. Differences
below compare energies ranked within each spin/irrep; they do not establish
state identities across bases.

| Difference, second minus first | DZ → TZ | TZ → aug-TZ |
|---|---:|---:|
| Ground energy | −0.825431 eV | −0.0248054 eV |
| Ensemble objective | −0.865111 eV | −0.0517186 eV |
| Maximum absolute ranked-root shift | 0.913800 eV | 0.100352 eV |
| Maximum absolute ranked-excitation shift | 0.0883686 eV | 0.0755466 eV |
| Ground z dipole | −0.00374061 a.u. | −0.0110101 a.u. |
| Minimum cross-basis active-subspace overlap | 0.985481 | 0.998189 |

Cross-basis overlaps use the AO cross-overlap integral at unchanged nuclear
coordinates. Their values are diagnostics; the same-basis restart overlap
gate is not applied to different bases. Independent fixed-orbital coverage,
native imports, competing starts and continuity across geometry remain
separate qualifications.

## Public evidence and reconstruction

The **170 payloads** contain both runs, the completed controller/source,
unchanged DZ seed, passing parent-import record and raw resource profiles.
The executable verifier reconstructs **80 raw final QC states**, fresh analysis,
both same-geometry basis comparisons and **10430 resource samples**. All payload
sizes/digests, byte-identical repackaging and a byte-for-byte public fetch are checked.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/stretched-basis-qc
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/stretched-basis-qc/stretched-basis-qc.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run --with pyscf==2.11.0 --with h5py==3.15.1 \
  python "$EVIDENCE_DIR/state-averaged-evidence/verify-stretched-basis-qc.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

The [QC recipe](../../calibration-sa11-stretched-basis-qc.json),
[parent stretched import](../stretched-import/README.md) and
[continuation](../../CONTINUATION.md) define lineage and remaining gates.
Successor PID **959637** advances 600-cycle fixed-orbital coverage on the
original checkpoints, followed by individually gated
[eight-root-per-sector native imports](../../calibration-sa11-stretched-basis-imports.json)
on CPUs 8–11. Each passing scan must satisfy all root/residual/spin flags and
space/root-count consistency; each import retains the independent physical
root/dipole and exact-decimal serialization gates.
