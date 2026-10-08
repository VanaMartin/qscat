# Compressed extraction and three-anchor orbital screens

Three finite controls run beside the large CAS11 sparse solve, using fresh
2-GiB containers and unchanged saved inputs. Their scientific execution exits
are all zero; the extraction gates reject and state identity stays unset.

## Compressed extraction

The [declared saved-data control](../../cas11-small-anchor-controls-contract.json)
fits both Pi sectors over requested-eV windows **1.5–6.5, 2–6 and 2.5–5.5**,
all below the first excited target threshold. Each window uses background
orders1–4, three starts and five held-point folds: **24 cases / 192 data fits /
120 held folds**. Actual fitting energies are `requested_eV * 0.0735 / 2`
Hartree. Exact synthetic and phase-branch invariance controls pass.

All cases are numerically valid and multistart checks pass, but all-declared
positions span **0.08993 eV**, widths **0.74645 eV** around a **2.31773-eV**
median (**32.206%**). These exceed the 0.05-eV / 5% gates. Constant-background
held-point errors reach **0.120 rad**, also rejecting their 0.05-rad gate;
higher-background held errors are smaller. Every declared background case
remains in the verdict. Both Pi sectors give identical results.

## Threshold and orbital screens

The matched CAS11/cc-pVDZ forty-state anchors at **R=1.9/2.1323/2.5 bohr**
retain native/QC target checks, Pi degeneracy and AO-metric orthogonality.
Adjacent physical core/active minimum singular values are **0.9132/0.9693**
and **0.8072/0.9319**; the end-to-end core/active minima are **0.6001/0.8249**.
The original core follow-up triggers remain recorded.

The separately declared [AO-transport follow-up](../../cas11-ao-transport-screen-contract.json)
uses the same ordered atom/basis/AO labels and symmetric metric square roots.
Adjacent AO-transport core/active minima are **0.9986/0.9873** and
**0.9983/0.9873**; end-to-end minima are **0.9940/0.9500**. The contrast shows
strong sensitivity of physical core overlaps to the moving atom-centred basis.
AO transport is an explicitly chosen coefficient-frame identification, not a
physical wavefunction overlap or a many-electron state assignment. Same-geometry
identity, reciprocity and subspace sign/rotation controls pass. Final state/pole
continuity remains unset.

| Control | Profiled wall / s | Kernel peak / GiB |
|---|---:|---:|
| Extraction | 0.605 | 0.0874 |
| Physical orbital/threshold screen | 0.402 | 0.0931 |
| AO-transport follow-up | 0.402 | 0.0683 |

## Public reconstruction

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/small-anchor-controls
REPLAY=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/small-anchor-controls/small-anchor-controls.tar.gz -C "$REPLAY"
EVIDENCE="$REPLAY/state-averaged-evidence"
uv run --with pyscf==2.11.0 --with h5py==3.15.1 python \
  "$EVIDENCE/verify-small-anchor-controls.py" "$EVIDENCE"
python "$EVIDENCE/pack-small-anchor-controls.py" "$EVIDENCE" "$REPLAY/repacked.tar.gz"
```

The portable verifier reconstructs **94 payloads / 10 resource samples**,
raw phase-point alignment, all fit/held-point/span gate decisions, target
threshold records, original checkpoint/source hashes and measured costs.
Cross-geometry AO integrals are recomputed. AO-transport values are checked
using SciPy metric square roots and Gram-matrix eigenvalues, independently of
the producer's NumPy eigendecomposition/SVD route. Public repackaging is
byte-identical; SHA256:
`194becd3530da71465e5d85defce4b081a3041be46c4adccd49972070d31b29a`.

The [separate compressed outer-grid control](../../cas11-compressed-outer-grid-contract.json)
tests energy-grid refinement from saved inner-region data. It cannot erase the
original all-background extraction rejection or assign a resonance pole.
