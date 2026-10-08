# Compressed CAS11 retained-data grid controls

The [finite outer-grid release](../../cas11-compressed-outer-grid-contract.json)
reuses the unchanged CAS11 compressed anchor at **R=1.9 bohr**, both Pi sectors,
and saved channels/inner-region amplitudes. All **20 native stages** complete
on **317/633-point** grids, requested-eV spacings **0.025/0.0125**, in
**634.56 seconds / 0.7989 GiB** with four ranks. These are outer-only costs.

Common original 0.05-eV-grid phases differ by at most **3.2e-6 rad modulo π**.
Complete grids, finite/nonnegative channel-summed cross sections and Pi symmetry
pass. All **48 complete-data window/background cases** retain three starts and
five held-point folds (**384 data fits / 240 held folds**). Matched finest-pair
position/full-width differences are at most **0.3651 / 2.0296 meV**, passing
the declared grid gates. The public report includes all three resolutions,
including the original saved-grid fits; changing grid spacing is distinct
from changing the fitted representation.

**Extraction remains rejected:** all-declared position spans are
**0.08886/0.08832 eV**, width spans **32.096/32.043%** of the median. These
remain outside the **0.05-eV / 5%** gates, consistent with the original
[saved-grid rejection](../small-anchor-controls/README.md). Constant-background
held errors also remain recorded. No background/window case is removed.

Native `RESON` reports **MAXFIT** on both fine grids; the finest grid has no
printed fitted candidates. Native fit acceptance and final pole identity stay
unset. Numerical grid stability does not qualify the one-candidate extraction
or select a production electronic model.

## Public reconstruction

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/compressed-outer-grid
REPLAY=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/compressed-outer-grid/compressed-outer-grid.tar.gz -C "$REPLAY"
EVIDENCE="$REPLAY/state-averaged-evidence"
PYTHONPATH="$EVIDENCE/outer-source" uv run --with pyscf==2.11.0 --with h5py==3.15.1 python \
  "$EVIDENCE/verify-compressed-outer-grid.py" "$EVIDENCE" --source "$EVIDENCE/outer-source"
uv run --with pyscf==2.11.0 --with h5py==3.15.1 python \
  "$EVIDENCE/verify-small-anchor-controls.py" "$EVIDENCE"
python "$EVIDENCE/pack-small-anchor-controls.py" "$EVIDENCE" "$REPLAY/repacked.tar.gz"
```

The companion verifies **257 payloads / 3170 profiles**: 3160 outer samples
plus the ten parent-screen samples. It contains the exact retained binary
inputs, native decks/logs/original exits, full phase/cross-section grids,
all multistart/held-fold fits, immutable source and the parent orbital screens.
The portable verifier reconstructs raw grids, fit residuals, gate decisions,
resource peaks and the full 0.05/0.025/0.0125-eV refinement sequence. Public
fetch/replay and byte-identical repackaging pass. SHA256:
`37437c02cc958971ce32ae3f5c7dbbc637dd1f53a0c609fc1102602f7b6da7c1`.
