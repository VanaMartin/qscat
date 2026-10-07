# CAS(10,10) native sparse boundary and omitted-spectrum controls

**The fixed-model 2048-root replay passes both Pi sectors.** Native PETSc
`mpisbaij`/SLEPc Krylov–Schur selects 2048 of 8350 scattering states, exports
the native MPI boundary data and reproduces the qualified dense reference's
complete 99-point grid and identical 1.8–3.3-requested-eV fits at background
orders 1–4. This qualifies this smaller-model implementation and its selected
spectrum; larger models require independent refinement/resource qualification.

| Selected poles | Maximum phase error / rad modulo π | Position/width gates, both sectors |
|---:|---:|---|
| 128 | 1.565795 | Reject at every fixed background |
| 512 | 1.469455 | Reject at every fixed background |
| 2048 | 2.0e-7 | Pass at all four backgrounds |

At 2048 poles, maximum eigenvalue error is **9.664e-13 Hartree**, absolute
residual **1.336e-10 Hartree**, and sign-aligned boundary-amplitude error
**1.003e-6**. Channels, thresholds and multipoles match the dense oracle
byte-exactly. Maximum position/full-width differences are **5.847e-8 /
2.381e-7 eV**. Native phases also match the same dense truncation exactly at
printed precision, separating iterative/export errors from omitted-spectrum error.

The separate dense-spectrum omission sequence is **128/512/2048/4096/6144/
8192/8350**, in both sectors. It confirms stability above 2048: phase error
falls to **1e-8 rad** at 4096 and vanishes at printed precision for higher
counts. Every higher-count fixed-window fit passes. The complete endpoint
reproduces the original dense boundary bytes and phase grid.

The omitted boundary-amplitude squared norm falls from **193.909** at 512 poles
to **5.7385e-4** at 2048, of a full norm **234.074**. This measures participation,
not an observable error bound; the propagated phase and fit differentials supply
the numerical gate. No partitioned or missing-state correction is applied.

## Measured replay costs

| Replay | Wall / seconds | Peak container / GiB |
|---|---:|---:|
| Original 128-root B1, including failed parser | 402.64 | 1.228 |
| 128-root B2 | 427.10 | 1.346 |
| 512-root B1+B2 | 1043.47 | 1.554 |
| 2048-root B1+B2 | 1340.74 | 1.861 |
| Fourteen dense-omission outer-only controls | 1898.71 | 0.374 |

These replays reuse electronic, integral and CONGEN data; prior full-pipeline
costs are separate. The 2048-root inner/export stages take **558.81 / 531.48 s**
for B1/B2. Both report 2548 Krylov vectors, 57 iterations and over 2048 converged
roots; exactly 2048 states are exported. The native selection remains the
unmodified largest-magnitude default; all selected roots are the lowest 2048
here because the entire active spectrum is negative.

## Fetch and independently reconstruct

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/sparse-boundary
mkdir -p /tmp/co-sparse-boundary
tar -xzf projects/ukrmol_co/qualification-evidence/sparse-boundary/sparse-boundary.tar.gz \
  -C /tmp/co-sparse-boundary
uv run python /tmp/co-sparse-boundary/state-averaged-evidence/verify-sparse-boundary.py \
  /tmp/co-sparse-boundary/state-averaged-evidence
```

The verifier checks **306 payloads / 25504 resource samples**, immutable owner
and producer-source hashes, original native exits, solver telemetry, all native
CI eigenvalues/coefficients, boundary/channel records, phase grids, native
candidates, identical fixed-window fits and all fourteen dense omissions. The
full dense native CI and boundary oracles travel with the archive, allowing
independent sign-aware reconstruction from original records.

The original B1 controller's analysis exit one and raw interleaved zero-error
RESON summary remain preserved. The repaired parser removes only explicitly
empty rank summaries. Capturing raw telemetry/phase units needed two additive
capture supplements; those scripts travel with the initial capture. Original
scientific outputs remain immutable. The earlier
[SWINTERF failure companion](../sparse-native-control/README.md) separately
retains the rejected partitioned path, pinned engine source/build and telemetry
sidecar provenance. Normalized tar/gzip repackaging is byte-identical.

The [CAS11 contract](../../sparse-scattering-cas11-contract.json) declares the
next finite dense-omission/native sequence. No CAS12 job is released by the
CAS10 result, and this implementation check does not establish active-space,
basis, ensemble/channel selection or geometry/pole continuity.
