# CAS(10,11) outer-propagation refinement

This companion checks radial propagation at fixed inner-region channels,
R-matrix amplitudes, forty target components and the 99-point requested
0.1–5.0-eV grid. Both Pi components use the same four refinement controls:

| Control | Gailitis matching radius (bohr) | Propagation subranges |
|---|---:|---:|
| Original baseline | 100 | 8 |
| h2 | 100 | 16 |
| h4 | 100 | 32 |
| r150 | 150 | 26 |
| r200 | 200 | 36 |

The Legendre order remains ten. The new radius controls use subintervals of
at most 5.125 bohr. All **forty native stages** complete, taking **1669.54
seconds / 1.868 GiB** together on four physical cores. Grids, cross-section
finiteness/nonnegativity/final-state sums, Pi agreement and the chosen
position/width/phase gates pass.

Identical printed refined phases are resolved with the **full-precision
unformatted K-matrices**, which are included for every control and baseline.
The pinned `outerio.f` writer stores each symmetric upper triangle by columns.
Reconstruction with `sum(arctan(eigvalsh(K)))` agrees with the native printed
phase sum modulo pi within **4.87e-8 rad**, its rounding precision. The largest
refinement-minus-baseline difference is **2.56845e-6 rad**. Differences between
successive refined controls stay below **2.28e-10 rad**. Linear-background
fits of the binary phases in fixed 1.6–3.5/2.0–3.1-input-eV intervals change
position/full width by at most **2.63e-8 / 4.35e-8 eV**. They use the actual
native Hartree energies, including the pinned 0.0735-Ryd/input-eV convention.

The source and runtime audit confirms the perturbations reached the active
path: `ASYM1` reads `BPROP` and builds all requested sectors through
`RPROP1_MPI`; `ASYM2` uses `CURLYR` and `GAILIT`; the energy loop then applies
`RPROPX` before forming K. Every requested sector is accounted for across all
four rank logs, and the reported matching radii/subrange counts match the
inputs. The four audited engine files exactly match the SHA256-pinned upstream
archive. This distinguishes a resolved numerical plateau from inactive knobs.

The reconstructed Cayley S-matrices pass unitarity/reciprocity checks; these
are algebraic checks implied by the stored symmetric K and do not independently
certify the electronic model. This refinement qualifies the radial propagation
for the retained model and observables. Active space, basis, continuum/channel
selection and background-fit sensitivity retain their separate gates.

The archive preserves the original successful controller exit, lease/resumption,
resources, inputs, raw outputs, all rank logs, binary K data, fitted records,
source snapshots and executable verifier. All **296 payloads** verify by size
and digest, and repackaging is byte-identical. The public client fetched and
verified the exact archive. The main supplement and other companions retain
their original snapshots.

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/outer-propagation
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/outer-propagation/outer-propagation.tar.gz \
  -C "$EVIDENCE_DIR"
PYTHONPATH=. uv run python "$EVIDENCE_DIR/state-averaged-evidence/verify-outer-propagation.py" \
  "$EVIDENCE_DIR/state-averaged-evidence"
```

See [the main evidence](../README.md),
[the qualification note](../../../../docs/physics/co-electronic-qualification.md#outer-only-propagation-refinement)
and [the continuation](../../CONTINUATION.md) for the remaining local controls.
