# CO evidence index

Public, content-addressed bundles preserve accepted calculations **and original
failures**. Each study directory contains a small `artifacts.json`, its original
manifest and a README stating the measured result, limits and exact replay command.
Fetch verifies sizes and SHA256 before accepting bytes. Extract source capsules
outside the checkout: their historical tests are not current source modules.

## Scientific reading order

| Evidence | What it establishes |
|---|---|
| [Near-equilibrium neutral](neutral-near-equilibrium/README.md) | 52 QZ/5Z well points; relative-basis/interpolation, not dissociation |
| [QZ atomic endpoints](atomic-qz-endpoints/README.md) | Four finite-method totals, all fragment sums and electron affinities |
| [Matched5Z endpoints](atomic-5z-endpoints/README.md) | Reconstructs QZ and 5Z; rejects finite-basis insensitivity, not the finite calculations |
| [Original5Z startup stop](atomic-5z-launch-stop/README.md) | Zero-quantum control failure, retained separately from corrected success |
| [Stretched DZ/TZ pole contrast](stretched-pole-counterpart/README.md) | Local finite poles; physical position/width sensitivity rejects |
| [Retained-target overlaps](stretched-target-overlap/README.md) | Physical target correspondence, not unique N+1 resonance identity |
| [Compressed CAS11/CAS12 poles](compressed-cas12-pole/README.md) | Active-space position insensitivity rejects |
| [Compressed CAS12/TZ QC](compressed-cas12-tz-qc/README.md), [A2 pilot](singlet-a2-pilot/README.md) | Failed prerequisites; no valid physical TZ contrast |
| [Phase identifiability](phase-identifiability-endpoint/README.md) | Known wrong widths can predict unused phases; corrected endpoint isolation and original defect preserved |
| [Direct-S screen](direct-scattering-screen/README.md) | Rejects the unqualified bare full-S fitting shortcut |
| [Aug-TZ properties](stretched-augtz-properties/README.md), [preparation](stretched-augtz-preparation/README.md), [export blocker](stretched-augtz-inner-preflight/README.md) | Reusable inputs; no inner/pole release across incompatible open-channel grid rows |

## Reusable numerical and implementation controls

- [Sparse native boundary](sparse-boundary/README.md),
  [CAS11 refinement](cas11-sparse-refinement/README.md),
  [stretched DZ](stretched-dz-inner/README.md) and
  [TZ exports](stretched-tz-inner/README.md).
- [Native conventions](native-exterior-method/README.md),
  [full reconstruction](native-full-reconstruction/README.md),
  [Numerov corners](numerov-corners/README.md),
  [native residual kernel](numerov-native-residual/README.md).
- [Physical shell](physical-pole/README.md),
  [coupled tail](coupled-tail/README.md),
  [dual transport](dual-transport-toy/README.md),
  [coupled Jost](coupled-jost/README.md),
  [original real-guard failure](stretched-real-guard/README.md) and
  [separate corrected-prescription controls](stretched-real-completion/README.md).
- [Physical ensemble-frame rejection](ensemble-frame/README.md),
  [native target reuse](native-target-interface/README.md) and
  [signed CSF correspondence](native-csf-map/README.md).

Other subdirectories retain narrower intermediate controls and stopped attempts.
Their READMEs describe those finite studies, not current successor authorization.
[Permanent conclusions](../../../docs/physics/co-results-and-rejected-routes.md)
and [production decisions](../PRODUCTION_DECISIONS.md) govern interpretation.

## Fetch and reconstruct

For example, reconstruct the matched atomic comparison without quantum calls:

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/atomic-qz-endpoints
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/atomic-5z-endpoints
uv run --no-sync python -m projects.ukrmol_co.verify_atomic_5z \
  projects/ukrmol_co/qualification-evidence/atomic-5z-endpoints/evidence.tar.gz \
  projects/ukrmol_co/qualification-evidence/atomic-5z-endpoints/collection.json \
  projects/ukrmol_co/qualification-evidence/atomic-qz-endpoints/evidence.tar.gz \
  projects/ukrmol_co/qualification-evidence/atomic-qz-endpoints/collection.json
```

The root bundle is the earlier 46-attempt electronic-qualification supplement,
including 160 fixed-orbital probes, checkpoints, source/licence/build evidence,
neutral failures and selected-root controls. Later companions do not overwrite it:

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/qualification-evidence.tar.gz \
  -C "$EVIDENCE_DIR"
EVIDENCE_ROOT="$EVIDENCE_DIR/state-averaged-evidence"
uv run python -m projects.ukrmol_co.collect "$EVIDENCE_ROOT" \
  --provenance "$EVIDENCE_ROOT/qualification-provenance.json" \
  --output "$EVIDENCE_DIR/reconstructed.json"
```

Exact reanalysis/repackaging uses the analysis source recorded in that manifest;
later analyzers record their own hashes. Raw saved-data reconstruction is not a
new integral/eigensolve. Native replay, numerical tolerances and absolute model
accuracy must be reported separately.

## Historical narrative

The [campaign-history bundle](campaign-history/README.md) preserves exact
superseded handoffs, progress figures and rendering inputs. It is optional audit
material: the clone's permanent notes already state methods, measurements,
assumptions, limitations and rejected routes. Small JSON regression inputs and
execution contracts remain tracked; untracked reports/scratch notes do not ship.
