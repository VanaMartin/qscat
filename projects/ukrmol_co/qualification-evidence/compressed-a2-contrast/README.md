# Compressed CAS11 five→nine outgoing-pole contrast

**Pass in both sectors at1.9bohr:** increasing retained A2 targets from five to
nine per spin lowers position by **10.876503meV** and full width by
**4.480322meV**. The declared budgets are50meV position and104.230719meV full
width, with frozen operational bands. This supports the finite compressed
pole-target contrast only; the stretched phase rejection and all electronic/
omitted-spectrum/continuity/full-scattering uncertainties remain.

| Own finite model | Position /physical eV | Full width /physical eV |
|---|---:|---:|
| CAS11 /five-A2 /250 channels /2048 roots | 3.771266190140 | 2.089094708262 |
| CAS11 /nine-A2 /298 channels /2048 roots, reused | 3.760389687393 | 2.084614386540 |

[Result/resources](../../progress-compressed-a2-contrast-result-20261025.json) ·
[Independent public reconstruction](../../verification-compressed-a2-contrast-public-20261025.json) ·
[Design](../../../../docs/physics/co-compressed-a2-pole-contrast.md) ·
[Execution contract](../../compressed-a2-contrast-execution-20261025.json)

Source/tests/exact owner freeze `dbc83e1` at23:17:52UTC precedes launch23:18:30UTC
and finish23:22:07UTC on25 October. A pre-execution note's23:20 caption was an
approximate prospective label, not a measured event; current captions are
corrected while the frozen capsule remains unchanged. All original scientific/
command/container/owner and nested pole exits are0. No timeout/OOM or retry.

The four own-model real-S guards pass; maximum entry/phase errors are
1.49918e-7 /4.23259e-7rad. Boundary roots, fixed32-node contour with nested16,
two derivative residues, visible rank-one Pi support and five-model B1/B2 pairing
all pass. Cross-model label-matched profile TV is.00084001, below the necessary
.25 identity screen; it is not continuous geometry tracking. Five leakage is
1.83348e-9 /1.83341e-9; same-model profile TV2.31364e-11 /norm mismatch2.17441e-10.

B1 uses six initial root evaluations; its tighter settings and all three B2
settings accept the existing candidate directly. Identical reported coordinates
are **not independent root iterates, exact degeneracy or zero uncertainty**.
Contour estimates remain distinct and within the same-root budgets.

Cost: **47 B1 +42 B2 =89/146 maps**,216.873008s container /215.532646s parent
wall,1.737797GiB cgroup peak,721.883760s cgroup CPU at final report save.
Parent30.346392CPU seconds excludes pole subprocesses/workers. Allocation:
eight physical cores0–7 /one BLAS thread /8GiB /3600s combined cap. Original
inputs verify123 before/after;42 source hashes also match frozen Git blobs.
Twenty owner/adapter/character controls plus eleven inherited pole controls per
sector pass in the pinned runtime;45 focused local tests preceded freeze.

Ten fresh public payloads authenticate168 result/source files and24 narrow saved
prerequisites. Independent native-Q contraction, opposite `I−YQ` solves, bordered
traces, contour moments, derivative residues, four own-native K-to-S guards and
all character/identity/full-width decisions reconstruct. Saved exterior maps
are authenticated, **not independently propagated**. Tiny replay drift is not an
absolute CO error bound. No nine repropagation, native solve or holdout analysis.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/compressed-a2-contrast
uv run --no-sync python projects/ukrmol_co/qualification-evidence/compressed-a2-contrast/reconstruct-saved.py \
  projects/ukrmol_co/qualification-evidence/compressed-a2-contrast/artifacts.json \
  verification-compressed-a2-contrast-fresh.json
```

**Decision:** retain nine; prioritize one matched retained-target CAS11/CAS12
model control using saved compressed target vectors. Five supplies a **local
matched comparator**, not a cheaper domain-wide production model. All114 holdout
S matrices and equilibrium pole exclusion stay protected; new50 remains parked.
No accurate electronic model, curve, fit or sweep is released.
