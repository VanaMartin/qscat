# Stretched DZ/TZ saved-input feasibility — passes for a finite64-state model

**Step1 of the queue is complete:** saved proof authenticates; choose all eight
available native roots in every spin/irrep,64 components per model. This avoids
repeating target solves and does not claim nine-A2 coverage, omitted-channel
bounds, electronic accuracy or production channel convergence.

[Public verification](../../verification-stretched-basis-feasibility-public-20261026.json) ·
[Reader](../../stretched_basis_feasibility.py) ·
[Route/objective](../../../../docs/physics/co-post-a2-model-route.md)

Nine freshly fetched public payloads verify217 archived files and177 original
frozen source/input entries. Actual DZ/TZ checkpoints, Molden files, native
CIDATA vectors, raw sector outputs, existing properties and full original source/
probe/qualification records are included. Large CONGEN CSF inputs are byte-pinned
in the collection rather than copied into this compact proof. Host collection
authenticates these selected inputs; no CSF transformation or native rerun occurs.

The reader independently decodes all64 roots per native file, integer phases,
complete record/set order, physical centers and CI-vector norms; checks raw
newly-solved-sector energies rather than inherited echoes; reconstructs property
counts and checks ordinary forty-state QC/objective flags. The TZ native binary
energies match the saved fixed-orbital CI within4.974e-13Ha; this tighter binary
comparison is not a replacement for the original rounded-output bound or an
absolute electronic-error estimate.128 saved widest-space vectors are norm-
checked and all48 saved probe records pass their original scalar gates.

**Independence limit:** original physical/spin action metrics and orbital-subspace
checks are authenticated original independent evidence. They are not freshly
recomputed here. No PySCF/package import, Hamiltonian action, CI/QC/native-engine
call, new property computation, outgoing map or holdout S analysis occurs.

The distinct stretched TZ-frame qualification remains exit0. Its original
seed-preservation failure/owner exit1 stays rejected. The original serialization
failure is preserved too. Neither case promotes the failed compressed TZ target.
A first packaging/reader attempt exits1 because tar retained fort.16/fort.24
symlinks without target members. Both that archive and failure record are public;
a fresh dereferenced package of the unchanged scientific bytes passes/exit0.

## Property inventory and minimum missing work

Native set order is singlet/triplet A1, B1, B2, A2.

| Inventory | DZ | TZ |
|---|---|---|
| Existing property counts |6,5,5,6,5,6,3,4|6,4,5,6,5,7,3,4|
| Existing states /multipoles |40 /827|40 /708|
| Missing states for selected64 |24|24|
| Original DENPROP alone |16742.48s (4.65h)|15842.32s (4.40h)|

Both **new64-state property tables are missing**. Reuse all existing native
vectors, target integrals and CSFs; no QC/target diagonalization/import repeat.
The next queue item is an exact property-only contract/owner freeze, not its
execution. Change only the existing `NTGTL` list to eight in all eight sets.
Validate own64 energy-ordered IDs, all required moments and ground dipole against
each frame. Never splice the old unequal tables into a claimed matched model.

Cost planning only: larger pair/state inventories can make these exports more
expensive than the measured40-state calls. Conditional6–16h per export,18h native /
19h owner ceiling,64GiB/eight physical cores is a proposed envelope, **not a
measured forecast, reserved job or threading-speedup claim**. Freeze exact input
pins, owner, floors and finite exit/stop rules separately before release. No
automatic second job after failure, import/preparation/scattering/pole or B2.

```sh
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence/stretched-basis-feasibility
uv run --no-sync python -B projects/ukrmol_co/qualification-evidence/stretched-basis-feasibility/reconstruct-saved.py \
  projects/ukrmol_co/qualification-evidence/stretched-basis-feasibility/saved-proof.tar.gz \
  projects/ukrmol_co/qualification-evidence/stretched-basis-feasibility/collection.json \
  verification-stretched-basis-feasibility-fresh.json
```

The pointer's Git commit records the saved-data collection checkout, not the
original producers. Each original frozen source/image/exit remains in the bundle;
the public reader is separately content-pinned and matches the tracked reader.
No model, curve, fit or sweep is qualified. Nine remains working reference;
2048/new50/equilibrium/114-holdout protections persist. No numerical job queued.
