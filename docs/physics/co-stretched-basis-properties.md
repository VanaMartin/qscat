# Two missing stretched CAS12 property inventories

## Decision served and non-release

The [saved-input review](co-post-a2-model-route.md) selects a finite stretched
DZ/TZ comparator with eight existing native roots per spin/irrep,64 components.
Both older property tables contain40 states, in different sector distributions,
and each omits24 selected states. **One new property inventory per basis is a
missing input, not a repeated successful identical computation.**

This stage freezes [source/tests/exact owner and contract](../../projects/ukrmol_co/stretched-basis-properties-contract-20261026.json).
It releases **no export**. Next is a separately reviewed DZ-only application;
TZ requires independent public verification of DZ and another explicit release.
No owner launches its peer. All preprocessing remains project-local; no new
physics method, solver, Rust kernel or qscat promotion is involved.

## Exact reuse and change

At R=2.5bohr, retain the existing optimized forty-component orbital objective,
two A1 core orbitals, ten active electrons and6/3/3/0 active orbital counts.
Compute properties on the existing **64 native CI vectors**, not a new64-state
orbital average. The original distinct-TZ-frame proof/seed-preservation failure
and compressed A2/TZ failures remain unchanged.

For each basis reuse exact CIDATA `fort.26`, target one-electron properties
`fort.16`, eight CONGEN expansions, Molden and checkpoint. Copy them into a
fresh writable work directory; original run and source mounts stay read-only.
Protect the old `fort.24` too, but do not copy it as a new output. Change only
DENPROP's `NTGTL` from the saved lists to `8,8,8,8,8,8,8,8`; preserve NFTSOR,
NTGTF/set order, units and every other deck byte. The two translated decks have
the identical prospective SHA256
`4fbc4ebbb87ebd2adf7ff1def36ffe9dc1356babba184c0f05974283ea0b77a4`.
Each basis still uses its own orbital/CI/integral data and ground dipole.

Use the previously validated property image and exact DENPROP binary, not a new
engine build or package installation. Thirty-six original file hashes/eighteen
per basis and five source prerequisites are pinned. No QC, orbital rotation,
target eigensolve/reimport, integral generation, CONGEN or scattering is requested.

The owner runs synthetic/mocked format/routing controls before the only native
call. Read all eight own CIDATA roots per sector, retaining their physical energy
convention; compare the new table with those energies, not another basis's roots.
Require exit0/completion,64 energy-ordered sequential IDs/eight per sector,
energy error≤1e-7Ha, finite unique dipole/quadrupole records covering all states
and both ranks, own ground z-dipole agreement≤1e-5au and unchanged original/copied
input hashes before/after. Reuse the established `property_inventory` check;
do not invent an unverified exact multipole count or claim a new completeness
theorem. Native output retains its original units; comparisons are in atomic
units with the ordinary Hermitian molecular convention, not the ECS c-product.

## Failure routing and bounded owner

Each owner/application can make **one** DENPROP call only. The persistent user
service bounds the entire owner to19h; the container client has18.5h and the
native process18h. Retain controls, copied-input/deck/output hashes, original
native/scientific/container/owner/service exits, log, source/image/archive pins
and `run_profiled` RAM/CPU/disk samples in the fresh capsule. No restart or
repair/retry follows a failed flag, property, input, dipole, control, cap or OOM.
Do not remove a failed container or reuse its output directory.

TZ accepts only a separate, hash-pinned **independent** DZ receipt within the
persistent root: correct contract digest,64 targets/eight per sector, all four
original scientific/command/container/owner exits0, public-reader exit0,
independent inventory verdict true and a property-output SHA256. The exact
receipt is a runtime output, so its path/hash must be frozen in the later TZ
launch record; it cannot be predicted now. Producer-only success, an unset
receipt or any failed/absent exit blocks TZ before a scientific container is
created. This is a serial dependency, not an automatic two-model queue.

Reserve eight distinct physical cores4–11/64GiB/1GiB shared memory only if later
released. Minimum available80GiB RAM and10GiB persistent disk; one library thread,
no network, no interference with unrelated workloads. Conditional6–16h/export
is planning only; old40-state DENPROP measured4.65/4.40h. No measured64-state
cost or eight-core speedup is claimed, and nothing is reserved at this freeze.
No numerical successor is released by a property pass; independently qualify the
actual outputs before considering preparation/B1 scattering or pole work.

## Freeze verification and remaining uncertainty

The [freeze-control evidence](../../projects/ukrmol_co/qualification-evidence/stretched-basis-properties-freeze/README.md)
records eight passing controls in the pinned image, twelve authenticated source
entries and all36 original inputs unchanged before/after. The original run root
was not mounted in the control container; all native/solver calls in controls
are mocked. Host owner imports with the standard library alone. Known-truth tests
reject changed units/counts/energies/missing moments and failed/unpinned DZ receipts;
they check one-call/no-retry routing, read-only/no-network mounts and service cap.
These are execution-readiness checks, not sweep or physical-accuracy progress.

The implementation/freeze allocation remains60minutes; no general workflow API
or scheduler was added. The finite64-state prescription is not nine-A2 or
channel/basis completeness. Prior stretched channel rejection, compressed B1
model-position rejection, compressed basis uncertainty, continuity, neutral
correlation and full-S embedding remain. Preserve2048 baseline, nine working
reference, parked new50, equilibrium exclusion and114 untouched holdout S matrices.
No model, curve, fit, sweep, native export or automatic successor is promoted.

## Separate DZ-only execution release

The later [DZ-only release](../../projects/ukrmol_co/stretched-basis-properties-dz-release-20261026.json)
uses the unchanged frozen source/contract, with a fresh committed-source capsule
and host eligibility review. It does not retroactively turn the freeze into a
scientific calculation. [Launch evidence](../../projects/ukrmol_co/qualification-evidence/stretched-basis-properties-dz-launch/README.md)
witnesses DENPROP active after all eight controls pass, with original/source/copy/
deck pins intact. Inventory, final cost and terminal exits remain pending.
No completed64-state property table or physical-model verdict is claimed by a
launch. Independently qualify/publish the actual DZ output before considering
the still separately gated TZ application; no automatic downstream work follows.

## Completed DZ input and independent public receipt

The [DZ table now qualifies](../../projects/ukrmol_co/qualification-evidence/stretched-basis-properties-dz/README.md):
64 states/eight per sector/2040 dipole/quadrupole records, own energy/dipole gates,
original exits0,4.552h/1.453GiB, no timeout/OOM. A separate standard-library reader
reconstructs the raw native energy/vector-norm and property checks plus resource
sample peaks from freshly downloaded evidence. It imports no producer code and
does not recompute Hamiltonian actions or integrals. This validates the missing
property input, not the target's physical accuracy or cross-basis continuity.

The [independent receipt](../../projects/ukrmol_co/qualification-evidence/stretched-basis-properties-dz-receipt/README.md)
pins contract/output/exits/public-reader0 and is itself publicly content-addressed.
The host fetches/hash-checks it before the [separate TZ-only release](../../projects/ukrmol_co/stretched-basis-properties-tz-release-20261031.json).
No new scientific method, orbital objective, preparation/scattering/pole job or
production model follows automatically. Preserve the original distinct-TZ-frame
qualification/seed failure and compressed A2/TZ stops. TZ uses its own saved
energies/dipole, not the newly successful DZ result as a numerical substitute.
