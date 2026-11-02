# Finite stretched CAS12 diffuse-flexibility control

## Physical question and scope

The qualified DZ/TZ finite poles differ by103.370meV in position and44.172meV in
full width, above50meV/5% physical budgets. Retained N-electron target subspaces
correspond predominantly, so an obvious lost/swapped retained manifold is not
observed. Neither model is a production electronic choice. Under the user's
continuation, select **one aug-TZ versus plain-TZ control atR=2.5bohr**, not a
QZ/CAS13/channel/root ladder. Diffuse radial flexibility can affect polarization
and a low-energy resonance even when valence target manifolds largely correspond.
The saved ground dipoles differ substantially; this motivates the control but
does not predict the pole shift or establish which model is more accurate.

Hold the equal40-component SA objective, active6/3/3/0 (12 active orbitals), two
frozen cores, eight roots per spin/irrep and64 retained targets fixed. Each model
uses its own qualified optimized frame, native B1/2048 export and matching
properties. Thus this measures **diffuse-basis-plus-frame sensitivity**, not an
isolated diffuse-function perturbation with frozen MOs. Units remain Hartree/bohr
internally; positions/full widths are reported in eV at the boundary.

## Reuse evidence and missing work

[Saved aug-TZ feasibility](../../projects/ukrmol_co/qualification-evidence/stretched-augtz-feasibility/README.md)
reconstructs the public48-probe supported80/160/240 coverage and the actual
native64-root import. Binary native energies differ from saved independent probes
by at most3.85122e-8Ha; the printed native-table difference3.84631e-8Ha differs by
expected print rounding. The actual imported checkpoint is separately pinned;
original seed-preservation error3.40484e-8Ha passes the original1e-7Ha gate. Saved
physical core/active subspaces pass. These are authenticated saved residual/frame
witnesses, not fresh Hamiltonian actions or independently regenerated AO integrals.

All eight native CSF hashes match the already qualified CAS12 plain-TZ definitions.
No QC, CI eigensolve, target import, integral or target CSF generation is repeated.
The actual native property file has40 states/812 unique moments, counts
5/3/6/7/6/7/3/3 in native order:24 of the common64 are missing. Reusing that file
as64, or padding missing moments with zero, is prohibited. Aug-TZ full-continuum
preparation/B1 eigenpair export and outgoing-pole qualification are not supplied
by the old target import.

## Staged execution and interpretation

Only the **first missing64-state DENPROP export** is released now, after exact
source/tests/contract/owner freeze and push. Change only NTGTL to eight in every
native sector; all other deck bytes, original orbitals/one-electron property
integrals/CSFs/native coefficients remain unchanged and read-only. The resulting
deck hash equals the already executed DZ/TZ64-state prescription. One native
DENPROP call maximum, no solver fallback/retry/automatic successor. Preserve old
40-state properties and all original exits, including the later compressed queue
failure and first symlink-packaging failure. A scientific/native/input/count/
dipole/integrity/timeout/OOM failure stops inconclusive, with no gate relaxation.

Gates reuse the established64-state property inventory: native completion/exit0,
64 sequential energy-ordered target IDs/eight per sector, own CIDATA energies
within1e-7Ha, finite unique valid rank1/2 moments covering all IDs and both ranks,
ground dipole within1e-5au, immutable copied/source inputs before and after.
Do not invent an exact common DZ/TZ/aug-TZ moment count. Independent public
reconstruction of the result is required before preparation is considered.

Later stages require separate decisions and freezes: own minimal full-continuum/
B1 preparation and2048-root native export; finite real/outgoing-map qualification
with own seed/box, numerical0.1meV position/full-width and1e-6rad phase budgets;
physical target correspondence and local residue/isolation review before a
finite-pole contrast. Nothing here preauthorizes those stages. There is no
automatic B2 rescue, geometry curve, fit, sweep or protected-holdout analysis.

The eventual finite comparison retains50meV position and5% full width with1meV
floor as physical candidate budgets. A rejection keeps simple basis-insensitivity
inference stopped and redirects electronic-model choice. A pass supports only
this local diffuse-flexibility contrast, not absolute basis/correlation/channel
convergence, unique N+1 identity or a production fitting target. No post-hoc
profile threshold or promotion because aug-TZ is larger is permitted.

## Resource and implementation envelope

The old aug-TZ40-state DENPROP stage measured15886.897s (4.413h); the **entire**
old import measured5.894h/41.114GiB, not a property-only memory measurement.
Matched64-state DZ/TZ exports measured4.552h/1.453GiB and4.421h/1.715GiB. The final
old aug-TZ DENPROP-duration sample window peaks near1.438GiB, an approximate
stage attribution, not a stage-specific kernel peak. Identical CSF dimensions
and matched64 prescription justify a conditional **4–7h** forecast, not an
aug-TZ64 measurement or eight-core speedup.

One property job: eight physical cores4–11, one library thread,8GiB/no swap/
network,1GiB shared memory, minimum16GiB available RAM/10GiB disk; native12h,
container12h30m, owner13h hard caps. Runtime no-native controls first, new persistent
capsule and service names, bounded read-only300s observations. At most60minutes
narrow implementation/freeze effort. Native later B1 scope would be a separate
hours/approximately30GiB-class allocation based on the measured DZ/TZ5.6–6.2h
exports; no full-route runtime is guaranteed before actual aug-TZ preparation.

This missing-input stage is necessary for the physical contrast, **not** itself
a measured scattering effect or production-readiness advance. No paid resources,
image rebuild, ComfyUI changes, held-out S analysis or qscat promotion is selected.

## Completed property stage and minimum preparation direction

[Independent public qualification](../../projects/ukrmol_co/qualification-evidence/stretched-augtz-properties/README.md)
supplies64 states/1928 valid unique moments; all eight original execution exits0.
Measured4.780938h/1.693340GiB lies inside the forecast, not a parallel speedup.
The first metadata-packaging failure is preserved beside corrected reading0;
maximum monitor gap300.614284s fails strict300s cadence, independently of native
and scientific success. No campaign numerical job/monitor remains active.

Next freeze only one own aug-TZ full-continuum integral/conditional algebraic B1
CSF reuse step. Actual native inventory/resolved deck/eight target-CSF hashes must
equal qualified DZ before copying signed CSFs and original CONGEN log; integral,
property, orbital and native-vector files remain aug-TZ-specific. On mismatch,
stop without fallback CONGEN. The existing release is DZ/TZ-only; a new exact
source/tests/receipt/owner freeze and push is needed before a call. Proposed
four cores/32GiB,30–180s conditional reuse forecast,2400s integral/3480s container/
3600s owner,30s observations and30minute freeze effort. Those are planning bounds
anchored to plain-TZ38.671s/4.174GiB, not measured aug-TZ preparation. No native
inner solve, pole or production work is released by this property qualification.
