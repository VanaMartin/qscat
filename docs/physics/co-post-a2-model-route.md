# Electronic-model route after the stopped A2 pilot

## Decision and scope

**14:56UTC completion:** [saved-input feasibility passes](../../projects/ukrmol_co/qualification-evidence/stretched-basis-feasibility/README.md)
for the finite all-eight-per-sector/64-state comparison. The saved distinct-TZ
proof is now public and its recorded gates reconstruct; no fresh Hamiltonian
action was run. Two64-state property exports are missing. Next is their exact
contract/owner freeze, not execution. Historical selection/release rationale follows.

**Redirect the next basis-model investigation to the already constructed
stretched CAS12 DZ/TZ targets at R=2.5bohr. Do not continue the compressed A2
driver/optimizer campaign.** This is a saved-evidence allocation decision,
not a basis-insensitivity result, new target qualification or numerical release.
The [compressed A2 failure](../../projects/ukrmol_co/qualification-evidence/singlet-a2-pilot/README.md)
and original failed TZ optimization stay rejected under their own contracts.
The86.774662meV compressed CAS11/CAS12 pole-position rejection also remains.

The decision-changing uncertainty is whether the finite CAS12 scattering model
is sensitive to DZ→TZ in a geometry where valid target construction already
exists. A stretched comparison cannot settle the missing compressed comparison,
prove accuracy, or repair the stretched channel-retention rejection. It avoids
making a new CI solver a prerequisite for every possible model test.

## The full ensemble must be a valid objective

For each atomic basis B and symmetry-preserving orbital rotation U, define

\[
F_B(U)=\frac1{40}\sum_{S\in\{0,1\}}\sum_{\Gamma\in\{A_1,B_1,B_2,A_2\}}
\sum_{j=1}^{5} E_{S\Gamma j}(U;B),
\]

where the E are normalized, orthogonal, physical Hamiltonian roots in the
specified total-spin/spatial-irrep subspace. Ten active electrons occupy
6/3/3/0 active orbitals; two inactive A1 cores stay doubly occupied. Weights are
per C2v component, not multiplied by triplet spin degeneracy. Energies are
Hartree, lengths bohr, and the orbital metric is the positive AO-overlap metric.
The molecular wavefunctions use the ordinary Hermitian inner product, not ECS.

An optimized target needs valid inner CI roots **and** orbital stationarity for
this objective. A spin penalty or alpha/beta-exchange symmetry is not a substitute
for total-spin/residual gates. Failed inner CI destroys the interpretation of a
small reported orbital gradient. Substituting frozen-frame A2 energies after a
failed optimization does not construct a stationary F_B.

The saved stretched candidates report all forty original/fresh flags, correct
weights/spins, objective arithmetic and final gradient/energy-step gates passing.
Their objective definition is unchanged across DZ/TZ; optimized frames change
with the atomic basis. This is a reoptimized basis-model contrast, not a pure
fixed-orbital perturbation, isolated orbital-count effect or state-identity proof.

## Alternatives considered

| Route | Disposition and reason |
|---|---|
| Increase A2 space/cycles, change seed, run eight roots, or relax residual gates | Not selected. Violates the stopped finite pilot and supplies no independent model direction. |
| Drop A2 root5/sector, renormalize35 weights, or reuse the failed TZ frame | Not selected. Changes the intended full ensemble or promotes an invalid optimization. A different ensemble would need a scientific definition, not numerical convenience. |
| Exact-spin/CSF inner solver inside a fresh forty-state CASSCF | Mathematically compatible with F_B, but not an available validated drop-in route in the reviewed evidence. Needs root/RDM/gradient and mixed-optimizer qualification before full QC; this would be a new solver/integration campaign. No implementation selected. |
| Project qualified DZ orbitals into TZ and use fixed-orbital native CI | A legitimate separately declared frozen-frame model could be built, but it does not optimize F_TZ or answer the selected relaxed-basis question. Not selected as a surrogate for the failed target. |
| CAS13, new50, more inner roots, or compressed CAS12 B2 | No stronger decision value established. B2 cannot undo the B1 rejection; new50 remains parked;2048 stays the numerical baseline. |
| Stretched CAS12 DZ/TZ, existing full-ensemble targets | **Selected for saved-input feasibility.** Same geometry/active counts/forty-state objective and existing64-root native records. No QC/target solve need be repeated. |
| Stretched DZ→aug-TZ | Not selected now. Adds diffuse augmentation to the zeta change and another branch; plain TZ already supplies the missing independent basis direction. Its existing results remain preserved. |

No package source or installed solver was run to make this choice.

## What the saved inventory establishes

[Public metadata](../../projects/ukrmol_co/qualification-evidence/stretched-basis-route-review/artifacts.json)
preserves seventeen saved records/pins, including exact config/target/diagnostics/
native-analyzer/resource JSON and the distinct TZ-frame qualification. Four public
payloads were freshly fetched/hash-checked, with objective/ranked-energy arithmetic
reconstructed. No PySCF import, eigensolve, Hamiltonian action, QC, native call,
outgoing map or holdout S analysis occurred. **This is metadata evidence, not a
new independent full-target verification.** Checkpoint/Molden entries are byte
pins, not the actual orbitals in this small metadata bundle.

| Saved model | DZ | TZ |
|---|---:|---:|
| Forty-state objective /Ha |−112.47157983595963|−112.50464389165757|
| Ground energy /Ha |−112.8878844899605|−112.91906900834246|
| Ground z-dipole /a.u. |−.180058755355614|−.1887579008757018|
| Final orbital gradient |5.68552e-8|3.72191e-8|
| Native root inventory |64|64|

The maximum difference of corresponding **ranked** excitations is.01085996877Ha.
This is a diagnostic, not a matched-state shift or a pole-budget test. No threshold
energies/dipoles are substituted for a scattering resonance comparison.

The TZ checkpoint must be the separately declared
[reoptimized target](../../projects/ukrmol_co/cas12-tz-reoptimized-target-qualification-contract.json),
SHA256 `8518dc191b8670af9ec87d7666b091eaf370187e7022684cd93716d240b40568`.
Its original seed-preservation failure/owner exit1 stays intact; the existing
distinct-frame record reports ordinary forty-root QC/native/full-action checks
passing. Unlike compressed TZ, the original optimization's CI flags already
pass. **These are not interchangeable cases.** Full public/raw reconstruction
of that distinct-frame proof is still a prerequisite, not newly claimed here.
The stretched DZ import already has
[public independent evidence](../../projects/ukrmol_co/qualification-evidence/dz-import-queue-failures/README.md).

## Minimum next deliverable and release boundary

The30minute saved-evidence review is complete. **Next: one saved-input feasibility
and retention/cost decision**, at most60minutes implementation/review effort,
zero new CI/QC/native/outgoing calls. Assemble the existing distinct-TZ-frame
proof, native root/property inventories and required original provenance into a
publicly replayable record; authenticate actual checkpoints/Molden/native vectors
against their saved proofs without recomputing successful targets. Stop if the
needed evidence is absent or inconsistent; do not recreate it by rerunning QC.

Declare a common retained target prescription before requesting properties or
scattering. Eight native roots per sector are not nine, and global-lowest40
properties are not automatically firstfive-per-sector. The stretched CAS11
five→nine phase rejection makes a local five-state comparison insufficient for
a production channel-convergence claim. A smaller matched retained model could
only be a labelled finite diagnostic; no pass would promote it or bind omitted
channels. Decide explicitly whether its rejection/pass changes a useful choice
or whether a larger common retained space is required. No automatic channel ladder.

Use this saved-only deliverable to identify exactly which property/preparation/
scattering inputs are missing and cost **those**, not new QC. Existing target
imports took5.93/6.53h and≈40.37/41.13GiB; these are sunk whole-import costs,
not forecasts for a property-only or scattering job. The compressed86352-state
inner solve's4.23h/22.49GiB on12 physical cores is a planning anchor, not a
measured stretched/TZ dimension or runtime. No multi-hour reservation or paid
provisioning is released by this review. Any later small/large/very-large jobs
retain the user's4/8/12-physical-core policy and need their own frozen caps.

Only after saved-input/retention feasibility passes may a separate exact execution
contract propose missing properties/preparation and eventually B1 scattering/
pole comparison. Preserve50meV position,5% reference FULL width with1meV floor,
0.05rad phase and the separate same-root budgets. A local pass would not settle
compressed basis convergence, active-space saturation, pole/state continuity,
neutral correlation, omitted spectrum, independent exterior or full-S embedding.

The failed compressed checkpoint is never an input to this redirect. Nine stays
the working reference, not a production model;2048/new50/equilibrium exclusion
and all114 untouched S holdouts remain protected. No geometry sweep, fit,
production model or successor numerical job is active or queued.

## Completed saved-input/retention/cost decision

The user asks to execute queue step1. Within the60minute cap, authenticate217
archived files/177 original frozen entries and freshly fetch nine public payloads.
Read both actual native CIDATA files:64 roots each, orthonormal to1.34e-14,
correct phases/centers/set order. The TZ binary energies agree with its existing
widest-space CI within4.974e-13Ha. All48 saved probe records and128 saved CI-vector
norms reconstruct; original action/subspace evidence is authenticated, not
recomputed. Original import/seed-preservation failure1 and distinct-frame
qualification0 remain separate; a packaging-only reader failure1 is preserved
beside dereferenced-package completion0.

**Select eight roots in every sector,64 total, in each model.** This uses the
maximum already available common native inventory without additional target
solves and avoids knowingly returning to the failed first-five stretched
prescription. Eight is still not nine and is not certified complete. The future
question is explicitly finite-model DZ/TZ sensitivity. A rejection would rule
out declaring this finite matched model basis-insensitive; a pass would support
only this local64-state contrast, not electronic accuracy or production promotion.
No outcome alone overturns stretched channel rejection or closes compressed
basis uncertainty. Ranked roots specify each model's truncation; no cross-basis
state-identity assertion is made.

Existing property counts are DZ `6,5,5,6,5,6,3,4` and TZ
`6,4,5,6,5,7,3,4` in native set order,40 states each,827/708 multipoles. Each
misses24 selected states. **Two64-state DENPROP exports are necessary missing
inputs**, not repeats of successful identical property calculations. Reuse exact
checkpoints/Molden, all64 CIDATA vectors, target integrals and eight CSF files per
frame; only NTGTL becomes `8,8,8,8,8,8,8,8`. Original DENPROP calls took4.65/4.40h.
The proposed conditional6–16h/export and18h native/19h owner ceilings are a
conservative planning envelope, not measured64-state costs or a threading gain.
Eight physical cores/64GiB follow the large-job policy. No resources are reserved.

**Queue now:** (1) saved-input feasibility complete; (2) freeze a property-only
contract/tests/exact owner with input pins, runtime floors and stops; (3) only on
separate release, compute/independently qualify the missing properties; (4) then
review minimum preparation and one B1 comparison. Step2 has at most60minutes
narrow implementation/freeze effort, zero scientific calls. No optimizer/CI
retuning, new target import, common-frame attribution project, root/channel
ladder or general solver API. A failed property prerequisite/job must block the
other model's launch and any successor. No automatic expensive run follows this
saved-data completion. All protected references/holdouts and scientific budgets
remain unchanged.
