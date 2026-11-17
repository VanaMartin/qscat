# CO: two fitting-data goals

## Goal table

| Goal | Data product | Reuse now | Missing output / completion test |
|---|---|---|---|
| **Neutral potential** | Own ground-state samples `{R, E0, method, state, zero, uncertainty}` → fitted `V0(R)` | Prepared homogeneous 5Z frozen-core RHF/CCSD(T) well and accepted molecular-zero local polynomial, 1.9–2.5 bohr; QZ kept separate | First local data/representation milestone complete; no global extrapolation, physical `De`, dissociation or factory-model export. Reproduction and original launch failure are in the [catalogue](experiments.md). |
| **Resonant potential** | Own threshold-relative `{R, Er, Gamma, estimator, model, identity, uncertainty}` → resonant effective-potential fit | Ten CAS11/cc-pVDZ/nine-A2/2048 paired geometries; eight development/two prospective validation rows, original failed scores preserved | Both frozen C2 tests pass, including public replay. Scoped local representation delivered alongside neutral; stop compute. No uniform/electronic/global accuracy or unique N+1 continuation claim. |

These are two inputs to a fitted potential, not a requirement to obtain somebody
else's potential. Local finite-model data, fit quality, electronic convergence
and physical prediction are different deliverables. A useful local fit does not
complete the global resonance/bound-anion or full-scattering objective.

This plan supersedes the earlier reference-first ordering and blanket bans on
all fitting/data collection. External author data, exact atomic affinity,
full-scattering embedding and a separate doorway/NRM Hamiltonian are **not
prerequisites** to these local data products. Original failures, accepted results
and validation exclusions are unchanged. The user approved these local routes.
Actual jobs/fits still need their separately frozen bounded releases; this
standing plan is not blanket execution or protected-data authority.

## Progress checkpoint — 15 November 2026

| Milestone | Current position |
|---|---|
| Local neutral data and representation | Complete: homogeneous5Z table and accepted molecular-zero local fit; QZ/MR kept separate. |
| Compatible resonance rows | Ten paired geometries accepted; original cubic/PCHIP scores retained. Eight development and two passing prospective C2 tests remain separate; no refit or successful-native repeat. |
| Reusable target inputs | Exact saved CI/checkpoints and retained native-energy inventory authenticated, including public replay. |
| Physical target-metric coverage | Complete:16 missing sector matrices pass numerical checks, eight qualified final-gap matrices reused. |
| Common-frame target/pole correspondence | Unresolved: weak compressed link and direct interior retained/buffer exchanges; neither automatic pole rejection nor a new testing release. |
| Resonance representation | Complete for scoped local use: both frozen C2 position/full-width tests pass, including public replay. Original candidate bytes unchanged; no uniform/electronic/global bound. |
| Electronic/global/dynamical accuracy | Unqualified; these later claims are not prerequisites for finishing the declared finite local dataset. |

The [canonical local-scope conclusion](../../docs/physics/co-results-and-rejected-routes.md#local-representations-accepted-scope)
states the completed representation decision and unchanged limitations, with
authenticated delivery/original-history pointers. Use the
[saved evaluator](README.md#using-the-local-representations), not a new acquisition
or validation-driven refit. All accepted rows, failed predictions and original
process failures remain preserved. The [current handoff](CONTINUATION.md) stops
compute; any different physical claim needs a separate scope decision.
Reference-first ordering concerns unqualified global claims, not these completed
local goals. No historical work package releases a new job.

## Computed neutral curve: immediate scope

Compute-first means `V0(R) = E0(R) - E_ref`, using our own energies. Use the
already-paid 5Z well as the primary local dataset; keep QZ separate. The 52
accepted QZ/5Z energies cover **26 distinct geometries**, not 52 geometries.
Keep absolute totals and each basis's own molecular reference at R=2.1323 bohr.
The existing 14-knot/12-midpoint tests measure representation, not independent
physical accuracy; these historical data are not globally blind.

**Energy-zero gate:** the current factory's EMO has a zero dissociation asymptote
and no additive offset in `fit_neutral`. Do not feed molecular-reference-zero
energies as though they were atom-relative. Use an explicitly offset/local
representation or table and defer physical `De` and the asymptote. Alternatively,
any computed fragment zero must match the molecular method/basis/core treatment;
MR atom totals cannot align the CCSD(T) well.

The first local representation now uses a fixed degree-6 Chebyshev polynomial
constrained to the molecular reference, with original knots and historical
midpoints kept separate. Its representation checks pass; no electronic accuracy
or asymptotic claim follows. Tables, coefficients, exact analysis and controls
are external records in the catalogue, not constants promoted to `qscat`.

The nine accepted all-electron aug-cc-pCVTZ CASSCF/SC-NEVPT2 endpoints are a
**separate** dataset: two atoms and seven molecular points, including R=3 and 4.
They preserve a finite bridge, not qualified dissociation. Their local-well
physical diagnostic failed; the R12 orbital optimization remains unqualified.
Do not splice CC and MR energies or repeat accepted endpoints. If extended-R
neutral coverage becomes essential, the bounded contingency is one missing MR
overlap point at R=2.5 under the accepted prescription, not another R12 retry.
It needs its own release and does not establish a global curve.

See [well evidence](qualification-evidence/neutral-near-equilibrium/README.md),
the latest neutral qualification in the [catalogue](experiments.md), and
[results and rejected routes](../../docs/physics/co-results-and-rejected-routes.md).

## Delivered finite resonance prescription

Use the existing operational outgoing definition, not new pole software:

- CAS11/cc-pVDZ: ten active electrons, two inactive A1 core pairs, active
  irreps `[5,3,3,0]`.
- Equal-weight **40-component orbital objective**; retained targets
  `[5,5,5,9]` per spin, 48 total. Nine-A2 is not a new 50-component optimization.
- Radius 18 bohr, l≤4, deletion 1e−6, uncorrected 2048 inner roots; retain all
  finite closed-channel/exterior couplings and actual target frames.
- Initially local, positive-energy scope within 1.9–2.5 bohr, only where the
  isolated branch can be followed. Threshold collision, lost isolation or
  ambiguous character stops/splits the segment.

This prescription maximizes compatible reuse; it is **not a converged physical
CO model**. Compressed active-space and separate stretched DZ/TZ sensitivity
tests failed their physical comparisons. Complex omitted-spectrum uncertainty
and multi-R identity remain unresolved. B1/B2 counterparts are identity checks,
not extra geometry nodes. Native real-axis fits, CAS12 poles and other bases
must remain separately labelled.

Keep the full [effective-target contract](../../docs/physics/co-effective-resonance-target.md),
[compressed reference](COMPRESSED_POLE_REFERENCE.md),
[paired character evidence](PAIRED_POLE_CHARACTER.md) and
[existing electronic prescription](a2-compressed-scattering-contract.json).
The latter is evidence, not renewed execution authority.

**Bounded alternative, not a parallel campaign:** a separately adopted native
RESON/Breit–Wigner dataset could use a single prospectively frozen same-model
recipe and its own window/background/grid adequacy gates. Audit its allowed
compressed development comparison first; stop if inadequate. Existing rejected
windows stay rejected. Such rows are estimator-defined, never relabelled outgoing
poles, and do not silently replace the outgoing target or lower its tolerances.

## Energy and state contract

Use Hartree and bohr internally; Gamma is the **full width**:

    z(R) = E_total,pole(R) - E_ground,native(R) = Er(R) - i Gamma(R)/2.
    Vres_hybrid(R) = V0_QC(R) + Er(R).

Keep individual native ground totals, total pole energies and incident energies.
Never use the state-average ensemble energy as a ground energy. A separately
correlated neutral plus native resonance is a **declared hybrid**, not an absolute
correlated anion calculation. UKRmol+ constructs target/scattering Hamiltonians
internally; `neutral.py` is a separate PySCF/Psi4 correlated-QC route.

Preserve the historical requested-eV grid conversion `E_inc=0.0735*E_requested/2`
separately from physical eV presentation/tolerances using `qscat.units`.
For a later dissociation claim, use matched `A0=E_C+E_O`, `AO=E_C+E_Ominus`,
`AC=E_Cminus+E_O`, so atom-relative anion limits are `-EA_O` and `-EA_C`.
State spin-free versus lowest-J conventions and core/charge/spin treatment.
`D0=De-ZPE`; a ground-fragment attachment threshold is `D0-EA`, not a fixed-R
resonance or an atom-to-fixed-molecular-reference offset. The
[atomic protocol](../../docs/physics/co-atomic-energy-references.md) governs those
later comparisons; they are not entrance gates to a molecular-zero local fit.

The current `ResonanceTarget`/fit API couples local data to `ea`, `R_inf` and
asymptotic setup. Store the local table first. A later local-fit interface must
explicitly defer unavailable global constraints, never invent constants.

## Beyond the local interval: crossing coverage

The accepted upper limit **R=2.5 bohr is a data/representation boundary, not a
physical endpoint** for resonance-to-bound continuation. That objective needs
compatible coverage on both sides of detachment, including a normalizable N+1
state below the matched neutral threshold with resolved binding uncertainty.
Re z crossing zero, an unresolved width, or intersecting independently zeroed
curves alone does not qualify binding. In the declared hybrid,
`Vres_hybrid−V0_QC=Er_native` by construction; a plotted crossing is not independent
electronic validation. Fixed-R binding also does not establish attachment dynamics.

Plan a small state-tracked crossing bracket before a dense extended curve;
do not defer this regime behind an unbounded local-accuracy campaign. Model
sensitivity and crossing coverage are separate decisions. Retain the frozen
1.9–2.5 evaluator without extrapolation or refit. The current outgoing map
rejects crossing its reference threshold interval: the physical bound-side
prescription and its controls must be qualified separately, not obtained by
removing that guard or clipping Gamma to zero. New geometries require their own
targets/thresholds and a bounded costed release. Extending the hybrid also needs
qualified neutral coverage; separate MR endpoints cannot be spliced into the CC
fit. Internuclear R is distinct from the electronic sphere radius a and exterior
endpoint; changing one does not automatically justify changing the others.

The separately bounded two-geometry pilot now supplies compatible finite-model
resonance-side and localized bound-side endpoint pairs; see the
[terminal handoff and evidence](CONTINUATION.md). That scoped pilot is complete
and acquisition stopped. These endpoints do not prove unique N+1 continuation,
locate a continuous crossing or extend the delivered evaluator. Further state/
model/reference/joining questions need a distinct bounded scope decision.

## Validity gates by stage

The operative checklist is [sweep readiness](SWEEP_READINESS.md):

1. **Accept a row:** immutable provenance, actual protocol/state/zero, final
   applicable residual/spin/density/orbital/import checks and extraction guards.
2. **Call it a curve:** compatible protocol and target/channel bookkeeping;
   physical target/subspace and pole/residue assignment in a common frame.
3. **Fit it:** actual sufficient nodes, bounded parameter count, declared zero,
   domain, numerical uncertainty and representation tests; freeze exclusions.
4. **Claim physical accuracy:** quantify electronic-model sensitivity and matched
   physical comparisons. A failed insensitivity test does not delete a valid
   finite row or prohibit all local collection.
5. **Claim global/dynamical predictions:** qualify fragments, dissociation,
   bound joining, tails and nuclear observables; full S additionally needs a
   separately specified channel/background model.

Existing numerical and physical budgets retain their original scopes; no
lower-accuracy production or exploratory tier is silently adopted.

## Completed local milestones and future decision boundary

The first-dataset milestones below are complete, not a current acquisition queue.
Before any new task, name the distinct physical decision it can change; more
interpolation points alone do not qualify electronic accuracy. Documentation
integration and the later separately released two-endpoint pilot are complete;
the terminal handoff, not historical proposals, governs any next scope decision.

| Order | Task | Completion test / bound |
|---|---|---|
| 1 | Prepare saved-data handoff — complete | Manifest-backed 5Z/QZ tables with totals/reference/status, separate nine-endpoint MR index and resonance-input ledger are published/authenticated. Compatible stretched CAS11 boundary inputs exist; no target/property/inner repeat is needed for that outgoing determination. |
| 2 | First missing-node release — complete | R=2.5 uses existing native boundary inputs, with no target/property/inner repeat. Both outgoing sectors, real guards and local counterpart tests pass; saved-data reconstruction preserves original source/input/runtime stops and numerical limits. |
| 3 | Deliver the scoped local representation | Complete: ten pairs, eight development/two prospective C2 PASSes and accepted neutral fit; public delivery verified. Stop new compute. Preserve R=2.1323/all114 S exclusions; no uniform/electronic/global claim. |
| 4 | Fit locally when each output suffices — complete | Local neutral/resonance representations and their evaluator are delivered; global factory/asymptotic integration remains deferred. Four actual compatible resonance training geometries are a floor, not proof of identifiability. Default `n_nodes=24` is not a 24-geometry acquisition requirement. |
| 5 | Freeze, then test reduced prediction — scoped design complete | Both prospective C2 tests pass with eight development nodes and no refit; original cubic/PCHIP rejections remain. The reserved equilibrium outgoing pole requires a separate later release and remains unread, as do all114 S holdouts. Selected-point passes are not a uniform error bound. |

Unknown whole-node costs must be resolved from saved inputs before launch.
Published outgoing-only timings exclude target/import/properties/inner work;
do not use them as new-geometry forecasts. Prefer four physical Sadaharu cores
for small jobs; eight/twelve require justification. No paid provisioning, image
rebuild or unrelated workload change by default. Active jobs target 600-s checks
and immediate terminal observation; no idle perpetual monitor.

Review at every deliverable and within one or two active working days. A
first-failure stop ends that release, not the scientific goal. New reusable code
is tested/committed/pushed before use; exact executed inputs/source/owners and
failures live externally under the [planning contract](../../docs/experiment-planning.md).
No successful electronic calculation is repeated merely to accumulate passes.

## Blind alleys and literature basis

Use the [no-repeat/reopening ledger](../../docs/physics/co-results-and-rejected-routes.md#no-repeat-and-reopening-conditions),
not historical “next” entries. External-author waiting, atomic precision ladders,
R12 restarts, full-S/doorway work and general exporters do not unblock the local
tables. Nine remains working, new50 parked, 2048 retained; equilibrium outgoing
validation and all 114 S holdouts stay protected. Cleanup/merge remain separate.

The committed reference notes support the relevant cautions:
[Dora 2016](../../reference/literature/dora-2016-epjd70-197.md), p. 6, Table 4,
documents model dependence and nonmonotonic improvement with larger CAS;
[Dora 2020](../../reference/literature/dora-2020-jpb53-195202.md), pp. 2–3 and
5–6, documents geometry/state-tracking difficulties and fit-sensitive widths.
Neither is a matched absolute-accuracy benchmark for our outgoing estimator.
Conventional eigenphase fits are not automatically outgoing poles; small atom
count alone does not diagnose a package bug. Additional reference-note work is
separate from this CO integration.

Detailed research, provenance and superseded plan snapshots belong in the
[external catalogue](experiments.md), not another chronological launch ledger.
