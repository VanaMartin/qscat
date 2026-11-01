# Finite stretched CAS12 atomic-basis control

## Decision served

Compare the already qualified cc-pVDZ and cc-pVTZ targets at R=2.5bohr, keeping
the forty-component equal-weight orbital objective, two frozen core orbitals,
ten active electrons and6/3/3/0 active orbital counts. Both native files contain
eight roots in every spin/irrep. The missing all-eight64-state property tables
now independently qualify: [DZ](../../projects/ukrmol_co/qualification-evidence/stretched-basis-properties-dz/README.md)
and [TZ](../../projects/ukrmol_co/qualification-evidence/stretched-basis-properties-tz/README.md).

The comparison changes the atomic basis **and its optimized orbital frame**;
it is not a perturbation on one frozen common orbital representation. Selecting
the same number of roots does not establish cross-basis state identity. The
finite64-state prescription maximizes available matched reuse but is not nine-A2,
omitted-channel completeness or production electronic convergence. Preserve the
stretched channel rejection and compressed CAS11/CAS12 pole-position rejection.
No successful QC/CI/import/property calculation needs repetition.

## Minimum new preparation

The saved target-only integrals omit the scattering continuum. Inventory finds
no complete full-continuum/B1 preparation on either selected frame. Generate only
one four-rank continuum integral set and one serial B1 N+1 CONGEN expansion per
basis under the [exact preparation contract](../../projects/ukrmol_co/stretched-basis-preparation-contract-20261031.json).
Copy the original eight target CSF files, CIDATA, Molden and checkpoint; copy the
new own64-state properties. Original inputs/source stay read-only. Do not rerun
target CONGEN, target diagonalization, orbital optimization or DENPROP.

Reuse the validated geometry-free CAS12 cards at radius18bohr, max_l4 and deletion
threshold1e-6. Only full orbital addresses resolve from the actual inventory;
the future SCATCI card selects all eight roots per sector. Preserve the CONGEN
target/L2 group definitions byte-for-byte. The SCATCI card is recorded, not run
by preparation. Derive the64-state target metadata from the own CIDATA energies;
retain the original forty-component orbital objective explicitly.

Require both independent public property receipts and exact input/source/engine
hashes before science, four-rank native integral inventory/error/completion gates,
nonempty integrals, and CONGEN's final projection/write-unit/raw-count/group gates.
Contracted dimension equals the actual L2 count plus the sum of eight retained
roots times each sector's actual radial continuum count; cap150000. The old
84942 L2/282 radial-sum anchor predicts87198, but no stretched dimension is forced.
Full-continuum creation alone does not certify native target-channel phase
handling, scattering interfaces or resonance poles.

## Bounds, serial release and progress checks

DZ preparation first; independently authenticate/publish its result before a
separately released TZ job. TZ requires a within-root hash-pinned independent
public DZ preparation success receipt. Each owner can make one integral and one
CONGEN call, never its peer or a numerical successor. Any failure stops without
repair/retry, retaining original exits, partial outputs and the failed capsule.

Four physical cores0–3/32GiB/one library thread/no network; minimum40GiB available
RAM/20GiB free persistent disk. Conditional3–20minutes per preparation is planning
only, based on the old compressed171s/10.70GiB result. Native integral40minutes,
CONGEN5minutes, container58minutes and durable owner60minutes are hard caps.
Source/test/owner freeze is bounded to60minutes implementation effort; no new
scheduler or solver framework. The shared dimensions reader gains an explicit
eight-root argument while retaining its five-root default and existing tests.

The user requests sustained work and regular checks: inspect each stage and at
most60seconds apart during preparation; any later long inner job needs at most
five-minute read-only status/resource checks. Resource sampling remains0.2seconds
memory/2seconds disk. Monitors may report status, never release successors or
hide failure. No recurring daily task, paid provisioning or unrelated workload
change is implied. User services survive SSH disconnect, not host shutdown or
all-session logout unless existing host configuration provides that persistence.

## What follows only after preparation qualifies

### Qualified DZ and separate conditional TZ reuse

[DZ preparation independently qualifies](../../projects/ukrmol_co/qualification-evidence/stretched-dz-preparation/README.md):
87198 contracted states,90.40s/8.288GiB, all original exits0. Full public raw
CSFs/integrals and independent inventories/counts/header authenticate. This does
not qualify CIRMAT channel phases, a scattering export or an electronic model.

The [separate TZ reuse release](../../projects/ukrmol_co/stretched-tz-preparation-reuse-release-20261031.json)
retains the executed DZ contract unchanged. CONGEN's occupation/spin/angular
algebra depends on numbered orbital addresses, its complete resolved deck and
the target-CSF definitions, not the radial Hamiltonian or optimized coefficients.
All eight original target-CSF hashes agree. After generating own TZ integrals,
require exact native inventory and full resolved deck equality, then copy only
the qualified B1 CSF definition and original DZ log, explicitly marked reused.
Identical algebra is not identical wavefunctions or cross-basis target identity.
No integral/CIDATA/property/Molden/checkpoint reuse across bases is allowed.
Mismatch stops without fallback CONGEN/retry. This avoids repeating a successful
identical algebraic calculation. Own TZ native/scientific exits and resource
costs remain separate from the original DZ CONGEN's exit and cost.

TZ [now also independently/publicly qualifies](../../projects/ukrmol_co/qualification-evidence/stretched-tz-preparation/README.md),
38.67s/4.174GiB, all exits0, same87198 contracted states. The exact conditional
reuse gates pass, and full public raw data/predecessor authenticate independently.

Review actual dimensions/cost and separately freeze a B1/2048 inner/export test,
then its counterpart in the other atomic basis. Preserve2048 as baseline; expand
only for a named discrepancy. Numerical boundary/eigenpair/own-K/export gates
must precede a coupled outgoing-pole comparison. Use each basis's own target
threshold/CF matrices and incident energy zero; compare physical eV position and
**full width**, not total molecular energies or raw native fit widths.

The later model question is sensitivity within50meV position,5% full width with
1meV floor and0.05rad phase. Pole/character/local-identity and numerical errors
must be qualified separately; an unstable or wrong-branch result is inconclusive,
not a basis verdict. One qualified B1 rejection is decision-changing; B2 is not
launched to rescue it. A pass alone does not promote either model or justify a
production sweep. No inner solve, pole search or B2/fit/sweep is released by this
preparation document. Keep nine working reference, parked new50, equilibrium pole
exclusion and all114 untouched S holdouts.

### Separate minimum own DZ B1/2048 allocation

Both preparations now qualify, so the [exact DZ-only inner/export contract](../../projects/ukrmol_co/stretched-dz-inner-contract-20261031.json)
selects one missing own native solve, after new source/tests/exact owner freeze
and push. Reuse all own DZ prepared inputs and the existing fully qualified
audited engine, with its same-image numerical-library/telemetry/installed-artifact
checks. Extend only explicit finite layout arguments: eight target roots per
sector,64 targets/400 boundary channels/2256 contracted continuum coefficients;
preserve five-root40/250 defaults. No new native solver or generic selection layer.

Exactly one twelve-rank B1/2048 SCATCI, then native RSOLVE/EIGENP/RESON for own99
real-K references at requested0.1–5.0eV; RSOLVE propagates to100bohr. The lowest
own DZ/TZ excited thresholds lie above5eV, so the six-open-channel reader remains
applicable only if actual emitted metadata confirms the closed gap. The native
RESON fit is diagnostic, not the coupled outgoing pole. No complex/Python map,
new target/property/preparation, B2 or root/channel ladder is released.

Twelve physical cores0–11/64GiB/8GiB shared backing;80GiB available RAM/20GiB disk
floors. Actual87198 dimension is0.98% above the compressed86352 case, but stored
continuum/channel counts are60% larger; do not extrapolate cost from dimension
alone. Previous compressed twelve-rank B1 measured4.227h/22.49GiB. Conditional
4–10h here,12h native/20min per outer stage/13h durable owner; implementation
bounded to60minutes. Check at most five minutes apart while running.

Original input/deck/engine/native/residual/export/finite-array/closed-gap gates
remain strict. A failed gate preserves original exits/partial evidence and stops
without retry. Independently/publicly qualify DZ export and measured cost before
reviewing a separately frozen TZ inner release. No automatic TZ/pole/B2/fit/sweep;
no CIRMAT-phase, omitted-spectrum, electronic-accuracy or state-identity claim.

### Qualified DZ export and separate TZ counterpart

The [own DZ B1/2048 export now independently/publicly qualifies](../../projects/ukrmol_co/qualification-evidence/stretched-dz-inner/README.md):
all exits0,6.22242h/28.662GiB;2048 roots/400 boundary channels/99 K references,
recorded residual1.27737e-10Ha, CI/boundary energies exact, K/native phase error
4.97976e-8rad moduloπ. This supplies an input, not a physical model verdict.
All94 monitor checks meet300s cadence. Independent reader failures remain
separate from successful native exits; exact final-allocation reconciliation
does not loosen numerical budgets or repeat science.

The [separate TZ-only counterpart](../../projects/ukrmol_co/stretched-tz-inner-contract-20261101.json)
requires that exact independent public DZ receipt, both preparation receipts,
own TZ input pins, frozen own64-target ID map and unchanged audited engine.
Keep the finite64/400/2048 model, geometry2.5bohr and forty-component orbital
objective. Atomic basis and its optimized frame change together; identical root
counts do not establish target-state identity. All own integrals/vectors/
properties remain basis-specific; only the already proved algebraic CSFs reuse.

No numerical method changes: add only the explicit serial predecessor gate and
own paths/map. Eight mocked runtime controls/source/tests/exact owner/frozen deck
must push first. Twelve physical cores/64GiB/8GiB backing,80GiB RAM/20GiB disk
floors, conditional5–10h/12h native/20min per outer stage/13h owner;45minutes
implementation. At most5min read-only checks. Preserve any failed gate/partial
output without repair/retry. Public independent TZ qualification precedes a
separately costed outgoing-pole/character/model-sensitivity proposal. No pole,
B2, fit, sweep, root/channel ladder, electronic-accuracy or omitted-spectrum
claim follows automatically; all protected references/holdouts remain untouched.

### Both own exports qualify; smallest pre-complex guard

[TZ now independently/publicly qualifies](../../projects/ukrmol_co/qualification-evidence/stretched-tz-inner/README.md),
5.561856h/27.132191GiB, all exits0; own2048 roots/64 targets/400 channels/2256
continuum coefficients,87198 dimension. Residual1.28373e-10Ha, CI/boundary energies
exact,99 K/native phase errors<=4.92617e-8rad. Independent public bytes reproduce
the own native arrays without a producer import or native repeat. One SSH timeout
leaves a495.443016s successful-observation gap; preserve that operational exception.

Before any complex search, [four real-axis guards](../../projects/ukrmol_co/stretched-real-guard-contract-20261101.json)
test each model against **its own** native real K. For each fixed development
index60/80, assemble the surface matrix Q from its2048 own poles and amplitudes at
total energy E0+E_Ry/2, integrate the same Riccati/Jost exterior to boundary18bohr,
and construct physical real S with the existing flux-normalized open amplitudes.
Compare to S(K) of that own native reference, not to the other basis. Use both
own incident thresholds/ground zero/CF matrices; physical energy arithmetic stays
native Rydberg in the radial equations and Hartree in the inner denominator.

Reuse the already validated method at120bohr tail parameter/40degrees/order10/
rtol1e-12. Require1e-6 absolute S,1e-6rad phase-moduloπ,1e-6 Wronskian/unitarity/
reciprocity and1e-7 flux errors,90s/RHS50000 per map. Exactly two development
maps per model, DZ before TZ; no legacy route-holdout access. Freeze source,
two exact independent public receipts, raw models, array hashes, tests and owner
before computation.4 physical cores/8GiB/600s science/900s owner; conditional
30–180s,60min implementation. This narrowly parameterizes an experiment, not a
new solver or production API. First failure stops without numerical rescue.

Retain actual Q,Y,A,S,native-reference matrices and all diagnostics/resources.
Passing is necessary own-model numerical consistency before complex use, not
electronic accuracy, basis insensitivity, CIRMAT-phase qualification or an
independent exterior oracle. No stretched coupled pole has been identified; a
separate physically justified seed/box and open-Pi residue/character/local-identity
proposal remains necessary after this finite guard. Do not promote native RESON
fits to outgoing poles or silently assume a compressed-pole identity at2.5bohr.
No complex search, B2, fit or sweep is released by the real-axis guard.
