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
