# Compressed CO A2 target inputs and matched scattering

**The compressed five-versus-nine-A2 control passes all declared observable
budgets. Independent local, host and public reconstruction confirm the result.**

## Completed scattering — 17 October reading

Owner 1365000 finishes **16 October, 16:55:18.733 UTC**, with all seven owner
steps, native stages, resource wrappers and comparison exits zero. Observer
1365474 ends normally after nine captures. Both services are inactive. The
[completion capture](progress-a2-compressed-scattering-20261017T0840.json)
reverifies 61 source / 59 input hashes and all 242 installed engine artifacts.

| Comparison | Maximum phase difference, modulo π | Position change | Full-width change |
|---|---:|---:|---:|
| Dense replay → fresh five-A2/2048 | 2.2e-6 rad | ≤0.007205 meV | ≤0.02509 meV |
| Five-A2 → nine-A2/2048 | **0.0115944 rad** | **13.241–18.210 meV lower** | **2.218–4.404 meV smaller** |
| Original scientific budget | 0.05 rad | 50 meV | 5%, with 1 meV floor |

Both sectors agree; the channel phase maximum is at **3.9 requested-eV**, with
no grid point exceeding budget. All four background orders pass in the declared
2–5 requested-eV window. Native residuals are at most **1.37e-10 Hartree**;
selected Pi splitting is zero. Complete exports confirm dimensions 27546/27818,
continuum 1410/1682 and boundary 250/298. Actual shared input hashes and each
model's own noncontiguous channel map/deck verify.

Measured five/nine costs are **1.651263 / 1.827292 hours**,
**6.011147 / 7.245762 GiB**. Dense outer replay is **4.248 minutes / 0.446 GiB**.
The [73-file reading](verification-a2-compressed-scattering-20261017.json)
reconstructs exports, telemetry, grids/Pi, maps/decks and both fit gates.
[Public replay](qualification-evidence/compressed-a2-scattering/README.md)
verifies eight payloads and 58 frozen source files. Host replay takes 0.815 s;
largest host/local fitted-parameter difference is 4.018e-9 eV. Four actual dense
inputs match the frozen complete 27546-pole source. Large shared native input
equality is checksum-audited remotely. Shared readers/fitter are not an
independent scattering engine.

## Three-anchor decision

| Anchor | Five→nine-A2 phase change | Outcome |
|---|---:|---|
| Compressed, 1.9 bohr | 0.0115944 rad | Independently/publicly reconstructed pass |
| Equilibrium, 2.1323 bohr | 0.0163665 rad | Independently/publicly reconstructed pass |
| Stretched, 2.5 bohr | 0.0679178 rad after named nine-A2 refinement | Reject |

This completes the missing anchor comparison. Five-A2 retention is not uniformly
insensitive over the domain; the stretched rejection persists. Nine-A2 is the
larger tested representation for subsequent matched controls, with full channel
convergence still unset. Passing compressed/equilibrium controls do not justify
extra root expansions. The established economical setting remains 2048 roots;
the retained refinement evidence is a measured sensitivity, not an omission bound.

The fresh equilibrium CAS11 baseline also enables a
[same-engine saved-output check](qualification-evidence/same-engine-active-space/README.md).
CAS11→CAS12 still rejects phase/position at **0.1212798 rad / 63.652–73.447 meV**;
CAS12/4096 preserves rejection. Old and fresh selected CAS11 phases are identical.
The scattering-engine/rank confound is removed, while target diagonalizer/import
and optimized-frame differences remain. The original strict deck-byte preflight
fails on two formatting-only lines; its failed/unset record is preserved beside
the separately declared supplementary exact-line audit and model rejection.

**Next decision:** validate the extraction estimator against known-pole phases
with smooth nonconstant backgrounds before selecting a new estimator or using
fit widths to rank expensive electronic models. Its design must freeze the
model-selection rule and unused-data tests before inspecting their outcomes;
all historical order/window rejections remain in force. In parallel planning,
use the now measured channel effects and costs to specify independent 40→50
orbital-average and basis controls, qualifying native target imports before any
matched scattering. Model/continuity/neutral-domain and restartable-pilot gates
still precede the production sweep and potential fitting.

Equilibrium five→nine-A2insensitivity passes at0.0163665rad,while stretched
insensitivity rejects at0.0679178rad after its named root refinement. The
compressed anchor completes the missing geometry-dependent channel comparison
before further three-anchor model/basis/ensemble choices.

## Qualified input records

The [target declaration](a2-compressed-target-preflight-contract.json) uses the
original compressedCAS11/cc-pVDZcheckpoint at1.9bohr,two frozen orbitals,ten
active electrons and forty-component state averaging. Explicit retention is
[5,5,5,9]per spin. Owner1360012completes **08:51:07.371UTC,16October**,with
all eight target solves,DENPROP,resource wrapper and owner exits0.
Observer1360347ends normally after3captures.

The [44-file reading](verification-a2-compressed-target-20261016.json) checks
48fullCIDATA vectors/Gram matrices,independent twelve-root energies/actions/
spins,original first-five preservation,Pi,915unique finite multipoles,dipole and
compressed own-model map. Singlet/triplet A2root9IDs are48/45.

| Check | Maximum raw result | Budget |
|---|---:|---:|
| Reference-energy difference | 8.669e-13Hartree | 1e-7Hartree |
| Native residual | 1.319e-10Hartree | 1e-9Hartree |
| Vector Gram error | 1.077e-14 | 1e-9 |
| Pi target-energy split | 5.684e-13Hartree | 1e-7Hartree |
| Property/binary energy difference | 4.818e-10Hartree | 1e-7Hartree |
| Ground dipole difference | 1.359e-11au | 1e-5au |

Cost is **51.178minutes /5.477GiB**,including47.068minutesDENPROP and
246.065seconds target solves. Checkpoint,Molden,integrals/CONGEN hashes remain
unchanged;61source/50input hashes and242installed artifacts reverify.
[Public replay](qualification-evidence/compressed-a2-target/README.md) passes
6payloads/3reader sources/44raw files. Host/macOS quoted metrics agree except
for8.882e-16in Gram error. Full native input reconstruction does not independently
validate every added transition moment or qualify full channel convergence.

## Historical matched-scattering declaration and launch state

The [separate scattering declaration](a2-compressed-scattering-contract.json)
freezes the same new48-state target representation,original compressed orbitals/
integrals/CONGEN,audited8-rankengine and2048scattering roots on both sides.
Only explicit A2retention5→9,contraction sequences and own-map interface change.
Expected dimensions27546/27818,continuum1410/1682,boundary250/298must be measured.

Common99-point0.1–5requested-eVgrid,radius18,maxl4,deletion1e-6,propagation100,
solver1e-12and native6GiBmemp remain fixed. Fits retain all four background
orders in the prospectively fixed **2–5requested-eVwindow**. Saved dense phases
show valid fit conditioning before native launch;the narrowed upper endpoint
relative to historical2–6preserves the established grid. This new matched
sensitivity does not replace the24original extraction cases or rescue their
**32.206%width-span rejection**. Actual fit energies are requested-eV×0.0735/2
Hartree. Budgets stay0.05rad moduloπ/0.05eVposition/5%fullwidth with0.001eVfloor.

At [15:14UTC](progress-a2-compressed-scattering-20261016T1514.json),owner1365000
and observer1365474are active. Runtime/input gates and dense replay pass.
Dense outer replay cost **4.248minutes /0.446GiB**;an
[eight-file independent read](verification-a2-compressed-dense-20261016.json)
checks complete99-point grids,Pi1.0e-8rad and nonempty/noMAXFIT native fits.
No dense Hamiltonian solve is repeated.

The [37-file independent baseline reconstruction](verification-a2-compressed-five-baseline-20261016.json)
confirms the fresh five-A2selected run: **1.651263h /6.011147GiB**,raw exports,
zero selected Pi split,own map/decks and all eight matched fits pass. Maximum
dense-reference phase difference is **2.2e-6rad**,position7.205e-6eV and
width2.509e-5eV. Supplementary fit-platform drift4.018e-9eV stays below1e-8eV.
The controller consequently starts the declared nine-A2model once.

Source is frozen at **`89a5275fdf5e55776dbfbe534793a974b493ed61`**,with61source/
59input hashes and242installed artifacts reverified. Allocation:8physical
CPUs4–11/eight MPIranks/one library thread/48GiBcgroup/8GiBshared;CPUs0–3
remain for small checks. Prelaunch persistent disk46.78GiB and available RAM
107GiB meet bounds. Forecast **2–4hours/model,4–8hours total plus replay/reading**
uses measured equilibrium and stretched costs. Hard stage8h/model12h/owner28h
limits are bounds,not predictions. Empty/MAXFIT/failed/bound fits stay unset;
false/unset stops without retries,relaxation or automatic successor.

Artifacts:`prepared/a2-compressed-scattering-source-20261016/`;services:
`co-a2-compressed-scattering-20261016.service`and
`co-a2-compressed-scattering-observer-20261016.service`. The observer records
immutable snapshots every30minutes after readable-owner/live-service PID match.

The completion and decisions above supersede this historical live state.
