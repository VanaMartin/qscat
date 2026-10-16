# Compressed CO A2 target inputs and matched scattering

**The compressed48-state target/property gate passes independent and public
reconstruction. The fresh five-A2/2048baseline passes its dense-observable gate;
the declared nine-A2/2048model is running on Sadaharu.**

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

## Matched scattering and current state

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

**Next:** checksum-fetch and reconstruct completed nine-A2exports,grids,Pi and
both-sector channel comparison. Decide from measured effects/cost;retain original
equilibrium/stretched outcomes and extraction failures. Model/basis/ensemble,
global continuity,absolute omission,extraction,neutral correlation/domain and
restartable-pilot qualification still gate the sweep and potential fitting.
