# Equilibrium CO nine-A2 input gate and matched scattering control

**The matched equilibrium five-versus-nine-A2 control passes its declared
observable budgets. Independent local,host and public raw reconstructions
confirm the pass.** The stretched control still rejects phase insensitivity.

## Completed scattering — 16 October

Owner1351889completes **02:59:15.260UTC**,with all native,resource,reader and
owner exits0. Observer1352372ends normally after9captures. Runtime stays frozen
at`a5954fb`;60source/61input hashes and242installed artifacts reverify.

| Comparison | Phase change moduloπ | Fitted position change | Full-width change |
|---|---:|---:|---:|
| Dense replay → fresh five-A2/2048 | 5.1e-6rad | ≤0.000767meV | ≤0.012905meV |
| Five-A2 → nine-A2/2048 | **0.0163665rad** | **9.427–10.924meV lower** | **2.672–3.645meV smaller** |
| Scientific budget | 0.05rad | 50meV | 5%with1meVfloor |

Both sectors give the same channel comparison,with the phase maximum at
**2.55requested-eV**,no grid points exceeding budget,and zero selected-model Pi
phase split. All orders1–4pass position/width comparisons in the frozen
1.8–3.3requested-eVwindow. Native selected residuals are at most1.301e-10Hartree;
complete exports confirm dimensions27546/27818,continuum1410/1682and boundary
250/298. Actual target/integral/CONGEN inputs match across models and actual
decks match their own noncontiguous channel maps.

Measured costs:five **1.656076h /6.407837GiB**,nine **1.817478h /6.896057GiB**;
dense outer replay **4.198minutes /0.446373GiB**. The
[73-file raw reading](verification-a2-equilibrium-scattering-20261016.json)
recomputes complete grids,CI/boundary exports,telemetry,Pi,decks/maps and both
fixed-window gates. [Public download/replay](qualification-evidence/equilibrium-a2-scattering/README.md)
passes8payload digests and58frozen reader sources. Host read costs1.167seconds;
maximum host/local fit difference8.293e-10eV. Four actual dense inputs match the
frozen complete27546-pole source;large shared input equality is checksum-audited
remotely. Shared readers/fitter are not an independent scattering engine.

**Decision:** equilibrium practical five→nine-A2insensitivity passes at2048.
It does not establish full channel convergence,absolute omitted-spectrum error,
global state/pole identity or extraction reliability. The stretched phase
rejection0.0679178rad survives its named refinement,and the saved equilibrium
CAS11→CAS12position/phase rejection remains. The
[compressed-anchor control](A2_COMPRESSED_CONTROL.md) measures the missing
geometry dependence before choosing further anchor models. No extra equilibrium
root expansion is justified by this passing comparison.

## Completed input gate

The [target declaration](a2-equilibrium-target-preflight-contract.json) retains
the original equilibrium CAS11/cc-pVDZ orbital checkpoint and forty-component
state average at2.1323bohr,two frozen core orbitals and ten active electrons.
Explicit target counts are[5,5,5,9]per spin,total48. Native owner1345746
completes at **15October22:15:40.865UTC**,with all eight target solves,DENPROP,
resource wrapper and owner exits0. Observer1346072ends normally after3captures.

| Check | Independent raw result | Budget |
|---|---:|---:|
| Complete native vectors | 48 across all8sets | Exact declared counts |
| Maximum reference-energy difference | 5.969e-13 Hartree | 1e-7 Hartree |
| Maximum native residual | 1.20921e-10 Hartree | 1e-9 Hartree |
| Maximum vector Gram error | 1.077e-14 | 1e-9 |
| Maximum Pi target-energy difference | 3.126e-13 Hartree | 1e-7 Hartree |
| DENPROP versus binary energy | 4.996e-10 Hartree | 1e-7 Hartree |
| Ground dipole difference | 7.984e-11 a.u. | 1e-5 a.u. |
| Unique finite dipole/quadrupole records | 930 | Valid own-state IDs/irreps |

The [44-file local read](verification-a2-equilibrium-target-20261016.json)
checks full binary vectors/headers/exhaustion,Gram matrices,independent
twelve-root actions/spins,first-five preservation,Pi,telemetry and property
inventory. Singlet/triplet A2root9are IDs48/45;each model uses its own
noncontiguous map. The original checkpoint,Molden,integrals andCONGEN hashes
remain unchanged;59source/50input hashes and242installed artifacts reverify.

Measured native cost is **51.817minutes / 5.425GiB**,including47.675minutes
forDENPROP and4.135minutes for target solves. The initial supplementary reader
rejects a two-atom center assumption;raw CIDATA includes carbon,oxygen and
the zero-charge scattering center. The preserved first reader and corrected
three-center check change no native outcome or numerical tolerance.

[Public download and reconstruction](qualification-evidence/equilibrium-a2-target/README.md)
passes seven payload digests,three frozen reader sources and all44raw files.
Host/macOS decisions match;the only quoted metric difference is1.443e-15in
maximum Gram error. The reader also reproduces the earlier stretched24-file
48-state/921-multipole gate. This certifies these native input records,not
independent accuracy of every transition moment or full channel convergence.

## Historical declaration and launch

The [scattering contract](a2-equilibrium-scattering-contract.json) freezes:

- The exact same equilibrium48-state target vectors/properties,original
  integrals/CONGEN,orbital frame and forty-component averaging in both models.
- Only retained A2channels change5→9per spin. Target contraction counts,
  sequences and own-map outer interface change accordingly. Expected contracted
  dimensions27546/27818,continuum1410/1682,boundary250/298must be measured.
- Both models use2048scattering roots,selected auditedO1engine,eight MPI ranks,
  CPUs4–11,one library thread,native6-GiBmemp,48-GiBcgroup and8-GiBshared memory.
- Same99-point0.1–5requested-eVgrid,radius18/maxl4/deletion1e-6/propagation100.
  Actual energies use requested-eV×0.0735/2Hartree. Fits use the existing
  **1.8–3.3requested-eV window/background orders1–4** on both sides.
- Budgets remain **0.05rad phase moduloπ /0.05eV position /5% full width with
  0.001eV floor**. Empty/MAXFIT/failed/bound fits stay unset;false/unset stops.

The owner first replays the saved **complete dense boundary** through the outer
region,without a new dense Hamiltonian solve. A fresh five-A2/2048model must
pass dense-observable checks on the same runtime/48-state input representation
before the expanded model starts. The older equilibrium selected output remains
supporting evidence;its release-engine/four-rank runtime differs.

Same-engine stretched five/nine models cost1.707/1.850hours and5.894/6.307GiB.
Conditional forecast is **2–4hours/model,4–8hours total plus replay/reading**;
stage8h/model12h/owner28hlimits are hard bounds. No larger-root run follows
automatically. Original stretched channel/extraction rejections remain preserved;
global continuity,basis/active-space/ensemble,absolute omission and neutral
qualification stay open. Production remains gated.

At [23:28UTC](progress-a2-equilibrium-scattering-20261015T2328.json),owner1351889
and observation-only observer1352372are active. Runtime and input gates pass;
dense outer replay is active. Source is frozen at`a5954fb`,with60source/61input
hashes and242installed engine artifacts reverified. The independent target
report is bound to actualfort24/fort26bytes;tests reject false/unset reports,
changed maps and replacement native files. **112focused tests pass.**

Services are`co-a2-equilibrium-scattering-20261016.service`and
`co-a2-equilibrium-scattering-observer-20261016.service`;artifacts are under
`prepared/a2-equilibrium-scattering-source-20261016/`. The observer records
immutable snapshots every30minutes,starting only after the readable-owner PID
matches the live service. Next inspect dense replay/five-A2completion and the
independent baseline verdict;then the declared owner may start nine-A2once.

**23:31UTC update:** dense replay completes all six native outer stages and
resource wrapper with exits0,in **4.198minutes /0.446GiB kernel peak**. The
[eight-file independent reading](verification-a2-equilibrium-dense-20261016.json)
verifies both99-point grids,Pi difference1.0e-8rad and nonempty/noMAXFIT native
fits. These are the original dense-model candidates,not an expanded-channel
result. The [23:31capture](progress-a2-equilibrium-scattering-20261015T2331.json)
records fresh five-A2/2048B1assembly at801.17%CPU/2.514GiB current container
memory. Dense-baseline and channel verdicts remain unset until comparison.
