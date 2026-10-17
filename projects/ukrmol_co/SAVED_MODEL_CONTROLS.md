# CO model/extraction controls — 15 October 2026

**Equilibrium CAS11→CAS12 rejects practical observable insensitivity, and nine
A2 channels do not resolve the stretched width-extraction sensitivity.** These
read-only controls use completed native outputs and the accepted low-root
guidance. The [contract](saved-model-controls-contract.json),
[quoted verification](verification-saved-model-controls-20261015.json) and
[public raw reconstruction](qualification-evidence/saved-model-controls/README.md)
retain the complete scope and original records.

## Equilibrium active-space sensitivity

**17 October supplement:** the fresh audited eight-rank CAS11 baseline exactly
reproduces the old selected phases. The
[same-engine reconstruction](qualification-evidence/same-engine-active-space/README.md)
preserves both model rejections below, removing the scattering-engine/rank
confound. Target diagonalizers/imports and optimized frames remain distinct.
The original strict deck-byte preflight fails on two formatting-only lines;
its failed/unset record is retained beside the supplementary exact-line audit.
The original study and results below remain unchanged.

Both Pi sectors reproduce these results, using the same1.8–3.3 requested-eV
fit window and background orders1–4:

| Comparison | Maximum phase change | Center shift range | Full-width change range |
|---|---:|---:|---:|
| CAS11/2048 → CAS12/2048 | 0.1212798 rad | 63.652–73.447 meV lower | 16.684–26.661 meV lower |
| CAS11/2048 → CAS12/4096 | 0.1212919 rad | 63.653–73.447 meV lower | 16.718–26.677 meV lower |

Maximum phase differences occur at2.55 requested-eV. Phase and all fitted
positions reject the unchanged0.05-rad/0.05-eV budgets; every full-width
comparison passes5% with0.001-eV floor. The selected CAS11 dense-reference
and CAS12 root-refinement checks still pass. This is strong practical model
sensitivity,with absolute omitted-spectrum error and full electronic convergence
remaining unqualified.

The current reader verifies124 checksum-fetched raw files,all original stage/
resource exits,complete binary CI/boundary exports,residuals,grids/Pi,and each
model's own noncontiguous DENPROP map/energies. Actual outer decks match byte
for byte. Common settings are equilibrium2.1323 bohr/cc-pVDZ,two frozen core
orbitals,ten active electrons,forty-component averaging,five retained roots
per spin/irrep,radius18,maxl4,deletion1e-6 and propagation100. The active space
changes[5,3,3,0]→[6,3,3,0],including its separately optimized target frame.

CAS12's original target-only config says deletion1e-7; its actual scattering
preparation says1e-6. The reader explicitly resolves that difference rather
than silently equating defaults. CAS11 selected scattering used4ranks/release
engine; CAS12 used8ranks/audited engine. Existing exact CAS12 engine reproduction
supports numerical consistency; this study includes no same-engine CAS11
reexecution. Original target solvers/import-action gates remain distinct.

## Extraction sensitivity after the channel expansion

The declared study retains **all96 window/order/sector/model cases**, with
**288 full-data starting-point fits and480 held-point folds**. Every case is
valid: optimization succeeds without a parameter bound,raw native fits are
nonempty/noMAXFIT,and original exits pass. No fit subset replaces the gate.

| Model | Full-width span / median | Maximum held-phase error | Extraction verdict |
|---|---:|---:|---|
| Equilibrium CAS11/five-A2 | 17.235% | 0.06870 rad | Reject |
| Equilibrium CAS12/five-A2 | 16.635% | 0.06819 rad | Reject |
| Stretched CAS11/five-A2 | 6.999% | 0.04113 rad | Reject |
| Stretched CAS11/nine-A2 | 6.991% | 0.04108 rad | Reject |

The width budget is5% of the median with0.001-eV floor; held-phase budget is
0.05rad. Position spans pass0.05eV. Equilibrium windows are1.8–3.3,1.5–3.5 and
1.2–3.8 requested-eV; stretched windows are0.4–1.6,0.5–1.5 and0.6–1.4. The
changed window/fold scope explains why these spreads differ from earlier
single-window/native summaries; earlier rejected records stay preserved.

The constant-background order1 dominates the measured sensitivity. For the
stretched nine-A2model,nonconstant orders2–4have per-order window width spans
≤0.153meV and maximum held-phase error≤0.002635rad. This is a **diagnosis**,
not permission to drop order1 and declare a passing extraction. A revised
estimator needs its own physical/analytic validation,predeclared model-selection
rule and unused-data check. Small phase-prediction errors do not establish
unique resonance parameters or resonance-pole identity.

## Independent reconstruction and historical next input gate

The pinned Sadaharu image and macOS reproduce all decisions with maximum
full-data fitted-parameter drift4.712e-9eV. The shared reader/fitter establishes
cross-runtime reconstruction consistency,not an independent scattering engine.
The host read costs1.716seconds including container startup. Public download
then reconstructs all124 raw inputs with57 verified frozen source files; all
seven published payload digests pass. **No native solve is repeated for this study.**

The subsequent native step was separately declared in
[the equilibrium nine-A2 input contract](a2-equilibrium-target-preflight-contract.json).
It qualifies48 native target vectors/properties at the unchanged original
equilibrium CAS11 orbital checkpoint,before a matched equilibrium channel
comparison. The same-model channel effect cannot be inferred from the stretched
rejection or the active-space shift.

At [21:28UTC](progress-a2-equilibrium-preflight-20261015T2128.json),owner1345746
and observation-only observer1346072 are active. **All eight target solves
exit0**, totaling248.123seconds; maximum native residual is1.20921e-10Hartree.
DENPROP is active;the import verdict and independent binary/property checks are
pending. Actual target decks match their prelaunch selection,all59source/
50input hashes and242installed engine artifacts reverify. Runtime source is
frozen at`2feb62d83d461df95fc376dc8971bd28d86bd066`.

Allocation is8physical CPUs4–11/eight MPI ranks/one library thread,
24-GiB cgroup/native6-GiB memp/8-GiB shared memory. Same-engine stretched
preflight cost51.686minutes/5.466GiB; equilibrium forecast45–90minutes is
conditional,with4-hour native/5-hour owner hard stops. Completion schedules
no scattering. Independently read the48-state vectors/Gram matrices/energies,
Pi,own-model multipole IDs and dipole before declaring the next comparison.
