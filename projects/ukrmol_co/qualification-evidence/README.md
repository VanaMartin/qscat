# CO electronic-qualification evidence

This supplement preserves **46 completed attempts in 18 batches**:
nine passing QC targets, ten passing neutral records and five neutral failures,
dense and five/eight-root SLEPc tight CAS(10,11) target imports, four rejected Davidson controls,
four passing small-model SLEPc controls, tight CAS(10,11) scattering and seven rejected CAS(10,11) ladder
entries, the matched tight CAS(10,10) scattering control and four numerical
CI-space repairs (two stretched-DZ passes and two projected-TZ failures).
It also includes all **160** CAS(10,10)/(10,11)/(10,12)
fixed-orbital CI coverage probes, final checkpoint comparisons and the
CAS(10,12) installed-library workspace audit and guarded utility checks,
including thirteen additional valid 16-/32-rank target-sector queries.

All three orbital starts recover the same fixed-ensemble solution: maximum
root/dipole differences are 7.69e-9 Hartree / 5.13e-8 a.u., with minimum final
active-subspace overlap singular value 0.999999999999924. Forty-six CI probes
converge completely; the eighth singlet B1/B2 root fails at trial-space size
40 and converges at 80/160. Every ensemble root converges in every probe.
The original nonzero diagnostic exit and failed controls remain preserved.
Both tight CAS(10,12) restart/RHF starts also pass: maximum root/dipole
differences are 4.17e-8 Hartree / 6.27e-8 a.u., with minimum final active
overlap singular value 0.999999999999496. All five ensemble roots pass every
CAS(10,12) coverage probe within 2.85e-13 Hartree. Two extra-root controls
fail at space 40 and converge at 80/160; a failed pre-probe import-path launch
and its fresh repaired retry are both retained.

The neutral pilot compares conventional frozen-core CCSD(T) at R=1.9/2.1323/2.5
bohr in aug-TZ/QZ bases. Its dipoles are CCSD lambda-density values. The
independent cc-pVDZ PySCF/Psi4 energy control agrees within 8.3e-11 Hartree;
the sentinel relative curves change by 76.04/68.47 meV for TZ-to-QZ and
17.24/15.82 meV for QZ-to-5Z. Four aug-TZ/QZ checks at R=3.0/4.0 fail external
RHF stability and stop before CCSD(T); that restricted route is blocked there.
These are qualification pilots, not a converged electronic model or neutral
curve. The tightened import passes all 40 averaged roots within 5.33e-10
Hartree and the independent dipole within 7.95e-11 a.u.; it took 97.75 minutes
with an 18.42-GiB kernel memory peak. The matched nine-job CAS(10,11) ladder
has two QC passes, three CI failures and four simultaneous geometry/basis
projection rejections. Its 50-component ensemble shifts the ground root by
0.45651 eV and the dipole by +125.09% relative to the 40-component objective.
The tight scattering run has complete 99-point grids, consistent Pi components,
two contracted dimensions of 27546 and a 4.423-hour wall / 25.20-GiB peak.
Its default candidate is at 2.520805 eV with full width 1.151002 eV. Position
and phase changes of 64.84 meV / 0.109809 rad from the earlier CAS(10,10) model
exceed the chosen gates; QC controls also differ, so the comparison does not
isolate the active-space change. Reference inputs and raw phase grids travel
with the archive. The matched tight CAS(10,10) rerun now passes in 21.38 minutes
at a 2.864-GiB peak; all forty imported roots/dipole and 48 independent coverage
probes pass. Tightening changes its earlier position/full width by only
−0.544/+0.395 micro-eV and phases by 1.30e-6 rad. Matching numerical controls
therefore leaves the active-space comparison failing at +64.8401 meV and
0.1098099 rad. The two passing stretched-DZ repairs still require their own
coverage/subspace and independent-import qualification; both projected TZ
repairs retain singlet B1/B2 CI failures after orbital convergence.

The later [near-equilibrium neutral companion](neutral-near-equilibrium/README.md)
independently verifies all 52 QZ/5Z points, including 46 new jobs and six unchanged
anchors. Cubic/PCHIP midpoint and 20-meV basis budgets pass; correlation treatment
and the four rejected stretched RHF references remain unresolved. Its public
archive verifies 959 payloads and 63867 profile samples.

The [native sparse control companion](sparse-native-control/README.md) preserves
the successful B1 128-root CAS(10,10) sparse eigensolve and its rejected legacy
SWINTERF continuation. Energies/residuals/continuum coefficients pass dense
differentials at a 1.228-GiB peak; phase and resonance verdicts remain unset for
that original attempt. The [ongoing replay guide](../SPARSE_SCATTERING.md)
describes native MPI boundary export and finite omitted-spectrum controls.

Forty-nine native saved-data replay attempts retain **48 successes and one
failure**: twelve background/detection controls, twenty-four two-component
window/background controls, and twelve partial B1 passes before a B2 unit-binding
failure. The fresh unit-fixed retry reads the template's actual `LUKMT`
(921/922), reproduces both default fits and checks raw saved-point fit grids.
Window clipping at fixed background changes position/width by at most
6.03 meV / 3.62%; background terms 1–4 change width by 12.65%, failing that gate.
The original exits/logs and corrected sources are retained; K-matrix/R-matrix
binaries stay on Sadaharu for execution. Twenty-four independent phase fits
reproduce the native fits, check printed residues and predict 120 interleaved
held-point sets; the poorer constant background remains evidence. The pinned
engine's 0.0735 Ryd/input-eV conversion (18.44-ppm offset from modern units) is
explicitly audited and retained.

Four outer-only pipelines refine the full grid from 99 to 197/393 points in
6.37 minutes at 0.719 GiB. Common-point phases agree at printed precision;
fixed-window linear-background position/full-width changes are at most
0.145/0.292 meV (0.0254% in width), within the chosen gates. Both fine automatic
native fits hit `MAXFIT=100`; their truncated records have unset acceptance
verdicts and independent fits use the full fixed intervals. Four original
driver failures remain preserved: two MPI-input EOFs, a pause-state race and
a completion-message mismatch after four successful native stages. Every
CPU-slot lease records resumption of its waiting owner.

The later [radial-propagation companion](outer-propagation/README.md) completes
eight pipelines and verifies forty native stages, all requested propagation
sectors and full-precision binary K matrices. Step/radius refinement changes
phases by at most 2.57e-6 rad, with successive refined differences below
2.28e-10 rad, passing the chosen numerical propagation gates.

The [fresh-basis companion](fresh-basis/README.md) preserves four completed
fresh TZ/aug-TZ QC starts: three passes and the original CAS(10,11) aug-TZ
fresh-CI rejection. Its 194 payloads and 160 raw final-state energies/spins
verify. Fresh CAS(10,12) TZ→aug-TZ lowers the ensemble objective by 1.80 eV
but raises the ground-state energy by 0.480 eV; competing-start and import
qualification remain open.

The subsequent [fresh-target coverage companion](fresh-coverage/README.md)
passes all 144 probes/936 eigenpair evaluations for those three targets under
the 600-cycle contract. All raw spectra/spins and 352 payload digests verify;
physical/spin-penalized residual maxima are 6.96e-10/9.99e-10 Hartree.

The [fresh aug-TZ repair companion](fresh-augtz-repairs/README.md) retains two
successful space-80/160 restarts of the rejected CAS(10,11) checkpoint.
All forty roots agree within 1.667e-11 Hartree, dipoles within 3.230e-12 a.u.
and core/active subspaces to roundoff. All 129 payloads and 80 raw final states
verify; the original rejected parent remains exact. Independent coverage,
competing starts and imports remain required.

The [fresh TZ import companion](fresh-tz-import/README.md) now checks all 64
CAS(10,11) roots within 3.196e-8 Hartree against the covered seed's independent
controls. The engine run passes in 54.38 minutes / 5.270 GiB. All 181 payloads
verify, including the original rounding-contract rejection, a CIDATA parser
recheck failure and the passing final-set/decimal recheck. Physical gates are
unchanged; competing-start/model qualification remains open.

The [repaired aug-TZ coverage companion](fresh-augtz-coverage/README.md) also
passes all 96 probes/624 eigenpair evaluations, with physical/penalized residual
maxima 6.62e-10/9.94e-10 Hartree. All 253 payloads verify. The 10.02-minute /
0.391-GiB scan releases a gated 64-root UKRmol import on Sadaharu.

That [aug-TZ import companion](fresh-augtz-import/README.md) now passes all 64
native roots against both repaired seeds within 5.969e-11/5.764e-11 Hartree.
Dipole and orbital-subspace gates pass; engine, supervisor and verifier all
exit zero. All 181 payloads and 96 raw reference spectra verify. Cost is
56.01 minutes / 5.268 GiB. Competing-start projections are now active.

The [competing-start companion](competing-starts/README.md) subsequently retains
two failed downward projections and two passing aug-TZ targets. Their space-80/160
pair passes, but both differ from the prior fresh branch: objective **−0.538388 eV**,
ground energy **−1.115756 eV**, z dipole **−0.06670992 a.u.** and minimum active
overlap **0.141085**. All 277 payloads and 80 raw final states verify; original
batch/controller exits stay one. A tested core+active projection interface repair
releases new-name downward retries after coverage of the new branch.

That [new-branch coverage](competing-augtz-coverage/README.md) now passes all
96 probes/624 eigenpair evaluations in 10.17 minutes / 0.380 GiB. Residual maxima
are 6.630e-10/9.977e-10 Hartree; all 335 payloads verify. The downward retry
queue starts after slot release, and the covered branch's native import is queued.

The [matched MPI/throughput companion](mpi-throughput/README.md) now passes all
seven CAS(10,10) replicas. One/two/four ranks take 3105.91/1830.62/1239.61 seconds;
four concurrent single-core jobs finish in 3444.36 seconds, **1.440×** the
throughput of repeated four-rank execution on the same allocation. All 1221
payloads and 98997 raw profile samples verify; numerical equivalence includes
the distinct initial-checkpoint lineage. Continuum work resumes after slot release.

The [CAS(10,12) import companion](cas12-import/README.md) subsequently passes all
40 native roots and the independent dipole on Sadaharu. Native roots agree with
current QC within **4.801e-11 Hartree**, covered seed spectra within **6.815e-9
Hartree**, and the separate RHF-start target within **3.485e-8 Hartree**; current
QC/import dipole error is **1.279e-11 a.u.** Complete cost is **6.201 hours /
40.405 GiB**. Serial DENPROP is **80.768%** of wall. All 335 payloads, 40 raw QC
states, 48 raw reference probes/312 eigenpairs and 111339 resource samples verify.
The original two extra-root coverage failures and prior controller retry failure
remain retained. The stretched CAS(10,11) queue starts on the released core group.

That [stretched import companion](stretched-import/README.md) now passes all
64 native roots against both covered DZ seeds within **3.598e-10 / 3.836e-10
Hartree**, with passing dipole/subspace gates. Cost is **54.25 minutes /
5.269 GiB**. Its 488 payloads preserve the original saved-table verifier failure
and passing exact-decimal recheck, and reconstruct 40 raw QC states, 96 reference
spectra/624 eigenpairs and 16229 resource samples. Same-geometry TZ/aug-TZ QC
subsequently passes in **14.57 / 20.33 minutes**, at **0.345 / 0.434 GiB**.
Its [170-payload companion](stretched-basis-qc/README.md) reconstructs 80 raw
QC states, same-geometry basis shifts and 10430 resource samples. Fixed-orbital
coverage and individually gated native imports continue on CPUs 8–11.

The [downward-projection companion](downward-projections/README.md) also passes
all four new-name equilibrium TZ starts and six same-basis restart comparisons.
Both aug-TZ lineages recover qualified fresh TZ within **4.156e-8 Hartree /
6.916e-8 a.u.**, with minimum active overlap **0.9999999999995466**. The 296
payloads reconstruct 160 raw QC states and 34132 resource samples. Original
projection-interface failures and distinct aug-TZ solutions remain retained;
the new branch's native import subsequently passes; electronic-model qualification
remains open.

The [continuum/competing-import companion](continuum-competing/README.md)
now passes both l=3/l=5 angular-cutoff comparisons and the competing aug-TZ
all-64-root native import. It verifies 1014 payloads, 36 native replays,
39 batch-source hashes and 219998 raw profile samples. Angular refinement
does not repair the 11.98–12.65% background-width sensitivity.

The [native scattering preflight](scattering-preflight/README.md) reproduces
CAS(10,11)'s dimension 27546 and measures CAS(10,12)'s dimension **86352**.
Valid aggregate workspace-array floors **222.391 / 222.555 GiB** block the full
current scattering solve on Sadaharu; the four-rank query retains its overflow
failure and no estimate. Its 814 payloads preserve the original scheduling
failure and four failed native preparations, and verify 245 source hashes and
1099 raw profile samples. Actual CAS(10,12) scattering runtime/scaling remain open.

Staged QC repairs, the neutral pilot and finite CAS(10,11) two-anchor follow-up
continue locally. The prior
eleven-/twelve-/eighteen-/thirty-seven-/thirty-eight-/forty-one-attempt supplements remain in
`manifest.json`'s `prior_snapshots` with its immutable URL and digest.

The Davidson controls request 5/8/16/32 roots and all fail required root imports
despite engine convergence. Independent diagonalization of the stored small
singlet-A1 Hamiltonian agrees with QC within 3.69e-13 Hartree. The SLEPc-enabled
source build passes 12/12 upstream checks; five-root and eight-root/tighter
small-model targets pass all 40 required roots within 4.98e-10 Hartree and
dipoles within 5.62e-11 a.u. These controls use about 0.12 GiB in 49.53/46.91
seconds. The full five-root CAS(10,11) SLEPc target passes all 40 roots within
5.34e-10 Hartree and dipole within 8.04e-11 a.u., at 52.55 minutes / 5.25 GiB.
Observed SCATCI/full-wall reductions are 89.17%/46.25%; these are not controlled
MPI scaling results. The source audit confirms dense PETSc Hamiltonian storage:
the largest CAS(10,12) matrix alone needs 37.41 GiB.
The eight-root/tighter CAS(10,11) import also passes all forty ensemble roots
within 5.33e-10 Hartree and the dipole within 7.90e-11 a.u., in 51.44 minutes
at a 5.27-GiB peak. Sixteen fixed-orbital probes at spaces 80/160 fully converge
and independently check all **64 requested roots** within 5.34e-10 Hartree,
including raw CASCI energies, root order and spins. Common dense/five-root
excitations match at printed precision. The scan takes 128.22 seconds at a
0.280-GiB peak; an initial pre-probe process-identity guard failure and its
repaired launch are retained. The later CAS(10,12) companion above establishes
the full forty-root/dipole target import locally.
The preceding 55-attempt target/scattering snapshot remains available in
[`../sa-evidence/`](../sa-evidence/README.md).

## Fetch and reconstruct

```bash
uv run qscat-run fetch projects/ukrmol_co/qualification-evidence
EVIDENCE_DIR=$(mktemp -d)
tar -xzf projects/ukrmol_co/qualification-evidence/qualification-evidence.tar.gz \
  -C "$EVIDENCE_DIR"
EVIDENCE_ROOT="$EVIDENCE_DIR/state-averaged-evidence"
uv run python -m projects.ukrmol_co.collect "$EVIDENCE_ROOT" \
  --provenance "$EVIDENCE_ROOT/qualification-provenance.json" \
  --output "$EVIDENCE_DIR/reconstructed.json"
```

Extract outside the checkout: historical source snapshots contain Python test
modules. The reconstructed JSON matches `qualification-results.json` in the
archive. `neutral-results.json` retains the fifteen-attempt neutral-only subset.
The file index records every payload's size and SHA256. Publication verification
checked 4635 payloads, 298 batch-source hashes, 17 embedded-image source hashes,
raw reanalysis of all 28 successful records, all 160 CI spectra/spins and
exact reconstruction of both aggregates. Repackaging is byte-identical.
All four Davidson rejections are reproduced; the stored Hamiltonian is
independently reconstructed and its eigenpair residuals checked.
All seven ladder failures, the four unstable neutral references and the
CAS(10,12) launch failure are checked against their retained diagnostics/logs.
The extra-root launch failure's original source/hash and exit are preserved;
the repaired worker's CPU-slot lease records verify resumption of its waiting
follow-on after the scan. Native scattering dimensions, raw phase comparisons,
all 49 replay attempts and saved-point fit grids are independently reconstructed.
The matched numerical configurations, generated diagonalizer inputs and raw
phase/candidate comparisons, independent/held-point fits, four outer grids and
their MAXFIT limits are also reconstructed. All 41 prior run records and sixteen
prior batch records remain exact.
The later [stretched-repair coverage companion](stretched-coverage/README.md)
retains another 96 probes and the failed strict coverage gate. It also measures
continuity of the matched CAS(10,10) orbitals. The main 46-attempt archive and
its prior snapshots remain frozen.
The [iteration-refinement companion](stretched-iteration/README.md) subsequently
repairs all nine failed stretched settings at 600 cycles and passes complete
96-probe scans with explicit Hamiltonian residuals. The original 200-cycle
outcomes remain preserved in their separate archive.
The full dense/SLEPc comparison and its storage-floor calculation are independently
reconstructed from raw stages, resources, roots, dipoles and pinned source.

The archive retains raw logs, configs, resource/stage records, source snapshots,
engine build/test provenance, upstream licences, RHF/CCSD/CASSCF checkpoints,
projected initial checkpoints and all workspace query commands. Valid eight-/
sixteen-rank triplet queries give array floors of 148.92/148.99 GiB; the
undersized four-rank query is explicitly rejected by the shipped utility.
The extended all-sector audit reaches 149.78/149.91 GiB on 16/32 ranks;
these are array floors before engine/library overhead.
Large UKRmol integral/channel/K-matrix artifacts remain on the host.

The archived `verify-qualification-evidence.py` reproduces the checks above:

```bash
PYTHONPATH=. uv run python "$EVIDENCE_ROOT/verify-qualification-evidence.py" \
  "$EVIDENCE_ROOT" \
  --original-archive projects/ukrmol_co/qualification-evidence/qualification-evidence.tar.gz \
  --repacked "$EVIDENCE_DIR/repacked.tar.gz"
```

Exact reanalysis/repackaging uses the analysis source commit recorded in the
manifest. Later analyzer versions deliberately record their own source hash.
The current utilities/packaging commit is
`28f5e1f6158abe775c8d2a70244504dacc46201a`, retaining the nested outer-engine
licence selection and independent phase-fit tools. `analysis-archive.py` retains the packaging
source; embedded image modules keep their actual original bytes and digests.
To replay the numerical recipes, copy archived `runs/` to a new host root
mounted at `/work` so retained checkpoint arguments resolve, and choose new
output/batch names. See [`../CONTINUATION.md`](../CONTINUATION.md) for the live
jobs and remaining gates, and
[the method contract](../../../docs/physics/co-electronic-qualification.md).
