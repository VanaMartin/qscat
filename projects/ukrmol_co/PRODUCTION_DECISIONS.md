# CO production decisions — 21 October 2026

The next milestone is a **matched electronic-model scattering comparison**.
The one-identity radial experiment is complete. Production readiness is measured
by decisions below, with each expensive successor separately costed from the
model and channel map it will actually use.

## First calculation: equilibrium orbital-ensemble control

The existing CAS11 fifty-component QC frame changes its ground energy by
**+0.456512 eV** and the common-forty objective by **+0.170147 eV**, relative to
the forty-component frame. These are diagnostics of different optimized frames,
not estimates of a resonance error. The independent fixed-orbital coverage already
passes48 probes on each of the CAS11/CAS12 frames. The missing step is native
import followed by matched scattering.

The [one-job recipe](calibration-sa11-fifty-import-20261021.json) and
[finite contract](fifty-import-contract-20261021.json) advance the equilibrium
CAS11 frame first:

- Same CAS(10,11), two frozen A1 spatial orbitals, cc-pVDZ and R=2.1323 bohr.
- Fifty equal-weight orbital-average components: `[10,5,5,5]` per spin.
- Native roots `[13,8,8,8]` per spin, **74 total**. All have independent
  fixed-orbital coverage. The DENPROP import setting uses50 states.
- Start from the existing converged fifty-component checkpoint. Preserve the
  ordinary200-cycle CI cap and original QC tolerances. Reoptimization must
  preserve every ensemble root within1e-7 Hartree, dipole within1e-5 a.u. and
  core/active singular values above1−1e-8.
- Verify all74 binary native eigenvalues/vectors and complete EOF, the native
  solver's requested/converged counts, all50 averaged-root imports and ground
  dipole. This qualifies the input frame, not ensemble convergence.

**Cost:** earlier same-basis64-root import3086.6s /5.266GiB on four ranks.
New forecast3600s expected /7200s upper /10800s hard cap; eight MPI ranks on
physical cores4–11, single-threaded native libraries and24GiB memory. The
existing PySCF target builder uses eight threads during QC. Host admission
requires32GiB available RAM /30GiB disk. Extra ranks and nonuniform74 roots have
not been benchmarked. All output, original exits and source/input hashes persist
under `/home/kooza/ukrmol/co-sa-20261006` on Sadaharu.

The subsequent scattering comparison must retain the **same physical channel
manifold** on both frames. Using50 DENPROP states here does not authorize a
50-versus40/48 channel comparison. Before a new inner solve: reconstruct explicit
property/dense-row maps, inspect physical subspace correspondence and retained/
omitted exchanges, and cost the resulting sparse2048-root solve. A2 root9 is not
covered by this74-root import. If a48-state nine-A2 comparison is selected, its
extra A2 input qualification needs a separately declared finite scope.

For orientation, the old equilibrium CAS11 full dense pipeline cost4.423h /
25.201GiB; this is **not** a forecast for a newly mapped sparse calculation.
The accepted CAS12/2048 baseline and2048→4096 control remain reusable.
No dense scattering Hamiltonian is authorized. No automatic channel/root ladder
follows a rejection.

**Completed21 October:** target and analyzer pass, but seed-frame preservation
rejects6/50 averaged roots, maximum3.38450e-7 Hartree against1e-7. Native energies
agree with the newly optimized frame to5.29468e-10 Hartree. Original verifier
exits1/1/owner1; overall import remains unset. Actual cost **2.899h /5.518GiB**,
dominated by135min QC and35.5min DENPROP, not the3.33min eigensolvers.
The [result/evidence](FIFTY_COMPONENT_IMPORT.md) preserves all153 original files.
The next candidate is independent coverage on the completed **new** fixed-orbital
frame, not another reoptimization or a weakened seed gate. Earlier96 probes on
two frames cost1280.58s; a one-frame forecast of10–20min requires its own finite
contract before submission. No successor is yet running.

**21 October,14:54UTC advance:** the [new fixed frame](FIFTY_FRAME_COVERAGE.md)
now passes all48 independent CI probes and all74 native binary root/vector/EOF
checks, exits0/0/0, in **5.144min /.27065GiB**. The old seed gate stays rejected.
The global-low50 property export contains only3 singlet-A2 states. An explicitly
selected40-state [property-only calculation](FIFTY_SELECTED_PROPERTIES.md) now
runs, forecast30min /upper45min /cap60min, using the qualified vectors and no new
QC/target diagonalization. Then establish physical cross-frame channel identity
and compatible scattering preparation/cost, rather than comparing unequal maps.

## Decision queue and release conditions

| Decision | Existing evidence | Next deliverable / release condition | Cost basis |
|---|---|---|---|
| Orbital ensemble | New frame independently qualifies48 probes /74 native roots; original seed rejection kept | Selected40 properties running; then physical cross-frame channel matching and a fixed-manifold scattering contrast | Import2.899h; new-frame coverage5.144min; properties forecast30min /cap60min; scattering cost after preparation/maps |
| Basis | Qualified CAS11 larger-basis target frames exist; ground-state accuracy alone is insufficient | One matched basis contrast on the selected ensemble/channel frame, phase first | Reuse imported targets; new inner solve requires its own measured forecast |
| Active space / channels | CAS11→CAS120.1213rad and64–73meV candidate-position shifts; stretched five→nine A2 changes0.0679rad | Choose a model only after independent ensemble/basis effects are measured; budget remaining uncertainty explicitly | Reuse existing2048 results; no automatic larger active-space or channel ladder |
| Resonance widths | Compressed/equilibrium required-window prediction rejects;94 known beyond-budget width counterexamples persist | A physically constrained extraction design with known-truth nonpolynomial successes, held-out phase predictions and explicit pole/width-identifiability veto | Existing calibration/native reading took13.9s /4.0s; new fit design remains uncosted; direct native complex poles need separate numerical qualification |
| Continuity | Retained/omitted exchanges and ambiguous coarse links remain | Neighboring-geometry physical subspace/resonance correspondence for the selected model | Bound the geometry pilot after the model is selected |
| Neutral reference |52 QZ/5Z points pass basis/interpolation over1.9–2.5bohr | Correlation-treatment uncertainty over this domain; widening to3–4bohr needs a stable reference | Reuse52 points; previous46 new points cost3.21h wall sum |
| Production launch | Model/extraction choices open | Frozen restartable small pilot, held-out geometries/energies, measured time/disk forecast | Cost from the selected model; no defensible total forecast yet |

The position/full-width/phase scientific budgets remain **0.05eV /5% with
0.001eV floor /0.05rad**. A failed extraction cannot supply a certified width
comparison; phases and target properties can still reject model insensitivity.
A passing pair is finite sensitivity evidence, not an absolute electronic error.

## Width work that can proceed from existing data

The discrepancy is named: compressed best inner-held phase errors0.00225936 /
0.00130390rad, equilibrium broad-window0.00108757rad, against0.001rad. Adding
background flexibility can reduce those errors while destroying identifiability.
The94 counterexamples establish that a local analytic/unitary phase factor does
not qualify as a causal finite-range CO Hamiltonian. No replacement width is
published merely because a more flexible function predicts phases.

The next extraction design must specify the physical background class and
threshold dependence, the meaning of its pole, its known-truth generators and
mandatory negative controls **before** a new CO application. Reuse the untouched
audit points, original rejected windows and counterexamples. A failure to find
another width is not proof of uniqueness. In the meantime, the model decision
uses matched real-axis phases without waiting for a new numerical-pole campaign.

## Numerical scope supporting these decisions

The radial full400 scientific record and three full captured operand sets are
exact. The original wrapper schema failure/exits1/1/1 remain preserved beside
the passing saved-only diagnosis. One radial identity does not qualify every
energy,800bohr or complex propagation. Original interior false/exits1/1 and
historical132/594 coverage remain. A further numerical pilot must name the
specific accuracy question affecting the proposed physical calculation and its
finite cost; filling every matrix or extending every radius is not a launch gate
by itself.
