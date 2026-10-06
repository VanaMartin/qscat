# UKRmol+ CO fixed-nuclei scattering experiment

**Status:** the execution pipeline works and the compact SEP continuum is
numerically stable at equilibrium, but the electronic model is not yet
accurate enough to provide production potential-fitting targets. The measured
calibration and throughput evidence below explains the remaining model issue.

This experiment preserves the first pilot and **40 calibration attempts: 30
validated calculations and 10 engine failures**, with 21 paired comparisons.
Here, “validated” means the output passes pipeline/internal-consistency checks,
not that its electronic model is physically converged.

## Convergence verdict

| Question | Finding | Decision |
|---|---|---|
| Does the engine/workflow reproduce reference calculations? | Five water target energies; 12 source-build serial/MPI checks; CO source/reference phases agree to `1e-7 rad` | Execution established |
| Is the compact equilibrium SEP continuum stable? | Tested angular/deletion changes shift position by less than 0.5 meV and width by less than 0.1%; background fits vary width by about 2.5% | Passes chosen numerical tolerances for this model/window |
| Is the electronic model converged? | SEP space/frozen-orbital changes shift position by 0.3–0.6 eV; ground-state-CASSCF CC target excitations and resonance fits remain model-sensitive | Further calibration required |
| Is the entire geometry curve qualified? | Compressed low feature is stable, high feature is not; stretched geometry has no usable automatic fit | Threshold diagnostics and model/continuity checks required |
| Can these data constrain a production potential? | No tested electronic configuration is qualified; a correlated neutral curve is still missing | Preserve as calibration evidence |

We have a reproducible execution platform and identified the dominant accuracy
problem. The resonance values have **not** settled onto an electronically
converged limit. The next useful refinement is balanced state-averaged target
construction, followed by active-space/basis/channel tests, rather than a denser
production geometry sweep of the current model.

Run a compact UKRmol+ static-exchange-plus-polarization (SEP) calculation
for electron scattering from CO. The first deck uses $R=2.1323$ bohr,
an RHF `aug-cc-pVDZ` target, two frozen core orbitals, five occupied valence
orbitals, and 20 virtual orbitals distributed as `[10, 4, 4, 2]` in UKRmol's
`[A1, B1, B2, A2]` order. Only the doublet `B1` and `B2` scattering symmetries
are calculated: the sectors containing the two components of $^2\Pi$ for CO
along the z axis. These C2v sectors also contain higher odd angular-projection
components, so their saved eigenphase sums are not isolated single-partial-wave
phases.
The continuum uses an 18-bohr sphere with Gaussian partial waves through
$l=4$ and a $10^{-7}$ orthogonalization deletion threshold. Double-precision
GBTOlib is selected by default. Psi4 uses full SCF integrals (`scf_type pk`),
with energy and density convergence thresholds of $10^{-10}$.
Their cross sections are **symmetry-resolved contributions**, not the full
elastic cross section of polar CO.

Geometry and propagation radii are in bohr; target and extracted resonance
energies are in Hartree. The external UKRmol electron-energy grid is explicitly
in eV: 491 points from 0.10 to 5.00 eV. Cross sections are in bohr² and
eigenphase sums in radians. Native `RESON` output uses Rydberg;
`analyze.py` converts it to Hartree at the input boundary.

## Measured pilot — 6 October 2026

The three-geometry, 20-virtual-orbital pilot passed the pipeline and symmetry
checks on Sadaharu with four MPI ranks. Native resonance fits, relative to
the neutral target threshold, were:

| R (bohr) | Virtual orbitals | Position (eV) | Width (eV) | Wall time (s) | Sampled memory peak (GiB) |
|---|---:|---:|---:|---:|---:|
| 1.9000 | 20 | 3.4351 | 2.3269 | 24.5 | 3.26 |
| 2.1323 | 20 | 1.9401 | 1.0496 | 24.7 | 3.26 |
| 2.5000 | 20 | 0.1522 | 0.0307 | 24.3 | 3.26 |
| 2.1323 | 30 | 1.3343 | 0.5347 | 27.9 | 3.66 |

The 20-orbital Hamiltonians have dimension 537 per scattering symmetry; the
30-orbital check has dimension 1158. The sampled allocated run-disk peaks
were about 62 MiB and 114 MiB respectively. Native `B1`/`B2` eigenphase sums
and fitted parameters agreed at text-output precision. The largest
Psi4/UKRmol target-energy difference was $5.4\times10^{-10}$ Hartree.
The upstream water example also reproduced all five reference target
energies at their printed precision.

**The virtual-space dependence is substantial:** adding ten virtual orbitals
moves the equilibrium resonance downward by approximately 0.61 eV and halves
its width. A convergence and balanced-model study is needed before these
data can constrain a fitted potential reliably.

The compressed geometry required a wider 0.1–8 eV grid, with 0.05 eV spacing,
for the automatic fit. At 2.5 bohr, the narrow resonance was checked on a
0.05–0.55 eV grid with 0.001 eV spacing. Its fitted position and width changed
by less than $10^{-5}$ eV from the original 0.01 eV grid. Complete inputs,
native fitted values, resource measurements and artifact directory names are
recorded in [`pilot-results.json`](pilot-results.json).

Two setup corrections mattered: the 13-bohr sphere failed the free-scattering
kinetic positivity test in both precision variants, while the 18-bohr sphere
passed in double precision; and full-integral PK SCF removed a target-energy
mismatch from Psi4's default density-fitting approximation. The exploratory
quadruple-precision calculation was canceled after the working double-precision
pilot completed; its intermediate files remain available.

## Accuracy calibration

The low-energy CO papers supply model comparisons, rather than a single
converged reference value:

| Calculation | Equilibrium position (eV) | Full width (eV) | Evidence |
|---|---:|---:|---|
| Laporta et al. (2012), six-electron SEP, STO/numerical continuum | 1.67 | 0.82 | [Reference note](../../reference/literature/laporta-2012-psst21-045005.md), preprint p. 4 |
| Dora et al. (2016), cc-pVTZ/CAS(10,10), 50-state CC | 1.73 | 0.84 | [Reference note](../../reference/literature/dora-2016-epjd70-197.md), p. 6, Table 4 |
| Dora and Tennyson (2020), cc-pV6Z, 27 physical states / 41 C2v components | 1.8744 | 1.2916 | [Reference note](../../reference/literature/dora-2020-jpb53-195202.md), p. 4, Table 2 |

Laporta's final nuclear cross sections use a 10% empirical width increase and
a 0.035 eV level shift (preprint p. 6). The larger 2020 basis was aimed at
higher-lying resonances; its lowest-resonance width differs substantially from
the 2016 model. Neither agreement with one paper nor increasing basis size
alone establishes accuracy.

Freezing four occupied A1 orbitals approximates Laporta's six-electron
polarization prescription. It reduces, but does not eliminate, virtual-space
sensitivity:

| Target / frozen occupied orbitals | Virtual space | Position (eV) | Width (eV) |
|---|---|---:|---:|
| aug-cc-pVDZ / 4 | [10,4,4,2] (20) | 2.1342 | 1.1771 |
| aug-cc-pVDZ / 4 | [14,6,6,4] (30) | 1.6472 | 0.7172 |
| aug-cc-pVDZ / 4 | [17,9,9,4] (all 39) | 1.3314 | 0.5673 |
| aug-cc-pVDZ / 2 | [17,9,9,4] (all 39) | 0.8980 | 0.3413 |
| aug-cc-pVTZ / 4 | [17,9,9,4] (39) | 1.6234 | 0.7684 |
| aug-cc-pVTZ / 4 | [24,14,14,8] (60) | 1.3292 | 0.5743 |

The almost matching DZ39/TZ60 values do not establish electronic convergence:
TZ60 is a truncated virtual space, and changing its size by 21 orbitals moves
the resonance by 0.29 eV. SEP improves the scattering correlation while keeping
an RHF target fixed; this imbalance is discussed by Dora et al. (2016), p. 5.

For DZ39 with four frozen orbitals, the equilibrium numerical changes are
small:

| Refinement | Position change (eV) | Width change | Maximum B1 phase change modulo pi, 0.1–5 eV |
|---|---:|---:|---:|
| l=3 versus l=4, deletion 1e-7 | 0.000057 | 0.040% | 0.0036 rad |
| deletion 1e-7 to 1e-6, l=4 | 0.000436 | 0.046% | 0.0190 rad |
| l=4 to l=5, both deletion 1e-6 | 0.000297 | 0.060% | 0.0089 rad |
| propagation radius 100 to 200 bohr | below native output precision | below native output precision | 0 |

The l=5 calculation needs deletion 1e-6 to pass the kinetic-positivity check.
The 15-bohr sphere and complete 85-orbital TZ virtual space fail that check at
both deletion thresholds in double precision. These failed refinements are
preserved; they are not evidence of convergence.

Changing native `RESON` background order from 1 to 4 moves the DZ39 fitted
position by about 1 meV and its width by about 2.5%; changing the automatic
detection threshold from 0.7 to 1.3 leaves the default fit unchanged. Replays
use the saved K-matrices, so they incur no new integral or scattering solve.

Compressed/stretched sentinel checks at 1.9/2.5 bohr compare l=3 and l=4.
The compressed 0.1–8 eV calculation returns two fitted features; selecting
the shape resonance requires inspecting the eigenphases. The broad feature
near 2.79 eV changes by 2.7 meV and 0.6% in width between l=3/l=4, with a
maximum phase change of 0.0173 rad over 0.1–4.5 eV. The additional fitted
feature at 5.90 eV appears only at l=4; phase differences reach 1.47 rad over
the full 0.1–8 eV grid. It is not a validated resonance target. At 2.5 bohr, all-39
six-electron SEP gives no usable automatic fit on 0.01–0.50 eV. This differs
from the original SEP20 pilot's narrow resonance, and does not by itself
establish that the anion is bound. A threshold/pole diagnostic is required.

### Acceptance criteria and target model

Provisional **testbed** numerical targets are 0.05 eV in resonance position,
5% in width with a 0.001 eV absolute floor, and 0.05 rad in eigenphase sums
modulo pi. These are chosen tolerances, not accuracy claims from the papers.
They must hold for energy-grid/background changes and continuum refinements
at compressed, equilibrium and stretched geometries. Atomic-basis, active-space
and channel-count changes must be checked separately. Near threshold, an
empty automatic fit needs a pole/bound-state check rather than a zero-width
replacement.

`CAS-A` enables contracted close-coupling diagnostics with a correlated
ten-electron target and explicit excited states. `--orbitals natural` currently
optimizes **ground-state CASSCF** orbitals. Those are distinct from the
multi-spin/multi-irrep state-averaged orbitals used in the papers. CAS(10,8)
ground-state orbitals place the first triplet excitation near 11.7 eV and
produce resonance positions 2.73/2.89 eV for DZ/TZ with 36 external orbitals.
That target representation is inadequate. The 2016 state-averaged models
place the lowest triplet Pi excitation around 6.3–6.5 eV (p. 4, Table 2;
spectroscopic and vertical excitation energies are distinct).

After increasing the workspaces, CAS(10,10)/40-state ground-state-CASSCF CC
completed with the following equilibrium fits:

| Basis / external orbitals | Position (eV) | Width (eV) | Lowest triplet Pi excitation (eV) |
|---|---:|---:|---:|
| aug-cc-pVDZ / 0 | 3.1452 | 1.5807 | 7.2925 |
| aug-cc-pVDZ / 10 | 2.7248 | 1.2239 | 7.2925 |
| aug-cc-pVTZ / 0 | 3.4501 | 1.8251 | 7.9820 |

Their basis/external-space shifts fail the provisional position criterion.
Ground-state orbital optimization remains different from the published
state-averaged target construction.

Larger-CAS generation requires raising all three `CONGEN` workspaces:
`NDIMX`, `CDIMX`, and `NODIMX`. `--congen-workspace` sets them in the copied
template, with the latter two one tenth of the first. These are input limits,
independent of the host's free RAM.

The production candidate is therefore a state-averaged CAS(10,10) or larger
target, compared in DZ/TZ bases with 40–50 C2v components, and an independently
correlated neutral curve. The existing ground-state-CASSCF runner is a
diagnostic toward that model. No production electronic configuration has yet
passed the acceptance criteria.

### State-averaged target checks

The PySCF multi-spin/multi-irrep continuation has completed independent target
import checks at equilibrium. [`sa-results.json`](sa-results.json) preserves
52 completed attempts: sixteen successful target-only pipelines, twelve successful
scattering pipelines, and 24 preserved setup/diagnostic failures. One rejected
CLI configuration failed before creating a run directory; its request and
batch log remain part of the aggregate and archive. Passing these
pipeline checks is distinct from electronic-model qualification.

| Target | Neutral energy (Hartree) | Lowest triplet Pi (eV) | Ground dipole z (a.u.) | Max QC/UKRmol energy difference (Hartree) |
|---|---:|---:|---:|---:|
| aug-DZ CAS(10,8), ground-only control | −112.885976288 | 11.7194 | 0.1405573 | 9.4e-11 |
| cc-pVDZ CAS(10,10), 40-component average | −112.894737477 | 6.3952 | 0.0920784 | 2.3e-9 |
| aug-cc-pVDZ CAS(10,10), 40-component average | −112.821499062 | 6.3497 | 0.1138490 | 5.1e-10 |
| cc-pVDZ CAS(10,8), 40-component average | −112.855375534 | 6.4858 | 0.2023335 | 8.0e-10 |
| cc-pVTZ CAS(10,10), 40-component average | −112.857271560 | 6.2744 | 0.0986170 | 6.9e-10 |
| cc-pVDZ CAS(10,10), tighter orbital optimization | −112.894737475 | 6.3952 | 0.0920732 | 2.6e-9 |
| aug-cc-pVDZ CAS(10,10), projected DZ start | −112.898457558 | 6.3070 | 0.0826292 | 6.1e-10 |
| cc-pVTZ CAS(10,10), projected DZ start | −112.925516422 | 6.3612 | 0.0878540 | 5.7e-10 |
| cc-pVTZ CAS(10,10), 50-component average, projected TZ start | −112.920679908 | 6.3503 | 0.1164896 | 7.9e-10 |

Every averaged root passed the 1e-7-Hartree import gate, and independent
ground dipoles agree within 1.9e-6 a.u. The DZ neutral energy, lowest triplet
excitation and dipole (about 0.234 D) are close to the published 40-state
model's printed values; see [Dora et al. 2016](../../reference/literature/dora-2016-epjd70-197.md),
p. 4, Table 2. This is a method plausibility check. The augmented-basis neutral
root is **0.07324 Hartree higher** (about 1.99 eV): common active-space identity
and starting-point stability must be checked before interpreting this as a
basis-refinement sequence. An improved triplet excitation alone is insufficient.

The eight-component truncated average converged but split a singlet Pi pair
by about 0.0325 eV and was rejected. The full 40-component equilibrium
ensembles preserved Pi degeneracy. Independent target-only pipelines took
360/405 seconds for DZ/aug-DZ, with approximately 1.82 GiB kernel peaks each.
The early import probes used an uncentered target; subsequent jobs use the
upstream center of mass. Basis/active-space/geometry refinements and scattering
checks are recorded in the continuation manifests as they complete.
Tightening the equilibrium orbital tolerances changes the lowest triplet
excitation by only 1.8e-6 eV and the dipole by 5.2e-6 a.u.; this is a two-point
stability check, not a complete optimization convergence sequence. At R=1.9
and 2.5 bohr, the default optimization missed the strict Pi gate by target
splittings of 2.56e-7 and 1.35e-7 Hartree. Fresh tighter-orbital/CI runs retain
the same gate. The first tighter retries stalled above the requested gradient
tolerance after 100 macroiterations; lowering the augmented-Hessian metric
cutoff is a separate retry. Four subsequent diagnostic-instrumentation failures
were preserved and repaired before new jobs were launched.

Projected DZ starts lower the aug-DZ/TZ ground roots by 2.0941/1.8570 eV,
with initial/final active-subspace overlap singular values of 0.9806–1.0000
and 0.9957–1.0000. For aug-DZ, the HF-selected start still has the lower
**ensemble** energy despite its higher **ground** energy; the projected TZ
start lowers both. The augmented-basis discontinuity is therefore a
starting-space/ensemble tradeoff, not an established basis-convergence trend.

### Initial state-averaged scattering and fit checks

Centered equilibrium CAS(10,10)/40-channel runs used the starting numerical
controls (18-bohr sphere, l=4, double precision, deletion 1e-6) and 99 energies
over 0.1–5.0 eV. The DZ native default fit is **2.4560/1.1196 eV**
(position/full width). The aug-DZ default produces two overlapping candidates,
**2.3849/1.1874** and **2.2397/1.4162 eV**. A candidate count is not a physical
resonance count; native automatic windows can cover the same feature.

The modulo-pi eigenphase difference between these two targets reaches
**0.5142 rad**, failing the provisional 0.05-rad gate. Twenty-four saved native
RESON replays vary one to four background terms and detection thresholds
0.7/1.0/1.3. Detection-threshold changes leave these candidates unchanged.
DZ background refinements give 2.4491–2.4614 eV and 1.1196–1.2547 eV, a
width spread exceeding the 5% criterion when the constant background is
included. Aug-DZ replay windows generate overlapping and sometimes
implausibly broad candidates (one exceeds 80 eV in width); inspect the raw
fit windows, residuals and goodness factors before assigning a feature.
These are extraction diagnostics, not certified pole parameters.

Projected-start aug-DZ and TZ give one candidate each:

| Target / starting subspace | Position (eV) | Full width (eV) |
|---|---:|---:|
| cc-pVDZ / HF-selected | 2.4560 | 1.1196 |
| aug-cc-pVDZ / projected DZ | 2.5294 | 1.2737 |
| cc-pVTZ / projected DZ | 2.5253 | 1.1928 |

The projected aug-DZ/TZ pair changes position by 0.0041 eV and phase by
0.0337 rad, but width by 6.36% relative to aug-DZ, failing the provisional
5% gate. The unprojected/projected aug-DZ phase change reaches 0.5977 rad.
Another 24 native replays give projected aug-DZ widths of 1.2737–1.4460 eV
and TZ widths of 1.1928–1.3475 eV for one to four background terms. Dropping
the constant-background case reduces each spread to about 1%, which is an
extraction stability observation rather than a justified choice of fit model.

Lowering the augmented-Hessian cutoff lets the R=2.5 target converge with
an orbital gradient of 8.5e-8 and Pi splitting of 1.12e-8 Hartree. The R=1.9
retry still stalls at 6.15e-7 with zero orbital rotation and is rejected.
Tightening the CI residual together with compatible CI/augmented-Hessian
cutoffs subsequently gives a **passing compressed target**: gradient
2.99e-8, Pi splitting 4.63e-8 Hartree, and maximum independent root-energy
error 6.58e-10 Hartree, in 8.08 minutes. Tightening only the residual while
keeping the default CI cutoff had failed CI convergence. Coupled-Newton
controls/retries still fail their gradient checks; they remain rejected
diagnostics. The successful route is the one-step optimizer.

The compressed R=1.9 scattering pipeline subsequently passed on 159 energies
over 0.1–8.0 eV in 20.30 minutes, with Pi phase agreement at printed precision
(1e-7 rad), a 2.86-GiB kernel memory peak and 594-MiB sampled disk peak. Its
default candidate is **3.5191/2.1026 eV**. Twelve native replays keep one
candidate, with positions 3.5109–3.5191 eV and widths 2.0529–2.4573 eV as
background order varies; detection-threshold changes do not alter the fits.
The width spread still fails the 5% extraction criterion. Native fitting uses
the 1.79–5.31 eV window and rejects a second detected feature near 6.65–6.75 eV.
These observations do not certify a pole or establish identity with the
equilibrium feature.

The matching stretched R=2.5 scattering request failed during tighter
CI-residual target optimization, after 3.29 minutes: all CI solvers converged,
but the orbital gradient stalled at 1.19e-6 with zero rotation. This is distinct
from its earlier passing target with looser CI controls. The augmented-Hessian
accuracy ladder resolves this new stall: both 1e-16 and 1e-20 eigensolver
tolerances give six-iteration convergence at an orbital gradient of 8.50e-8
and pass all independent checks in 7.32/7.40 minutes. Their maximum root-energy
difference is 1.17e-11 Hartree and dipole difference 4.35e-12 a.u.; maximum
QC/UKRmol root error is 5.16e-10 Hartree. Pi splittings are 7.18e-12/1.61e-11
Hartree. This is a two-point optimizer-accuracy check following a failed looser
solve, not electronic-model convergence. A fresh stretched scattering
calculation uses the passing 1e-16 checkpoint.

That stretched R=2.5 scattering retry passed in **18.70 minutes**, with a
2.85-GiB kernel memory peak and 595-MiB sampled disk peak. On 300 energies over
0.01–3.0 eV, its default candidate is **0.9739/0.2818 eV**. Twelve native
replays retain one candidate in the 0.755–1.205 eV fit window: positions span
0.97367–0.97395 eV and widths 0.28178–0.28915 eV. The 2.61% maximum width
change relative to the default meets the chosen background-order tolerance;
detection changes leave it unchanged. These are fit candidates, not certified
poles or an electronic-model qualification.

The stretched **0.020/0.010/0.005-eV** spacing sequence gives:

| Spacing (eV) | Points | Position (eV) | Full width (eV) | Calculation wall (min) |
|---:|---:|---:|---:|---:|
| 0.020 | 150 | 0.9738685 | 0.2817588 | 20.46 |
| 0.010 | 300 | 0.9738707 | 0.2817833 | 18.70 |
| 0.005 | 599 | 0.9738722 | 0.2817949 | 22.16 |

The finest two fits differ by 1.50e-6 eV in position and 0.0041% in width;
phases agree at shared printed energies. The coarse-to-fine span is
3.63e-6/3.61e-5 eV in position/width. Automatic fit windows vary from
0.762–1.219 to 0.751–1.199 eV across the sequence; an independently varied
window remains a separate check. Another 24 background/detection replays on
the coarse/fine grids keep maximum width changes below 2.66%. Measured cost
depends on concurrent workloads; the coarser run is not a timing baseline for
the middle run.

The equilibrium DZ CAS(10,10)/40-channel **l=3/4/5** sequence gives
positions 2.45894/2.45597/2.45576 eV and widths 1.12623/1.11964/1.12248 eV.
Across the sequence, the largest pairwise changes are 3.18 meV in position,
0.585% in width and 0.02615 rad in phase. All meet the chosen gates over
0.1–5.0 eV. Width and phase changes are nonmonotone; this establishes tested
angular-cutoff stability within those gates, not a monotone extrapolation.
Radius, deletion, propagation and extraction remain model-specific controls.
The runs took 14.00/17.10/22.80 minutes, with 2.62/2.84/3.21-GiB kernel peaks.

The initial centered DZ/aug-DZ jobs finished together in **20.25 minutes** on eight
physical cores: 18.82/20.24 minutes per four-rank job, approximately
2.86/2.89 GiB kernel peaks and 593 MiB sampled run-disk peaks. Projected-target
scattering takes 21.88/22.59 minutes, with 2.89/3.45 GiB kernel peaks and
593/594 MiB sampled run-disk peaks. The 40-to-50-channel refinements test
channel retention separately; the new model still needs continuum/grid/fitting
qualification.

The DZ 40-to-50-channel check retains the 40-component orbital ensemble and
gives 2.4471/1.1178 eV in 22.86 minutes. Relative to the 40-channel run,
position changes by 0.0089 eV, width by 0.17%, and phase by 0.0167 rad, meeting
the provisional criteria for this pair. This is a two-point channel-retention
check, not a full channel-convergence sequence. Its Hamiltonian dimension is
8690 per Pi component, with a 2.96-GiB kernel peak and 665-MiB sampled disk peak.

Separately, increasing the projected TZ **orbital ensemble** to 50 components
raises the ground root by 0.1316 eV and changes its dipole by 32.6%; the lowest
triplet Pi is 6.3503 eV. The 50-component target passed all independent gates
in 6.56 minutes. These ensemble and channel-count changes measure different
effects and should not be pooled into a single convergence sequence.

### Larger active spaces and CI root coverage

CAS(10,11) finished all engine stages in 102.22 minutes, with an 18.39-GiB
kernel peak, but its fifth triplet-A1 QC root was 0.129757 eV above UKRmol's
fifth root. All original QC convergence/spin flags had passed. Six fresh
fixed-orbital CASCI probes (five/eight requested roots and Davidson spaces
40/80/160) recover UKRmol's lowest five roots within 2.75e-10 Hartree; common
energies agree across probes within 4.27e-14 Hartree. The original purported
fifth root matches the fresh sixth root. The target remains rejected: the
probe validates the fixed-orbital spectrum, not minimization of the intended
lowest-root orbital ensemble.

`--target-ci-fresh-start` tests the root-selection issue during reoptimization.
A fresh-CI audit now checks final root-selection stability before the expensive
UKRmol pipeline. The single-component ground control passes. An initial mixed
solver override failed because PySCF copied the first solver's instance kernel
into the mixer; installing the overrides after mixer construction fixes that
interface. The failed revision is preserved. A two-spin A1-only control then
passed its averaged-root import check but failed computed Pi degeneracy, so
it is not a passing physical target. The balanced DZ CAS(10,8) 40-component
control passes in 2.45 minutes: fresh/optimized roots agree within 8.53e-14
Hartree, Pi splitting is 2.85e-13 Hartree, and maximum UKRmol import error is
5.44e-10 Hartree. The CAS(10,11) restart retains all gates.

CAS(10,12) QC passed its original gates, but MPI-SCATCI stopped before target
diagonalization: its 43194-dimensional singlet-A1 matrix needed a 3.74-GB
local block, exceeding the template's 2.5-GiB internal budget per process.
The 32-GiB container limit does not enlarge that budget. `--scatci-memory-gib`
sets `memp` in both target and scattering templates; its default remains 2.5.
The checkpoint-based retry uses 6 GiB per process in an 80-GiB container,
allowing headroom for the later density-property stage. Choose and measure
the process/container budgets together rather than inferring them from host RAM.

## Measured cost and parallel execution

For equilibrium DZ39/four-frozen SEP on Sadaharu, with 491 energies:

| MPI ranks | Calculation wall time (s) | Aggregate CPU time (s) | Kernel memory peak (GiB) |
|---:|---:|---:|---:|
| 1 | 54.7 | 54.4 | 1.56 |
| 2 | 41.5 | 81.6 | 2.44 |
| 4 | 33.0 | 127.1 | 4.16 |

These measurements use the same four-physical-core CPU group (12–15);
other calibration calculations were running on disjoint cores. All three
rank counts produced identical printed fitted parameters. Four independent
one-rank jobs pinned individually to those four cores completed together in
**63.1 seconds**, about 2.1 times the throughput of one four-rank job on the
same allocation. Each used about 1.6 GiB; the combined envelope is roughly
6–7 GiB. The sampled DZ39 run-disk peak is about 180–188 MiB per job.

For **25–40 geometries at this SEP size**, four one-rank workers imply about
7–11 minutes by rounding up to full 63-second batches, versus 14–22 minutes
with one four-rank worker. A planning allowance of **10–15 minutes** covers
ordinary geometry variation and reruns. Using 12 one-rank workers projects
roughly 3–5 minutes, but that twelve-way throughput has not been measured.
This estimates an exploratory SEP campaign, not a converged electronic model.

Costs change strongly with the electronic model. Equilibrium TZ39/TZ60 SEP
needed about 85/109 seconds on four ranks and sampled 6.7/8.4 GiB peaks.
CAS(10,8)/40-state ground-state CASSCF CC needed about 145 seconds (DZ) or
267 seconds (TZ) for only 99 energies; outer-region propagation dominated.
An 80-state HF-orbital CC diagnostic took about 31 minutes for 491 energies
and found no usable automatic resonance fit. Channel count can therefore
cost much more than adding SEP energy points. The final balanced-model campaign
time must be remeasured after its target and channel space are selected.

The CAS(10,10)/40-state trials needed **13.4–15.1 minutes per geometry**
for 99 energies; kernel memory peaks were **3.2–5.5 GiB** and sampled run-disk
peaks **610–661 MiB**. Their three-job batch finished in **15.1 minutes** on
12 physical cores. At that measured size, 25–40 geometries with three workers
would take about **2.3–3.5 hours for 99 energies**, before allowance for
geometry variation. Increasing to 300–800 energies gives a stage-scaled
estimate using

`T(N) = T(99) - T_rsolve(99) + T_rsolve(99) * N / 99`.

Measured `rsolve` totals span 79–156 seconds; Hamiltonian solution and target
density processing dominate the remainder. This yields **16–33 minutes per
geometry**, or approximately **2.4–7.8 hours for 25–40 geometries on three
workers**, depending on model and energy count. This is a conditional planning
estimate for the measured CAS size: it assumes comparable geometries, the same
energy interval/channel structure and linear propagation cost. State averaging,
larger active spaces, additional channels, refinement runs and fitting add
cost. The accuracy-qualified configuration still needs its own timing.

For the state-averaged CAS(10,10) examples, measured 99-point propagation
totals range from 76.6 seconds (DZ40) to 290.6 seconds (DZ50); target and
scattering diagonalizations dominate the remaining cost. The same stage-scaled
formula predicts **21–57 minutes per geometry** at 300–800 energies, or
**3.2–13.3 hours for 25–40 geometries on three four-core workers** for these
measured 40/50-channel cases. This assumes the same geometry cost and comparable
threshold structure; it excludes refinement/retry time and is not a budget
for a larger, electronically qualified active space. The state-averaged model
has not yet had its own one/two/four-rank and concurrent-geometry scaling sweep.

Even these CC diagnostics used far less than Sadaharu's approximately 123.5
GiB RAM. Three 4-rank jobs use 12 of its 16 physical cores; this is a practical
starting layout for larger CC models. Per-job limits of 16–24 GiB leave ample
headroom for the completed small models; kernel and sampled peaks are recorded
where available. `/home` scratch and retained artifacts are the storage budget
to track, not just each `scratch/` subdirectory.

## Reproduce on a Docker host

The shared [Docker guide](../../docker/ukrmol-plus/README.md) describes engine
sources/checksums, all six build targets, dependency additions, SSH deployment
and bind-mount ownership. Build from the repository root on the host:

```bash
docker build --build-arg BUILD_JOBS=4 -f docker/ukrmol-plus/Dockerfile \
  -t qmodeling/ukrmol-co:source .
```

Choose a writable persistent `RUN_ROOT` under `/home`, following the Docker
guide's UID/GID setup. Run and analyze a fresh equilibrium pilot:

```bash
IMAGE=qmodeling/ukrmol-co:source
docker run --rm --cpuset-cpus=0-3 --memory=16g --memory-swap=16g \
  --shm-size=1g -v "$RUN_ROOT:/work" \
  --entrypoint /opt/ukrmolp/entrypoint.sh "$IMAGE" bash -c \
  'python3 -m projects.ukrmol_co.run --workdir /work/runs/co-equilibrium-001 && \
   python3 -m projects.ukrmol_co.analyze /work/runs/co-equilibrium-001'
```

The default is the initial two-frozen/20-virtual SEP pilot. For the more
numerically checked **diagnostic** six-electron SEP model, add
`--frozen-orbitals 4 --virtual-orbitals 17 9 9 4 --deletion-threshold 1e-6`.
Neither is a qualified production configuration.

### Multi-spin/multi-irrep target orbitals

`--orbitals state-averaged --model CAS-A` uses the pinned PySCF backend in the
CO image. Supply the active space explicitly. `--target-only` first compares
its QC roots against independent UKRmol target diagonalizations:

```bash
docker run --rm --cpuset-cpus=0-3 --memory=16g --memory-swap=16g \
  --shm-size=1g -v "$RUN_ROOT:/work" \
  --entrypoint /opt/ukrmolp/entrypoint.sh "$IMAGE" bash -c \
  'python3 -m projects.ukrmol_co.run --workdir /work/runs/co-sa-target-001 \
     --model CAS-A --orbitals state-averaged --active-orbitals 4 3 3 0 \
     --virtual-orbitals 0 0 0 0 --target-only --ranks 4 && \
   python3 -m projects.ukrmol_co.analyze /work/runs/co-sa-target-001'
```

The default orbital ensemble is five equally weighted roots per spin/irrep,
40 C2v components in total. `--sa-singlet-roots` and `--sa-triplet-roots` each
take four counts in A1/B1/B2/A2 order; `--target-roots` must cover every averaged
root. For nonuniform target expansions, `--target-singlet-roots` and
`--target-triplet-roots` override the four computed-root counts per spin.
For example, `10 5 5 5` for each spin computes 50 C2v components; the default
`--target-states-used` retains all of them. The computed counts must cover the
corresponding orbital-ensemble counts, while scattering-channel retention
can be varied separately. Changing computed roots alone preserves the orbital
ensemble; changing `--sa-*` counts reoptimizes the common orbitals.
Equal B1/B2 counts are necessary but do not guarantee a rotationally balanced
ensemble: a truncated average may omit other degenerate partners. A failed
Pi-degeneracy check is preserved as a failed calculation.

The target builder rejects unconverged or wrong-spin roots, checks the AO
orbital metric and Pi degeneracy, and exports core/active/external blocks in
Molden order. The analyzer checks every averaged root's absolute energy to
1e-7 Hartree and verifies the imported Molden checksum. `target.json` records
the common-orbital ensemble, individual energies, spin diagnostics, **ground
root** dipole, MO inventory and source hash. Logs retain the upstream slot name
`target.psi4.out`; successful records explicitly identify PySCF as the backend.
See the [method and verification contract](../../docs/physics/co-state-averaged-target.md).

Optimization controls are `--target-energy-tolerance` (default 1e-9 Hartree),
`--target-gradient-tolerance` (1e-5), `--target-ci-tolerance` (1e-10),
`--target-max-cycles` (100), and
`--target-memory-mb` (8000). The QC backend uses `--ranks` local threads before
the MPI target stages start. These controls require refinement checks alongside
the electronic model. Passing the import checks does not qualify a fitting target.
`--target-ah-lindep` (1e-14) and `--target-ah-start-tolerance` (2.5) control the
augmented-Hessian orbital optimizer. Inspect `target-diagnostics.json` and the
QC log when iterations stall: tighter energy/gradient thresholds alone need
not overcome a dropped trial vector. Diagnostics are retained on failure.
`--target-ci-residual-tolerance` optionally specifies a separate CI residual
norm, and `--target-optimizer newton` selects the coupled orbital/CI algorithm
instead of the default `one-step`. The Newton history records a combined
orbital/CI gradient; compare it with the correct requested tolerance.
`--target-ci-lindep` (1e-14) controls removal of CI residual trial vectors;
tight residual tolerances need a compatible squared-norm cutoff. The cutoff
retry recipe is in `calibration-sa-ci-cutoff.json`.
The saved Newton controls do not meet their gradient tolerances. Use that
option for optimizer diagnostics; the successful continuation uses `one-step`.
`--target-ah-tolerance` (1e-12) sets the inner augmented-Hessian eigensolver
accuracy, separately from the outer gradient threshold. Use
`--target-verbosity 6` to retain trial-space diagnostics when it stalls;
`calibration-sa-ah-accuracy.json` supplies a tight stretched-geometry ladder.
The final-orbital fresh-CI audit records root-energy differences and spin
diagnostics in `target-diagnostics.json`. `--target-ci-fresh-start` additionally
discards warm guesses during orbital optimization; both routes retain the
independent import gates. Probe a suspect sector without reoptimizing orbitals
inside the CO image:

```bash
python3 -m projects.ukrmol_co.ci_probe \
  /work/runs/previous/output/CO/geom1/co.casscf.chk \
  --output /work/diagnostics/ci-coverage-new \
  --spin 2 --irrep A1 --roots 5 8 --spaces 40 80 160 --ranks 4
```

`--target-initial-checkpoint /work/runs/previous/output/CO/geom1/co.casscf.chk`
projects a prior target's core/active orbitals into the new AO basis. It requires
the same active-space size and irrep counts, preserves a copy in the new run,
and records its SHA256. Compare this start against the default lowest-HF-MO
selection when adding diffuse basis functions or following geometry: a lowered
ensemble energy can still trade away ground-state quality. Orbital/root identity
needs explicit overlap checks beyond comparing scalar energies.

The packaged image includes the CA bundle and Python's `SSL_CERT_FILE` setting
needed to download scripts into a fresh cache. Earlier runs used a populated
cache, which did not exercise that deployment path.

`run.py` downloads [UKRmol-scripts 1.0](https://zenodo.org/records/7851856),
verifies its published archive checksum, and copies the library and input
templates into the run. Its local library copy records stage timings and stops
on failed commands; upstream otherwise prints an error and continues.
Every run requires a **new directory**, preserving earlier inputs and results.
It selects matching executable **and shared-library** precision variants;
the upstream entrypoint defaults to double-precision libraries.

All calculation files live under the mounted `/work`, backed by `RUN_ROOT`.
Both `TMPDIR` and `PSI_SCRATCH` are placed inside
each run's `scratch/`. UKRmol also writes integrals and `fort.*` files directly
in the geometry working directory; these are included in total run disk usage.

### Source-engine evidence

The source engine is gated by the Dockerfile's `build` → `test` → `runtime`
stages, and the `co-source` target adds this experiment. `co-pilot` retains the
reference engine. The source runtime carries `/opt/ukrmolp/source-provenance.json`
and upstream logs in `/opt/ukrmolp/source-tests/`.

The source image was built successfully on Sadaharu; all **12** selected
upstream checks passed. The CO DZ39/four-frozen/one-rank comparison against
the reference-binary image gave identical printed resonance position/width
and a maximum phase difference of $10^{-7}$ rad over 0.1–5 eV. The measured
source image ID is
`sha256:599eeee919f56a27e6a0a4e214c5fc80a2f5780752f07b654b21fe1339f61ae6`.
Its source/test provenance is recorded in
[`source-build-results.json`](source-build-results.json).

The reorganized Docker targets were independently rebuilt and passed the same
12 checks. A fresh-cache, one-rank DZ39 calculation passed after adding the CA
bundle/Python certificate-path fix, and again matched the reference position/
width at printed precision with a maximum phase change of `1e-7 rad`.
[`packaging-results.json`](packaging-results.json) records this build's image
and Dockerfile/helper hashes and the 61.7-second calculation timing.

## Finite remote batches and fit replays

Each job manifest is a list of named jobs and runner arguments.
`batch.py` runs them through Docker on the host, snapshots the project source
read-only, records its SHA256 hashes and image ID, and uses disjoint CPU slots:

```bash
python3 -m projects.ukrmol_co.batch projects/ukrmol_co/calibration-throughput.json \
  --root "$RUN_ROOT" --image qmodeling/ukrmol-co:source \
  --cpu-groups 0 1 2 3 --memory 4g
```

Manifest names and run names must be new for each invocation. Results/logs
are under `batches/<manifest-stem>/`, calculations under `runs/<job-name>/`.
This finite job interface supports SSH execution now; an HTTP service and AWS
deployment have not yet been implemented.

Inside the image, replay an existing fit with:

```bash
python3 -m projects.ukrmol_co.refit /work/runs/co-eq-f4-sep39 --nback 3
```

The local QSCAT workspace collector creates the calibration JSON snapshot
from copied lightweight artifacts (`runs/` and `batches/`):

```bash
uv run python -m projects.ukrmol_co.collect /path/to/copied-artifacts \
  --pairs projects/ukrmol_co/calibration-pairs.json \
  --provenance projects/ukrmol_co/campaign-provenance.json \
  --output projects/ukrmol_co/calibration-results.json
```

Use the saved provenance only when rebuilding the 6 October campaign; supply
your own metadata for a new campaign. The collector does not assume a host.
Grid-refinement comparisons use shared printed native energies and record the
compared point count when grids differ. They do not interpolate through a
feature absent from a coarser grid; fewer than two common points is an error.

### Experiment files

| Files | Role |
|---|---|
| `co.pl`, `run.py` | Molecular deck, checksum-pinned scripts, strict PK SCF, workspace/precision selection and profiling |
| `analyze.py`, `test_analyze.py` | Pipeline/symmetry checks, native fits and a CC column-block regression |
| `batch.py` | Finite fresh-container execution with CPU slots and source/image provenance |
| `refit.py` | RESON replays from retained channels, amplitudes and K-matrices |
| `pilot-jobs.json`, `pilot-results.json` | Original pilot recipe and quoted results |
| `calibration-*.json` | Calibration job recipes, paired comparisons and the full results snapshot |
| `collect.py`, `campaign-provenance.json` | Rebuild that snapshot from saved lightweight evidence |
| `target.py`, `state_average.pm`, `sa-results.json` | Common state-averaged orbitals, explicit upstream adapter and continuation evidence |
| `ci_probe.py` | Fixed-orbital CI root-count/trial-space coverage diagnostic |
| `archive.py`, [`sa-evidence/`](sa-evidence/README.md) | Deterministic continuation archive, including orbital checkpoints and source snapshots |
| `source-build-results.json` | Original source-engine pins, upstream gate and CO differential evidence |
| `packaging-results.json` | Reorganized targets, cold-cache acquisition fix and independent CO differential check |
| [`evidence/`](evidence/README.md) | Public checksum-verified lightweight archive and snapshot reconstruction instructions |
| [`CONTINUATION.md`](CONTINUATION.md) | Specific host/artifact locations, campaign replay map and next electronic-model qualification gate |

The runner evolved during calibration. Manifests preserve the requested
arguments; per-run `config` and generated inputs, per-batch source hashes and
available source snapshots preserve the historical execution. The first four
batches have hashes but no source snapshot; their decks/templates are retained.
In particular, replaying early CAS manifests
with today's raised workspaces will not recreate the old CONGEN limit failures.
Timing results also depend on CPU slots, concurrency and image, not just the
electronic configuration.

## Artifacts and checks

Each run preserves:

- `config.json`, `co.pl`, the upstream driver and its instrumented library;
- `run.log` and all generated UKRmol/Psi4 inputs and outputs;
- `stages.tsv`: task, program, spin, symmetry, exact wall seconds, raw wait status;
- `resources.csv` and `resources.json`: sampled aggregate container memory and
  allocated run/scratch disk usage, plus kernel cgroup peak memory and aggregate
  CPU time in newer runs;
- `output/CO/geom1/eigenph.doublet.{B1,B2}` and
  `xsec.doublet.{B1,B2}`;
- channels, R-matrix amplitudes, K-matrices, and integrals for further analysis.

The analyzer checks the complete energy grids, finite nonnegative cross
sections, `B1`/`B2` degeneracy, and consistency of the RHF target energies.
For ground-state CASSCF it instead checks the correlated target energy against
Psi4's final CASSCF energy. CC cross-section columns are joined across native
output blocks, and the total is checked against the final-state sum.
Eigenphase agreement is checked modulo $\pi$ with a 0.002-radian absolute
tolerance because the native text output has limited precision. Native
resonance fits are recorded only when positive finite positions and widths
are present. An empty fit list means the automatic fit found no usable result.

Memory is sampled every 0.2 s and disk every 2 s, so reported peaks are sampled
lower bounds. Container memory includes filesystem page cache and any other
processes in that container. Newer runs also include `memory.peak`, the kernel's
peak for the container lifetime. Batch workers use one fresh container per
calculation, so each peak belongs to one job. Files remain available after the
calculation exits.

A successful pilot establishes the execution pipeline and internal consistency.
Before using its resonance curve as a potential-fitting target, vary the
virtual space, continuum angular cutoff, sphere radius and deletion threshold,
and inspect eigenphases and resonance-fit stability. The compact RHF/SEP
neutral energy is not a correlated neutral potential curve.
