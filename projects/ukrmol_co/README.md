# UKRmol+ CO fixed-nuclei scattering experiment

**Status:** the source-built engine and common state-averaged target/import
pipeline work. Equilibrium angular-cutoff and stretched energy-grid checks
meet the chosen tolerances for the tested CAS(10,10) model. Electronic-model
and resonance-extraction qualification remain open, so production
potential-fitting targets are not yet qualified.

The completed evidence preserves the first pilot and two calibration snapshots:

| Campaign | Attempts | Validated pipelines | Preserved failures | Comparisons | Public evidence |
|---|---:|---:|---:|---:|---|
| Original SEP/ground-state-CASSCF calibration | 40 | 30 | 10 | 21 | [`evidence/`](evidence/README.md) |
| Multi-spin/multi-irrep state-averaged continuation | 55 | 29 (17 target-only, 12 scattering) | 26 | 12 | [`sa-evidence/`](sa-evidence/README.md) |

The continuation also includes 96 native RESON replays and six fixed-orbital
CI coverage probes. “Validated” means solver/import and pipeline consistency,
not electronic-model convergence. The CAS(10,11) fresh-CI restart passes all
40 averaged-root imports and the ground-dipole check, but active-space changes
remain substantial. The CAS(10,12) memory retry passed QC/fresh-CI checks and
singlet-A1 diagonalization, then exceeded its triplet-A1 memory budget.
[`CONTINUATION.md`](CONTINUATION.md) records their host state and next gates.

## Convergence verdict

| Question | Finding | Decision |
|---|---|---|
| Does the engine/workflow reproduce reference calculations? | Five water target energies; 12 source-build serial/MPI checks; CO source/reference phases agree to `1e-7 rad` | Execution established |
| Is the compact equilibrium SEP continuum stable? | Tested angular/deletion changes shift position by less than 0.5 meV and width by less than 0.1%; background fits vary width by about 2.5% | Passes chosen numerical tolerances for this model/window |
| Are the tested SA numerical controls stable? | Equilibrium l=3/4/5 stays within 3.18 meV, 0.585% width and 0.02615 rad; stretched 0.020/0.010/0.005-eV grids give a finest-pair width change of 0.0041% | Passes chosen angular/grid gates for these models |
| Is the electronic model converged? | Basis, starting subspace, active space and orbital ensemble change target properties and scattering; projected aug-DZ/TZ widths differ by 6.36% | Further qualification required |
| Is the entire geometry curve qualified? | SA compressed/stretched scattering passes pipeline checks; compressed fit widths remain background-sensitive and feature identity across geometry is unresolved | Extraction, continuity and threshold diagnostics required |
| Can these data constrain a production potential? | No tested electronic configuration or correlated neutral curve is qualified | Preserve as calibration evidence |

We have a reproducible execution platform and identified the dominant accuracy
problem. The resonance values have **not** settled onto an electronically
converged limit. The next useful refinement is to finish the larger-active-space
import checks and compare basis, ensemble and channel refinements before
choosing a production geometry sweep.

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
55 completed attempts: seventeen successful target-only pipelines, twelve successful
scattering pipelines, and 26 preserved setup/diagnostic failures. One rejected
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
| cc-pVDZ CAS(10,11), fresh-CI reoptimization | −112.924284941 | 6.3851 | 0.0277316 | 5.2e-10 |

Every averaged root passed the 1e-7-Hartree import gate, and independent
ground dipoles agree within 1.9e-6 a.u. The default CAS(10,10) DZ neutral energy,
lowest triplet excitation and dipole (about 0.234 D) are close to the published 40-state
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
fifth root matches the fresh sixth root. The original target remains rejected: the
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
5.44e-10 Hartree.

The CAS(10,11) fresh-CI restart completed in **127.80 minutes**, with an
18.41-GiB kernel memory peak and a 303.29-MiB sampled run-disk peak. Independent
reanalysis passes all 40 averaged-root imports (maximum error 5.20e-10 Hartree),
computed Pi partners and the ground dipole (difference 8.10e-11 a.u.). Fresh-CI
energies agree within 1.14e-13 Hartree; the QC Pi splitting is 5.29e-11 Hartree
and MO orthogonality error 2.89e-15. Its nine-iteration orbital optimization
reaches gradient 4.77e-6 at the requested 1e-5 tolerance. This resolves the
observed missing-root mismatch for this restart; globally preferable orbitals
and electronic-model convergence remain unestablished.

The balanced 40-component DZ targets give the following active-space comparison:

| Active space | Ensemble energy (Hartree) | Ground energy (Hartree) | Lowest triplet Pi (eV) | Ground dipole z (a.u.) | Requested gradient tolerance |
|---|---:|---:|---:|---:|---:|
| CAS(10,8), fresh-CI control | −112.303519400 | −112.855375650 | 6.4858 | 0.2023325 | 1e-7 |
| CAS(10,10), tight equilibrium target | −112.344874513 | −112.894737475 | 6.3952 | 0.0920732 | 1e-6 |
| CAS(10,11), fresh-CI restart | −112.370969233 | −112.924284941 | 6.3851 | 0.0277316 | 1e-5 |

CAS(10,10)-to-(10,11) lowers the ground root by 0.80403 eV and changes the
dipole by −69.88%, despite only a −10.10-meV triplet-Pi excitation change.
Initial/final CAS(10,11) active-overlap singular values span 0.94139–1.00000.
The differing optimizer tolerances and unresolved orbital identities prevent
treating this as a certified convergence sequence. Tighten and compare starts
before selecting the electronic model; stable excitation energies alone do not
qualify its neutral curve, dipole or scattering.

CAS(10,12) QC passed its original gates, but MPI-SCATCI stopped before target
diagonalization: its 43194-dimensional singlet-A1 matrix needed a 3.74-GB
local block, exceeding the template's 2.5-GiB internal budget per process.
The 32-GiB container limit does not enlarge that budget. `--scatci-memory-gib`
sets `memp` in both target and scattering templates; its default remains 2.5.
The checkpoint-based retry used 6 GiB per process in an 80-GiB container. QC
and the fresh-CI audit passed; singlet-A1 diagonalization completed in 95.86
minutes. Triplet-A1 then failed its internal budget: its dimension is 70674
and each of four ranks requested a 9,994,717,728-byte (9.31-GiB) matrix block.
The run exited with code 25 after 98.80 minutes, with a 56.10-GiB kernel
memory peak and a 144.79-MiB sampled run-disk peak. It remains an engine
failure, without full target-import/dipole validation. Choose and measure
the process/container budgets together, including diagonalization workspace;
the matrix allocation alone is not a total-memory estimate.

## Measured cost and parallel execution

### Electronic qualification tools

`--qc-only` runs the state-averaged optimizer and fresh-CI/spin/Pi/MO checks
without UKRmol target diagonalization. These records are `qc_validated`;
`--target-only` retains the independent UKRmol root/dipole checks. The
CAS(10,11) tight-start recipes are in `calibration-sa11-qc-tight.json`, with
the selected import recipe in `calibration-sa11-tight-import.json`.
All three tightened starts pass in 6.08/13.13/61.89 minutes. Their ensemble
objectives agree within 5.69e-14 Hartree, individual roots within 7.70e-9
Hartree, dipoles within 5.14e-8 a.u., and minimum final active-subspace overlap
singular values exceed 0.99999999999992. Forty-six of 48 fixed-orbital
root-count/trial-space probes converge completely; the eighth singlet B1/B2
root fails at space 40 and converges at 80/160. All five ensemble roots converge
in every probe. The completed refinements support stability of this fixed
orbital model. Its tightened independent UKRmol import now passes all 40 roots
within 5.33e-10 Hartree and the ground dipole within 7.95e-11 a.u., with a final
gradient of 2.92e-8. It took 97.75 minutes at an 18.42-GiB kernel memory peak.
The subsequent equilibrium scattering pilot completes in 4.423 hours at a
25.20-GiB peak. Both contracted scattering dimensions are 27546, against raw
CONGEN counts of 344124. Its default candidate is 2.520805 eV with full width
1.151002 eV. Relative to the earlier CAS(10,10) model, position/phase changes
of 64.84 meV / 0.109809 rad exceed the chosen gates; width changes by 2.80%.
`calibration-sa10-tight-scattering.json` now completes with matching tight
numerical settings and CAS(10,10)'s own checkpoint in 21.38 minutes at a
2.864-GiB peak. All forty imported roots/dipole and 48 lowest-root probes pass.
Tightening changes the earlier CAS(10,10) position/full width by just
−0.544/+0.395 micro-eV and phases by 1.30e-6 rad. The matched active-space
comparison still fails position/phase gates (+64.8401 meV / 0.1098099 rad).
`qualification-results.json` and
[`qualification-evidence/`](qualification-evidence/README.md) preserve the
completed target/neutral supplement and its verified public archive.

The qualification supplement also preserves four rejected serial-Davidson
targets and two passing small-model SLEPc targets. The latter reproduce all
40 required roots within 4.98e-10 Hartree and dipoles within 5.62e-11 a.u.,
including extra-root/tighter-tolerance checks, using about 0.12 GiB in
49.53/46.91 seconds. The full five-root CAS(10,11) SLEPc import now passes all
40 roots within 5.34e-10 Hartree and dipole within 8.04e-11 a.u., in 52.55
minutes at a 5.25-GiB kernel peak. Its dense-reference excitation energies
match at printed precision. The observed SCATCI/full-wall reductions are
89.17%/46.25%; DENPROP still takes 39.42 minutes. The eight-root/tighter control
now also passes all forty ensemble-root/dipole imports in 51.44 minutes at a
5.27-GiB peak. Sixteen fully converged fixed-orbital CI probes at spaces 80/160
independently check all 64 requested roots within 5.34e-10 Hartree. Common dense
and five-root excitations match at printed precision. This completes the
CAS(10,11) extra-root solver gate; the local CAS(10,12) import has started after
its recorded controls and CPU-slot release. The verified supplement now preserves
46 attempts in eighteen batches and 160 raw CI probes. Both small numerical
80/160-vector controls pass all forty imported roots/dipoles and have identical
printed excitations. The completed four-entry failed-checkpoint repair batch
has two QC passes and two retained CI failures.
Both stretched-DZ CAS(10,11) 80/160-vector repairs now pass QC; both projected
TZ retries still fail singlet B1/B2 CI convergence after orbital convergence.
Their paired roots/objective/dipole/subspaces agree, but three five-root
space-40 probes fail an ensemble-root convergence flag. Both repairs reject
the strict all-probe coverage gate despite all 64 space-80/160 probes passing.
The [coverage companion](qualification-evidence/stretched-coverage/README.md)
preserves those 96 probes and confirms matched CAS(10,10) orbital continuity.
The [iteration follow-on](qualification-evidence/stretched-iteration/README.md)
repairs all nine failed settings at 600 cycles and passes full 48-probe scans
of both checkpoints with explicit physical/spin-penalized residuals below
7.07e-10/9.94e-10 Hartree. Independent stretched import and same-geometry
TZ/aug-TZ QC are now queued under that refined gate.
Forty-nine saved-data native replay attempts preserve 48 successes and one B2
unit-binding failure with a fresh corrected retry. Window clipping at fixed
background meets the chosen gates; background terms 1–4 change width by 12.65%,
so resonance-extraction qualification remains open.
`phase_fit.py` and `phase_diagnostics.py` independently reproduce all 24 native
window fits within 0.48/0.25 micro-eV in position/full width and audit printed
residues against the saved grid. Analytic unitary-S tests and interleaved
held-point fits show that the constant background describes the data less well
than terms 2–4. The pinned engine uses 0.0735 Ryd per requested eV; its 18.44-ppm
energy-grid convention is explicitly retained at this boundary. See
[the phase-fit contract](../../docs/physics/co-electronic-qualification.md#independent-phase-fit-diagnostic).
The fresh equilibrium TZ/aug-TZ QC batch completes with three passes and the
original CAS(10,11) aug-TZ fresh-CI rejection, retained in the
[fresh-basis companion](qualification-evidence/fresh-basis/README.md). Fresh
CAS(10,12) TZ→aug-TZ lowers the averaged objective by 1.80 eV while raising
the ground-state energy by 0.480 eV; root coverage, competing starts and
independent imports are still required. Finite coverage/residual scans and
same-basis fresh-lineage CI-space repairs are now supervised on Sadaharu.
The [fresh-target coverage companion](qualification-evidence/fresh-coverage/README.md)
now passes all 144 probes/936 eigenpair evaluations, with physical/penalized
residual maxima 6.96e-10/9.99e-10 Hartree, in 32.36 minutes / 0.728 GiB.
The covered 64-root fresh CAS(10,11) TZ import now passes within 3.196e-8
Hartree against independent seed-orbital controls, in 54.38 minutes / 5.270 GiB.
The [import companion](qualification-evidence/fresh-tz-import/README.md) retains
the passing engine run, two verifier failures and the passing final-set/decimal
recheck. Competing-start and electronic-model gates remain open.
Both fresh-lineage CAS(10,11) aug-TZ space-80/160 QC repairs now pass in
10.69/10.64 minutes at 0.419/0.477 GiB. All forty roots/dipoles agree within
1.667e-11 Hartree / 3.230e-12 a.u. and core/active subspaces match to roundoff.
The [repair companion](qualification-evidence/fresh-augtz-repairs/README.md)
retains both successes and the exact rejected parent. Ninety-six independent
coverage/residual probes now pass, with 624 eigenpair evaluations, in
10.02 minutes / 0.391 GiB; the
[coverage companion](qualification-evidence/fresh-augtz-coverage/README.md)
preserves the raw spectra and residual gates. The
[64-root aug-TZ import](qualification-evidence/fresh-augtz-import/README.md)
now passes against both repaired seeds within 5.969e-11/5.764e-11 Hartree,
with dipole/subspace gates accepted, in 56.01 minutes / 5.268 GiB. Four
same-geometry TZ↔aug-TZ space-80/160 projections now test competing starts;
electronic-model qualification remains required.
Outer-only 99→197→393-point checks also pass: identical common-point phases and
fixed-window position/full-width changes at most 0.145/0.292 meV (0.0254% in
width), at 6.37 minutes / 0.719 GiB. Native fine-grid fits hit `MAXFIT=100`;
their truncation is retained and the independent diagnostic uses full intervals.
Radial propagation also passes 8→16→32-subrange and 100→150→200-bohr matching
radius controls. A full-precision binary K audit verifies every requested sector
and limits baseline phase changes to 2.57e-6 rad, with successive refined
differences below 2.28e-10 rad. Binary-phase fits change position/full width by
at most 2.63e-8/4.35e-8 eV. The forty-stage diagnostic takes 27.83 minutes /
1.868 GiB; its [public companion](qualification-evidence/outer-propagation/README.md)
preserves all binary outputs, rank logs and the reconstruction.
The finite continuation also queues `calibration-sa11-tight-continuum.json`:
l=3/l=5 controls against the completed tight l=4 scattering baseline, using the
qualified selected-root target solver and the existing scattering diagonalizer.
Their baseline/extra-root gates pass, and a finite successor now advances them
on CPUs 12–15 while CAS(10,12) import continues on 8–11. Each passing run receives
saved-K-matrix background/detection
replays; [the live handoff](CONTINUATION.md#local-continuation--7-october-2026)
records the exact queues, resource gates and new fresh-TZ/CI observations.
This target route stores a dense PETSc Hamiltonian; its largest CAS(10,12)
matrix floor is 37.41 GiB. The prepared local limit is 64 GiB, with 16-GiB
internal budgets and an 80-GiB available-host-memory gate.

Both tight CAS(10,12) QC starts pass: restart/RHF walls are 21.17/152.18
minutes. Objectives differ by 2.84e-14 Hartree, roots by at most 4.17e-8
Hartree and dipoles by 6.27e-8 a.u.; the minimum active-subspace overlap
singular value is 0.999999999999496. All five ensemble roots converge in all
48 fixed-orbital probes; extra singlet A1/A2 roots fail at space 40 and
converge at 80/160. Independent import remains required.

The first nine-job CAS(10,11) ladder has two QC passes and seven failures:
three CI-convergence failures and four simultaneous geometry/basis projection
rejections. The 50-component DZ ensemble raises the ground root by 0.45651
eV and changes the dipole by +0.0346869 a.u. (+125.09%). This is an orbital
objective change, distinct from channel retention. Staged same-geometry basis
projections and CI-space refinements are recorded in `CONTINUATION.md`.

The neutral pilot uses `python -m projects.ukrmol_co.neutral` or the existing
batch launcher with `--runner neutral`. Its aug-TZ/aug-QZ recipes are in
`calibration-neutral-pilot.json`; `neutral-results.json` and
`neutral-provenance.json` preserve its completed attempts, including the initial
failed Psi4 Python import and the passing executable-based control. Psi4 has
its own Python environment in the pinned toolchain.

The independent cc-pVDZ control agrees in RHF/CCSD(T) energy within
1.5e-13/8.3e-11 Hartree. The original seven neutral records pass convergence, reference
stability and CCSD-density checks; this is not a qualified neutral curve.
The six basis-pilot calculations finished in 90.93 seconds on two disjoint
four-core workers. Aug-TZ/QZ jobs take about 5–7/38–41 seconds, with kernel
memory peaks near 0.50/3.21 GiB.

| R (bohr) | Relative CCSD(T) energy, aug-TZ (eV) | Relative CCSD(T) energy, aug-QZ (eV) | CCSD dipole z, aug-QZ (a.u.) | aug-QZ amplitude norm / sqrt(10) | aug-QZ largest amplitude singular value |
|---|---:|---:|---:|---:|---:|
| 1.9000 | 1.35893 | 1.28289 | 0.2055393 | 0.0123970 | 0.0247273 |
| 2.1323 | 0 | 0 | 0.0501085 | 0.0178684 | 0.0372068 |
| 2.5000 | 1.32736 | 1.39583 | −0.1925559 | 0.0316658 | 0.0675234 |

Each relative curve is referenced to its own R=2.1323 energy. Basis changes
shift the sentinel relative energies by 76.04/68.47 meV; the stretched
amplitude diagnostics increase. Basis and correlation-treatment checks remain
open, and the dipoles are CCSD lambda-density values, not CCSD(T) derivatives.
`calibration-neutral-5z.json` completes the same three aug-cc-pV5Z sentinels
in 18.97 minutes on one four-core worker. Relative 5Z energies at R=1.9/2.5
are 1.26565/1.41165 eV, each referenced to its own equilibrium energy.
QZ-to-5Z shifts are −17.24/+15.82 meV. The all-memory peaks reach the 24-GiB
container limit, including page cache; sampled anonymous peaks are about
17.14 GiB and run-disk peaks 12.81–14.24 GiB.
`calibration-neutral-stretched.json` tests R=3.0/4.0 in aug-TZ/QZ. All four
attempts fail external RHF stability, preserving converged RHF diagnostics and
original exits; a different correlation treatment is required for that region.
The current neutral aggregate retains fifteen attempts: ten passing records,
the initial failed import and four rejected unstable stretched references.
See [the qualification method](../../docs/physics/co-electronic-qualification.md)
for the numerical contract and the CAS(10,12) workspace audit. Its valid
eight-/sixteen-rank array floors are about 149 GiB, exceeding the current
host RAM; the four-rank triplet query is undersized.

The prepared larger-memory campaign is documented in
[`LARGE_HOST.md`](LARGE_HOST.md), with matched QC recipes, ten CAS(10,12) target
imports and the prerequisites for the first scattering pilot. The selected
Frankfurt `r8a.16xlarge` has 64 physical x86-64 cores and 512 GiB RAM. Its initial
full-campaign forecast is 168 instance-hours, with a provisional 96–216-hour planning range;
the documented forecast separates measured stages from unmeasured MPI and
larger-space costs. Paid provisioning is deferred: qualification continues on
Sadaharu, and the first external-host experiment must fit approximately $200,
at least 80% below that earlier forecast. Selected-root target diagonalization
is being investigated against dense controls before claiming local feasibility.

### Scattering model costs

The tight **CAS(10,10)/cc-pVDZ/40-channel** matched benchmark now passes all
seven replicas, with both contracted scattering dimensions 8350 and 99 energies:

| MPI ranks | Runner wall (s) | Aggregate CPU (s) | Kernel memory peak (GiB) | Full-job speedup |
|---:|---:|---:|---:|---:|
| 1 | 3105.91 | 3105.86 | 2.849 | 1.000× |
| 2 | 1830.62 | 3435.09 | 2.743 | 1.697× |
| 4 | 1239.61 | 4268.46 | 2.822 | 2.506× |

The sequential MPI trials use CPUs 12–15; four one-rank replicas pinned
individually to those same physical cores finish together in **3444.36 seconds /
57.41 minutes**, **1.440×** the throughput of repeating the measured four-rank
job. Their summed per-job memory peaks give an **11.342-GiB upper envelope**.
All forty roots/dipoles and orbital subspaces match the baseline; phases differ
by at most 2e-7 rad and position/full width by 0.137/0.041 micro-eV. The
[benchmark companion](qualification-evidence/mpi-throughput/README.md) retains
all seven runs, frozen sources, raw profiles and the completed CPU-slot lease.
These are single-shot timings with recorded concurrent host work; larger models
need their own measurement. The continuum supervisor resumes after this trial.

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
has not yet completed its own one/two/four-rank and concurrent-geometry scaling
sweep. The matched CAS(10,10) scaling recipe is now queued on the same four-core
allocation, followed by four individually pinned one-rank replicas; see
[`calibration-sa10-tight-mpi-scaling.json`](calibration-sa10-tight-mpi-scaling.json)
and [the live continuation](CONTINUATION.md#local-continuation--7-october-2026).

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

The default `--target-diagonalizer auto` uses the upstream target dispatch.
Experimental `--target-diagonalizer slepc` requires an image built with
`--build-arg WITH_SLEPC=ON`; it selects distributed Krylov–Schur for the target
while retaining the existing scattering solver. `--target-diagonalizer-tolerance`
and `--target-diagonalizer-max-cycles` default to 1e-12 and 500.
The raw-log analyzer checks the solver identity and completed sectors before
the usual root/dipole import gates. `davidson-serial` preserves the rejected
controls; it failed the required spectrum even with additional roots and a
tighter tolerance. See [the selected-root contract](../../docs/physics/co-electronic-qualification.md#selected-root-target-experiment).

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
`--target-ci-max-space` separately refines the PySCF CI trial-vector space;
its default remains `max(40,8*ensemble roots)` per sector. The trial-space
override must exceed every ensemble root count and is recorded by sector.
The matched CAS(10,11) ladder's stretched DZ triplet-A1 and equilibrium TZ
singlet-B1/B2 CI failures motivate 80/160-vector controls and fresh retries.
These change a numerical solver limit, while retaining the orbital objective
and its import/convergence gates. Recipes and queued checks are recorded in
[`CONTINUATION.md`](CONTINUATION.md).
For geometry/basis changes together, stage through a passing DZ checkpoint at
the requested geometry, then change basis at fixed geometry. The pinned PySCF
projection helper rejects a simultaneous change; the original four rejections
remain calibration evidence.
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
