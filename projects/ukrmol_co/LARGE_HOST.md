# CO large-memory qualification campaign

**Paid provisioning is deferred.** Continue qualification on Sadaharu and keep
the first paid experiment within **approximately $200**, reducing the earlier
full-campaign compute forecast by at least 80%. The retained hardware candidate
is EC2 `r8a.16xlarge` in `eu-central-1`; its **168-hour / $1036 full-campaign
estimate below is historical planning, not the current execution budget**.
The queued recipes are inputs, not passing calculations.

## Cost-constrained continuation

Complete the tight CAS(10,12) QC starts, matched basis/geometry/ensemble QC
ladder, fixed-orbital root coverage and all feasible CAS(10,11) imports and
scattering on Sadaharu. Also finish local neutral refinements and saved-data
fit replays. These results determine which single larger calculation has the
most value; the ten-entry dense target manifest is a retained refinement queue.

Only five roots per sector are needed by the 40-component target contract;
the dense path allocates all eigenvectors. Four initial serial-Davidson
controls (5/8/16/32 roots) used about 0.12 GiB and 46–48 seconds each, but all
failed required root imports despite native convergence reports. Their maximum
energy errors were 0.13931/0.07915/0.05290/0.02823 Hartree; retain the failed
controls. Independent diagonalization of the first stored singlet-A1
Hamiltonian agrees with QC within 3.69e-13 Hartree, isolating the discrepancy
to the selected-root solve rather than that Hamiltonian.

The SLEPc-enabled source build now passes all twelve upstream serial/MPI
checks, with pinned toolchain PETSc/SLEPc library digests. Its explicit
distributed Krylov–Schur small-model controls both pass all 40 roots within
4.98e-10 Hartree and dipoles within 5.62e-11 a.u., using about 0.12 GiB in
49.53/46.91 seconds. Its full five-root CAS(10,11) dense-reference comparison
passes all 40 roots within 5.34e-10 Hartree and dipole within 8.04e-11 a.u.,
in 52.55 minutes at 5.25 GiB. Observed SCATCI/full-wall reductions are
89.17%/46.25%; DENPROP dominates the remaining time. Its eight-root/tighter
refinement now passes in 51.44 minutes at 5.27 GiB; sixteen independent CI
probes check all 64 requested roots within 5.34e-10 Hartree. See
[the solver contract](../../docs/physics/co-electronic-qualification.md#selected-root-target-experiment).
A selected-root target solver does not remove the all-spectrum requirement
of the existing scattering calculation.

The extra-root gate and CAS(10,12) restart coverage release a 64-GiB local target
import. It now completes with all forty native roots and the ground dipole
passing: **6.201 hours / 40.405 GiB**, current-QC native-root error **4.801e-11
Hartree**, covered-seed error **6.815e-9 Hartree**, independent RHF-start error
**3.485e-8 Hartree**, current-QC dipole error **1.279e-11 a.u.** The
[CAS(10,12) companion](qualification-evidence/cas12-import/README.md) verifies
all 335 payloads and raw spectra/profiles. The restart remains an experimental
model input; competing starts remain part of electronic qualification.
Separately, two matched CAS(10,11) QC entries fail CI convergence after orbital
convergence. Small-model 80/160-vector trial-space controls now pass, and fresh
reoptimizations of those rejected checkpoints have started; their original
failures remain evidence.

The selected-root target route retains a dense PETSc Hamiltonian. The largest
CAS(10,12) matrix floor is 37.41 GiB (about 9.35 GiB per rank on four ranks).
Its local recipe now uses 16-GiB internal budgets; the 64-GiB container trial
requires 80 GiB available host RAM and 20 GiB free scratch before launch.
The completed job measures a **40.405-GiB** complete-job kernel peak. Internal
and container limits remain separate, above the matrix floor.

Both CAS(10,12) tight starts and the restart's ensemble-root coverage pass;
the common final active subspaces agree to a minimum overlap singular value
of 0.999999999999496. A finite staged CAS(10,12) QC ladder is running locally.
The nine-job CAS(10,11) ladder finishes with two passes, three CI failures and
four rejected simultaneous geometry/basis projections. New basis starts use
passing same-geometry DZ checkpoints. The retained original manifests remain
historical recipes; newly derived staged manifests record their own hashes.

The neutral TZ/QZ/5Z relative-energy sequence now shifts by 17.24/15.82 meV
for QZ-to-5Z at R=1.9/2.5. At R=3.0/4.0, all four aug-TZ/QZ references fail
external RHF stability. The restricted neutral route is blocked there pending
a different validated correlation treatment; increasing its basis alone does
not resolve that reference defect.

The matched local CAS(10,10) scattering benchmark now measures **3105.91 /
1830.62 / 1239.61 seconds** on one/two/four ranks within the same four physical
cores. Four concurrent one-rank jobs finish in **3444.36 seconds**, **1.440×**
the throughput of repeating the measured four-rank job, with **19.90% less
aggregate CPU**. Their summed per-job peaks bound the memory envelope at
**11.342 GiB**. All seven replicas pass numerical equivalence; the
[public companion](qualification-evidence/mpi-throughput/README.md) retains
inputs, host-workload samples, profiles and the reconstruction. Use this measured
layout for independently qualified calculations at that compact model size.
CAS(10,12) now has a measured four-rank target memory/stage baseline: eight
SCATCI sectors **46.81 minutes**, serial DENPROP **5.009 hours / 80.768%** of
wall. Its full scattering memory/scaling remains unmeasured; the 64-core forecast
and approximately $200 first-experiment cap remain conditional. Increasing MPI
ranks alone does not accelerate the measured serial density bottleneck.

Competing larger-basis starts also expose a genuine electronic-model gate:
two reproducible TZ-projected aug-TZ targets have an objective **0.538388 eV
lower** than the earlier fresh aug-TZ branch, ground shift **−1.115756 eV**,
dipole shift **−0.06670992 a.u.** and minimum active overlap **0.141085**.
The [competing-start companion](qualification-evidence/competing-starts/README.md)
retains the failed restart verdicts and two projection-interface exceptions.
New-branch coverage passes and tested downward-projection retries continue locally;
numerically passing import alone does not prequalify a unique orbital model for
the paid experiment.

At the recorded Frankfurt rates, a **24-hour `r8a.16xlarge` window costs
$148.01 in compute**; a **48-hour `r8a.8xlarge` window also costs $148.01**.
The 256-GiB instance can accommodate one audited dense target worker, while
the 512-GiB instance provides more room for the scattering pilot. Storage and
transfer need to fit the remaining allowance. Neither window is yet a measured
completion forecast. A first paid experiment should cover one prequalified
equilibrium calculation, with its required numerical checks, after local work
establishes that it fits the time/memory budget. Other hardware options can use
the same x86-64 Docker recipes and resource audits.

## Scientific purpose

CAS(10,12) distributes ten valence electrons among twelve spatial orbitals,
with the two core orbitals doubly occupied. The requested C2v active counts
are `[6,3,3,0]` in A1/B1/B2/A2 order. It adds another sigma-like A1 orbital to
CAS(10,11), retaining the paired Pi spaces and the 40-component orbital ensemble.
That enlarges target correlation and the short-range scattering configuration
space. It tests whether the smaller active spaces omit important flexibility.

The existing DZ CAS(10,10)-to-(10,11) comparison lowers the ground root by
0.80403 eV. The provisional CAS(10,11)-to-(10,12) QC comparison lowers it by
another 0.08043 eV and changes the ground dipole from 0.0277316 to 0.0175114 a.u.,
about 36.9%. These values use the older optimizer tolerances; tightened,
matched-start and orbital-identity checks must precede a convergence claim.
The relative dipole change is large partly because the dipole is small; retain
absolute changes as well. None of these target changes is itself a measured
change in resonance position or width.

We need an active-space convergence test. CAS(10,12) is the next controlled
test, not a predetermined production model or a guarantee of convergence.
QC already evaluates this active space using iterative CI on Sadaharu. The
large-memory requirement remains for the all-spectrum dense route and any
selected-model scattering calculation. The qualified selected-root target
route now fits Sadaharu at its measured 40.405-GiB peak.

## Architecture and execution layout

[AWS's instance specifications](https://docs.aws.amazon.com/ec2/latest/instancetypes/mo.html#mo_summary)
identify R8a as **AMD x86_64**. `r8a.16xlarge` uses the AMD EPYC 9R45
(fifth-generation EPYC/Turin) and has **64 vCPUs, 64 physical cores, one thread
per core and 512 GiB RAM**. The
[R8a product page](https://aws.amazon.com/ec2/instance-types/r8a/) confirms this
mapping and AVX-512 support. Select an amd64/x86_64 Linux AMI, such as Ubuntu
24.04 LTS. This matches the pinned Linux x86-64 compiler/MPI/Psi4 toolchain;
use the existing Docker build and engine gates on the new host.

The audited maximum target array floors are **149.78 GiB on 16 ranks** and
**149.91 GiB on 32 ranks**. These include both full matrices, workspaces and
eigenvalues, and exclude other live engine arrays and MPI/library overhead.
Increasing ranks distributes memory but does not remove its aggregate floor.
The four-rank triplet query is invalid; preserve its overflow rejection.

Start target jobs with two disjoint **32-physical-core** slots and **224-GiB
container limits**; this leaves 64 GiB outside their combined limits. Recipes
set `--scatci-memory-gib 16` for the internal matrix budget independently of
the container limit. Confirm actual process grids and library queries on EC2,
and measure a 16-/32-rank sector comparison before fixing the throughput layout.
Inspect CPU/NUMA topology rather than assuming the logical ID ordering.

Scattering has an additional electron and continuum configurations; its memory
budget must be checked from actual CONGEN dimensions. Initially allow one
scattering worker with a provisional 384-GiB container limit, subject to a valid
workspace query and headroom check. Two concurrent target jobs fitting in RAM
does not establish that two scattering jobs fit.

The tight CAS(10,11) scattering run on Sadaharu completes in 4.423 hours at
a 25.20-GiB kernel peak. B1/B2 solve 27546-dimensional **contracted** Hamiltonians
in 73.12/74.92 minutes. Each raw CONGEN count is 344124; it is not the dense
diagonalizer's dimension. The earlier B1-only 24.18-GiB sampled observation is
retained in the diagnostic history. Preserve both counts and query the actual contracted
CAS(10,12) dimension before applying a dense-memory or cubic-time forecast.
The matched tight CAS(10,10) control completes locally in 21.38 minutes at
2.864 GiB, with all forty imported roots/dipole and 48 lowest-root probes passing.
Numerical tightening changes its position by only −0.544 micro-eV and phases
by 1.30e-6 rad. Matched CAS(10,11)/(10,10) position/phase changes of
64.8401 meV / 0.1098099 rad exceed the chosen gates. Saved-data fit-window
controls pass at fixed background, but background
terms 1–4 change width by 12.65%; production fitting remains unqualified.

Use persistent EBS storage mounted under `/home`, with approximately 500 GiB
for scratch, retained attempts and transferred evidence, plus space for Docker
and the source build on the root volume. Follow the shared
[Docker deployment guide](../../docker/ukrmol-plus/README.md) for bind-mount
ownership, engine checks and source/image provenance. Use On-Demand for these
initial long diagonalizations: the current pipeline has no demonstrated
restart of an interrupted dense eigensolver.

## Prepared queue and dependencies

| Recipe or work | Count | Purpose / prerequisite |
|---|---:|---|
| `calibration-sa12-qc-tight.json` | 2 QC jobs | Both equilibrium DZ starts and final subspace/root/dipole comparison pass |
| `calibration-large-host-qc.json` | 18 QC jobs | Matched CAS(10,11)/(10,12): eight new geometry/basis combinations each, plus an equilibrium 50-component ensemble for each |
| `calibration-tz-qc-fresh-starts.json` | 4 QC jobs | Equilibrium cc-TZ/aug-TZ RHF starts in both active spaces; queued behind the CAS(10,11) ladder on Sadaharu |
| Fixed-orbital coverage/start checks | As needed | Check all eight sectors; compare first ensemble roots, spins, ensemble objectives, dipoles and final active subspaces before promoting a checkpoint |
| `calibration-large-host-targets.json` | 10 dense jobs | Nine CAS(10,12) 40-component targets: three geometries x three bases; one DZ equilibrium 50-component orbital-ensemble target |
| Equilibrium CAS(10,12) scattering pilot | 1 initial pilot | Full target/import pass, valid scattering dimensions/workspaces; measure phases, candidates and cost against CAS(10,11)/(10,10) |

The geometry sentinels are **1.9, 2.1323 and 2.5 bohr**. Bases are
**cc-pVDZ, cc-pVTZ and aug-cc-pVTZ**: compare cardinal refinement first, then
diffuse augmentation at the same TZ cardinality. All 40-component targets use
five equally weighted roots in each singlet/triplet C2v sector. The 50-component
ensemble uses `[10,5,5,5]` roots in each spin and is a separate orbital-objective
test, not merely a channel-retention test.

The first six dense entries prioritize DZ geometry, equilibrium basis and
equilibrium ensemble dependence. The final four complete the TZ/aug-TZ
compressed/stretched combinations. Run only recipes whose exact initial
checkpoints have passed QC and coverage checks. If a competing start finds a
different preferred ensemble solution, revise the corresponding checkpoint
input rather than treating the currently named restart as automatically best.
Add fresh RHF controls at equilibrium in both TZ bases before declaring basis
start stability. Matched CAS(10,11) full imports fit Sadaharu and can be prepared
there, leaving the cloud host's RAM for CAS(10,12).

The existing same-active-space projection restriction applies: never feed a
CAS(10,11) checkpoint into a CAS(10,12) recipe. Retain the tightened restart and
all checkpoint hashes in the transferred evidence root so `/work/runs/...`
paths resolve without modifying historical directories.

The QC continuation can use three four-core slots on Sadaharu after the seed
checks pass and competing workloads are reassessed. Give each worker 32 GiB.
The neutral basis/correlation refinements and existing K-matrix fit-window/
background replays also belong in this pre-provisioning preparation.
`calibration-neutral-5z.json` prepares the same three frozen-core CCSD(T)
sentinels in aug-cc-pV5Z to complete the TZ/QZ/5Z basis sequence locally.
Its worker is queued behind the exact CAS(10,12) RHF-start container on CPUs
12–15 and has completed. The restart's fixed-orbital coverage scan also
passes after a preserved pre-probe import-path failure. CPUs 0–3 now run
the staged CAS(10,12) QC ladder. CPUs 12–15 have completed the extra-root
SLEPc control, all-root CI scan and saved-data fit replays. A short matched tight
CAS(10,10) scattering control and its coverage/independent-phase diagnostic now
complete there. Outer energy-grid and stretched-DZ repair coverage checks use
short gated leases; the finite worker then awaits the whole fresh-TZ batch
before coverage/repair/continuum jobs. CPUs 8–11 have started the gated local
CAS(10,12) SLEPc import after v3's controls pass, and that import now completes.
The gated stretched CAS(10,11) import/basis worker advances on released 8–11.
Both stretched-DZ CAS(10,11)
80/160-vector repairs pass QC; both projected equilibrium TZ retries still fail
singlet B1/B2 CI convergence after orbital convergence. Larger trial space does
not by itself resolve basis/orbital continuity.

The dense recipe is directly compatible with `batch.py`. First make a fresh
pilot manifest containing its equilibrium DZ entry. After its first singlet
and triplet memory/workspace checks, launch a second worker from the remaining
passing QC checkpoints while the pilot continues. Include this staggered
startup in the setup allowance below. Scattering still requires the pilot's
final all-root/import/dipole pass. A full two-worker batch on a checked host
has the form:

```bash
python3 -m projects.ukrmol_co.batch projects/ukrmol_co/calibration-large-host-targets.json \
  --root "$RUN_ROOT" --image "$IMAGE" --cpu-groups 0-31 32-63 --memory 224g
```

Resolve the two CPU groups against the observed topology. The launcher is a
finite queue, not a dependency scheduler: it does not wait for QC gates, stop
remaining jobs when a model fails, or adapt memory/ranks automatically. Use
filtered, newly named manifests for staged execution and all retries.

## Retained full-campaign time estimate

The measured CAS(10,12) **singlet-A1 stage alone** took **5751.68 seconds
(1.60 hours)** on four Sadaharu cores at dimension 43194. Its subsequent triplet
allocation failed, so 98.80 minutes is not a completed target runtime. There
are four singlet sectors near dimension 42000 and four triplet sectors near
71000. A cubic dense-eigensolver estimate gives

`T_diag(4) = sum[5751.68 * (N_sector / 43194)^3] = 34.20 hours`.

Ideal rank scaling would reduce this eigensolver-only estimate to 8.55 hours
at 16 ranks or 4.28 hours at 32 ranks. These are extrapolations, not measured
EC2 timings. Hamiltonian construction, communication, imperfect scaling,
shared memory bandwidth, QC restarts and DENPROP add time. The previous
CAS(10,11) import spent 37.81 minutes in DENPROP; its larger-space growth is
unmeasured. Use **12–30 hours per complete CAS(10,12) target** as the initial
allowance on a 16–32-core worker. Recompute from the first EC2 singlet/triplet
stages and final density stage, rather than applying a vendor speedup claim.

| Retained full-campaign work | Elapsed instance-hour allowance |
|---|---:|
| Ten dense target/import jobs, two workers | 60–150 h (five waves x 12–30 h) |
| One equilibrium scattering pilot, including its repeated target stages | 18–48 h; broad provisional allowance pending actual dimensions |
| Build/engine gate, staggered startup, sector scaling, collection and transfer | 8–12 h |
| Total arithmetic range | 86–210 h |
| Rounded planning range | **96–216 h (4–9 days)** |
| Earlier full-campaign budget, now deferred | **168 h (7 days)** |

This historical budget covers the target qualification queue and one measured scattering
pilot. Additional sentinel scattering, continuum/extraction refinements and the
eventual 25–40-geometry production campaign need a renewed forecast from that
pilot. Their decks depend on model selection and scattering workspace checks.
Keep the instance busy with prepared, independently useful jobs; stop it when
that queue is exhausted while retaining evidence on EBS.

As checked on 6 October 2026, the official Frankfurt Linux On-Demand price is
**$6.16704/hour**. Compute-only budgets are **$592.04 for 96 h**, **$1036.06 for
168 h**, and **$1332.08 for 216 h**. EBS, transfer and taxes are additional.
[`large-host-pricing.json`](large-host-pricing.json) records the official feed,
retrieval time, decoded response digest and rate codes. Recheck the rate and
regional instance offering when provisioning.

## Resource evidence and readiness

[`large-host-resources.json`](large-host-resources.json) records the audited
dimensions, maximum queries and forecast assumptions. The thirteen new valid
installed-library queries, source/licence copies, binary digest and execution
commands are retained in `diagnostics/cas12-large-host-workspaces/` under the
state-averaged evidence root. The combinatorial singlet/triplet-A1 counts match
the two independent existing CONGEN outputs. The thirteen queries and tight
CAS(10,11) import are included in the eighteen-attempt qualification supplement;
The current 46-attempt supplement additionally preserves the five/eight-root
CAS(10,11) SLEPc imports, independent all-64-root CI scan, tight CAS(10,11)
scattering, 49 native replay attempts including one failure, small 80/160-vector
controls, both tight CAS(10,12) starts and restart coverage, the nine-job
CAS(10,11) ladder, 5Z neutral sentinels and four rejected unstable stretched
references, the matched tight CAS(10,10) control and coverage, four CI-space
repairs, independent/held-point phase fits, finer outer grids and original
driver failures. Its later CAS(10,12) companion establishes forty-root/dipole
import and measured target memory/stages locally. Staged QC and continuum
controls remain active; all-spectrum scattering still needs its own audit.

Before proposing the cost-constrained paid experiment, finish the tight CAS(10,12) starts,
matched QC ladder and root/subspace checks; finalize the selected checkpoints
and their hashes; finish matching CAS(10,11) controls; verify the transferred
bundle and staged manifests. Use the established selected-root target feasibility
on Sadaharu and narrow the external-host request to the remaining bottleneck.
Rebuild and test on the chosen external host when it becomes available, then
measure its first heavy sector and update the hours.
All failed attempts remain evidence. Passing this queue establishes import
and model-comparison evidence, not automatic production qualification.
