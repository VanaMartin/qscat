# CO experiment: resources and continuation

## Live qualification handoff

Fresh CAS(10,11) TZ and both repaired aug-TZ targets now pass QC, independent
600-cycle fixed-orbital coverage and all-64-root UKRmol import/dipole checks.
The original rejected fresh aug-TZ record remains exact. Subsequent competing
starts find a distinct, reproducible aug-TZ branch whose objective is 0.538388 eV
lower; its coverage and repaired downward projections are the next local gates.
The matched CAS(10,10)
one-/two-/four-rank and concurrent benchmarks also pass all seven replicas.
CAS(10,12)'s tight QC starts and coverage pass, but full target-density/dipole
processing remains active; its compressed DZ CI failure still gates dependent
repairs. The neutral stretched RHF reference remains rejected. The next decisions
depend on these finite workers:

| Physical CPUs | Active work | Follow-on dependency |
|---|---|---|
| 0–3 | Staged CAS(10,12) QC ladder, 32-GiB containers | Passing DZ seeds release basis pairs; a new finite queue repairs eligible CI failures and checks both trial-space seeds before releasing missing basis pairs |
| 4–7 | Both new CAS(10,11) aug-TZ branch coverage/residual scans, 8-GiB container | Passing scans release four fresh-name core+active downward-projection retries, then resume the waiting staged CAS(10,11) owner |
| 8–11 | CAS(10,12) target DENPROP/dipole processing, 64-GiB import container | Final independent import qualification releases the recorded stretched import/basis and other dependent workers |
| 12–15 | Original CAS(10,11) l=3/l=5 continuum controls, 48-GiB containers | Their scientific gates pass; a finite successor advances them after the idle predecessor's recorded scheduling retirement, then replays saved-data fits |

The experiment roots remain `/home/kooza/ukrmol/co-sa-20261006` and
`/home/kooza/ukrmol/co-neutral-20261006`. The approximately $200 first paid
experiment ceiling remains deferred; all launches here use Sadaharu.
The new queues are detailed under [Local continuation — 7 October](#local-continuation--7-october-2026).

## Saved campaign — 6 October 2026

The calculation host is **Sadaharu**, reached with `ssh sadaharu` as `kooza`.
It has an AMD Ryzen 9 9950X3D (16 physical cores, 32 threads) and approximately
123.5 GiB RAM. In this campaign CPU IDs 0–15 each select a distinct physical
core. Recheck topology and competing workloads before a new batch.

Persistent files are under `/home/kooza/ukrmol/co-pilot`, mounted as `/work`:

- `runs/<name>/`: inputs, all UKRmol/Psi4 outputs, integrals, scratch, timing,
  resources and analysis; failed directories are retained.
- `batches/<manifest-stem>/`: manifests, commands, original exit codes, wall
  times and, for later batches, read-only source snapshots. The first four
  batches predate snapshotting; their hashes and generated run inputs remain.
- `upstream/`: UKRmol-scripts 1.0 archive and extracted sources.
- `qmodeling/`: the original deployed experiment source.
- `source-build-v5.log`: the original successful source build log.
- `pr-checkout/`: the packaged Docker build context; `packaged-build-coldcache.log`
  records its independent build/test verification.

The existing interactive `ukrmol-co-pilot` container was used for early setup;
finite fresh-container jobs are the reproducible execution interface. The
other observed container, `comfyui`, is unrelated to this experiment.

Historical images:

| Image | Recorded ID |
|---|---|
| Reference engine / original CO pilot | `sha256:e25463b7d4a3fa4ba268b2a621e9e9e895893d8826d5e22dc042c8f6f35b9bd8` |
| Original source engine | `sha256:599eeee919f56a27e6a0a4e214c5fc80a2f5780752f07b654b21fe1339f61ae6` |

These are local image IDs, not published registry digests. Rebuild elsewhere
using `docker/ukrmol-plus/Dockerfile` and record the new ID. Source checksums
and comparisons, not equality of rebuilt image bytes, are the evidence.

## Replaying the campaign

The committed JSON recipes name the jobs used for the quoted findings:

| Manifest | Question |
|---|---|
| `pilot-jobs.json` | Original three geometries and equilibrium SEP30 check |
| `calibration-freezing.json` | Two/four frozen A1 orbitals and DZ virtual-space growth |
| `calibration-basis-continuum.json` | Basis/virtual space, l cutoff, sphere, deletion and propagation |
| `calibration-continuum-retry.json` | Numerical-failure retries with deletion 1e-6 |
| `calibration-mpi.json` | One/two/four-rank timing at the same equilibrium model |
| `calibration-throughput.json` | Four independent one-rank jobs |
| `calibration-sentinels.json` | Compressed/stretched l=3/4 comparisons |
| `calibration-cc-hf.json` | CAS(10,8) HF-orbital CC40/CC80 diagnostics |
| `calibration-cc-casscf.json`, `calibration-cc-tz.json` | Ground-state-CASSCF CAS diagnostics and basis comparison |
| `calibration-cc-no-virtuals.json` | CAS(10,8) with no external virtuals |
| `calibration-cc-workspace.json` | Intermediate CAS(10,10) NDIMX-only workspace retry; failed |
| `calibration-cc-workspaces.json` | Successful CAS(10,10) runs with all workspaces raised |
| `calibration-source.json` | Source-engine differential check |
| `calibration-pairs.json` | 21 phase/fit comparisons; this is not a job manifest |

Copy the checkout to a new host directory, or use the packaged checkout on
Sadaharu. Give each repeat campaign a fresh root or fresh manifest/job names.
For example, on the host from the checkout root:

```bash
python3 -m projects.ukrmol_co.batch projects/ukrmol_co/calibration-cc-workspaces.json \
  --root /home/kooza/ukrmol/co-replay-001 \
  --image qmodeling/ukrmol-co:source --cpu-groups 0-3 4-7 8-11 --memory 24g
```

Prepare the new root with host-login and container-user write access. The
measured three-worker CAS(10,10) batch used 12 physical cores and finished in
15.1 minutes; the published cost forecast in the README is conditional on that
model size, not a budget for an accuracy-qualified state-averaged target.

The final runner sets `NDIMX=200000`, `CDIMX=NODIMX=20000`. Early manifests
refer to older snapshots with smaller or only partly corrected input limits.
Use the historical snapshot and templates to reproduce such failures; current
defaults intentionally execute the repaired configuration. Some original
nonzero batch codes were analyzer errors repaired after the engine succeeded;
`initial_batch_exit_code` retains this distinction in `calibration-results.json`.

Early batch source snapshots are absent for `calibration-freezing`,
`calibration-basis-continuum`, `calibration-continuum-retry` and
`calibration-cc-hf`. Their exact molecular decks, instrumented libraries,
templates, generated program inputs and outputs are included in the evidence;
do not claim that every historical Python harness revision is archived.

For a new aggregate, give `collect.py` fresh campaign metadata. For a historical
snapshot, use `campaign-provenance.json` as shown in the README. Native fits
remain in Hartree; the paired comparisons label their differences in eV.

## Current state-averaged campaign

### State-averaged continuation — 6 October 2026

Work is on `feat-co-state-averaging`. The new persistent root is
`/home/kooza/ukrmol/co-sa-20261006`, with its deployed sources in `checkout/`.
The image tag is `qmodeling/ukrmol-co:state-averaged`; the initial image ID is
`sha256:e662c6e21ebc082b0e1c0d3a3e7ecdb8c60bbe4bec9eaf8b12a2a2972ddb8e56`.
Each batch binds its own source snapshot, so the recorded source hashes also
matter. The PySCF/h5py additions follow the cached, previously tested engine.
The refreshed embedded-source image is
`sha256:24aeadcd0b27242064d248ffb76262e3d8c1598a1ad4c419a6a4ea7e6ed053e0`;
its build log is `build-final.log`. The package/module import smoke check passed.
The latest standalone source layer is tagged
`qmodeling/ukrmol-co:state-averaged-ci`, image ID
`sha256:b251e85f9fe7c1909bfd2e11f516756793874111e09b9cdf16f02221b04b3b19`.
`build-ci-coverage.log` records its cached engine gate and refreshed
embedded source. PySCF/h5py/target/archive/ci_probe imports passed; all 13 Python/Perl
source digests matched the 55-attempt publication checkout. Recorded jobs name their original
image IDs and immutable source snapshots.

### Completed snapshot

`sa-results.json` contains **55 completed attempts in 24 batches**:

| Outcome | Count |
|---|---:|
| Validated target-only pipelines | 17 |
| Validated scattering pipelines | 12 |
| Engine/runner failures | 23 |
| Failed target-import analysis | 2 |
| Pre-run setup rejection | 1 |
| Paired comparisons | 12 |
| Native RESON replays | 96 |
| Fixed-orbital CI coverage probes | 6 |

“Validated” means solver/import and pipeline consistency. **No electronic
model is yet qualified for production potential fitting.** Numerical findings,
historical failures and portable reproduction commands are consolidated in
[`README.md`](README.md); each completed batch retains its exact manifest,
source snapshot, image ID, commands and original exit codes. Early target-only
import probes used C at the origin; later jobs match the upstream center of
mass. Projected starts retain copied binary checkpoints and their hashes.

### Findings that determine the next calculation

- **Execution and target import work.** The engine passed 12 upstream checks;
  common-orbital targets have all-root energy and ground-dipole import checks.
  Tight CAS(10,10) sentinel targets pass at R=1.9/2.5 bohr, with gradients
  2.99e-8/8.50e-8 and root-import errors below 6.59e-10 Hartree. Compatible
  CI/augmented-Hessian cutoffs and eigensolver accuracy resolved their stalls.
  Recorded Newton controls still fail their gradient gates; use `one-step`.
- **Tested numerical controls meet the chosen gates.** Equilibrium DZ
  l=3/4/5 changes position by at most 3.18 meV, width by 0.585%, and phase by
  0.02615 rad. The stretched 0.020/0.010/0.005-eV grid sequence has finest-pair
  differences of 1.50e-6 eV in position and 0.0041% in width. Angular width/phase
  changes are nonmonotone; these are model-conditional stability results.
- **Fit qualification varies with geometry.** The compressed candidate is
  3.5191/2.1026 eV; background-order widths span 2.0529–2.4573 eV and fail
  the 5% criterion. The stretched candidate is 0.9738707/0.2817833 eV, with
  maximum background width change 2.61%. Neither establishes feature identity
  across geometry or certifies a pole/bound-state interpretation.
- **Electronic-model dependence remains substantial.** Projected aug-DZ/TZ
  widths differ by 6.36%. Changing the projected TZ orbital ensemble from 40
  to 50 components raises the ground root by 0.1316 eV and changes its dipole
  by 32.6%. The separate DZ 40-to-50 retained-channel check, with its orbital
  ensemble fixed, passes the pairwise gates. Do not combine these refinements.
- **CAS(10,11) exposed a missing lowest CI root.** The original engine run
  completed in 102.22 minutes with an 18.39-GiB kernel peak, but its purported
  fifth triplet-A1 QC root matches the fresh sixth root: import error 0.129757
  eV. Six fixed-orbital probes recover UKRmol's lowest five roots within
  2.75e-10 Hartree; their common roots agree within 4.27e-14 Hartree.
  `diagnostics/cas11-ci-coverage/` preserves the checkpoint, exact source and
  execution record. The fresh-CI reoptimization now passes all 40 root imports
  within 5.20e-10 Hartree and the dipole within 8.10e-11 a.u. in 127.80 minutes,
  at an 18.41-GiB kernel peak. Fresh-CI energies agree within 1.14e-13 Hartree.
  Its final gradient is 4.77e-6 at the requested 1e-5 tolerance. Tightened
  optimizer/start checks remain open before selecting the electronic model.
  Relative to CAS(10,10), the ground root drops 0.80403 eV and the dipole
  falls 69.88%, while the lowest triplet-Pi excitation changes only −10.10 meV.
- **Fresh-CI controls are established.** `calibration-sa-fresh-ci.json`
  preserves a passing ground control and a failed mixed-solver override.
  Installing overrides after mixer construction fixes that interface. The
  balanced CAS(10,8) 40-component control passes in 2.45 minutes, with fresh-CI
  root difference 8.53e-14 Hartree and import error below 5.44e-10 Hartree.
  A two-spin A1-only control failed computed Pi degeneracy and remains rejected.
- **CAS(10,12) needs an explicit engine memory budget.** Its original QC
  target passed, but its 43194-dimensional singlet-A1 matrix exhausted the
  2.5-GiB SCATCI per-process limit. A larger container alone does not change
  that limit. The 6-GiB-per-process retry completed singlet-A1 diagonalization,
  then failed the 70674-dimensional triplet-A1 matrix allocation: 9.31 GiB
  per rank. Total wall time was 98.80 minutes and kernel memory peaked at
  56.10 GiB in its 80-GiB container. Preserve this as a second engine-memory
  failure; full root-import/dipole checks remain incomplete.

### Latest larger-active-space jobs

Both retries have final batch records. CAS(10,11) passes target/import checks;
CAS(10,12) completed with engine exit code 25. The fresh-CI batch exits nonzero
because its preceding two-spin A1-only control fails computed Pi degeneracy,
despite engine success and agreement of its two averaged A1 roots.

| Manifest / run | CPUs | Container limit | Latest observed progress |
|---|---|---:|---|
| `calibration-sa-fresh-ci-retry.json` / `co-eq-ccdz-sa11-40-target-fresh-ci-mixfix` | 0–3 | 32 GiB | Completed success: QC, fresh CI, every averaged-root import and dipole passed |
| `calibration-sa-active-memory.json` / `co-eq-ccdz-sa12-40-target-memp6` | 8–11 | 80 GiB | Completed failure: singlet-A1 passed; triplet-A1 exceeded the 6-GiB per-process budget |

QC values below are in the remote runs' `target.json`. CAS(10,11) is now
import-qualified; CAS(10,12) remains incomplete. Both completed batches,
including the failed A1-only control, are included in the archive.
Both runs use their requested default 1e-5 orbital-gradient tolerance.

| QC quantity | CAS(10,11), fresh-CI restart | CAS(10,12), memory retry |
|---|---:|---:|
| Ground-root energy (Hartree) | −112.924284941 | −112.927240475 |
| Ground dipole z (a.u.) | 0.02773156 | 0.01751140 |
| Final macroiteration gradient | 4.77e-6 | 9.99e-7 |
| Maximum fresh-CI energy difference (Hartree) | 1.14e-13 | 8.88e-11 |
| Maximum Pi splitting (Hartree) | 5.29e-11 | 6.70e-9 |

Those batches released their CPU allocations. The following qualification
continuation uses new run names and immutable source snapshots.

### Tight-target and neutral qualification continuation

`calibration-sa11-qc-tight.json` runs CAS(10,11) at energy/gradient tolerances
1e-11/1e-7, with a fixed 40-component ensemble and fresh CI. The successful
checkpoint, original rejected-checkpoint and RHF starts pass QC reanalysis in
6.08/13.13/61.89 minutes, at gradients 3.22e-8/6.22e-8/8.30e-8. Across the
three starts, maximum root/dipole differences are 7.69e-9 Hartree / 5.13e-8
a.u.; the final active-subspace minimum overlap singular value is
0.999999999999924. Their common ensemble objective is −112.370969232576
Hartree. This supports start stability for the fixed ensemble/active space,
not global orbital optimality or electronic-model convergence. The completed
QC batch released CPUs 0–7.

The tightened restart checkpoint is
`runs/co-eq-ccdz-sa11-qc-tight-restart/co.casscf.chk`, SHA256
`64131e79fab6e1a55c1489646555f090cb5186530a66babfe39fd46bbbfcb402`.
`calibration-sa11-tight-import.json` starts
`co-eq-ccdz-sa11-tight-import` from this checkpoint. The full independent
UKRmol target check completed on CPUs 8–11 in a 32-GiB container in 97.75
minutes (97.76 minutes including launch/analysis). Independent raw reanalysis
passes all 40 averaged roots within 5.33e-10 Hartree and the ground dipole
within 7.95e-11 a.u.; its kernel memory peak is 18.42 GiB. The final gradient
is 2.92e-8, and fresh-CI energies agree within 1.14e-13 Hartree. The imported
checkpoint is `runs/co-eq-ccdz-sa11-tight-import/output/CO/geom1/co.casscf.chk`,
SHA256 `3480f6127549d963596e67c3592550f533f7a6e30db0d6343d336a7d9e794554`.

`calibration-sa11-tight-scattering.json` now runs
`co-eq-ccdz-sa11-tight-cc40` from that imported checkpoint, on CPUs 8–11 in a
48-GiB container, with 6-GiB internal SCATCI matrix budgets. It retains the
tight QC controls, 40 channels, 18-bohr/l=4 continuum, deletion 1e-6 and
99 native energies over 0.1–5.0 eV. Launch output is
`launch-sa11-tight-scattering.log`; the attached supervisor is
`sh_112ad1252001dZHiXMnz7jPTuw`. Final phases, dimensions and runtime are pending.

Its first doublet-B1 block has now completed: raw CONGEN reports **344124
uncontracted CSFs**, while SCATCI solves a **27546-dimensional contracted
Hamiltonian**. The block takes 4387.34 seconds (73.12 minutes); its printed
ScaLAPACK workspace is `lwork=380327961, liwork=192840`. The sampled run
memory peak to this observation is 24.18 GiB, a lower bound on the eventual
whole-run peak. `diagnostics/cas11-scattering-first-block/` retains the completed
B1 raw logs, hashes and observation contract. B2/whole-run analysis remains
active; this partial observation is excluded from the public completed snapshot.
Use the contracted dimension for its dense-array audit rather than the raw
CONGEN count. CAS(10,12) scattering still needs its own measured dimension.

`diagnostics/cas11-tight-ci-coverage/` completed 48 fixed-orbital probes in
116.91 seconds on CPUs 4–7 and 12–15: five/eight roots at trial-space sizes
40/80/160 in all eight sectors. Forty-six probes converge completely; the
eighth singlet B1/B2 roots fail at space 40 and converge at 80/160. Every
ensemble root converges in every probe. Converged first-five energies agree
with the target within 1.28e-13 Hartree, with maximum common-root spread
1.71e-13 Hartree. Raw CASCI logs and spins were independently checked. The
original nonzero campaign exit remains preserved. Final checkpoint subspace
comparison of all three starts is in `diagnostics/cas11-tight-start-comparison/`.

The CAS(10,12) workspace audit is complete:
`diagnostics/cas12-scalapack-workspace-audit/` preserves raw installed-library
queries and the upstream source/licence; `diagnostics/cas12-scalapack-guard/`
verifies the shipped utility, including explicit rejection (exit 5) of the
undersized four-rank triplet query. Valid eight-/sixteen-rank queries require
148.92/148.99 GiB for the main arrays alone. This exceeds host RAM; do not
repeat the current dense calculation by increasing only `memp`. See
[the qualification method](../../docs/physics/co-electronic-qualification.md)
for the included/excluded arrays and overflow signature.

The independent neutral pilot is complete under
`/home/kooza/ukrmol/co-neutral-20261006`. Three batches retain eight attempts:
seven passing RHF/frozen-core CCSD(T) records and the initial failed Psi4
Python import. The executable-based cc-pVDZ PySCF/Psi4 control agrees within
8.3e-11 Hartree. Six aug-TZ/QZ calculations at R=1.9/2.1323/2.5 took 90.93
seconds on two four-core workers. `neutral-results.json` and
`neutral-provenance.json` retain the aggregate and method/reference convention.
The basis change shifts the sentinel relative energies by 76.04/68.47 meV;
the R=2.5 amplitude diagnostics increase, so further basis and correlation
checks remain necessary before selecting a full neutral curve.

The completed supplement is published through
[`qualification-evidence/`](qualification-evidence/README.md).
`qualification-results.json` contains 46 attempts in eighteen completed batches:
nine passing QC targets, ten passing neutral records and five neutral failures,
dense and five/eight-root SLEPc tight CAS(10,11) all-root/dipole imports, four rejected Davidson
controls, four passing small-model SLEPc controls, tight CAS(10,11) scattering
and seven rejected CAS(10,11) ladder entries, the matched tight CAS(10,10)
scattering control, and four 80/160-vector repairs (two passes/two CI failures).
Its public archive retains 160 CI probes and twenty-one diagnostic directories,
including the independent all-64-root check, 49 native fit-replay attempts,
24 independent/120 held-point fits, four finer outer-only grids with native
MAXFIT diagnostics, and preserved failed launches/retries.
Verification checked 4635 payload digests, 298 batch-source hashes, 17 image-source
hashes, all 28 successes, all four Davidson rejections, both aggregates, the
ladder/reference/repair/driver failures, matched scattering comparison, outer
grids and dense/SLEPc resource/storage comparison;
repackaging is byte-identical. Staged QC, continuum controls and
the full CAS(10,12) import remain pending. The earlier supplements remain
available through the manifest's prior-snapshot pointers.

The immutable current supplement is
`https://data.qscat.org/ukrmol-co-electronic-qualification-2026-10-06/qualification-evidence.tar.a4c0738b6d01.gz`,
SHA256 `a4c0738b6d0168c49df3600c109a4f8baab67b5dcb7f81f8130fea9cc52d586a`,
29,897,245 bytes. The public client fetched and verified those bytes, and the
extracted verifier reproduced all checks and byte-identical repackaging.
Publication bundle, pointers, aggregate/source provenance, verification report
and four new/updated verifier attachments are retained on Sadaharu in
`prepared/publication-46-20261007/`. The 41 earlier run records and sixteen batch
records remain exact; the 41-attempt archive remains an immutable prior snapshot.

The refreshed standalone image is `qmodeling/ukrmol-co:electronic-qualification`,
ID `sha256:b5d5f7fc4d638d19eaeb65399b5acf20f56f4fe4526c63761b8945a0b01c5d08`.
`build-electronic-qualification.log` and `electronic-qualification-image.json`
record the cached engine gate, runtime imports, pinned dependency versions and
all 17 embedded Python/Perl/Fortran source digests. Numerical jobs retain their
earlier image IDs and immutable source snapshots.

Repository handoff: analysis/packaging code is pinned at
`28f5e1f6158abe775c8d2a70244504dacc46201a`. The committed-main index is complete
and upstream-current at `0f9768e4e89f2accb6f3ff5c2d3dc921719325eb`; the branch
adds state-averaged/fresh-CI work plus the QC-only, correlated-neutral and
workspace qualification tools, selected-root configuration/solver checks,
optional SLEPc source build, CI trial-space control, independent phase fits and
matched scattering/outer-grid controls with their evidence. Source anchors and local
changes were resolved against that baseline before editing and publication.

### Published evidence and collection

Completed evidence is published through [`sa-evidence/`](sa-evidence/README.md).
The 55-attempt snapshot has 29 validated pipelines, 23 engine failures, two
failed target-import analyses and one pre-run setup rejection. The public
client verified all 6,616 payload digests and exactly reconstructed twelve
comparisons. All 29 successful runs passed raw-output reanalysis, six fixed-CI
probes match the independent UKRmol roots, and repackaging is byte-identical.
The prior 37/41/52/53-attempt archives were publicly fetched with matching
digests. Every batch in that historical snapshot is complete and included. The
newer qualification continuation above has separate completed evidence and live jobs.
The immutable current archive is
`https://data.qscat.org/ukrmol-co-state-averaged-2026-10-06/state-averaged-evidence.tar.99d79a4bfdf2.gz`,
SHA256 `99d79a4bfdf2433f4c90e3f549b8f7eb84ba3af16c7b1e6f70b87759f3fd7a54`.

Collect fresh completed artifacts outside the checkout, then use:

```bash
uv run python -m projects.ukrmol_co.collect "$COPIED_ROOT" \
  --provenance projects/ukrmol_co/sa-provenance.json \
  --pairs projects/ukrmol_co/calibration-sa-pairs.json \
  --output projects/ukrmol_co/sa-results.json
```

The collector ignores batches without a final `result.json`, retains original
batch exit codes, and includes `target.json` even when a later engine stage
fails. Check native ground dipoles with the updated analyzer; its DENPROP
comparison was added after the initial successful import batch.
The collector also retains pre-run setup failures without fabricating resource
measurements. Reanalyze successful copied runs before collecting a current
snapshot: a fresh copy from the host may contain older analyzer summaries.
Reanalyze after the final sync and before collection/packaging; another broad
sync can overwrite locally refreshed summaries with older host records.

### Remaining gates

The larger-host preparation is recorded in [`LARGE_HOST.md`](LARGE_HOST.md).
The retained hardware candidate is Frankfurt `r8a.16xlarge` (AMD x86_64,
64 physical cores, 512 GiB). Its earlier **168-hour** full-campaign estimate
and provisional **96–216-hour** range cover ten imports plus one scattering pilot. The
official 6 October Linux On-Demand price is $6.16704/hour, recorded in
`large-host-pricing.json`; this is approximately $1036 for seven days of compute.
**Paid provisioning is deferred:** the first external-host experiment must be
reduced by at least 80%, to approximately $200, while qualification continues
on Sadaharu. Keep the full queue and estimate as documentation. Investigate
selected-root UKRmol target diagonalization against the existing dense
controls; only residual memory/time bottlenecks should motivate a paid job.
The same audited x86-64 execution bundle can be used on other supplied hardware.

The nine CAS(10,11) entries from `calibration-large-host-qc.json` completed
on CPUs 4–7 in one 32-GiB worker in 133.62 minutes, with two QC passes and
seven preserved failures. Their derived manifest and selection/hash
record are under `prepared/calibration-large-host-cas11-qc{,.provenance}.json`;
batch records are `batches/calibration-large-host-cas11-qc/`. Launch output is
`launch-large-host-cas11-qc.log`; supervisor `sh_112adc835001cNIQlav3d6Gp81`.
The complementary CAS(10,12) ladder now uses the completed start/coverage gates
and staged projection procedure below.
`calibration-tz-qc-fresh-starts.json` adds four equilibrium RHF-start controls
for CAS(10,11)/(10,12) in cc-pVTZ and aug-cc-pVTZ. They are queued behind the
whole nine-job CAS(10,11) ladder on CPUs 4–7, with 32-GiB limits and a frozen
source at `prepared/tz-fresh-starts-source/`. Supervisor
`sh_112e5e5b6001xMjUai178wgBYx` uses a Linux process-exit notification for
ladder PID 898757, rather than waiting for only its current container.
The next launch log is `launch-tz-qc-fresh-starts.log`. These starts preserve
the same 40-component objective and need root/dipole/subspace comparison
against the projected starts before the basis sequence can be qualified.

`calibration-sa12-qc-tight.json` runs the equilibrium restart/RHF-start
QC checks on CPUs 0–3 and 12–15, two 32-GiB containers, using the
`electronic-qualification` image. Launch output is `launch-sa12-qc-tight.log`;
the attached supervisor is `sh_112a1235b001TlaKHQwphj8Oi5`. These checks precede
the prepared eighteen-job matched QC ladder and ten-job large-host target queue.
The restart completed in 21.17 minutes, with gradient 5.64e-8, ensemble
energy −112.3847060392213 Hartree, ground energy −112.92724051225153 Hartree
and dipole z 0.0175118735858473 a.u. The RHF start passes in 152.18 minutes
at gradient 4.20e-8. Their objectives differ by 2.84e-14 Hartree, roots by at
most 4.17e-8 Hartree and dipoles by 6.27e-8 a.u. The final active-subspace
minimum overlap singular value is 0.999999999999496; the comparison is in
`diagnostics/cas12-tight-start-comparison/`. The restart checkpoint SHA256 is
`9451c7baed762ab674e9d36111cc220202972703ac52346e1ce13dcdc5ac55a0`;
the RHF-start digest is
`05acb847b242c5a93da0fec9015a281bcff056915fe4355f255453820b91677a`.
The new `diagnostics/cas12-large-host-workspaces/` audit passed thirteen queries
covering all target-sector dimensions on 16-/32-rank grids, including both
32-rank orientations for the largest sector. Maximum array floors are
149.78/149.91 GiB, before engine overhead. `large-host-resources.json` retains
their summary and the stage-scaled forecasting assumptions. The extended
diagnostic is included in the eighteen-attempt supplement; live QC attempts are excluded.

### Local selected-root and neutral continuation

`calibration-target-davidson-small.json` and
`calibration-target-davidson-coverage.json` completed four rejected CAS(10,8)
targets on CPUs 0–3 with 8-GiB limits. All eight sectors report native Davidson
success, yet 5/8/16/32 requested roots give maximum required-root errors of
0.13931/0.07915/0.05290/0.02823 Hartree. The 32-root control also tightens
`crite` to 1e-13. Keep the original nonzero analysis/batch exits and raw logs.
`diagnostics/target-davidson-controls/` retains the small stored singlet-A1
Hamiltonian and its independent NumPy/SciPy reconstruction: lowest five roots
agree with QC within 3.69e-13 Hartree, residuals below 1.61e-13 Hartree.
Larger serial-Davidson jobs are gated out.

The SLEPc-enabled source build is `qmodeling/ukrmol-co:selected-roots`, image
ID `sha256:6b3b0ededa85494076229b111efe969985407023e404bff888f17cafb4630666`.
`build-selected-roots.log` records a fresh pinned engine build and 12/12
upstream checks. Runtime linking/import checks pass. Source provenance records
PETSc/SLEPc library digests from the same pinned toolchain; this image is a
separate configuration from the dense source images.
`calibration-target-slepc-small.json` completed on CPUs 0–3, one 8-GiB slot,
with five roots, then eight roots and a tighter tolerance. Its launch log is
`launch-target-slepc-small.log`. Both controls pass all forty averaged roots
within 4.98e-10 Hartree and the dipoles within 5.62e-11 a.u.; eight raw sectors
identify Krylov–Schur in each run. Walls are 49.53/46.91 seconds, peaks
0.1122/0.1120 GiB. Common first-five roots agree to the printed precision.
The raw-log analyzer rejects solver fallback/incomplete sectors before applying
the usual all-averaged-root/dipole gates. `calibration-target-slepc-cas11.json`
completed on CPUs 0–3 with a 16-GiB limit, supervisor
`sh_112d036710018xEwempAilQec9`, log `launch-target-slepc-cas11.log`.
Its full independent import passes all 40 roots within 5.34e-10 Hartree and
the ground dipole within 8.04e-11 a.u. Complete wall/peak are 3152.81 seconds /
5.25 GiB, versus the dense reference's 5865.15 seconds / 18.42 GiB. SCATCI
totals are 351.70 versus 3248.81 seconds (89.17% reduction), but DENPROP takes
2365.46 seconds, limiting the observed full-wall reduction to 46.25%.
Concurrent workloads and CPU placement differ; this is not MPI scaling evidence.
`diagnostics/cas11-selected-root-comparison/` independently reconstructs the
comparison, verifies raw logs and retains the pinned target source showing
dense PETSc Hamiltonian storage.

`calibration-neutral-5z.json` completed aug-cc-pV5Z CCSD(T) at the existing
three sentinel geometries for a TZ/QZ/5Z basis sequence. Its finite worker
used CPUs 12–15 with a 24-GiB limit in the neutral evidence root,
using the selected-root image by ID and a frozen source copy at
`prepared/neutral-5z-source/`. Supervisor `sh_112e0abf3001S8k4V3u0CJgThH`
blocks on the exact RHF-start container's exit, then runs all three jobs;
it does not depend on that orbital start passing. Its log will be
`/home/kooza/ukrmol/co-neutral-20261006/launch-neutral-5z.log`.
All three records pass in 18.97 minutes. Their own-reference relative energies
at R=1.9/2.5 are 1.26565/1.41165 eV; QZ-to-5Z shifts are −17.24/+15.82 meV.
Memory peaks reach the 24-GiB cgroup limit including page cache, with about
17.14-GiB anonymous peaks and 12.81–14.24-GiB run-disk peaks. The subsequent
`calibration-neutral-stretched.json` batch tests R=3.0/4.0 in aug-TZ/QZ on the
same CPUs, with 16-GiB limits. All four references converge RHF but fail external
stability in 15.74 seconds total; retain the original exits and diagnostics.
Supervisor `sh_11364c06f0016E0feSb5dQSaXI` and `launch-neutral-stretched.log`
identify that completed rejection batch. The restricted neutral route is
blocked there pending a different validated correlation treatment.

The CAS(10,12) restart's original coverage launch used CPUs
0–3 with an 8-GiB limit, supervisor `sh_112e25f9e001jXYdP8JvH67gD6`.
`prepared/cas12-coverage-followup.py` waits for the exact CAS(10,11) SLEPc
container to exit, independently rechecks its all-root/dipole import and dense
reference agreement, then scans five/eight roots and trial spaces 40/80/160
in all eight sectors. It retains individual exits and gates every ensemble
root; failed extra-root controls remain evidence. Resources and source hashes
are preserved under `diagnostics/cas12-tight-ci-coverage/`; launch output is
`launch-cas12-tight-ci-coverage.log`. This does not select the restart over the
RHF start or qualify an independent CAS(10,12) import. Its child process failed
before any probes because the diagnostic working directory could not import
`projects`; the original exit/resources/source are preserved. The repaired
launch explicitly binds the frozen package and sets `PYTHONPATH=/opt/qmodeling`
for all children. Supervisor `sh_1130398bb001eNc4f7Em2WfdMf`, script
`prepared/cas12-coverage-pathfix.py`, log
`launch-cas12-tight-ci-coverage-pathfix.log` and
`diagnostics/cas12-tight-ci-coverage-pathfix/` retain the fresh retry.
It completes all 48 probes in 705.79 seconds: all five ensemble roots pass
within 2.85e-13 Hartree, 46 probes converge fully, and the extra eighth singlet
A1 and seventh/eighth singlet A2 roots fail at space 40 and converge at 80/160.
Independent raw CASCI energies and every recorded spin have been verified.

### CI convergence refinement and gated local CAS(10,12) trial

The first three CAS(10,11) ladder observations, recorded before the whole
batch completed, are saved at `prepared/cas11-qc-initial-observations.json`:

| Job | Outcome | Wall (min) | Kernel memory peak (GiB) |
|---|---|---:|---:|
| Compressed DZ, R=1.9 | QC passes: gradient 5.67e-8, fresh-CI agreement 1.14e-13 Hartree | 12.52 | 0.301 |
| Stretched DZ, R=2.5 | Orbital convergence passes; triplet A1 CI fails | 16.47 | 0.303 |
| Equilibrium TZ | Orbital convergence passes; singlet B1/B2 CI fails | 23.31 | 0.327 |

The compressed ground/ensemble energies are −112.86483724227091 /
−112.21029375789468 Hartree, dipole z 0.18617598314225828 a.u., Pi splitting
6.96e-13 Hartree and MO orthogonality error 2.55e-15. It remains a QC-only
record pending fixed-orbital coverage and independent import. The two failures
are preserved with original exits and checkpoints; neither is a passing target.
The completed batch additionally passes the equilibrium 50-component DZ
ensemble (52.74 minutes), fails equilibrium aug-TZ singlet B1/B2 CI (28.49
minutes), and rejects all four combined geometry/basis projections before
optimization. The pinned PySCF helper does not support changing both at once.
The 50-component ground energy/dipole are −112.90750858551057 Hartree /
0.0624155060260635 a.u.; compared with the 40-component ensemble, the ground
root rises 0.45651 eV and the dipole changes +0.0346869 a.u. (+125.09%). This
is an orbital-ensemble dependence, not a channel-retention result.

`--target-ci-max-space` now exposes the PySCF trial-vector limit independently
of ensemble roots/weights. The default is still `max(40,8*roots)`; limits 80
and 160 are experimental refinements, with recorded actual sector values.
`calibration-ci-trial-space-control.json` first tests both limits against the
passing small-model dense spectrum/dipole. Only passing controls release
`calibration-cas11-ci-space-retries.json`, four fresh reoptimizations from the
two rejected final checkpoints. Existing convergence/spin/Pi/MO gates remain
required. The live frozen source is `prepared/ci-space-and-cas12-source-v3/`.

Supervisor `sh_11366e4ab001aML66KbDnQqQ6e` queues this finite continuation
behind the exact tight-scattering container on CPUs 8–11, with 32-GiB limits
for controls/retries and 64 GiB for the gated CAS(10,12) import.
Its script is `prepared/ci-space-and-cas12-followup-v3.py`, log
`launch-ci-space-and-cas12-followup-v3.log`; per-step commands, original exits,
source hashes and stop reasons are recorded in the frozen source's
`execution.json`. After the CI-space controls/retries, it independently
rechecks the five-root CAS(10,11) SLEPc/dense comparison and waits for
`calibration-target-slepc-cas11-coverage.json`: eight UKRmol roots per sector
with `crite=1e-13`, retaining the same 40-component orbital ensemble. This
independent control ran concurrently on the freed neutral CPUs 12–15,
with a 16-GiB limit, supervisor `sh_11366bf60001sviokgYA0SUSuq`. Its log is
`launch-cas11-extra-root-parallel.log`; the recorded batch PID/command are in
`prepared/cas11-extra-root-parallel-launch.json`. Pre-launch supervisor revisions
v1/v2 were superseded before running any experiment; their records and exits
remain under the original frozen source directories.
Sixteen fixed-orbital PySCF probes then check **all 64 requested roots**, using
trial spaces 80/160, and retain evidence under
`diagnostics/cas11-slepc-extra-root-coverage/`.
The completed early independent scan is recorded below; the original v3
supervisor retains its own scheduled verification before the larger import.

Only that all-root/dipole/extra-root gate, plus the completed 48-probe
CAS(10,12) restart coverage gate and unchanged checkpoint hash, can release
`calibration-target-slepc-cas12.json`. The supervisor also requires at least
80 GiB available RAM and 20 GiB free scratch before that launch. The target
uses 16-GiB internal budgets because its largest dense PETSc matrix alone
requires 37.41 GiB (9.35 GiB per rank on four ranks). This is an
independent **tight-restart import trial**; it does not select that start over
the RHF start or establish basis/model/scattering convergence. No paid host
is provisioned by this continuation.

### Staged projection ladders

Supervisor `sh_113359889001lD3ePn6nKY7Yfi` runs
`prepared/cas12-qc-staged-followup.py` on CPUs 0–3, with frozen source
`prepared/cas12-qc-staged-source/` and log `launch-cas12-qc-staged-followup.log`.
It checks the completed restart/RHF subspaces and 48-probe restart coverage,
then runs five base QC entries: compressed/stretched DZ, equilibrium TZ/aug-TZ
and the equilibrium 50-component DZ ensemble. Its base batch is
`calibration-large-host-cas12-qc-base`. Each passing sentinel DZ checkpoint
releases two same-geometry basis projections, with derived manifests
`calibration-large-host-cas12-qc-r1900-basis` and
`calibration-large-host-cas12-qc-r2500-basis`. A rejected seed stops only its
dependent geometry; all source/checkpoint/manifest hashes and exits are retained.

Supervisor `sh_113674883001zx4Z2Cxt7WwKT5` queues
`prepared/cas11-qc-staged-followup-v2.py` on CPUs 4–7, using
`prepared/cas11-qc-staged-source-v2/` and log
`launch-cas11-qc-staged-followup-v2.log`. It waits for the whole fresh-TZ batch
and live v3 continuation, requires both small CI-space controls to pass, then
runs two aug-TZ 80/160-vector retries and four staged basis starts. The compressed
DZ seed is the passing host-preparation checkpoint; stretched starts require a
passing 80/160-vector DZ retry. Derived manifests and exact checkpoint hashes
are recorded in the frozen source and each batch. All such records remain QC
qualification work until coverage/import/start/basis checks pass.

### Local continuation — 7 October 2026

The latest resource inspection finds approximately 100 GiB available RAM and
159 GiB free on `/home`. All four physical-core groups have active work again,
including the matched tight CAS(10,10) control on CPUs 12–15. SMT partners
remain outside the experiment's CPU allocation. The following **completed
individual attempts belong to unfinished batches**, so they are observations
outside the published 41-attempt supplement:

| Run | Outcome | Wall (min) | Kernel memory peak (GiB) |
|---|---|---:|---:|
| `co-eq-cctz-sa11-qc-tight-hf-start` | QC passes; final gradient 6.81e-8, fresh-CI root agreement 9.95e-14 Hartree | 39.75 | 0.323 |
| `co-eq-aug-tz-sa11-qc-tight-hf-start` | Orbital/optimization CI flags pass; fresh singlet-A1 CI audit fails | 22.80 | 0.428 |
| `co-r1900-ccdz-sa12-qc-hostprep` | Orbital flag passes; singlet-A1 CI fails at trial space 40 | 59.20 | 0.716 |

`prepared/live-qc-observations-20261007.json` retains their original exits,
resource records, final iterations and input/log/checkpoint hashes. The passing
fresh TZ target has ensemble/ground energies −112.38809524704945 /
−112.90185648989254 Hartree and dipole z 0.143618428763436 a.u. Its Pi splitting
is 2.22e-11 Hartree and MO metric error 8.13e-15, but the smallest initial/final
active-overlap singular value is 0.171594: scalar convergence does not establish
active-space continuity or the preferable orbital solution. Its checkpoint SHA256
is `c41b08888d05483c6845f894c94c568ce9fadb7adc831fdefc57382c32d2a73a`.
The rejected fresh aug-TZ audit has energy agreement at 1.42e-13 Hartree;
the unconverged CI flag still rejects it. A small energy difference cannot
replace the CI convergence gate.

Two additional finite supervisors have been launched, each using a fresh source
copy and the pinned selected-root image. They block on Linux PID-exit
notifications and complete predecessor records, retaining generated manifests,
checkpoint/source hashes, commands, original exits and explicit stop reasons:

| Queue | CPU group | Frozen source / execution record | Launch log / attached supervisor |
|---|---|---|---|
| Fresh-start coverage, repairs and continuum | 12–15 | `prepared/local-coverage-scattering-source/` | `launch-local-coverage-scattering.log`; `sh_1137e33f1001jM36cG30kdTCeG` |
| CAS(10,12) base-seed repairs and staged projections | 0–3 | `prepared/local-cas12-repairs-source/` | `launch-local-cas12-repairs.log`; `sh_1137e33f8001btj40rUlpPyXW8` |

Each source directory contains `sadaharu-followup.py`, `compare-checkpoints.py`,
`dependencies.json` and `execution.json`. Recorded host PIDs are 918330/918331;
the supervisor SHA256 is
`1185613f615a8221174cd8e386ff39fab0eb608e66e2e610ddc0d99eb4ca0cc2`.
These are finite campaign drivers retained with the experimental evidence.

The CPU-12–15 queue requires the whole extra-root SLEPc and fresh-TZ batches
to release their allocations. It reanalyzes every passing fresh TZ checkpoint
and performs 48 independent fixed-orbital probes: five/eight roots at spaces
40/80/160 in all eight sectors. A failed extra root remains evidence; every
ensemble root must converge with the correct spin and agree within 1e-7 Hartree
in every probe. Diagnostics are `diagnostics/lowest-roots-<run>/`.

After the live v3 supervisor finishes and both small trial-space controls pass,
eligible fresh-start CI failures receive **freshly named 80/160-vector
reoptimizations**, including a failed lowest-root coverage check on an otherwise
passing fresh start. Only an orbital-converged checkpoint with a CI-convergence
failure or measured coverage failure is retried. Other failures are recorded
for review. Passing repairs receive the same 48-probe scan, and paired passing
repairs have their roots, ensemble objective, dipoles and core/active subspaces
compared. The derived batch is `calibration-tz-fresh-ci-space-repairs`.

The CPU-0–3 queue waits for the whole existing CAS(10,12) staged supervisor and
the live v3 controls/import supervisor. It makes 80/160-vector reoptimizations
of eligible CI failures in `calibration-large-host-cas12-qc-base`, preserving
the original requested orbital ensemble. Its derived batch is
`calibration-cas12-base-ci-space-repairs`. A failed sentinel DZ seed releases
the previously blocked same-geometry TZ/aug-TZ pair only when **both repaired
DZ seeds pass QC and the 48-probe coverage scan**, then agree within
1e-7 Hartree in roots/objective, 1e-5 a.u. in dipole and 1e-8 in core/active
subspace singular-value distance from unity. The first passing trial-space seed
is then a provisional projection input, not a selected production orbital model.
Derived basis batches are `calibration-cas12-r1900-basis-repaired` and
`calibration-cas12-r2500-basis-repaired`.

The checkpoint-comparison helper was exercised against the already-qualified
CAS(10,11) tight restart/RHF pair: it reproduces the 7.69e-9-Hartree root,
5.13e-8-a.u. dipole and 0.999999999999924 active-subspace agreement. A
different-basis negative control is rejected. Commands, original exits, source
hash and logs remain in
`diagnostics/local-followup-checkpoint-comparison-smoke/`.

The eight-root/tighter CAS(10,11) SLEPc batch completes successfully in
3086.61 seconds (51.44 minutes) at a 5.27-GiB kernel peak. All eight raw sectors
report Krylov–Schur with eight requested eigenpairs, and all forty ensemble
root/dipole imports pass within 5.33e-10 Hartree / 7.90e-11 a.u. Common
excitation energies agree exactly at printed precision with both dense and
five-root SLEPc references; ground dipole differences are −1.54e-11 /
+8.00e-13 a.u. SCATCI/DENPROP totals are 299.95/2409.43 seconds. These timings
retain their distinct CPU placement and concurrent-workload conditions.

`diagnostics/cas11-slepc-extra-root-coverage-early/` uses the released CPU slot
to complete sixteen independent fixed-orbital probes: all eight requested
roots at spaces 80/160 in every sector. All probes fully converge and check all
64 requested roots within 5.34e-10 Hartree, with verified raw CASCI energies,
root order and spins. Its profiled wall/peak are 128.22 seconds / 0.280 GiB.
The final imported checkpoint SHA256 is
`658ae6e4758c641524983372bfeb30aa8b33be09f9fc54609755d2ecbde8a4d7`.

The first early-scan controller stops before any probes because its command
guard compared space-separated arguments against Linux's NUL-separated
`/proc/<pid>/cmdline`. The repaired controller normalizes the separator and
checks that the queued slot owner is paused and no calculation occupies the
CPU group. It resumes that waiting worker on exit. Original source/hash,
exit 1, zero-probe failure record, repaired command, exact pause/resume times
and resource records are retained in the completed diagnostic. Attached
repaired supervisor: `sh_11399be0e001NHvOpBDG2vNqYQ`; script:
`prepared/launch-early-slepc-coverage-pathfix.py`.
The 41-attempt public supplement includes this complete gate and its failed
launch. The original v3 queue still rechecks its own all-64-root gate before
the local CAS(10,12) import; the larger import remains pending.

The tight CAS(10,11) l=4 scattering baseline completes in **15924.02 seconds
(4.423 hours)** on CPUs 8–11 at a **25.20-GiB** kernel peak (23.14-GiB sampled
anonymous peak). Both B1/B2 raw CONGEN counts are **344124**; both contracted
SCATCI dimensions are **27546**, with `lwork=380327961, liwork=192840`.
Scattering SCATCI stages take 4387.34/4495.11 seconds; RSOLVE takes
56.46/59.67 seconds. All forty imported roots/dipole pass within
5.33e-10 Hartree / 7.98e-11 a.u.; the 99-point grid is complete and Pi
phase/cross-section differences are 1.00e-8 rad / 1.00e-7 bohr².
The default fitted position/full width is **2.520805/1.151002 eV**.

Relative to `co-eq-ccdz-sa10-cc40-continuum-l4`, position increases by
**64.84 meV**, width by 31.36 meV (**2.80%**), and the maximum phase difference
modulo pi is **0.109809 rad**. Position/phase gates fail; the width gate passes.
The CAS(10,11) run also tightens QC/fresh-CI controls, so this recorded-model
comparison does not isolate the active-space change. Reference inputs, raw
phase grids and exact QC-control differences are retained in
`diagnostics/cas11-tight-scattering-replays/comparison-inputs/cas10-l4/`.

Twelve saved-K-matrix replays in that diagnostic use background terms 1–4 and
detection thresholds 0.7/1.0/1.3. All succeed and reproduce the default fit.
Positions span 2.516114–2.527453 eV; full widths span 1.151002–1.296633 eV.
Maximum width change from default is **12.65%**, failing the background-width
gate. Terms 2–4 have a smaller spread; that does not remove the retained
term-1 sensitivity. The worker takes 1.20 seconds at a 0.113-GiB peak.

Twenty-four window/background replays also pass in
`diagnostics/cas11-tight-scattering-windows-unitfix/`: both B1/B2, terms 1–4,
detection 1.0, automatic intervals and clips to 1.8–3.3 / 2.0–3.15 eV.
The pinned engine's `RESONC` intersects its detected interval with
`ABVTHR`/`BELTHR`; `GETETA=.false.` uses only saved energy points. Printed
fit grids contain 39/30/23 points, approximately 1.6–3.5 / 1.8–3.25 / 2.0–3.1 eV.
Exact native grids, conversion, source/manual/licence and inputs are retained.
B1/B2 fit differences remain below 1e-6 eV. At fixed background, clipping
changes position by at most 6.03 meV and width by 3.62%; for terms 2–4 the width
changes are 0.63/0.63/0.45%. These window tests pass the chosen gates.

The original `diagnostics/cas11-tight-scattering-windows/` retains twelve
passing B1 replays and a failed first B2 replay: the driver bound unit 921,
while B2 requests 922. Native exit 2 reports EOF on `fort.922`; the controller
stops with exit 1. The fresh retry reads each template's `LUKMT` and initializes
each standalone report before B2 appends. It reproduces both default fits and
passes raw grid/fit checks. Failed/repaired walls are 0.40/0.80 seconds.
Original exits, frozen scripts and resumed CPU-slot owners are preserved.
The supplement retains **49 native replay attempts: 48 successes and one
failure**, including the twelve partial B1 successes. K-matrix/R-matrix binaries
remain on Sadaharu; the archive contains inputs, logs, fits and raw grids.

Both small CI-space controls pass in 38.50/34.69 seconds at 0.112-GiB peaks.
Each imports all forty roots within 4.98e-10 Hartree and dipole within
4.50e-11/3.07e-11 a.u., preserving the forty-component ensemble. Their printed
excitations agree; imported dipoles differ by 1.70e-11 a.u. Relative to the
earlier small dense model, the largest excitation change is 2.10e-8 Hartree.
The v3 supervisor records these gates and completes the stretched-DZ/TZ
repairs in `calibration-cas11-ci-space-retries`. Both stretched-DZ repairs
pass QC at spaces 80/160 in 389.55/363.26 seconds (0.327/0.390-GiB peaks).
Both projected equilibrium TZ retries still fail singlet B1/B2 CI convergence,
despite orbital convergence, after 807.23/806.98 seconds. The four-job batch
takes 2369.01 seconds. All original failures/checkpoints remain retained.
V3's own sixteen all-64-root CI checks also pass within 5.333e-10 Hartree,
and it records 106.58 GiB available RAM and the unchanged CAS(10,12) seed before
starting `co-eq-ccdz-sa12-slepc-target` on CPUs 8–11 in its gated 64-GiB container.

The CPU-12–15 group completes
[`calibration-sa10-tight-scattering.json`](calibration-sa10-tight-scattering.json),
`co-eq-ccdz-sa10-tight-cc40-matched`. Every numerical/continuum/channel argument
matches the completed CAS(10,11) l=4 recipe, changing only active counts to
`[4,3,3,0]` and using its own same-active-space final checkpoint. It uses the
same dense image ID `sha256:b5d5f7fc4d638d19eaeb65399b5acf20f56f4fe4526c63761b8945a0b01c5d08`.
The earlier CAS(10,10) baseline cost 1025.95 seconds / 2.84 GiB, making this
finite numerical-control repair reasonable locally. The new job has a 16-GiB
container, 6-GiB internal budgets, and pre-launch gates of 24 GiB available RAM
and 20 GiB free scratch. Source, checkpoint hash, command and slot lease are in
`prepared/cas10-matched-source/execution.json`; log:
`launch-cas10-tight-matched.log`; attached supervisor:
`sh_113b141ed001biI6gyZSACSdmN`. Its controller pauses only the still-waiting
CPU-12–15 follow-on, verifies the slot is free, and resumes it after the single
job. The successful job takes **1282.76 seconds / 2.864 GiB**, with all forty
imported roots within 4.69e-10 Hartree and dipole within 3.28e-11 a.u. Its default
candidate is **2.455964890/1.119642705 eV**. Relative to the loose CAS(10,10)
control, position/full width change by −0.544/+0.395 micro-eV and phases by
1.30e-6 rad. The matched CAS(10,10)→CAS(10,11) position/phase differences remain
**+64.8401 meV / 0.1098099 rad**, failing those gates; width changes by about
+2.80%. The new automatic-diagonalizer metadata defaults are inactive; the
eight generated target SCATCI inputs are byte-identical. The controller records
original exit zero and resumption of the waiting CPU-slot owner.

The additional finite follow-on completes after the matched controller's exact
PID exit: `prepared/matched-scattering-followup-source/`, PID **924661**,
attached supervisor `sh_113c32aec001RdAueCErEvvBzJ`. Its predecessor is PID
923482, `launch-cas10-matched.py`; its source commit is
`6d0a02dc7153ff9e58e43321c3e4fbd843d1f471` and controller SHA256
`e7b50637f2b954fbb96581c922864757227eb4b9eeb1e3515eb0eebefb1c4857`.
It leases CPUs 12–15 only if the existing coverage/scattering worker is still
waiting, has no child process, and has no calculation in the slot. Otherwise it
waits for that worker's exact exit and finished execution record. An 8-GiB
container requires 16 GiB available RAM and 20 GiB free scratch. All pause/resume
events, original exits, commands, source/unit hashes and dependencies are
recorded; the current frozen source does not alter the running predecessors.

The follow-on reanalyzes the successful matched CAS(10,10) run and performs 48
fixed-orbital probes (five/eight roots, spaces 40/80/160 in all eight sectors).
Every ensemble root must converge with the right spin and agree within
1e-7 Hartree. It compares matched CAS(10,10) with its earlier loose-control
baseline and the tight CAS(10,11) model using raw phases and native candidates.
All five ensemble roots pass every probe within 9.95e-14 Hartree, with correct
spins. Forty-six probes fully converge; the extra eighth singlet B1/B2 roots
fail at space 40 and converge at 80/160. It also runs the independent CAS(10,11)
phase-residual/held-point diagnostic in the pinned Linux image, with 24 native
fits, 72 separated starts and 120 held-point fits. This combined worker takes
**89.29 seconds / 0.215 GiB** and resumes the waiting owner. Evidence is under
`diagnostics/matched-scattering-and-independent-phase/`; the new results are
published in the verified 46-attempt supplement, preserving the 41-attempt snapshot.

The independent phase fitter has already passed 14 analytic unitary-S/mesh/
branch/held-point tests locally; all 76 project fast tests pass. Its local raw
reanalysis reproduces all 24 native window fits within 4.74e-7/2.50e-7 eV in
position/full width. Three separated starts agree to 1.55e-9 eV / 5.90e-10
relative width. Constant/linear/quadratic/cubic background held-point RMS
residuals for the automatic B1 interval are 0.016926/0.001147/0.000510/0.000128
rad. This measures the constant background's poorer predictive fit, while
retaining its width sensitivity and the existing gates. The source audit
verifies that native RSOLVE/EIGENP use 0.0735 Ryd per requested eV, an 18.44-ppm
offset from modern units; independent fits use the actual native Hartree grid.
Fitted Rydberg candidates were already converted with the modern constant.
The completed local report and actual Darwin/Python/NumPy/SciPy provenance are
retained in `prepared/independent-phase-local-analysis-20261007/`; the separate
Linux execution now completes and its fits/held-point results reproduce locally.
See [the phase-fit contract](../../docs/physics/co-electronic-qualification.md#independent-phase-fit-diagnostic).

An outer-only energy-grid diagnostic reuses the retained CAS(10,11) channel and
R-matrix-amplitude files, refining requested steps 0.05→0.025→0.0125 eV over the
same 0.1–5.0-eV interval (99/197/393 points). Its frozen controller is
`prepared/cas11-outer-grid-completionfix-source/cas11-outer-grid-completionfix.py`,
PID **929921**, attached supervisor `sh_113d31634001UT5b07TRAiNLCA`. It leases the
still-waiting CPU-12–15 group, with 4-GiB container / 16-GiB available-RAM /
20-GiB scratch gates, and preserves native stage exits, grids, cross sections,
fits and input hashes in `diagnostics/cas11-tight-outer-energy-grids-completionfix/`.
The worker completes in **381.99 seconds / 0.719 GiB**, with all twenty native
stages succeeding. Both 197-/393-point phase/cross-section grids pass completion,
finiteness, nonnegative/final-state-sum and Pi checks. The 99 common-point phases
match the baseline exactly at printed precision. Fixed-window linear-background
fits shift position/full width by at most **0.14459/0.29127 meV** (0.0254% in
width) across all three resolutions, passing the chosen energy-grid gates.
Its independent linear-background fits hold labeled-input-eV intervals
1.6–3.5 and 2.0–3.1 fixed across the three grids, using the actual native Hartree
abscissa. Native `RESONC` caps saved fit grids at **MAXFIT=100**; truncated fits
are explicitly marked and do not receive native position/width pass verdicts.

Four earlier driver attempts remain preserved: serial RSOLVE under MPI gives
rank stdin EOF; switching to MPI_RSOLVE without its required `inp` binding gives
EOF on that file; an immediate SIGSTOP-state assertion races signal delivery
before container launch; and an exact completion-string check rejects a
successful I_XSECS stage (native text is `Task has been successfully completed`).
The latter retains successful RSOLVE/EIGENP/T_MATRX/I_XSECS stages and its full
197-point B1 grid. Fresh sources repair the executable/input protocol, use a
bounded two-second pause acknowledgement and recognize both native completion
phrases. Every lease resumes the waiting worker in its exit handler.

The finite `prepared/stretched-repair-coverage-source/stretched-repair-coverage.py`
worker (PID **930259**, `sh_113d4663a001ogdQOs6BFS7aG7`) waits for PID 929921's
exact exit and completion record before borrowing that group under an 8-GiB
container. It scans
both passing stretched-DZ repairs with 48 probes each, compares roots/objective/
dipoles/core/active subspaces under the existing gates, and measures the matched
CAS(10,10) initial/final subspace overlaps. Evidence is destined for
`diagnostics/cas11-stretched-ci-repair-coverage/`. These are qualification
diagnostics; further independent import, basis and continuum checks remain open.
This worker subsequently completes in **559.24 seconds / 0.340 GiB**, with
original exit **1**, and resumes the waiting owner. Both strict all-probe gates
reject the repaired checkpoints: three five-root space-40 probes leave the
fifth ensemble root unconverged (triplet A1 on both and singlet A2 on the
space-160 checkpoint). Of 96 probes, 87 fully converge; all 64 probes at spaces
80/160 fully converge. All first-five spectra agree within 9.95e-14 Hartree,
including the failed-flag probes. Pair roots/objective/dipole agree within
2.63e-12 Hartree / 7.11e-14 Hartree / 4.81e-11 a.u., and the minimum active
overlap is 0.9999999999999984. That pair agreement does not erase the three
required convergence-flag failures. The matched CAS(10,10) continuity check
passes: minimum initial/final active/core overlaps are
0.9999999998858343 / 0.9999999999996424.

This completed diagnostic is published separately in
[`qualification-evidence/stretched-coverage/`](qualification-evidence/stretched-coverage/README.md),
preserving the frozen 46-attempt archive. Its 199 payloads and 96 spectra/spins
verify, with byte-identical repackaging. Immutable archive:
`https://data.qscat.org/ukrmol-co-stretched-coverage-2026-10-07/stretched-repair-coverage.tar.0d8f79e3de6a.gz`,
SHA256 `0d8f79e3de6a270a870dd5efc32d1452f46b46605d1a6cf89ce8f8280e933e63`.
The post-publication inspection confirms all four disjoint CPU groups are
active: staged CAS(10,12) stretched-DZ QC on 0–3, fresh CAS(10,12) aug-TZ QC on
4–7, the gated CAS(10,12) target import on 8–11, and this coverage/subspace
diagnostic on 12–15. Available host RAM is then 68 GiB and free `/home` scratch
158 GiB; the larger import is within its recorded 64-GiB container limit.

The residual/iteration follow-on completes on the released group:
`prepared/stretched-residual-iteration-source/stretched-residual-iteration.py`,
PID **932074**, attached `sh_113eac4b4001swWK2NNSGbZ5a5`. It repeats the nine
failed CI settings at the original 200 and refined 600 iterations, holding
orbitals, roots, trial space and all tolerances fixed. It explicitly measures
normalized physical and spin-penalized Hamiltonian residuals, spin and CI-vector
orthogonality. An analytic Hubbard-dimer energy/residual check and a perturbed
vector control gate that measurement. Only if all extended retries pass does it
launch full 48-probe, 600-iteration scans on both checkpoints. The original
200-cycle rejection remains in its public companion; this new diagnostic is
`diagnostics/cas11-stretched-ci-residual-iteration/`. It uses an 8-GiB container,
16-GiB available-RAM and 20-GiB scratch gates and records owner resumption.
Its original exit is **zero**, wall **789.89 seconds** and kernel peak
**0.285 GiB**. All nine original settings reproduce their failed flags at
200 cycles and pass at 600. Both full 48-probe scans then pass all requested
roots, flags, spins, orthogonality and physical/spin-penalized residual gates;
the latter maxima are **7.07e-10 / 9.94e-10 Hartree**, respectively. Maximum
ensemble-energy difference is 1.28e-13 Hartree. The analytic Hubbard-dimer
ground state has a 1.74e-15 physical residual; the perturbed vector gives
0.0504. All 114 probe attempts (nine 200/600 pairs plus 96 full-scan probes)
are retained, with 624 eigenpairs in the full scans. Raw CASCI energies, spins
and whole-probe convergence messages reconstruct from the profiled `run.log`.

The [iteration companion](qualification-evidence/stretched-iteration/README.md)
is published, publicly fetched and verified: 281 payload digests, complete raw
spectra/spins and byte-identical repackaging. Its immutable archive is
`https://data.qscat.org/ukrmol-co-stretched-ci-iteration-2026-10-07/stretched-ci-iteration.tar.c8927f5b7572.gz`,
SHA256 `c8927f5b75726424a02dff624f9f294d23851b3dc8b4152a3ef93d49259bfd3f`.
Publication bundle and verifier are retained in
`prepared/publication-stretched-iteration-20261007/`. The original 200-cycle
coverage companion remains exact. Passing the refined **600-cycle** numerical
contract releases the next import/basis controls, without qualifying a model.

The finite `prepared/stretched-import-basis-source/stretched-import-basis.py`
supervisor, PID **933694**, attached `sh_113fdbc74001wgM6TF4p34y6Xq`, waits for
v3 PID 916425's exact completion, then uses its free **8–11** group. It requires
both refined coverage scans, measured residual gates, pair agreement and the
unchanged stretched-DZ seed. The first recipe,
[`calibration-sa11-stretched-import.json`](calibration-sa11-stretched-import.json),
requests eight UKRmol roots per sector under the qualified tighter SLEPc
settings, with forty-component orbital averaging and forty retained channels.
Its all-forty-root/dipole import is followed by independent comparison of all
**64** native roots against the eight-root/space-160/600-cycle CASCI controls,
and core/active subspace comparison. It uses a 16-GiB container / 6-GiB internal
budget, with 24 GiB available RAM and 20 GiB free scratch required.

After that import passes,
[`calibration-sa11-stretched-basis-qc.json`](calibration-sa11-stretched-basis-qc.json)
projects the passing R=2.5-bohr DZ seed into same-geometry TZ/aug-TZ QC with
space 80. An orbital-converged CI failure permits one freshly named space-160
restart; all original failures remain preserved. These use 32-GiB containers
with 40 GiB available RAM / 20 GiB scratch gates. Passing basis records still
need coverage, competing-start and independent-import qualification. Commands,
sources, checkpoint/manifest hashes, original exits and resources are retained
in `diagnostics/cas11-stretched-import-basis/` and the derived batch snapshots.
CAS(10,12)'s outcome is an execution dependency here, not a scientific gate on
the already covered CAS(10,11) stretched seed.

The short retained-data propagation worker follows that exact completion:
`prepared/cas11-outer-propagation-source/cas11-outer-propagation.py`, PID
**932465**, attached `sh_113ed2b42001o6JKlUVZbZgMDZ`. The baseline has eight
propagation subranges to 100 bohr, with ten Legendre polynomials per subrange
and Gailitis matching. New controls use 16/32 subranges at 100 bohr, then
150/200-bohr matching radii with 26/36 subranges (maximum interval length
5.125 bohr). Both components retain the same 99-point energy grid and inner
amplitudes. Native and full fixed-window independent fits are retained in
`diagnostics/cas11-tight-outer-propagation/`. The pinned source's `MAXI/MAXF`
select saved T-matrix columns, not the coupled-channel space; changing them
alone would not be a channel-convergence test.

This worker now completes with original exit **zero**, **1669.54 seconds /
1.868 GiB**, and resumes the waiting owner. All forty native stages and both
components' grid/cross-section/Pi checks pass. All four refined controls have
identical printed phases, prompting a full-precision binary K audit. The audit
verifies every requested sector across all four rank logs and follows the
active pinned `BPROP`→`RPROP1_MPI`→`CURLYR`/`GAILIT`→`RPROPX` path. All four
audited engine files match the checksum-pinned upstream archive.

Binary K eigenvalue phase sums reconstruct native phases within 4.87e-8 rad.
Maximum refinement-minus-baseline phase difference is **2.56845e-6 rad modulo
pi**, while successive refined differences stay below **2.28e-10 rad**.
Fixed-window linear-background fits using full-precision binary phases change
position/full width by at most **2.63e-8 / 4.35e-8 eV**. Radial propagation
passes its chosen numerical gates for this fixed inner-region model.
The [propagation companion](qualification-evidence/outer-propagation/README.md)
preserves all binary K outputs and baseline, raw rank logs and the executable
verifier. Its 296 payloads verify, with public byte-for-byte fetch and
byte-identical repackaging. Archive:
`https://data.qscat.org/ukrmol-co-outer-propagation-2026-10-07/outer-propagation.tar.09ee4c814d9a.gz`,
SHA256 `09ee4c814d9af4148e2e54094a325322ae57a240406c88d69ffd44305b487ea5`.
Publication bundle/verifier persist at
`prepared/publication-outer-propagation-20261007/`.

The finite state-averaged throughput probe follows PID 932465:
`prepared/sa10-mpi-throughput-source/sa10-mpi-throughput.py`, PID **932604**, attached
`sh_113ef4b36001u3ews9SJlngFTy`. Its tracked input is
[`calibration-sa10-tight-mpi-scaling.json`](calibration-sa10-tight-mpi-scaling.json).
It runs the matched CAS(10,10) model at 1/2/4 ranks on CPUs 12–15, then four
one-rank jobs on individual cores in that same group. The derived concurrent
manifest is frozen with the driver. It uses the original matched dense image
and its own final checkpoint; all model/QC controls match. It checks roots,
dipoles, core/active subspaces, phases and candidates before accepting throughput
comparisons, retaining per-stage wall times, CPU time, kernel peaks, disk peaks
and concurrent-workload records. Serial jobs use 16-GiB limits / 24-GiB available
RAM; four concurrent jobs use 8-GiB limits each / 48-GiB available RAM, with
20-GiB scratch gates. Evidence will be `diagnostics/cas10-tight-mpi-throughput/`.
All three follow-ons borrow the original owner only while it is still waiting
without children; otherwise they wait for its exact recorded completion.

The subsequent resource inspection finds the first one-rank matched scaling
job running on CPUs **12–15**, alongside staged CAS(10,12) TZ QC on **0–3**,
fresh CAS(10,12) aug-TZ QC on **4–7** and the CAS(10,12) target import on **8–11**.
The host has **106 GiB available RAM** and **158 GiB free `/home`**. The remaining
scaling/concurrency entries, stretched import/basis follow-on and original
coverage/continuum queue retain their dependency gates.

The fresh-start supervisor now completes with original exit **one**. Its four
entries take **16030.77 seconds / 4.453 hours** on CPUs 4–7:

| Fresh start | Original exit | Wall (minutes) | Kernel peak bytes |
|---|---:|---:|---:|
| CAS(10,11) cc-pVTZ | 0 | 39.75 | 346562560 |
| CAS(10,11) aug-cc-pVTZ | 1 | 22.80 | 459587584 |
| CAS(10,12) cc-pVTZ | 0 | 112.87 | 799301632 |
| CAS(10,12) aug-cc-pVTZ | 0 | 91.74 | 875302912 |

The rejected aug-TZ CAS(10,11) entry converges orbitals and optimized CI, but
its fresh singlet-A1 CI flag is false. Fresh root energies match within
1.42e-13 Hartree; the flag still rejects export. Both CAS(10,12) fresh basis
starts pass all QC gates, including fresh lowest-root checks, spins, Pi
degeneracy and MO orthogonality. Their minimum initial-to-final active overlaps
are **0.152714 / 0.036790**, indicating considerable orbital reorganization.
For fresh CAS(10,12), TZ→aug-TZ lowers the averaged objective **1.795712 eV**,
but raises ground energy **0.479725 eV** and changes z dipole **−0.041475 a.u.**
The objective is a forty-component average, not the individual ground energy.
Competing starts and state/active-subspace continuity remain required.

The [fresh-basis companion](qualification-evidence/fresh-basis/README.md)
preserves all four attempts: 194 payloads and 160 raw final-state energies/spins
verify, with public byte-for-byte fetch and byte-identical repackaging.
Its archive is
`https://data.qscat.org/ukrmol-co-fresh-basis-qc-2026-10-07/fresh-basis-qc.tar.636f4d07416a.gz`,
SHA256 `636f4d07416a10ea5419fbaab7451a001df29650d7c0d9d5193b4534ce231c80`.
The publication/verifier persist at `prepared/publication-fresh-basis-qc-20261007/`.

The released **4–7** slot now runs
`prepared/fresh-basis-residual-coverage-source/fresh-basis-residual-coverage.py`,
PID **935553**, attached `sh_1142209c4001gvQ4Z7yd9u79m3`. It safely leases the
still-waiting staged CAS(10,11) owner PID **916479**, with acknowledged SIGSTOP,
no owner children and no active container in the group. Its 8-GiB container
launch records **107.17 GiB available RAM / 156.76 GiB free scratch**, above
16-/20-GiB resource gates. The three passing fresh targets receive 48
five/eight-root, space-40/80/160 probes each under the validated **600-cycle**
numerical contract, keeping tolerances, spin penalty and orbitals fixed.
The Hubbard-dimer/perturbed-vector residual controls run first. Every probe
measures physical and spin-penalized residuals, spins, CI-vector orthogonality,
ensemble energies and same-root-count agreement across trial spaces. All
144 probes/936 eigenpairs are scheduled, with evidence written to
`diagnostics/fresh-basis-ci-residual-coverage/`; qualification is pending.

That scan now completes with original exit **zero**, wall **1941.87 seconds /
32.36 minutes** and kernel peak **781508608 bytes / 0.728 GiB**. All three
targets pass every convergence flag, spin, physical/spin-penalized residual and
CI-vector orthogonality gate in all 144 probes/936 eigenpair evaluations.
Maximum physical/penalized residuals are **6.96e-10 / 9.99e-10 Hartree**;
maximum ensemble-energy, trial-space and common-first-five root-count differences
are **2.84e-13 / 3.13e-13 / 3.98e-13 Hartree**. The analytic Hubbard-dimer and
perturbed-vector controls pass. All raw spectra/spins/whole-probe convergence
blocks reconstruct from `run.log`; recorded residuals and frozen measurement
source are retained. The worker resumes PID 916479, and the fresh aug-TZ
repair immediately acquires that waiting owner's slot.

The [fresh-coverage companion](qualification-evidence/fresh-coverage/README.md)
publishes all 352 verified payloads, with public byte-for-byte fetch and
byte-identical repackaging. It includes unchanged parent QC checkpoint inputs.
Archive:
`https://data.qscat.org/ukrmol-co-fresh-basis-coverage-2026-10-07/fresh-basis-coverage.tar.2058e2e4c269.gz`,
SHA256 `2058e2e4c269a98540e6da9af7b34604ea5f277c12cf9fad50fa10b778ac2016`.
Publication/verifier persist at `prepared/publication-fresh-basis-coverage-20261007/`.

The finite fresh-lineage repair follows that exact worker completion:
`prepared/fresh-augtz-ci-repairs-source/fresh-augtz-ci-repairs.py`, PID
**935990**, attached `sh_11424888f001DdVAMZaDpGFCs1`. Its tracked input is
[`calibration-sa11-fresh-augtz-ci-repairs.json`](calibration-sa11-fresh-augtz-ci-repairs.json).
It separately restarts the rejected fresh aug-TZ checkpoint at spaces **80/160**,
with all physical settings, ensemble counts and tight tolerances unchanged.
The default 200 CI cycles remain unchanged in these QC repairs. The completed
fresh-target residual scan is an execution dependency, not a scientific gate
for the independent aug-TZ repair. It borrows the same waiting owner only while
the group is child/container-free; otherwise it waits for the owner to complete.
32-GiB containers require 40 GiB available RAM / 20 GiB scratch. Original
exits, frozen source/manifest/seed hashes and owner resumption are recorded.
Passing repair records still need coverage, competing starts and imports.

The first fresh aug-TZ space-80 repair is now active on **4–7**. The finite
`prepared/fresh-tz-import-source/fresh-tz-import.py` supervisor, PID **937849**,
attached `sh_11445fce8001MeKw10vN1gjYmq`, follows that repair supervisor's exact
completion. Its input is
[`calibration-sa11-fresh-tz-import.json`](calibration-sa11-fresh-tz-import.json).
It requires the passing fresh CAS(10,11) TZ scan, all measured residual gates
and unchanged seed SHA256
`c41b08888d05483c6845f894c94c568ce9fadb7adc831fdefc57382c32d2a73a`.
The same-basis import uses five-root/four-sector singlet/triplet orbital
averaging, eight computed roots per sector, forty retained channels, CI space
80 and SLEPc tolerance **1e-13 / 1000 cycles**. QC's 200-cycle default is
unchanged. Target-only execution uses a **16-GiB** container / **6-GiB** internal
budget, with **24 GiB available RAM / 20 GiB scratch** required. It leases the
still-waiting staged owner only while child/container-free, otherwise waiting
for that owner's completion.

The follow-on compares all **64 native roots** against eight-root/space-160/
600-cycle controls, checks the all-forty-root/dipole import, compares the
reoptimized QC roots/dipole directly with the seed, and requires matching
core/active subspaces. It records original exits, raw sectors, timings and
source/manifest/checkpoint hashes in `diagnostics/cas11-fresh-tz-independent-import/`.
This is a covered fresh-orbital import trial; competing-start and model gates
remain open. The aug-TZ repair's success is not a scientific prerequisite for
the independent, already covered TZ seed.

The matched one-rank CAS(10,10) scattering replica completes successfully in
**3105.91 seconds / 51.77 minutes**, with **3058868224 bytes / 2.849 GiB** kernel
peak and **3105.86 CPU seconds**. Its forty imported roots/dipole agree within
**5.04e-10 Hartree / 3.98e-11 a.u.**, both contracted scattering dimensions are
8350, and both complete 99-point pipelines pass their symmetry/grid checks.
The two-rank replica is now active on **12–15**. Full matched phase/candidate,
core/active-subspace and concurrency comparisons await the finite benchmark's
remaining entries; this single replica does not establish the best layout.

The latest target-import stage audit finds **all eight CAS(10,12) SCATCI sectors
completed successfully**; density/dipole processing remains active. Full import
qualification awaits that stage and the original supervisor's analysis.

Finally, the CPU-12–15 queue runs
[`calibration-sa11-tight-continuum.json`](calibration-sa11-tight-continuum.json):
two equilibrium CAS(10,11)/40-channel calculations at **l=3 and l=5**, using
the exact tight imported DZ checkpoint and the unchanged 40-component ensemble,
18-bohr sphere, deletion 1e-6 and 99-point 0.1–5.0-eV grid. They use the
qualified SLEPc target route with the existing all-spectrum scattering solver.
The completed l=4 baseline, all-64-root SLEPc verification and unchanged seed
hash are prerequisites; each new run's imported roots, threshold and dipole
are checked against the baseline. The baseline plus passing controls receive
twelve native saved-K-matrix replays each (background orders 1–4 and detection
thresholds 0.7/1.0/1.3). Phase/fit comparisons and actual contracted dimensions
still need analysis before declaring the larger model's numerical stability.

QC batches use 32-GiB containers, with at least 40 GiB available host RAM and
20 GiB free scratch required at launch. The two scattering jobs run serially
in 48-GiB containers with 6-GiB internal budgets, requiring at least 64 GiB
available RAM and 20 GiB scratch. They wait for the live v3 CAS(10,12) import
supervisor to finish, so this new large-memory work does not overlap that trial.
Every queue stops at a missing dependency or failed qualification/resource gate;
individual rejected QC attempts retain their original exits and permit the
other independent finite entries to be inspected.

The continuation is currently at a **pending-results dependency**: all four
physical-core groups have active work and finite follow-ons,
and further launches are handled by the
recorded finite queues. The restricted neutral route at R=3.0/4.0 still needs
a separately designed and validated correlation treatment. Production fitting,
the full CAS(10,12) import and actual CAS(10,12) scattering resources remain
unqualified.

### Fresh aug-TZ repair completion and remaining local work

Both fresh-lineage CAS(10,11) aug-TZ repairs complete successfully, original
batch exit **zero**, batch wall **1280.74 seconds / 21.35 minutes**:

| CI space | Original exit | Wall (seconds) | Kernel peak bytes |
|---|---:|---:|---:|
| 80 | 0 | 641.47 | 449634304 |
| 160 | 0 | 638.23 | 512544768 |

Every optimized/fresh CI flag and the usual gradient/spin/Pi/MO gates pass.
The maximum forty-root difference is **1.667e-11 Hartree**, averaged-objective
difference **−2.842e-13 Hartree** and z dipole difference **−3.230e-12 a.u.**
Minimum core/active-subspace overlaps are
**0.9999999999999999 / 0.9999999999999992**. Both share the original rejected
fresh checkpoint SHA256
`39005317ac4d60d3a91bc07107b2d8864c419a381103f983faf21387c5aa8020`;
their maximum recorded-root changes from that parent are below 5.45e-8 Hartree.
The original rejected parent remains exit one. Neither numerical agreement nor
successful later repair changes that rejection.

The pair comparator takes a small read-only, 2-GiB-bounded analysis container
on **4–7**, briefly overlapping the newly released TZ import. Its exact command,
source and overlap scope are retained in `diagnostics/fresh-augtz-ci-space-pair/`.
These measured orbital-subspace values pass the same numerical restart gate as
earlier same-basis controls. They do not qualify competing projected starts.

The [fresh aug-TZ repair companion](qualification-evidence/fresh-augtz-repairs/README.md)
publishes both successes, the exact rejected parent and the pair comparison.
All **129 payloads / 80 raw final states** verify, with public byte-for-byte
fetch and byte-identical repackaging. Archive:
`https://data.qscat.org/ukrmol-co-fresh-augtz-repairs-2026-10-07/fresh-augtz-repairs.tar.e98e41526e5e.gz`,
SHA256 `e98e41526e5e22bb18ec0446c1b788e3147e3fa29aad27643db5eefaf2ec974e`.
Publication/verifier persist at `prepared/publication-fresh-augtz-repairs-20261007/`.

The fresh TZ import PID **937849** is now active on **4–7**. It follows the
completed repair supervisor and leases waiting staged owner PID **916479**.
The finite `prepared/fresh-augtz-residual-coverage-source/fresh-augtz-residual-coverage.py`
worker, PID **939061**, attached `sh_1145a0cf2001IhceJzrfFJSBuv`, waits for that
exact import completion and its slot release. It requires the passing numerical
pair comparison and unchanged repaired seeds:

- Space 80: `744e08d2a9b5771a10e17ced9a114884e4e6af0c6a82acd3d052fc98cc8ef1df`.
- Space 160: `45f9bf646e9781cdf149d3de46d2d28967b5b57a21fc9f1be0b0f55f8d9bea4d`.

The two checkpoints receive **96 probes / 624 eigenpair evaluations**, five/eight
roots and spaces 40/80/160 with 600 CI cycles, retaining the tight tolerances,
analytic residual controls and all-root gates. Evidence is written to
`diagnostics/fresh-augtz-ci-residual-coverage/`. An **8-GiB** container requires
16 GiB available RAM / 20 GiB scratch. It uses the same acknowledged,
child/container-free waiting-owner lease or waits for that owner's completion.
The TZ import is an execution dependency; its scientific outcome does not
determine the separate augmented-TZ coverage verdict.

The two-rank matched CAS(10,10) replica has also finished; the four-rank entry
is now active on **12–15**. Matched phase/candidate, resources and concurrent
throughput verdicts await the finite benchmark's remaining entries.
CAS(10,12) DZ target density/dipole processing remains active on **8–11**,
with a measured **single-core DENPROP** process continuing to consume CPU.
Projected CAS(10,12) TZ QC remains active on **0–3**, with its averaged objective
near −112.43189855 Hartree and orbital gradients still around 1e-5 to 1e-6,
above the requested 1e-7 threshold. The original finite 150-cycle limit applies.
At this inspection the host has **104 GiB available RAM / 156 GiB free `/home`**.
These are active numerical/cost investigations; no paid host is yet required.

### Fresh TZ import completion and verifier repair

The fresh CAS(10,11) TZ engine/batch completes with original exit **zero**,
wall **3262.90 seconds / 54.38 minutes** and kernel peak **5658492928 bytes /
5.270 GiB**. All eight SLEPc sectors converge with eight roots each. Import-run
QC vs UKRmol passes all forty ensemble roots within **5.253e-10 Hartree** and
the ground dipole within **3.366e-11 a.u.**. Stage totals are QC adapter
537.91 seconds, integral preparation 0.34, CONGEN 2.86, SCATCI 371.32 and
DENPROP **2350.20 seconds**. This trial shows that target-density processing
can dominate after the selected-root eigensolver becomes inexpensive.

Its final supervisor verification originally exits **one**. A consistency
assertion uses `atol=1e-10` for the native ten-decimal energies against the
saved table; `lib/ukrmollib.pm::save_target_energies` prints `%16.9f`, so half
the last stored unit is **5e-10 Hartree**. The first read-only recheck also exits
one: later native SCATCI logs echo earlier CIDATA spectra, so parsing the first
`EIGEN-ENERGIES` block includes inherited sets. The repaired verifier selects
the unique spectrum after `CI data will be stored as set number`, checks the
sector number and compares exact decimal tokens under that format contract.
Analytic controls reject corrupted tables, inherited-only spectra and missing
final sets. Neither failure's source/log/exit is overwritten; physical root
and dipole gates remain **1e-7 Hartree / 1e-5 a.u.**.

The passing recheck verifies all **64 native roots** against eight-root/
space-160/600-cycle controls within **3.196e-8 Hartree**. The original fresh
seed is slightly reoptimized in the import run: forty QC roots change by at
most **3.109e-8 Hartree**, while native dipole changes **−5.119e-8 a.u.**
The all-64 comparison includes that reoptimization. Minimum core/active
subspace overlaps are **0.9999999999999983 / 0.9999999999997468**. The read-only
single-core verifier takes **0.654 seconds**, with brief overlap of the active
augmented-TZ coverage explicitly retained in its execution record.

The [fresh TZ import companion](qualification-evidence/fresh-tz-import/README.md)
preserves the successful engine run, two verifier failures and passing recheck.
All **181 payloads**, 64 final-set native roots and 48 raw reference spectra/spins
verify, with public byte-for-byte fetch and byte-identical repackaging. Archive:
`https://data.qscat.org/ukrmol-co-fresh-tz-import-2026-10-07/fresh-tz-import.tar.e91bbc4ac41e.gz`,
SHA256 `e91bbc4ac41e6c49c228fe33bc5634cd3d0d19759e0c8afac1ab2237e4105d49`.
Publication/verifier persist at `prepared/publication-fresh-tz-import-20261007/`.
Failure/recheck evidence is in `diagnostics/cas11-fresh-tz-independent-import/`,
`diagnostics/cas11-fresh-tz-import-rounding-recheck/` and its `-v2` companion.
Competing-start and electronic-model qualification remain open.

The independent augmented-TZ residual worker PID **939061** now runs on **4–7**.
The finite `prepared/fresh-augtz-import-source/fresh-augtz-import.py` supervisor,
PID **942825**, attached `sh_1148b7b36001uIeqQdq1hlO2N4`, follows its exact
completion and requires **both** 48-probe scans, the numerical pair gate and
the unchanged repaired seed hashes to pass. Its tracked input is
[`calibration-sa11-fresh-augtz-import.json`](calibration-sa11-fresh-augtz-import.json).
It uses the space-80 seed, eight target roots per sector, forty retained channels,
SLEPc tolerance 1e-13 / 1000 cycles and the same tight QC/physical settings.
Target-only execution uses a **16-GiB** container / **6-GiB** internal budget,
with **24 GiB available RAM / 20 GiB scratch** required. It borrows the waiting
staged owner only under an acknowledged child/container-free lease, otherwise
waiting for that owner to finish. The repaired final-set/decimal verifier is
frozen as `verify-fresh-target-import.py`; all 64 roots, ground dipole and
core/active subspaces must pass before completion. Evidence is written to
`diagnostics/cas11-fresh-augtz-independent-import/`. This qualification trial
does not settle competing starts or electronic-model convergence.

The parameterized `verify-fresh-target-import.py` also receives a known-TZ
differential control in `diagnostics/fresh-target-verifier-control/`: original
exit **zero**, **0.293 seconds**, reproducing every final-set/decimal/root/dipole/
subspace result of the passing TZ recheck. This brief one-thread analysis runs
inside the active CAS(10,12) import cgroup and records its command/source hash.

The augmented-TZ coverage worker now completes with original exit **zero**,
wall **600.90 seconds / 10.02 minutes** and kernel peak **419938304 bytes /
0.391 GiB**. Both checkpoints pass all **96 probes / 624 eigenpair evaluations**,
including every CI convergence, spin, physical/spin-penalized residual and
orthogonality gate. Maximum physical/penalized residuals are **6.62e-10 /
9.94e-10 Hartree**; maximum ensemble, trial-space and common-first-five
root-count energy differences are **1.85e-13 / 1.99e-13 / 1.99e-13 Hartree**.
The analytic Hubbard-dimer and perturbed-vector controls pass. Every raw
spectrum/spin/whole-probe convergence block reconstructs.

The [augmented-TZ coverage companion](qualification-evidence/fresh-augtz-coverage/README.md)
publishes all **253 verified payloads**, with public byte-for-byte fetch and
byte-identical repackaging. Archive:
`https://data.qscat.org/ukrmol-co-fresh-augtz-coverage-2026-10-07/fresh-augtz-coverage.tar.3f8817e75644.gz`,
SHA256 `3f8817e756449ae26054a69bee936dfa549afd5691c295a45a78144218521d1b`.
Publication/verifier persist at `prepared/publication-fresh-augtz-coverage-20261007/`.
The coverage worker resumes PID 916479; the import supervisor PID **942825**
immediately acquires the waiting owner's acknowledged lease. The 64-root
aug-TZ import is now active on **4–7**, with **106.60 GiB available RAM /
151.46 GiB free scratch** recorded before its 16-GiB container launch.

All three matched MPI scaling entries now complete with original exits zero.
Container wall times are **3106.46 / 1831.14 / 1240.12 seconds** for one/two/four
ranks on the same physical-core group; the batch takes **6177.72 seconds**.
Four concurrent one-rank replicas are now active on CPUs **12 / 13 / 14 / 15**.
Final numerical-equivalence and throughput comparisons await their completion.
CAS(10,12) density/dipole processing and projected TZ QC continue on **8–11**
and **0–3**, respectively. Further local work remains reasonable.

### Matched MPI and throughput completion

The matched benchmark supervisor PID **932604** now completes with original
exit **zero**, resuming CPU-slot owner PID **918330** after the acknowledged
lease. All seven runs and both batches have original exits zero. The resumed
coverage/continuum queue successfully reanalyzes the fresh TZ QC record and
advances to its original finite coverage/continuum prerequisites on **12–15**.

| MPI ranks | Runner wall (seconds) | Aggregate CPU (seconds) | Kernel memory peak (GiB) | Full-job speedup |
|---:|---:|---:|---:|---:|
| 1 | 3105.91 | 3105.86 | 2.849 | 1.000× |
| 2 | 1830.62 | 3435.09 | 2.743 | 1.697× |
| 4 | 1239.61 | 4268.46 | 2.822 | 2.506× |

The scattering SCATCI stage alone speeds up **1.980× / 3.764×** at two/four
ranks. All three sequential rank entries use CPUs **12–15**, with 16-GiB
limits and one thread per numerical library; batch wall is **6177.72 seconds**.
Four one-rank replicas pinned individually to **12 / 13 / 14 / 15**, with
8-GiB limits each, finish together in **3444.36 seconds / 57.41 minutes**.
Their runner walls are 3442.00/3412.68/3380.81/3443.62 seconds.
`4 × 1239.61 / 3444.36 = 1.43958`: **43.96% higher throughput / 30.54% less
four-job elapsed time** against repeating the measured four-rank runner.
The four-job serial reference is extrapolated, not a separate measured batch.
Aggregate CPU is **13675.56 seconds**, **19.90% less** than four repetitions
of the measured four-rank aggregate CPU. Summed per-job kernel peaks bound
the concurrent envelope at **12178423808 bytes / 11.342 GiB**, with a
**5368893440-byte / 5.000-GiB** summed run-disk peak envelope. Profiles use
per-run elapsed clocks, so these are upper envelopes rather than synchronized
instantaneous measurements. The exact concurrent host workloads are recorded
before both batches; timings are single-shot measurements on one geometry.

All seven replicas pass target/import, grid, nonnegative cross-section, final-state
sum, Pi and matched numerical gates. Both contracted scattering dimensions
are **8350**. Their seed is the baseline's final checkpoint SHA256
`97d8bc1166e34d0d98d826e71d98bc6413cfc644572ef19229500726f5ac31d4`;
the baseline itself started from an earlier continuum checkpoint. That lineage
change passes explicit checks: forty QC roots/dipoles differ by at most
**6.119e-9 Hartree / 2.986e-8 a.u.**, minimum active-subspace overlap
**0.9999999999998141**, and core subspaces satisfy the same numerical gate.
Import roots/dipoles agree with each replica's QC within **5.044e-10 Hartree /
4.295e-11 a.u.**. Maximum baseline phase, native position and full-width
changes are **2e-7 rad / 0.137 micro-eV / 0.041 micro-eV**. This establishes
execution-layout equivalence at that model size; larger spaces/channels need
their own measurement, and electronic-model/extraction gates remain open.

The [MPI/throughput companion](qualification-evidence/mpi-throughput/README.md)
publishes all seven runs, both frozen batches, exact published baseline inputs,
native outputs, profiles and source/lease records. All **1221 payloads** verify;
disposable-copy reanalysis of all eight records and **98997 raw profile samples**
reconstruct the numerical and resource checks. Public fetch is byte-for-byte and
repackaging byte-identical. Archive:
`https://data.qscat.org/ukrmol-co-matched-mpi-throughput-2026-10-07/mpi-throughput.tar.619491c1b75a.gz`,
SHA256 `619491c1b75afba0f468a393ab2d4f0aeff3cb1df86992da27fd1cf7b85cae02`.
Publication/verifier persist at `prepared/publication-mpi-throughput-20261007-v2/`.
The first publication is retained in the companion manifest; only verifier
formatting and its analysis-source digest change, with identical numerical evidence.
Tracked recipes are `calibration-sa10-tight-mpi-scaling.json` and
`calibration-sa10-tight-throughput.json`; completed directories remain evidence.
The approximately $200 external bottleneck cap remains unchanged. Continue
aug-TZ and CAS(10,12) imports, staged QC and the resumed finite continuum queue
on Sadaharu before choosing that experiment.

### Repaired fresh aug-TZ import completion

The fresh-lineage aug-TZ import supervisor PID **942825** completes with original
engine/batch, supervisor and final-verifier exits all **zero**. Wall is
**3360.53 seconds / 56.01 minutes**, kernel peak **5656363008 bytes / 5.268 GiB**.
Stage totals are QC adapter **643.11 seconds**, integral preparation 1.25,
CONGEN 2.81, eight SCATCI sectors **394.14**, and DENPROP **2319.03**.
All **64 native roots** agree with independent eight-root/space-160/600-cycle
controls on the space-80 seed within **5.969e-11 Hartree**; the separate repaired
space-160 seed agrees within **5.764e-11 Hartree**. All forty ensemble import
roots agree with import-run QC within **4.672e-10 Hartree**, and DENPROP dipole
within **5.169e-12 a.u.**. Seed reoptimization changes a QC root by at most
**3.621e-11 Hartree**, and native dipole changes **+7.474e-11 a.u.** from that
seed. Minimum core/active subspace overlaps are
**0.9999999999999984 / 0.9999999999999990**. Physical root/dipole gates remain
1e-7 Hartree / 1e-5 a.u.; the final-set parser and nine-decimal saved-table
contract pass their previously validated controls.

The [aug-TZ import companion](qualification-evidence/fresh-augtz-import/README.md)
publishes **181 verified payloads**, both exact repaired seed inputs, all
64 native roots and all **96 raw reference spectra / 624 eigenpair evaluations**.
Public fetch is byte-for-byte and repackaging byte-identical. Archive:
`https://data.qscat.org/ukrmol-co-fresh-augtz-import-2026-10-07/fresh-augtz-import.tar.3568cc77de32.gz`,
SHA256 `3568cc77de32e415fcac6b403649ab7002636ee93b076f8fdcf7b2f0f75d0437`.
Publication/verifier persist at `prepared/publication-fresh-augtz-import-20261007/`.
The rejected fresh parent stays exact in the earlier repair companion.
Competing starts and electronic-model convergence remain open.

### Released-core continuation

The completed aug-TZ importer resumes waiting staged owner PID **916479**.
The finite `prepared/fresh-basis-competing-starts-source/fresh-basis-competing-starts.py`
worker, PID **951086**, attached `sh_114cd009700149SQwzGgPv3IBf`, then safely
leases **4–7** under an acknowledged child/container-free stop. Its recipe is
[`calibration-sa11-fresh-basis-competing-starts.json`](calibration-sa11-fresh-basis-competing-starts.json):
four same-geometry TZ↔aug-TZ projections at spaces **80/160**. Both qualified
seeds retain their exact checksums, and both prior 64-root import verdicts must
pass before launch. Ten active electrons, two inactive core orbitals, forty
equal-weight ensemble components, tight tolerances, **200 CI cycles / 150 orbital
cycles** and default root counts remain fixed. Each run uses a **16-GiB** container,
with **24 GiB available RAM / 20 GiB scratch** required. Passing targets receive
up to **six** comparisons: two same-basis numerical pairs and four comparisons
against the corresponding qualified fresh target, including roots, averaged
objective, dipole and core/active subspaces. Individual failed QC or comparison
exits stay recorded. The worker writes evidence to
`diagnostics/cas11-fresh-basis-competing-starts/` and resumes the staged owner.

After MPI completion, the original coverage/continuum owner PID **918330**
passes its three fresh-target ensemble coverage checks, then waits for the
whole CAS(10,12) v3 supervisor before advancing. The fixed-DZ continuum work's
own launch gates already pass: both small CI-space imports, the tight scattering
baseline, 64-root SLEPc coverage and the unchanged qualified seed. A finite
successor, PID **950692**, attached `sh_114cab6e0001xf6xkkOHRS4TQT`, now advances
the original l=3/l=5 recipes on **12–15**, independently of that long density stage.
Source/record: `prepared/local-continuum-ready-source/continuum-ready.py`.

The predecessor is stopped and acknowledged, checked for no children or active
slot containers, then receives **SIGINT** and resumes to run its final-record
handler. Its original SSH exit is **130**, a recorded intentional scheduling
stop of an idle supervisor; no numerical calculation is interrupted.
`predecessor-retirement.json`, `predecessor-completion.json` and the raw
`predecessor-console.log` preserve the before/final execution and reason.
The successor verifies every scientific gate before the transition. Negative
scheduling controls reject an active child and an occupied CPU slot without
sending SIGINT; their results are frozen in `scheduling-controls.json`.

The successor retains the original `calibration-sa11-tight-continuum` batch/job
names, **48-GiB** limits and **64 GiB available RAM / 20 GiB scratch** gates.
Each passing continuum target is compared with the baseline threshold, all
common excitations and dipole. Twelve native B1 background/detection replays
follow for the baseline and each passing l=3/l=5 record; exact passing existing
replays can be reused by digest. The earlier fresh-basis QC/repair/coverage/import
records remain preserved in their published companions.

All four physical-core groups now have active work: projected CAS(10,12) TZ QC
on **0–3**, cross-basis CAS(10,11) QC on **4–7**, CAS(10,12) density processing
on **8–11**, and the l=3 continuum entry on **12–15**. Host resources are
**106 GiB available RAM / 150 GiB free `/home`** at this inspection. Further
local work remains reasonable; the approximately $200 first external bottleneck
cap stays deferred until these finite outcomes identify an actual blocker.

### Competing-start completion and projection repair

The four-entry competing-start supervisor PID **951086** completes with original
batch/controller exits **one**, and resumes staged owner PID **916479**. Batch
wall is **3444.76 seconds**. Both aug-TZ→TZ projections fail in **0.403 seconds**
before QC with `RuntimeError: Too many orbitals in mo_init (try passing only the
occupied orbitals)`. Both TZ→aug-TZ entries pass QC in **1737.00 / 1705.00 seconds**,
with peaks **455372800 / 515391488 bytes**. Their CI-space-80/160 pair passes:
forty-root difference **5.941e-12 Hartree**, objective difference **−1.137e-13
Hartree**, z dipole difference **+1.901e-11 a.u.**, minimum active overlap
**0.9999999999999984**. Both comparisons against the earlier fresh aug-TZ target
exit one: objective **−0.538388 eV**, ground **−1.115756 eV**, maximum root shift
**2.126321 eV**, z dipole **−0.06670992 a.u.**, minimum core/active overlaps
**0.99983956 / 0.14108460**. This is a separately converged competing orbital
solution; earlier fixed-orbital coverage/import does not establish a unique
optimum. Retain both branches and qualify the new one before using it.

The [competing-start companion](qualification-evidence/competing-starts/README.md)
publishes all four attempts, original **1/1/0/0** exits, all three comparison
records with original **0/1/1** exits, source/lease records and exact qualified
seed inputs. All **277 payloads** and **80 raw final states** verify; every
root/objective/dipole comparison reconstructs and all three orbital-subspace
comparisons independently recompute from checkpoints. Public fetch is
byte-for-byte and repackaging byte-identical. Archive:
`https://data.qscat.org/ukrmol-co-competing-starts-2026-10-07/competing-starts.tar.b5193060808b.gz`,
SHA256 `b5193060808b6d5eabe0b75cfb2d6d7e0edf0c1cfabaf785a51fcacae2ceafeb`.
Publication/verifier persist at `prepared/publication-competing-starts-20261007/`.

`target.project_initial_orbitals` repairs the actual downward interface: when
the source MO column count exceeds the destination, pass only inactive-core
plus active columns to PySCF, which builds the destination virtual complement.
The exact failed **92→60-column** checkpoint pair now projects its **13**
core+active columns with active irreps `[5,3,3,0]` and **1.932e-14** orthogonality
error. Existing same-basis and upward routes retain the original full-input API;
simultaneous geometry/basis changes still reject. Initial-orbital provenance
records source and supplied column counts. Four genuine PySCF regressions and
all **80 CO tests with PySCF 2.11.0** pass. Patched source, exact-checkpoint control
and regression log are included separately from the unchanged original run sources.

The new aug-TZ checkpoints are space-80
`2f0cc1e4cbd31ede3012cf19a9584ddefa88b9240b990ec908d5f2c100bc78fb` and space-160
`95ba854810510cb75580bc1cdc7eb058a3bd06ff983fdbb91bd99a77d24b9d7e`.
The coverage worker PID **953666**, attached `sh_115094bcd001FbtnuSGswyXvFq`,
leases waiting owner **916479** on **4–7** after the acknowledged child/container-free
stop. Source: `prepared/competing-augtz-residual-coverage-source/competing-augtz-residual-coverage.py`.
It checks both passing original targets and their numerical pair, exact seed
hashes, the analytic Hubbard-dimer/residual negative control, then all **96**
fixed-orbital probes: 5/8 roots × spaces 40/80/160 × eight sectors × two seeds,
with **600 cycles** and unchanged 1e-12 energy / 1e-9 residual tolerances.
Raw spectra, spins, solver flags, physical/penalized residuals and resources go
to `diagnostics/competing-augtz-ci-residual-coverage/`. Limit is **8 GiB**, with
16 GiB available RAM / 20 GiB scratch required; launch records **100.77 GiB /
148.13 GiB**. Ordinary target CI remains 200 cycles.

The downward-retry supervisor PID **954343**, attached `sh_1150eecdc0016n7grzLOHvIeDn`,
waits for that exact coverage worker and both passing scan verdicts. It then
leases the same waiting staged owner for
[`calibration-sa11-tz-projection-coreactive-retries.json`](calibration-sa11-tz-projection-coreactive-retries.json):
four new-name TZ QC projections, spaces 80/160 from each original-fresh and
new-competing aug-TZ branch. Limits/tolerances remain **16 GiB / 200 CI cycles /
150 orbital cycles**; source/seed hashes and original import gates stay explicit.
Passing runs receive up to six pair/fresh-TZ comparisons. All exits are retained,
and the worker resumes its owner in its final handler. Source:
`prepared/tz-projection-coreactive-retries-source/tz-projection-coreactive-retries.py`.
CAS(10,12) QC/density, the independent l=3/l=5 continuum queue and later stretched
import/basis dependencies continue. Available host RAM is approximately **100 GiB**;
the approximately $200 external experiment remains deferred.

### Scientific decisions still pending

1. **Qualify the larger-active-space targets.** The tightened CAS(10,11)
   import passes; retain its completed start and root-count refinements and
   finish the running equilibrium scattering and matched QC controls.
   CAS(10,12) needs a larger-RAM dense execution host or
   separately gated eigensolver configuration; the installed-library audit rules
   out the current dense layouts on this host. Compare
   against the smaller active spaces before adopting a target for scattering.
2. **Qualify the common orbitals.** The `state-averaged` backend constructs
   multi-spin/multi-irrep CAS targets. The existing `natural` option optimizes
   only the ground singlet. Verify ensemble completeness, excited energies,
   dipoles, orbital/state identities and continuity as R changes.
3. **Check electronic convergence.** Compare DZ/TZ atomic bases, active/external
   spaces and about 40–50 C2v target components at R=1.9, 2.1323 and 2.5 bohr.
   Published values are plausibility checks with model-dependent spread, not
   error bars or a requirement to force agreement.
4. **Recheck continuum and extraction.** Start from radius 18 bohr, l=4, double
   precision and deletion 1e-6, then vary those controls and fit window/grid/
   background for the new electronic model. Use the provisional 0.05 eV,
   5% width (0.001 eV floor) and 0.05 rad modulo-π criteria from the README.
   Diagnose poles/bound states at threshold and track the same physical feature.
5. **Qualify the correlated neutral curve.** Extend the completed CCSD(T) pilot's
   basis sequence and test the stretched correlation treatment before computing
   a full curve. Keep the relative scattering resonance and the neutral curve's
   reference convention explicit.
6. **Remeasure cost, then cover geometry/energy.** Qualify the electronic model
   before the tentative 25–40-geometry, 300–800-energy campaign. Multi-day
   execution is acceptable. Use adaptive windows where fits narrow, and measure
   throughput/memory again at the selected model/channel size.
7. **Connect to potential fitting.** `projects/potential_factory/target.py`
   supplies neutral/resonance target holders; eigenphase fitting and the polar
   tail still need development. Fixed-nuclei data do not uniquely identify a
   local radial potential. HTTP and AWS execution are subsequent capabilities.
