# CO experiment: resources and continuation

## Live qualification handoff

**Local low-memory milestone:** the fixed CAS(10,10) native sparse/MPI
2048-root replay now qualifies in both Pi sectors at **1.861 GiB**, with
complete-grid phase error **2.0e-7 rad** and passing identical fixed-window
position/full-width comparisons. Its independent dense-omission refinements
remain stable through the complete spectrum. Public reconstruction passes.
CAS11 dense omission also passes at 2048 poles; its native finite sequence
continues locally. Both stretched CAS12 DZ repair checkpoints now pass their
independently reconstructed 600-cycle coverage scans. Paid provisioning remains
deferred while these gated Sadaharu workers run.

Fresh CAS(10,11) TZ and both repaired aug-TZ targets now pass QC, independent
600-cycle fixed-orbital coverage and all-64-root UKRmol import/dipole checks.
The original rejected fresh aug-TZ record remains exact. Subsequent competing
starts find a distinct, reproducible aug-TZ branch whose objective is 0.538388 eV
lower; both new-branch coverage scans and all four repaired downward projections
pass. Its independent native import now passes all 64 roots; both aug-TZ
branches remain distinct.
The matched CAS(10,10)
one-/two-/four-rank and concurrent benchmarks also pass all seven replicas.
CAS(10,12)'s tight QC starts, ensemble coverage and full forty-root/dipole import
now pass locally in 6.201 hours / 40.405 GiB. Native scattering preparation
measures dimension 86352 and valid aggregate array floors above 222 GiB,
blocking the audited dense full-spectrum scattering solve on Sadaharu. Its staged
DZ/basis failures feed the eligible numerical-repair queue.
Stretched CAS(10,11) DZ also passes its 64-root/dipole import after a
retained-data decimal-format recheck; its original verifier failure is preserved
and same-geometry TZ/aug-TZ QC also passes. Both independent basis coverage
scans pass; the stretched aug-TZ 64-root import also passes. Stretched TZ's
original extra-root gate rejects triplet-A1 root eight, so a fixed-orbital
diagnostic on the imported checkpoint now separates native solver error from
seed-orbital drift: all native roots match independent CI at those actual
orbitals within 5.01e-11 Hartree. The original seed-gate failure is retained.
A 26-point-per-basis near-equilibrium neutral pilot has completed and its
held-point/basis diagnostics independently reconstruct from the public companion.
The stretched RHF reference remains rejected. The next decisions
depend on these finite workers:

| Physical CPUs | Active work | Follow-on dependency |
|---|---|---|
| 0–3 | Released: CAS11 dense-omission stage completes with exit zero | All six pole counts pass both-sector phase/fixed-window gates; native qualification continues on 4–7 |
| 4–7 | CAS11 sparse owner PID 1020516, 24-GiB native containers | Native 2048-root control passes both sectors; 4096-root refinement active, then declared 8192/16384 controls |
| 8–11 | Electronic successors PID 1022822, stretched CAS12 aug-TZ QC active | All 96 fifty-component coverage probes and stretched TZ QC pass; DZ 64-root import follows with its separate resource/import gates |
| 12–15 | Released: anchor owner PID 988205 completes with controller exit zero | Both compressed and stretched scattering anchors complete; stretched pipeline independently reconstructs, while extraction/continuity qualification remains open |

The experiment roots remain `/home/kooza/ukrmol/co-sa-20261006` and
`/home/kooza/ukrmol/co-neutral-20261006`. The approximately $200 first paid
experiment ceiling remains deferred; all launches here use Sadaharu.
**Preferred next solver task:** qualify native sparse/iterative scattering
locally before considering a large-host allocation. See
[the task list](#sparseiterative-scattering-qualification--preferred-before-a-large-host)
and the [reproducer contract](../../docs/physics/co-electronic-qualification.md#sparseiterative-scattering-qualification-contract).
The new queues are detailed under [Local continuation — 7 October](#local-continuation--7-october-2026).

### Hourly observation — 8 October 2026

The user requests hourly observation and gated continuation for approximately
seven hours. Watch supervisor **95764**, attached shell
`sh_118e88aa6001JoIVbq5z7laOtT`, starts at **2026-10-08 00:27:35 UTC** and
checks immediately, then hourly through **07:27:35 UTC**. It keeps the local
machine awake while the bounded watcher runs; existing scientific queues keep
their original finite contracts and independently advance eligible stages.

Each check records owner steps/original exits, physical-core/container ownership,
available RAM and artifact disk. Completed components receive fresh additive
captures and independent reconstruction of raw records, source/checkpoint
hashes, resource profiles and applicable scientific gates. Failed or incomplete
gates remain explicit blockers. No running source/output directory is reused,
no failure is erased and no CAS12 scattering or paid provisioning is released
by the watcher.

The persistent host journal, raw captures and mirrored `STATUS.md`/owner record
live under `prepared/hourly-observation-20261008/`. Local execution, captures,
verification reports and consolidated handoff live under
`/private/var/folders/4k/bqw_pvrj3_z97dhlkzpnqs3r0000gn/T/opencode/co-hourly-observation-20261008/`.
Frozen watcher/collector/verifier sources and their digests travel with those
records. Hour-zero reconstruction already passes the **96 fifty-component
coverage probes** and the **stretched CAS11 anchor** pipeline checks.

The latest completed CAS11 **2048-root** native control has both-sector
observable verdicts true: maximum phase error **5.1e-6 rad**, energy error
**1.194e-12 Hartree**, residual **1.330e-10 Hartree**, and fixed-window
position/full-width differences **7.669e-7 / 1.291e-5 eV**. Two-sector replay
cost is **8836.11 seconds / 7.693 GiB**. Higher native spectral refinement is
still required before a larger-model release.

Both fifty-component coverage scans pass all **48 probes per checkpoint**;
ensemble-root errors are below **1.848e-13 Hartree**. Scan cost is **1280.58
seconds / 0.6583 GiB**. The stretched CAS12 TZ QC successor also passes its
orbital and original/fresh CI flags in **7160.31 seconds / 0.8181 GiB**.
Native fifty-component import and independent larger-basis coverage remain open.

Stretched CAS11 anchor owner **988205** completes all five steps with controller
exit zero at `finished_unix=1791403926.7180297`. The stretched calculation costs
**12239.26 seconds / 25.197 GiB**; target energy consistency is **5.460e-10
Hartree** and Pi phase splitting **1.0e-8 rad**. The native candidate is
**0.983905 / 0.292887 eV** position/full width, near the upper edge of its
0.01–1.0-requested-eV pilot grid. Its fit qualification remains unset pending
window/background and near-threshold identity/continuity controls. Raw native
outputs do not show a MAXFIT diagnostic. The independent pipeline reconstruction
checks **407 payloads / 60977 profiles**; public companion publication remains
separate follow-on work.

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

### Competing aug-TZ branch coverage completion

Coverage supervisor PID **953666** finishes with original profiled/container
exit **zero**, releasing its acknowledged lease and resuming staged owner
**916479**. Both scans pass all **96 probes / 624 eigenpair evaluations**, with
maximum physical/penalized residuals **6.630e-10 / 9.977e-10 Hartree**. Ensemble,
same-root-count space and common-first-five root-count differences are at most
**2.274e-13 / 1.706e-13 / 1.564e-13 Hartree**, CI orthogonality **4.771e-14**.
Wall is **610.35 seconds / 10.17 minutes**, kernel peak **407740416 bytes /
0.380 GiB**, aggregate CPU **1953.75 seconds**. Original analytic controls pass;
ordinary target CI remains 200 cycles. Both seed hashes remain unchanged.

The [coverage companion](qualification-evidence/competing-augtz-coverage/README.md)
publishes **335 verified payloads**, with every raw spectrum/spin/convergence
block reconstructed, recorded action residuals checked and comparisons recomputed.
Public fetch is byte-for-byte and repackaging byte-identical. Archive:
`https://data.qscat.org/ukrmol-co-competing-augtz-coverage-2026-10-07/competing-augtz-coverage.tar.b0e24bad20f5.gz`,
SHA256 `b0e24bad20f5e97f67bc7693184a4b9ed1da6b903b5b748fb6e579dfd32732ae`.
Publication/verifier persist at `prepared/publication-competing-augtz-coverage-20261007/`.

Retry supervisor **954343** immediately acquires the acknowledged waiting-owner
lease on **4–7**, launching its four new-name downward projections after both
coverage verdicts pass. Launch resource gate records **100.47 GiB available RAM /
148.11 GiB free scratch**; each QC container remains **16 GiB**.

The finite import supervisor PID **954948**, attached `sh_11515817e001Una5gIFHyTtyl6`,
now waits for that exact worker's completion. Source:
`prepared/competing-augtz-import-source/competing-augtz-import.py`; recipe:
[`calibration-sa11-competing-augtz-import.json`](calibration-sa11-competing-augtz-import.json).
It imports the new space-80 seed with eight roots/sector, forty equal-weight
ensemble components and forty configured retained states. SLEPc tolerance/cap
remain **1e-13 / 1000**, internal budget **6 GiB**, container **16 GiB**, with
24 GiB available RAM / 20 GiB scratch required. Both full coverage scans and the
numerical pair, exact seed hashes and all recorded physical/penalized residuals
must pass. Final verification uses the validated final-CIDATA-set/decimal parser,
all 64 roots, dipole and core/active subspaces against the covered reference.

The scheduling dependency is completion of the downward worker, retaining its
original exit even when comparisons fail. It does not redefine that scientific
verdict: this independent import has its own explicit QC/pair/coverage gates.
The previous fresh-lineage import and new-branch disagreement remain evidence.
The import worker leases/resumes the same staged owner in its final handler;
CAS(10,12) and fixed-DZ continuum tasks remain on their disjoint core groups.
Further local work remains reasonable, with no paid provisioning.

### CAS(10,12) local target-import completion

Supervisor PID **916425**, attached `sh_113d599540014G0jyyj7Y6Sndi`, completes with
original SSH/controller exit **zero**. Its original finite step exits are
**0 / 1 / 0**: passing small CI-space controls, preserved CAS(10,11) retry failure,
and passing `calibration-target-slepc-cas12`. That batch/engine exits zero in
**22324.47 seconds / 6.201 hours**, aggregate CPU **9.437 hours**, kernel peak
**43384504320 bytes / 40.405 GiB** in the 64-GiB container. The exact restart
seed remains `9451c7baed762ab674e9d36111cc220202972703ac52346e1ce13dcdc5ac55a0`.

All **40 native target roots** pass against current QC within **4.801e-11
Hartree**, covered seed first-five spectra within **6.815e-9 Hartree**, and the
independent RHF-start target within **3.485e-8 Hartree**. Saved-table/current-QC
error is **5.465e-10 Hartree**, with every native-to-saved exact-decimal rounding
error ≤5e-10 Hartree. Ground DENPROP dipole agrees with current QC within
**1.279e-11 a.u.**, the seed within **3.810e-9 a.u.** and the RHF start within
**5.881e-8 a.u.** Minimum initial/final core/active overlaps are
**0.9999999999999997 / 0.9999999999999988**; RHF-start/final active overlap is
**0.9999999999995508**. Tight QC/fresh-CI, spin, MO/Pi and physical import gates
remain unchanged. The computed native root inventory is five per eight sectors;
the original extra-root/space-40 coverage failures remain separate failures.

Stage totals: QC **1472.31 seconds**, integral preparation 0.095, eight CONGEN
stages 12.22, eight SCATCI sectors **2808.79 seconds / 46.81 minutes**, serial
DENPROP **18030.97 seconds / 5.009 hours**. DENPROP is **80.768%** of profiled
wall, so more MPI ranks alone do not remove the measured density bottleneck.
The largest dense PETSc target matrix floor stays **37.41 GiB**, now with an
observed complete-job peak of **40.405 GiB**. Numerical target import fits locally;
the all-spectrum scattering path still needs its actual contracted dimensions,
workspaces and memory/scaling measurements.

The [CAS(10,12) companion](qualification-evidence/cas12-import/README.md) verifies
**335 payloads**, all native spectra/dipole, **40 raw QC states**, **48 raw
reference probes / 312 eigenpairs**, both seed/RHF-start subspace comparisons
and **111339 raw resource samples**. All 174 copied prior QC/start/coverage
inputs match the immutable main archive; two original extra-root failures and
the original v3 controller's failed CAS(10,11) batch remain preserved. Public
fetch is byte-for-byte and repackaging byte-identical. Archive:
`https://data.qscat.org/ukrmol-co-cas12-import-2026-10-07/cas12-import.tar.0173df04a769.gz`,
SHA256 `0173df04a7694f1d13a230bc1bc8a6c4f495c78ea23e780d08f9b94acdd0d195`.
Publication/verifier persist at `prepared/publication-cas12-import-20261007/`.

The stretched import/basis supervisor PID **933694**, attached
`sh_113fdbc74001wgM6TF4p34y6Xq`, now advances on released **8–11**. Its
`co-r2500-ccdz-sa11-slepc-roots8-ci-space80-target` container is active; the
subsequent same-geometry TZ/aug-TZ QC and conditional space-160 repairs retain
their independent QC/coverage/import gates. Staged CAS(10,12) TZ QC remains
active on **0–3**, downward competing-start retries on **4–7**, and l=3 continuum
on **12–15**. Available RAM is approximately **99 GiB** and free `/home` **148
GiB** at this transition. Further finite Sadaharu work remains reasonable;
the approximately $200 paid experiment stays deferred.

### Queued native CAS(10,12) scattering-dimension preflight

The passing forty-root/dipole target releases a further local resource audit.
Preflight supervisor PID **956517**, attached `sh_1153462ce001bUWk5LqtsP2iS8`,
waits for the exact competing aug-TZ import worker **954948** to finish, then
leases/resumes waiting staged owner **916479** on **4–7**. If that owner advances
or has active children/containers, the guard rejects the lease and preserves the
original error. Source/record:
`prepared/cas12-scattering-preflight-source/cas12-scattering-preflight.py`,
SHA256 `472e27521d048a02e0d67c00f5bd56ea0bed460129ea346a3aec7f8a7a1a6daf`.
The tracked [preflight contract](scattering-preflight-contract.json) preserves
the inputs, resource limits and controls. Container limit is **16 GiB**, with
**24 GiB available RAM / 20 GiB free scratch** required.

The finite worker first regenerates native integral/CONGEN preparation for the
completed tight CAS(10,11) baseline, requiring the inferred contracted dimension
to reproduce native SCATCI's **27546**, raw CONGEN **344124**, uncontracted L2
**26136**, contracted continuum **1410**, in both Pi sectors. It then prepares
the same fixed-DZ/radius-18/l=4/deletion-1e-6 scattering model from the newly
qualified CAS(10,12) target. Each preflight makes a fresh disposable copy of
retained target inputs/binaries and checks the original target/CI data hashes.
The program filter executes only continuum-containing integral preparation and
the B1/B2 native CONGEN stages; the disposable main driver writes native SCATCI
inputs and exits before any Hamiltonian construction/diagonalization. QC,
target eigenpairs and DENPROP are reused from their passing records.

The contracted dimension is inferred independently from the native CONGEN L2
CSF range plus `sum(NUMTGT*NOTGT)` in the generated scattering input. The
completed CAS(10,11) native output confirms that formula at dimension 27546;
negative controls reject a missing L2 group and a changed target inventory.
Their raw reference files, results and source provenance are frozen with the
worker. Limits are `NDIMX=10000000`, `CDIMX=NODIMX=1000000`, independent of RAM.
The inferred dimension remains a preparation result, not a completed SCATCI
dimension measurement or a physical scattering result.

Finally the already audited installed-library utility queries the inferred
CAS(10,12) dimension on **2×2 / 4×4 / 4×8** process grids, oversubscribed for
these lightweight queries only. Binary SHA256 stays
`ef38d6762972549112458999cdc1be0c2f2cb4c4bb76951c8bdb4a8ae1907c2b`.
Overflow/undersized queries retain their nonzero exits; only valid queries emit
array floors. Outputs/resources/source/original exits go to
`diagnostics/cas12-scattering-dimension-preflight/`. This finite local audit
reduces the remaining external-host uncertainty without repeating the five-hour
density stage. Paid provisioning remains deferred.

### Stretched import verifier failure, recheck and basis continuation

Original supervisor PID **933694**, attached `sh_113fdbc74001wgM6TF4p34y6Xq`,
finishes with SSH/controller exit **one** after a successful engine/batch exit
**zero**. Its post-import verifier fails at the saved/native comparison in
`stretched-import-basis.py`: `atol=1e-10` rejects seven singlet-A1 values with
maximum difference **4.99994712e-10 Hartree**. The saved table uses `%16.9f`.
The original driver, execution record, exception and raw native results remain
at their original paths; expensive QC, target eigensolves and DENPROP are reused.

New finite successor PID **957385**, attached `sh_115500536001s9jfjjnzUAzJe9`,
uses the freed **8–11** slot. Source:
`prepared/stretched-import-recheck-basis-source/stretched-import-recheck-basis.py`,
SHA256 `e6253cf789948d20a76e36f7a0160a08c47d1680d3ba404e8e0b3094e24e689a`.
It first checks the completed original owner, unchanged source/seed hashes,
engine/batch passes and slot availability, then runs the previously validated
`verify-fresh-target-import.py`, SHA256
`5d6f6bbd24172e9a9098aef9ce79a6a6456c2e3464eab123ac62d7616f341766`.
The unique-final-CIDATA parser and exact-decimal half-unit **5e-10-Hartree**
format bound pass all saved/native pairs. Physical root/dipole gates remain
**1e-7 Hartree / 1e-5 a.u.** Recheck exit is **zero**, recorded separately in
`diagnostics/cas11-stretched-import-decimal-recheck/`.

All **64 native roots** pass against the covered space-80/160 stretched-DZ
seeds within **3.598e-10 / 3.836e-10 Hartree**; forty native ensemble roots match
current QC within **4.939e-11 Hartree**. Current-QC dipole error is
**1.691e-11 a.u.**, seed dipole errors **8.856e-11 / 1.367e-10 a.u.**, minimum
active-subspace overlaps **0.9999999999999990 / 0.9999999999999987**.
Complete import wall is **3254.85 seconds / 54.25 minutes**, kernel peak
**5657108480 bytes / 5.269 GiB**. QC totals **371.96 seconds**, eight SCATCI
sectors **360.39 seconds**, DENPROP **2519.30 seconds**.

The [stretched-import companion](qualification-evidence/stretched-import/README.md)
verifies **488 payloads**, **40 raw QC states**, both seeds' **96 spectra / 624
eigenpair evaluations**, both checkpoint-based subspace comparisons and
**16229 raw resource samples**. It preserves the original failure and checks
the refined coverage's archived input bytes. Public fetch is byte-for-byte and
repackaging byte-identical.
Archive:
`https://data.qscat.org/ukrmol-co-stretched-import-2026-10-07/stretched-import.tar.ce39b22258f5.gz`,
SHA256 `ce39b22258f5ae8551354e3105af2ea56d6fbafc00d55688f681767f5c60fcb8`.
The publication bundle and executable verifier persist at
`prepared/publication-stretched-import-20261007/`.

The successor now runs `co-r2500-cctz-sa11-qc-dz80-start-ci-space80`, followed
by the same-geometry aug-TZ QC entry. Both retain the ordinary **200-cycle CI**
cap and original ensemble/tolerances. Only an orbital-converged CI failure
permits a newly named space-160 retry. Containers use **32 GiB**, with **40 GiB
available RAM / 20 GiB free scratch** required. Launch headroom is
**107512430592 available RAM bytes / 157826818048 free disk bytes**. Passing
basis attempts still require coverage, independent import and competing-start
qualification. The other three CPU groups retain their existing owners.

### Completed stretched basis QC and gated coverage/import successor

Successor PID **957385**, attached `sh_115500536001s9jfjjnzUAzJe9`, finishes
with original SSH/controller exit **zero**. Both
`calibration-sa11-stretched-basis-qc` entries pass without repairs. At **R=2.5
bohr**, cc-pVTZ QC takes **874.04 seconds / 14.57 minutes / 0.345 GiB** and
aug-cc-pVTZ **1219.58 seconds / 20.33 minutes / 0.434 GiB**. Final orbital
gradients are **7.077e-8 / 5.791e-8**, below the 1e-7 gate; ordinary/fresh-CI,
MO/spin/Pi checks all pass. Sequential batch wall is **2094.73 seconds**.

DZ→TZ and TZ→aug-TZ changes, respectively: ground energy **−0.825431 /
−0.0248054 eV**, averaged objective **−0.865111 / −0.0517186 eV**, maximum
ranked-root shift **0.913800 / 0.100352 eV**, maximum ranked-excitation shift
**0.0883686 / 0.0755466 eV**, z dipole **−0.00374061 / −0.0110101 a.u.**
Minimum cross-basis active overlaps are **0.985481 / 0.998189**, calculated
with the AO cross-overlap at fixed nuclear coordinates. These are basis
diagnostics, not same-basis numerical agreement or identified states.

The [stretched-basis QC companion](qualification-evidence/stretched-basis-qc/README.md)
verifies **170 payloads**, **80 raw final QC states**, both cross-basis
comparisons and **10430 resource samples**. Archive:
`https://data.qscat.org/ukrmol-co-stretched-basis-qc-2026-10-07/stretched-basis-qc.tar.3350694250b4.gz`,
SHA256 `3350694250b4fbe4f480e21db93d93da723ee0b943dcdbaee87ee595a32e5841`.
Public fetch is byte-for-byte and repackaging byte-identical. Original runs,
parent checkpoint and source hashes are preserved; publication
and verifier persist at `prepared/publication-stretched-basis-qc-20261007/`.

Finite successor PID **959637**, attached `sh_115796107001U2Jw3GtBtzP4H4`,
owns **8–11**. Source:
`prepared/stretched-basis-coverage-imports-source/stretched-basis-coverage-imports.py`,
SHA256 `2834b523a9ca8e2aef6be7a40597d33a189215b9604c5dae0ffaf6788b4ec74d`.
It checks the completed predecessor/source hashes and unchanged seeds:

| Basis | QC checkpoint SHA256 |
|---|---|
| TZ | `7aeab21d97c9c72a8e4e3f213238585a31650acaa3bcb67402eed28a921e1560` |
| aug-TZ | `1a1305c1df8f11b8aece22a531ae24b6fbe505433ac41d0a983356aaab36ec48` |

Each original checkpoint receives **48** independent probes: both spins,
four irreps, five/eight roots, trial spaces 40/80/160 at **600 cycles**. The
residual/analytic-control functions are AST-identical to the previously
validated controls; helper SHA256
`095e81ca0b57ded9a5eb6c6f082860864fdd641bd72eda247cd696ea6e2f9cfd`.
The Hubbard analytic eigenstate and perturbed-vector residual controls pass.
Every root needs convergence/spin/orthogonality, physical and spin-penalized
residuals ≤1e-9 Hartree, ensemble agreement and same-root-space/common-first-five
agreement. A failed scan preserves its original nonzero diagnostic exit and
only blocks that target's import. Ordinary target CI caps stay 200 cycles.

The [native import recipe](calibration-sa11-stretched-basis-imports.json),
SHA256 `8340337245b42ee962c230cc2568555e4371f595bff876c32a4feb361acc130d`,
requests eight roots per sector, forty-component averaging and forty retained
states with tight SLEPc settings. Passing scans release individual 16-GiB
imports, followed by all-64-root/dipole/subspace checks and the existing
unique-final-CIDATA/exact-decimal parser. Coverage uses an 8-GiB container /
16 GiB available-RAM floor; native imports require **24 GiB available RAM**;
all stages require **20 GiB free scratch** and an idle exclusive slot.
Launch headroom is **114908250112 available RAM bytes / 157269692416 free disk
bytes**. Diagnostic path: `diagnostics/cas11-stretched-basis-ci-residual-coverage/`.
Competing-start and across-geometry model qualification remain separate gates.

### Completed equilibrium downward-projection qualification

PID **954343**, attached `sh_1150eecdc0016n7grzLOHvIeDn`, completes original
SSH/controller exit **zero**. All four new-name core+active TZ starts and all
six spectrum/objective/dipole/subspace comparisons pass. Source aug-TZ matrices
are **92×92**, destination TZ matrices **60×60**; the thirteen-column
core+active projection repair resolves the original interface failure.

| aug-TZ seed lineage | CI space | QC seconds | Kernel peak GiB |
|---|---:|---:|---:|
| Repaired fresh RHF | 80 | 2085.26 | 0.352 |
| Repaired fresh RHF | 160 | 2049.96 | 0.409 |
| Upward fresh TZ | 80 | 1386.34 | 0.353 |
| Upward fresh TZ | 160 | 1326.20 | 0.409 |

Sequential batch wall is **6849.90 seconds / 1.903 hours**. Across the two
space-pair and four fresh-TZ comparisons, maximum root difference is
**4.156e-8 Hartree**, dipole difference **6.916e-8 a.u.**, objective difference
**8.527e-14 Hartree**, minimum active overlap **0.9999999999995466**. Both
aug-TZ lineages recover the existing qualified fresh-TZ solution; their two
aug-TZ solutions remain distinct and neither is declared a unique optimum.

The [downward-projection companion](qualification-evidence/downward-projections/README.md)
verifies **296 payloads**, **160 raw final QC states**, all six checkpoint-based
restart comparisons and **34132 resource samples**, with byte-for-byte public
fetch and byte-identical repackaging. Original interface failures
remain in the earlier companion. Archive:
`https://data.qscat.org/ukrmol-co-downward-projections-2026-10-07/downward-projections.tar.ed05f231aa41.gz`,
SHA256 `ed05f231aa412a22c47261a39a2866eb98e3f9c19d89876859c1a698fcd4f8d1`.
Publication/verifier persist at `prepared/publication-downward-projections-20261007/`.

The downward supervisor resumes staged owner PID **916479**. It advances on
**4–7** to `co-eq-augtz-sa11-qc-ci-space80` before competing-import PID
**954948** can lease the slot, so that import now follows the existing exact
staged-owner completion contract. Preflight PID **956517** continues to wait
for PID 954948's completion, independently of its scientific verdict; preflight's
own target gates remain explicit. No active calculation is retired for this
scheduling transition. CPUs **0–3** retain staged CAS(10,12) work and **12–15**
now run the queued l=5 continuum control.

### Stretched import verdicts, released CPU slot and neutral pilot

PID **959637** completes with original controller exit **one** and step
exits **[0,0,1,0,0]**. Both seed-orbital coverage scans pass all **96 probes /
624 eigenpair evaluations**. Both native engine/batches pass; aug-TZ's
independent all-64-root verifier passes. TZ's original verifier rejects
**triplet A1, root eight**, whose difference from the covered seed is
**1.56396e-7 Hartree** against the unchanged **1e-7** gate. Forty required
roots and the current-QC dipole pass. This failure occurs at the physical
root comparison, not the saved-table decimal-format gate.

Finite successor PID **965433**, attached `sh_115db4f6f001Vy3QiLPA5amPuc`,
owns **8–11** after the original controller completes and its process exits.
Source: `prepared/import-checkpoint-neutral-source/import-checkpoint-neutral.py`,
SHA256 `ef5306180287a187ded272b219c0dc0e06dd570feac744aa86118650ca9084c6`.
It validates the predecessor's source hashes and preserves the original exit
one before using the released slot. The imported-orbital diagnostic requests
eight roots in all eight sectors at spaces **80/160**, **600 cycles**, with
the unchanged validated residual helper. It checks native spectra against
the actual import checkpoint, separately records checkpoint-minus-seed
shifts, and preserves original seed-gate verdicts. Diagnostic:
`diagnostics/cas11-stretched-import-checkpoint-ci/`, **8-GiB** container,
**16 GiB available RAM / 20 GiB free scratch** floors.

The imported-orbital diagnostic now passes all **32 probes / 256 eigenpair
evaluations**. All 64 native roots agree with independent CI at the actual
TZ/aug-TZ checkpoints within **4.984e-11 / 5.009e-11 Hartree**. The seed-to-import
orbital change therefore explains the original TZ discrepancy; its strict
seed-gate rejection is retained. Aug-TZ matches its seed within **4.755e-8
Hartree**, passing the original gate. Native-import costs are **2986.70 /
3121.02 seconds**, with **5.265 / 5.261 GiB** peaks. Seed coverage takes
**524.56 seconds / 0.377 GiB**; the imported-orbital diagnostic takes
**307.42 seconds / 0.448 GiB**.

The [stretched-basis import companion](qualification-evidence/stretched-basis-imports/README.md)
reconstructs all **128 raw probes / 880 eigenpair evaluations**, native
spectra/serialization/dipoles, original exits, **478 retained source hashes**
and **34607 resource samples**. Its **1296 payloads** pass byte-for-byte public
fetch and byte-identical repackaging. Archive:
`https://data.qscat.org/ukrmol-co-stretched-basis-imports-2026-10-07/stretched-basis-imports.tar.44ddb0fd8d23.gz`,
SHA256 `44ddb0fd8d23776b0715eb87da416c5b19dc67349110e91933051da140d02e52`.
Publication, frozen evidence and executable reconstruction remain under
`prepared/publication-stretched-basis-imports-20261007/`.

The same finite supervisor then executes the independent
[near-equilibrium neutral recipe](calibration-neutral-near-equilibrium.json),
regardless of that scattering diagnostic's scientific verdict. It reuses
the original QZ/5Z points at **1.9/2.1323/2.5 bohr** and adds **23 points per
basis** on a **0.025-bohr** grid: **26 geometries per basis / 46 new jobs**.
The [contract](neutral-near-equilibrium-contract.json) specifies frozen-core
spherical RHF/CCSD(T), CCSD lambda-density dipoles, stable references and
per-basis relative-energy zero. Predeclared held-point comparisons use cubic
not-a-knot and PCHIP without extrapolation, with **1 meV / 0.001 a.u.**
diagnostic budgets; QZ→5Z relative-energy budget is **20 meV**. All points
retain convergence/density/amplitude checks. This is a labelled pilot over
the already stable reference domain, with correlation treatment still open.
One sequential four-core worker uses **24-GiB containers**, **20000 MB**
PySCF memory, **32 GiB available RAM / 40 GiB free scratch** floors.
Batch/analysis/anchor hashes remain under the neutral root in
`neutral-pilot-inputs/`; each fresh run remains under `runs/`.
Supervisor **965433** subsequently completes with original controller exit zero
and `finished_unix=1791380172.6187947`. All **46 new jobs**, imported-checkpoint
CI and the held-point/basis analysis exit zero. Independent public reconstruction
verifies all **52 points**, **959 payloads** and **63867 profiles**. Cubic/PCHIP
energy errors are **0.11386/0.95428 meV** at QZ and **0.11330/0.91453 meV** at
5Z; maximum held dipole error is **7.21e-6 a.u.**. All midpoint budgets pass;
QZ→5Z relative-energy maximum **17.2433 meV** passes the 20-meV pilot budget.
New-job wall sum is **11557.49 seconds**, kernel container peak **24 GiB** and
sampled anonymous peak **17.138 GiB**. The
[public companion](qualification-evidence/neutral-near-equilibrium/README.md)
retains raw logs, checkpoints/source provenance, byte-exact anchors and the four
stretched RHF rejections. Archive SHA256:
`d8a394e9ad94d1085218b9fc874539d6e03f194c832f81df06c48c912e33f411`.
The near-equilibrium basis/interpolation pilot passes; the correlation treatment
and stretched/dissociation domain remain open. CPUs **8–11** are released.

Staged owner PID **916479** also completes its three batches with exits
**[1,1,0]**: both equilibrium projected aug-TZ 80/160 retries still reject
singlet B1/B2 CI/spin flags; both R=1.9 TZ/aug-TZ staged starts reject
singlet A2 CI/spin flags despite orbital convergence; both R=2.5 basis
starts pass. Original equilibrium objectives agreeing to roundoff does not
override the failed flags. The separately completed fresh/downward targets
remain distinct, accepted lineages. PID **954948** subsequently completes its
covered competing aug-TZ native import on **4–7** with controller/verifier exits zero.
CAS(10,12)'s projected TZ attempt preserves its **29192.74-second** failure:
150 macroiterations end at gradient **1.05004e-5**, above **1e-7**, with
singlet-A2 CI rejected. Its aug-TZ attempt subsequently rejects triplet A2 CI
despite orbital convergence in **8676.26 seconds / 0.824 GiB**. The full staged
CAS(10,12) batch completes with exits **[1,1,1,1,0]**: compressed DZ rejects
singlet A1 CI in **3551.82 seconds**; stretched DZ rejects singlet A2 and triplet
B1/B2 CI in **7795.03 seconds**; projected TZ rejects as above; projected aug-TZ
rejects as above; fifty-component equilibrium DZ QC passes in **2450.85 seconds /
0.878 GiB**. The [staged-target companion](qualification-evidence/staged-targets/README.md)
now reconstructs all four original rejections, both CAS11/CAS12 fifty-component QC
pilots and their tight forty-component equilibrium references: **545 payloads /
281438 profiles**. Common-forty objective changes at the new orbitals are
**+0.170147 / +0.004204 eV** for CAS11/CAS12; reported ground-energy changes
are **+0.456512 / +0.007930 eV**. These are QC diagnostics, with independent
lowest-root coverage and native fifty-component import still unset. Archive SHA256:
`f06fbe3d0d71400f684905afd74a967e6c018b1e279008efb5bb206844f7064f`.
Eligible repair supervisor PID **918331** completes
at `finished_unix=1791385782.7340915`. The compressed space-80/160 pair still
rejects QC; both stretched DZ repairs pass QC but fail their independent
ensemble-coverage gates. Neither DZ pair releases a larger-basis successor.
The [original repair companion](qualification-evidence/cas12-base-repairs/README.md)
now reconstructs **437 payloads / 28844 profiles** and all **96** original
200-cycle coverage probes. Compressed QC retains a singlet-A1 CI rejection
despite orbital convergence. Stretched QC walls/peaks are **1427.43/1350.82 s**
and **0.786/0.968 GiB**; the pair agrees in ground energy within **1.973e-11
Hartree**. Original coverage fifth roots fail only at space 40 (singlet A2,
triplet B1/B2), and converge at larger spaces. Extra eighth-root flags also
reject. Archive SHA256:
`d4abeb85a47f27b0a6b16d9491fb4bdb35453d58db20086b84624d9052d51e3b`.
Its unresolved projected-TZ failure remains excluded from automatic retry.

The separately owned fixed-orbital residual successor PID **1017316** uses the
already validated **600-cycle** scan on physical CPUs **8–11**, released by the
completed neutral owner. Orbitals, active/frozen spaces, tolerance, spin penalty,
five/eight requested roots and spaces **40/80/160** remain unchanged. All **96
probes** execute before any qualification decision; physical/penalized residuals,
spin and root-coverage gates remain mandatory. Its source/owner and diagnostics
persist under `prepared/cas12-stretched-repair-600-source-v2/` and
`diagnostics/cas12-stretched-dz-repair-ci-residual-coverage-600-v2/`.
The first wrapper PID **1016190** rejects before any scientific probe because
`hashlib` was absent from its extracted driver imports. Original source, logs,
resources and exit one are retained separately and copied into the fresh
successor's inputs. No larger-basis/import successor is automatically released.

That successor now completes with original exit zero at
`finished_unix=1791398960.561587`, in **2240.11 seconds / 0.6252 GiB**. Both
unchanged CAS12 stretched checkpoints pass all **96** probes. Independent
reconstruction gives maximum physical/penalized residuals **7.638e-10 /
9.997e-10 Hartree**, spin-squared error **8.216e-15**, ensemble energy error
**2.701e-13 Hartree**, and trial-space/common-root-count differences below
**1.706e-13 Hartree**. The [residual companion](qualification-evidence/cas12-residual/README.md)
verifies **627 payloads / 11138 scan profiles** after public fetch, including
original failed 200-cycle probes and the failed wrapper. Archive SHA256:
`94b6757579a30ad63e0c86a806360a890abf08dd1f3d1fab9cc91e21a73b3447`.

Electronic successor PID **1022822** explicitly consumes that independently
reconstructed coverage gate and owns physical **8–11**. A finite fifty-component
equilibrium residual scan runs first on the existing CAS11/CAS12 checkpoints:
**10/13 A1 roots**, **5/8 roots** in the other sectors, both spins, spaces
**40/80/160**, **600 cycles**, all **96** probes and unchanged residual/spin/
orthogonality gates. A scientific coverage rejection retains its unset native
fifty-component-import verdict; it does not invalidate independently qualified
CAS12 stretched seeds.

The separately eligible [CAS12 stretched basis QC pair](calibration-sa12-stretched-basis-qc.json)
then projects the exact space-80 checkpoint
`7a5dd2f91ea3bf96f1548d5fd49c4730b53deb35e979f6867c84b412948b82da`
to TZ/aug-TZ, preserving CAS(10,12), the forty-component ensemble, tight
tolerances and the 150-macroiteration/200-CI-cycle QC contract. Each runs on
8–11 in a 32-GiB container with a 48-GiB available-RAM floor. Original numerical
rejections persist; no automatic larger-basis target import follows.

The [DZ 64-root import](calibration-sa12-stretched-import.json) also consumes
that exact checkpoint, using dense native target storage with selected-root
SLEPc, 16-GiB SCATCI workspace and a **64-GiB container / 80-GiB available-RAM**
floor. Its independent verifier retains all 64 roots, forty-state dipole and
core/active-subspace gates; the previously validated CAS11 verifier changes only
the twelve-orbital active-subspace endpoint. Source/proofs/checkpoint hashes,
explicit commands, resource floors and original exits persist under
`prepared/electronic-successors-source/`. No scattering successor is automatic.

The completed l=3 scattering control independently reanalyzes against l=4:
**+3.9602 meV** position, **+0.7434%** full width, **0.0166575 rad** maximum
phase difference modulo pi over the common 99-point grid. Both Pi sectors
pass. Independent fixed-window fits at **1.8–3.3 / 2.0–3.15 requested eV**
and backgrounds 1–4 stay within **3.982 meV / 0.7214%** between angular
cutoffs. Control cost is **10621.10 seconds / 2.950 hours / 24.329 GiB**.
The executable reanalysis and input-digest report persist under
`prepared/continuum-l3-reanalysis-20261007/`.
The l=5 calculation and queued native background replays subsequently finish
with original controller exit zero.

### Completed continuum/import evidence and measured scattering-memory blocker

Both B1/B2 l=5 comparisons against l=4 pass: **+0.07987 meV** position,
**+0.2655%** full width and **0.0162174 rad** phase modulo pi. Cost is
**14353.38 seconds / 26.300 GiB**. Common-window independent fits stay within
**0.3930 meV / 0.2279%**. All **36** native angular/background/detection replays
reconstruct. Background-width spans remain **12.60% / 12.65% / 11.98%** at
l=3/4/5, exceeding the unchanged 5% extraction gate.

The competing aug-TZ import independently passes all 64 roots within
**3.8811e-8 Hartree**, with native/seed dipole error **1.8756e-7 a.u.**;
required forty-root/current-QC error is **4.9431e-10 Hartree**. Cost is
**3259.47 seconds / 5.260 GiB**. Both orbital branches are retained.
The [continuum/import companion](qualification-evidence/continuum-competing/README.md)
verifies **1014 payloads**, **39 source hashes**, **219998 profile samples**,
byte-identical repackaging and public fetch. Archive SHA256:
`70c739dfebc709128872a722deed2c238bedfc941f715b76666814fa3ce37c5f`.

Preflight PID **956517** retains the original pre-native scheduling failure:
`ValueError('Staged owner advanced; preflight slot requires a new explicit transition')`.
Fresh released-slot preparation PID **982643** rejects a copied `moints` file;
PID **983294** generates new integrals but exposes the CONGEN workspace floor;
PID **983855** exposes missing retained-sector initialization from skipped target
stages; PID **985572** restores native CIDATA order but retains the independent
`NBMX` workspace overflow. The final new-name `cas12-scattering-preflight-nbmx-source`
successor completes with exit zero, preserving every predecessor's bytes/exits.
It verifies completed owner records and target/checkpoint/CI hashes before launch,
uses `LNDO=NBMX=100000000`, and initializes/checks forty-state inventory without
executing DENPROP. All preparations use the original 16-GiB container and
24-GiB RAM / 20-GiB scratch floors.

The native CAS(10,11) control reproduces **344124 / 27546** raw/contracted
configurations, including **26136 L² + 1410 continuum** configurations.
CAS(10,12) measures **990990 / 86352**, including **84942 L² + 1410 continuum**.
Preparation costs are **46.92 / 172.99 seconds**, at **2.984 / 12.730 GiB**.
One dense CAS(10,12) matrix needs **55.556 GiB**; valid installed-library grids
4×4/4×8 need **222.391 / 222.555 GiB** aggregate array floors. The 2×2 query
retains exit five and no estimate. Actual scattering runtime/peak/scaling remain
unmeasured. This proves the audited dense full-spectrum solve exceeds Sadaharu RAM.
The [preflight companion](qualification-evidence/scattering-preflight/README.md)
verifies **814 payloads**, **245 source hashes**, **1099 profile samples**,
byte-identical repackaging and public fetch. Archive SHA256:
`10e95b5389901990ad02824e75aa2f5fdfb68065ab8cd2b253c03ff6631c6470`.
Publication/reconstruction bundles persist in the two `prepared/publication-*-20261007/`
directories named by those companions.

### Finite two-anchor CAS(10,11) follow-up

Successor PID **988205** explicitly consumes completed continuum owner PID
**950692** on CPUs **12–15**, using image
`sha256:6b3b0ededa85494076229b111efe969985407023e404bff888f17cafb4630666`.
Its frozen driver SHA256 is
`136bc3f31796328ba44f0e25c8c320a87a0d9320011deb3b0d1b917440402e2e`.
Source/CLI/hash checks pass before deployment. Execution persists at
`prepared/cas11-anchor-followup-source/execution.json`.

Compressed R=1.9-bohr DZ fixed-orbital coverage uses five/eight roots, spaces
40/80/160 and the validated 600-cycle contract, with unchanged residual/spin
gates. That scan passes and releases the
[compressed 64-root import](calibration-sa11-compressed-import.json).
Only a passing all-root/dipole verifier releases its scattering entry.
The separately qualified R=2.5 DZ 64-root checkpoint has exact SHA256
`48e45536a2a795136db6d19706100ea56b6fba71980343064ea8d8a5f65874d7`;
its scattering eligibility is independent of a compressed import failure.
The [two-entry scattering recipe](calibration-sa11-anchor-scattering.json)
uses the same l=4/radius-18/deletion-1e-6/five-root-per-sector model as equilibrium,
on **0.1–8.0 requested eV / 0.05-eV spacing** at R=1.9 and
**0.01–1.0 requested eV / 0.005-eV spacing** at R=2.5. Each full pipeline uses
a 48-GiB container with 64-GiB available-RAM and 20-GiB scratch floors.
This is a geometry/extraction pilot: empty automatic candidates remain unset,
and bound-state/pole/state-continuity qualification still precedes a full sweep.

The [compressed-anchor companion](qualification-evidence/compressed-anchor/README.md)
now independently reconstructs **1115 payloads / 75913 run-profile samples**.
All **48** fixed-orbital probes pass; physical/penalized residual maxima are
**7.034e-10 / 9.716e-10 Hartree**. All **64** native roots agree with the
independent CI coverage oracle within **9.050e-9 Hartree**. Import cost is
**3035.78 s / 5.262 GiB**. Compressed scattering completes both 27546-dimensional
Pi sectors in **11439.93 s / 25.277 GiB**, with **1e-9 rad** Pi phase splitting
and fitted candidate **3.826384 / 2.322242 eV** position/full width. Raw outputs,
source/checkpoint provenance and a historical parent-owner snapshot travel with
the public archive, SHA256
`eafc854bef6069b3d4fc2b804ed1317395a740464cf5b40ced641abcae78e974`.
The separately owned stretched calculation remains active; this compressed
pipeline pass does not qualify electronic selection, extraction or continuity.

### Sparse/iterative scattering qualification — preferred before a large host

**CAS(10,10) 2048-root implementation qualifies; finite CAS(10,11) qualification launched.** Investigate the pinned engine's native
sparse/iterative contracted-scattering route as the preferred response to the
CAS(10,12) dense-workspace blocker. Paid provisioning remains deferred. The
222.4-GiB array floor describes the audited dense ScaLAPACK path, not every possible
CAS(10,12) solver. The successful target SLEPc runs used dense Hamiltonian storage
and selected eigenpairs; they are not a completed sparse-scattering demonstration.

The pinned archive now resolves the native dispatcher, sparse constructor,
selected-state export and downstream SWINTERF inputs. The
[executable replay guide](SPARSE_SCATTERING.md) and
[finite control contract](sparse-scattering-contract.json) describe the first
CAS(10,10) sequence, **128/512/2048** scattering eigenpairs at fixed 40-state
target/model inputs. The first named-input owner PID **1006374** completes its
B1 solve but rejects the downstream SWINTERF partitioned path. Its frozen source
and execution record persist under
`prepared/sparse-scattering-control-named-input-source/`. The
[first public companion](qualification-evidence/sparse-native-control/README.md)
verifies **212 payloads / 1464 profiles** and dense root/residual/continuum checks
at **292.93 seconds / 1.228 GiB**; the original outer failure remains exact.
Archive SHA256:
`d27863143d93817f2ca740a338601f109d8fe3c79a22f5f1ac571f81622067b2`.

MPI-boundary owner PID **1006959** then passes all four B1 native stages at
128 poles; its analysis rejects a zero-error native rank-summary block interleaved
with fitted values. The guarded parser fix removes only an explicitly empty
summary. Original logs/source/exits remain in
`prepared/sparse-scattering-control-mpi-boundary-source/` and
`diagnostics/cas10-native-sparse-roots128-mpi-boundary-20261007/`.
Dense static channel/threshold/multipole records match byte-exactly, boundary
energies within **5.97e-13 Hartree**, and amplitudes within **3.46e-8**. The
complete phase grid differs by **1.56 rad modulo π**, failing the gate.
The route has no partitioned correction: omitted eigenstates contribute zero
to its truncated spectral pole sum.

Successor PID **1008121** uses physical CPUs **4–7**, four ranks and 16-GiB
containers. Its frozen source/owner record is
`prepared/sparse-scattering-control-boundary-differential-source/`.
Its finite queue contains the 128-root **B2** control, then **512/2048** controls
in both sectors, retaining channels/thresholds/multipoles and boundary-differential
checks. A scientific phase rejection permits the next declared refinement;
native/inner-region failures stop the queue. No CAS(10,11)/(10,12) successor
is released by these incomplete small-model controls.

Independent dense-spectrum omission owner PID **1011041** uses released physical
CPUs **0–3**, a 2-GiB container and the original dense boundary/channel data.
Its finite **128/512/2048/4096/6144/8192/8350** sequence runs both Pi sectors
through unchanged RSOLVE/EIGENP grids and independent **1.8–3.3 requested-eV**
fits at background orders 1–4. The full-spectrum 8350 endpoint must reproduce
the dense boundary bytes and phase grid. This directly measures the approximation
required of selected-root scattering, independently of iterative-solver residuals.
Source, pinned native energy-conversion files and owner records persist under
`prepared/dense-spectrum-omission-source/`; diagnostic output under
`diagnostics/cas10-dense-spectrum-omission-20261007/`.

The dense-omission owner subsequently completes with original exit zero in
**1898.71 seconds / 0.374 GiB**. Both Pi sectors give the same phase errors:
**1.56579 / 1.46946 rad** at 128/512 poles, **2.0e-7 rad** at 2048 poles and
**1.0e-8 rad** at 4096 poles; higher counts reproduce the full printed phase
grid. All four identical fixed-window fits reject at 128/512 and pass at
2048/4096/6144/8192/8350. At 2048 poles, position/full-width differentials are
below **1e-6 eV**. The full-count endpoint reproduces the original dense
boundary bytes and phase grid. The omitted boundary-amplitude squared norm
falls from **193.909** at 512 poles to **5.7385e-4** at 2048, of a full-spectrum
norm **234.074**. This is a participation diagnostic, not an observable bound;
the complete propagated phases/fits supply the numerical gate.

Native owner **1008121** completes with original controller exit zero at
`finished_unix=1791397854.027752`. Its **2048-root** control passes both sectors:
maximum energy error **9.664e-13 Hartree**, residual **1.336e-10 Hartree**,
boundary-amplitude error **1.003e-6**, phase error **2.0e-7 rad**, and identical
fixed-window position/full-width differences **5.847e-8 / 2.381e-7 eV**. Static
channels/thresholds/multipoles match byte-exactly. Native phases match the
corresponding dense 2048-pole truncation at printed precision. Wall/peak are
**1340.74 seconds / 1.861 GiB** for both sectors; inner/export stages take
**558.81 / 531.48 seconds**. Both use 2548 Krylov vectors and 57 iterations.
This qualifies the selected-spectrum implementation for this fixed CAS10 model;
it does not qualify electronic selection or a larger scattering model.

The [public boundary companion](qualification-evidence/sparse-boundary/README.md)
independently reconstructs **306 payloads / 25504 profiles**, all six native
sector controls and all fourteen dense omissions. The complete dense native CI
and boundary oracles, raw telemetry, original parser failure, source/owner
snapshots and additive capture supplements travel with the archive. Public
fetch/reconstruction and byte-identical repackaging pass. Archive SHA256:
`8f31e6e1af3da152b5845175efcc7b4af8066c209d80e754e2af7b0848e68f06`.

CAS11 successor PID **1020516** explicitly consumes those completed gates and
the released physical **0–3 / 4–7** slots. Its
[predeclared contract](sparse-scattering-cas11-contract.json) uses the qualified
`co-eq-ccdz-sa11-tight-cc40` dense oracle, dimension **27546**, unchanged forty
neutral states and both Pi sectors. A **2048/4096/8192/16384/24576/27546**
dense-omission sequence runs first in a 2-GiB container on 0–3. If at least one
declared native cutoff passes and remains stable at all higher dense cutoffs,
the **2048/4096/8192/16384** native sequence runs on 4–7 with four ranks,
**24-GiB** containers, **48-GiB available-RAM / 30-GiB free-disk** floors and
the unchanged native solver/export/gates. Otherwise it records a specific
spectral-coverage blocker without launching the native sequence.

Fresh physical/SMT ownership and resource checks precede every container.
Frozen source, original commands/exits, CAS10 independent reconstruction,
reference input hashes and owner record persist under
`prepared/cas11-sparse-qualification-source/`. Native controls receive identical
fixed-window/background comparisons after each solve. Scientific phase/fit
rejections permit the next declared refinement; engine/eigenpair/export failures
halt the queue. No CAS12 scattering launch is automatic.

The CAS11 dense-omission stage subsequently completes in **1456.15 seconds /
0.6691 GiB**, with original exit zero. Both sectors pass all six declared counts.
At 2048 poles, maximum phase error is **5.1e-6 rad**, and fixed-window
position/full-width changes are below **7.669e-7 / 1.291e-5 eV**. Higher counts
remain stable; the complete endpoint reproduces the dense boundary bytes and
printed phase grid. The predeclared native **2048-root** control is active on
4–7; this dense-oracle finding does not substitute for its native differential.

Original scheduler/build-environment attempts and the MPI-stdin native abort
remain in their original directories. The successor passes the named input file
to every MPI rank. Read-only `EPSSolve` telemetry will establish actual matrix
storage, nonzero/allocation counts, selection and residuals; source inspection
alone does not establish those measurements.

- [x] **Audit the native path and downstream contract.** Resolve the pinned
  contracted-Hamiltonian sparse initialization, PETSc storage and SLEPc dispatcher;
  record actual nonzero counts, allocation, eigenpair-selection policy, `nstat`
  semantics and boundary-amplitude/export behavior. Sparse preallocation guesses
  are not measured sparsity. Preserve matching source/library/image digests.
- [x] **Establish a small-model differential control.** Reuse a qualified
  CAS(10,10) dense reference with identical target/orbital/continuum inputs.
  Check eigenpair residuals, symmetry, boundary amplitudes and downstream execution;
  reject silent dense fallback or incomplete eigensolver output.
- [ ] **Qualify the CAS(10,11) scattering approximation.** Predeclare a finite
  increasing-eigenpair sequence and energy coverage after the path audit. Compare
  B1/B2 full common-grid phases and fixed-window position/full width against the
  completed 27546-dimensional dense reference, with identical extraction settings.
  Require both reference agreement and stability on further spectral refinement;
  diagnose omitted-state/background contributions rather than fitting away the error.
- [ ] **Measure the complete resource envelope.** Record wall/CPU time, aggregate
  container peaks, matrix nonzeros, eigenvector/workspace storage and persistent
  scratch for every passing and failed control. Declare owner/CPU/RAM/disk gates
  before launching; sparse storage alone does not prove the full calculation fits.
- [ ] **Gate a local CAS(10,12) pilot on those results.** Keep its qualified
  forty-state target and physical model fixed. Attempt only when the measured
  resource envelope fits Sadaharu with headroom; retain failure/blocker evidence
  if spectral convergence or memory is inadequate.
- [ ] **Publish a reproducible low-memory recipe and verdict.** Include executable
  configurations, root-count/refinement tables, source pins, raw phases, boundary
  data, profiles and an independent verifier. State the validated energy/model
  domain and distinguish numerical scattering convergence from electronic-model
  convergence and the existing background-width sensitivity.

The [qualification contract](../../docs/physics/co-electronic-qualification.md#sparseiterative-scattering-qualification-contract)
defines unchanged observable gates and the evidence needed by reproducers. A
larger host becomes the fallback only after this local route has a documented
feasibility or convergence verdict.

Repository verification for this handoff passes **90 tests / 1 skipped**, project
Ruff, all three new CAS12 recipe entries against the live runner's argument
parser, **205 relative documentation links**, five public archive hashes/fetches/
portable reconstructions, and `git diff --check`. The main-only search index is
complete and upstream-current at
`0f9768e4e89f2accb6f3ff5c2d3dc921719325eb`; relevant branch-only source and
qualification additions are resolved against current files and their local
delta. Frozen producer snapshots remain at their original executable hashes.

### Scientific decisions still pending

1. **Qualify the larger-active-space targets.** The tightened CAS(10,11)
   import passes; retain its completed start and root-count refinements and
   finish the running equilibrium scattering and matched QC controls.
    CAS(10,12)'s separately gated selected-root configuration now passes its full
    forty-root/dipole import locally; the installed-library audit still rules out
    current all-spectrum dense target layouts. The measured 86352-dimensional
    dense scattering path requires array floors above 222 GiB. Qualify the native
    sparse/iterative route above before considering a larger-host fallback, and compare
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
