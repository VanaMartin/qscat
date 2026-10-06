# CO experiment: resources and continuation

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
UKRmol target check is running on CPUs 8–11 in a 32-GiB container; its earlier
model cost was 127.80 minutes, not a guaranteed duration for this retry.
Wait for its final batch/analyzer records before launching the equilibrium
scattering pilot. The active SSH job has an attached background supervisor.

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
`qualification-results.json` contains eleven attempts in four completed batches:
three passing QC targets, seven passing neutral records and one failed neutral
runner. Its public archive retains 48 CI probes and the four completed diagnostic
directories. Verification checked 728 payload digests, 62 batch-source hashes,
17 image-source hashes, all ten successful records and both aggregates;
repackaging is byte-identical. The active target-import job is excluded.

The refreshed standalone image is `qmodeling/ukrmol-co:electronic-qualification`,
ID `sha256:b5d5f7fc4d638d19eaeb65399b5acf20f56f4fe4526c63761b8945a0b01c5d08`.
`build-electronic-qualification.log` and `electronic-qualification-image.json`
record the cached engine gate, runtime imports, pinned dependency versions and
all 17 embedded Python/Perl/Fortran source digests. Numerical jobs retain their
earlier image IDs and immutable source snapshots.

Repository handoff: analysis/packaging code is pinned at
`443b27fbc9bc021e812a64f8c49fececc56e9540`. The committed-main index is complete
and upstream-current at `0f9768e4e89f2accb6f3ff5c2d3dc921719325eb`; the branch
adds state-averaged/fresh-CI work plus the QC-only, correlated-neutral and
workspace qualification tools and their evidence. Source anchors and local
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

1. **Qualify the larger-active-space targets.** Finish the tightened CAS(10,11)
   import check, retaining the completed start and root-count refinements.
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
