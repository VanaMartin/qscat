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

## Next qualification gate

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

`sa-results.json` records completed batches; `sa-provenance.json` supplies
their host/root metadata. The first seven attempts include three successful
target-only import checks and four preserved setup/diagnostic failures. The
eight-component truncated average failed Pi degeneracy; the full 40-component
cc-pVDZ and aug-cc-pVDZ ensembles passed all-root UKRmol energy comparisons.
Their lowest triplet Pi excitations are 6.3952 and 6.3497 eV, respectively.
The augmented-basis ground root is 0.07324 Hartree higher, so comparable active
spaces and starting-point stability remain open checks. Early target-only
import probes put C at the origin; centered continuation jobs match the
upstream center of mass and sphere origin.
The completed target ladder adds three successful pipelines (CAS(10,8), TZ,
and tightened equilibrium DZ) and two failed sentinel degeneracy checks.
Their Pi splittings were 2.56e-7/1.35e-7 Hartree at R=1.9/2.5 bohr. The snapshot
currently contains 22 attempts: ten validated pipelines and 12 preserved
failures, including four diagnostic-callback setup failures.

Launched manifests:

- `calibration-sa-target-ladder.json`: completed CAS(10,8), TZ, R=1.9/2.5 and
  tighter optimization controls, using CPU slots 0–3 and 12–15.
- `calibration-sa-scattering.json`: completed centered DZ/aug-DZ CAS(10,10),
  40 channels, 99 energies, using slots 4–7 and 8–11. Twenty-four native
  RESON replays are saved in the two runs' `refits/` directories.

- `calibration-sa-projected.json`: completed projected DZ-checkpoint starts in
  aug-DZ/TZ, followed by tighter R=1.9/2.5 retries, on slots 0–3 and 12–15.
  The sentinel retries use 1e-11-Hartree orbital energy, 1e-7 orbital gradient
  and 1e-12 CI convergence tolerances, retaining the 1e-7-Hartree Pi gate.

The projected starts compare against the lowest-HF-orbital starts. Preserve
the source checkpoint copy and hash as part of each new run. The source
checkpoint remains in the full remote workspace; the lightweight copied
evidence does not include binary checkpoints.
Both projected QC stages and native target checks passed. Their ground energies
are −112.898457558 Hartree
(aug-DZ) and −112.925516422 Hartree (TZ). Initial/final active-space overlap
singular values range from 0.9806 to 1.0000 and 0.9957 to 1.0000, respectively.
For aug-DZ, the default start has the lower **ensemble** energy despite the
projected start's lower **ground** energy. This is an ensemble/active-space
tradeoff, not a proof that the projected result is globally preferable. For TZ,
the projected start lowers both ensemble and ground energies.

The tight sentinel retries did not converge after 100 macroiterations. Their
orbital gradients stalled near 6.8e-7/1.5e-7 with zero steps, above the requested
1e-7 threshold. `calibration-sa-hessian.json` and
`calibration-sa-projected-scattering.json` subsequently exposed a diagnostic
callback error during microiterations; all four failures remain preserved.
`calibration-sa-retries.json` runs the corrected callback on four disjoint
slots: tighter-Hessian target checks at R=1.9/2.5 and projected aug-DZ/TZ
scattering. It uses the same strict Pi and convergence gates.

Additional running manifests are `calibration-sa-channels.json` (DZ 50
computed/retained channels with the 40-component orbital ensemble fixed,
then projected TZ 50-component averaging) and `calibration-sa-active.json`
(DZ CAS(10,11)/(10,12) target-only checks). They use the released slots 4–7
and 0–3, respectively. The latter uses 32-GiB containers and larger CONGEN
workspaces. Projected scattering remains on slots 8–11 and 12–15.

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

### Remaining gates

1. **Qualify the common orbitals.** The `state-averaged` backend constructs
    multi-spin/multi-irrep CAS targets. The existing `natural` option optimizes
    only the ground singlet. Verify ensemble completeness, excited energies,
   dipoles, orbital/state identities and continuity as R changes.
2. **Check electronic convergence.** Compare DZ/TZ atomic bases, active/external
   spaces and about 40–50 C2v target components at R=1.9, 2.1323 and 2.5 bohr.
   Published values are plausibility checks with model-dependent spread, not
   error bars or a requirement to force agreement.
3. **Recheck continuum and extraction.** Start from radius 18 bohr, l=4, double
   precision and deletion 1e-6, then vary those controls and fit window/grid/
   background for the new electronic model. Use the provisional 0.05 eV,
   5% width (0.001 eV floor) and 0.05 rad modulo-π criteria from the README.
   Diagnose poles/bound states at threshold and track the same physical feature.
4. **Calculate a correlated neutral curve.** RHF pilot energies are not adequate
   nuclear-dynamics inputs. Keep the relative scattering resonance and the
   neutral curve's reference convention explicit.
5. **Remeasure cost, then cover geometry/energy.** Qualify the electronic model
   before the tentative 25–40-geometry, 300–800-energy campaign. Multi-day
   execution is acceptable. Use adaptive windows where fits narrow, and measure
   throughput/memory again at the selected model/channel size.
6. **Connect to potential fitting.** `projects/potential_factory/target.py`
   supplies neutral/resonance target holders; eigenphase fitting and the polar
   tail still need development. Fixed-nuclei data do not uniquely identify a
   local radial potential. HTTP and AWS execution are subsequent capabilities.
