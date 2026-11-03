# Running and reproducing the CO experiment

## Boundaries

These are diagnostic execution interfaces, not a production electronic-model
prescription. Use new run and manifest names: the runner refuses to overwrite
old calculations. Historical frozen contracts must run with their archived
sources; current code is not a byte-identical replacement for an old capsule.

## Docker pilot

The [Docker guide](../../docker/ukrmol-plus/README.md) specifies engine pins,
precision variants, serial/MPI tests, persistent scratch and mount ownership.
From the repository root on the Docker host:

```sh
docker build --build-arg BUILD_JOBS=4 -f docker/ukrmol-plus/Dockerfile \
  -t qmodeling/ukrmol-co:source .
```

Choose a writable persistent `RUN_ROOT`, then run the default two-frozen,
20-virtual-orbital equilibrium SEP diagnostic:

```sh
IMAGE=qmodeling/ukrmol-co:source
docker run --rm --cpuset-cpus=0-3 --memory=16g --memory-swap=16g \
  --shm-size=1g -v "$RUN_ROOT:/work" \
  --entrypoint /opt/ukrmolp/entrypoint.sh "$IMAGE" bash -c \
  'python3 -m projects.ukrmol_co.run --workdir /work/runs/co-equilibrium-001 &&
   python3 -m projects.ukrmol_co.analyze /work/runs/co-equilibrium-001'
```

The six-electron SEP diagnostic adds `--frozen-orbitals 4
--virtual-orbitals 17 9 9 4 --deletion-threshold 1e-6`. Neither is a production
configuration. The sphere is 18 bohr; angular/deletion checks must be repeated
when the physical model changes. The analyzer converts native RESON Rydbergs
at its input boundary. B1/B2 are symmetry-resolved contributions, not the total
elastic cross section of polar CO.

## State-averaged target

```sh
docker run --rm --cpuset-cpus=0-3 --memory=16g --memory-swap=16g \
  --shm-size=1g -v "$RUN_ROOT:/work" \
  --entrypoint /opt/ukrmolp/entrypoint.sh "$IMAGE" bash -c \
  'python3 -m projects.ukrmol_co.run --workdir /work/runs/co-sa-target-001 \
     --model CAS-A --orbitals state-averaged --active-orbitals 4 3 3 0 \
     --virtual-orbitals 0 0 0 0 --target-only --ranks 4 &&
   python3 -m projects.ukrmol_co.analyze /work/runs/co-sa-target-001'
```

The default ensemble averages five equal-weight roots in each spin/irrep:
40 C2v components, not 40 distinct physical electronic states. Irrep order is
A1/B1/B2/A2. `--sa-singlet-roots` and `--sa-triplet-roots` change the orbital
objective. `--target-singlet-roots` and `--target-triplet-roots` change computed
root counts; they must cover the ensemble. Scattering retention is separate.
Increasing computed roots alone must not silently change the ensemble.

`--orbitals natural` means ground-state CASSCF, not the multi-state average.
`--qc-only` stops before native target construction. `--target-only` includes
independent UKRmol roots/dipoles. `target.json` records the backend, ensemble,
ground-root dipole, spin diagnostics, MO inventory and source hashes. The
historical QC log slot `target.psi4.out` can contain an explicitly identified
PySCF calculation. Pi degeneracy, every averaged root and the Molden digest
must pass; optimizer convergence alone is insufficient.

The current parser (`python3 -m projects.ukrmol_co.run --help`) is the authority
for flags/defaults. Important independent controls include:

- Orbital energy/gradient thresholds, CI tolerance/residual, cycles and memory.
- `--target-ci-max-space` for trial-space refinement without changing the
  ensemble; it must exceed all ensemble root counts. Inspect fresh-CI flags,
  physical/spin-penalized residuals and spins, not energy agreement alone.
- `--target-initial-checkpoint` projects core/active orbitals, preserves the
  source checkpoint and its digest, and requires matching active-space counts.
  Stage simultaneous geometry/basis changes through a passing same-geometry
  checkpoint; the pinned PySCF helper rejects changing both at once.
- `--target-diagonalizer slepc` requires the SLEPc-enabled image. Native target
  storage and scattering storage are different questions. `davidson-serial`
  retains rejected spectrum controls and is not an equivalent validated choice.
- `--congen-workspace` raises NDIMX, CDIMX and NODIMX input limits; free host RAM
  does not override those native workspaces.

Method details and measured limits live in
[state-averaged targets](../../docs/physics/co-state-averaged-target.md) and
[electronic qualification](../../docs/physics/co-electronic-qualification.md).

## Batch, replay and collection

```sh
python3 -m projects.ukrmol_co.batch projects/ukrmol_co/calibration-throughput.json \
  --root "$RUN_ROOT" --image qmodeling/ukrmol-co:source \
  --cpu-groups 0 1 2 3 --memory 4g
```

`batch.py` takes a finite manifest, snapshots source read-only, records hashes
and image identity, and uses disjoint CPU slots. Outputs live under
`batches/<manifest-stem>/` and `runs/<job-name>/`. Some archived manifests have
their own stricter owners/resources: this example does not override them.
The QC stage uses local threads before native MPI stages; record physical-core
and SMT ownership rather than assuming ranks describe all CPU use.

Inside the image, replay a native fit without new integrals or inner solves:

```sh
python3 -m projects.ukrmol_co.refit /work/runs/co-eq-f4-sep39 --nback 3
```

Collect copied lightweight evidence locally:

```sh
uv run python -m projects.ukrmol_co.collect /path/to/copied-artifacts \
  --pairs projects/ukrmol_co/calibration-pairs.json \
  --provenance projects/ukrmol_co/campaign-provenance.json \
  --output /path/to/new-results.json
```

Use historical provenance only to reconstruct that historical campaign. New
calculations need their own metadata. Grid comparisons use shared native energies,
not interpolation through an unresolved feature; fewer than two common points
is an error. Extraction windows/backgrounds and uncertainty checks remain part
of the scientific contract even when a native fit returns a value.

## Evidence and resource interpretation

Run evidence includes config/decks, instrumented upstream driver/library,
logs/generated native inputs and outputs, `stages.tsv`, sampled resources,
channels/amplitudes/K matrices and target checkpoints. TMPDIR and PSI_SCRATCH
are inside each persistent run's scratch directory; native integrals and fort
files may also live in the geometry directory. Keep them in total disk accounting.

Memory is sampled every 0.2 s and allocated disk every 2 s. Sampled peaks are
lower bounds and container memory includes page cache. Newer records include
the kernel lifetime `memory.peak`; sequential atoms within one container do not
have independent kernel peaks. Fresh containers isolate per-calculation peaks.
Exit statuses, timeout/OOM, source integrity and observed monitoring gaps must
remain separate from numerical and physical verdicts.

Fetch and independently reconstruct accepted and failed calculations through
the [evidence index](qualification-evidence/README.md). Extract historical
archives **outside the checkout**, since source capsules contain test modules.
Original engine licences/build pins and exact calculation source are retained.
The first four calibration batches have hashes/decks but no complete source
snapshot; today's runner cannot recreate their earlier workspace failures.
