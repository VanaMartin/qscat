# Native sparse/iterative scattering controls

This is an experimental lower-memory route for the same CO scattering model.
It is being tested locally before considering large-host allocation. The
[qualification contract](../../docs/physics/co-electronic-qualification.md#sparseiterative-scattering-qualification-contract)
requires a dense differential, spectral refinement and measured resources.
See [the live handoff](CONTINUATION.md#sparseiterative-scattering-qualification--preferred-before-a-large-host)
for execution status; launching a control is not a passing scattering verdict.

## What the pinned source establishes

The checksum-pinned UKRmol-in 3.3.0 archive, SHA256
`1b81e423b7e2a5a530c336cfc7f7ab50d6f2705abf771c71d6119fc065d9cac2`,
contains the relevant sources under `source/mpi-ci-diag/`:

- `Modules/Dispatcher/Dispatcher_module.F90`,
  `DispatchMatrixAndDiagonalizer`: explicit `igh=-1` chooses the iterative
  backend; a SLEPc-enabled build uses `SLEPCMatrix` and `SLEPCDiagonalizer`.
  Automatic dispatch selects a dense solver above 20% of the full dimension.
- `Modules/Hamiltonian/Contracted_Hamiltonian_module.f90`,
  `build_contracted_hamiltonian`: requests `MAT_SPARSE`. In the SLEPc matrix
  implementation this reaches `MatCreateSBAIJ`, with block size one. Its
  preallocation guess is not a sparsity or peak-memory measurement.
- `Modules/Options/Options_module.f90`, `compute_expansions`: `nstat=0`
  expands to the complete contracted dimension. Positive `nstat` requests a
  selected scattering spectrum, independently of the retained neutral states.
- `Modules/Diagonalizers/SLEPCDiagonalizer_module.F90`, `diagonalize_slepc`:
  selects Krylov–Schur and sets the eigenpair count and tolerance. The active
  code does not set `EPSWhich` or call `EPSSetFromOptions`; commented selection
  code and generic PETSc command-line flags do not establish an available
  spectral-selection interface. Measure the library default and match the
  emitted states to the dense reference. Largest magnitude is not generally
  equivalent to lowest energy.
- `Modules/Hamiltonian/SolutionHandler_module.F90`: writes the requested
  eigenpairs to native CIDATA. `vecstore=1` retains continuum coefficients.
  `vecstore=3` additionally keeps vectors in memory for the native MPI boundary
  exporter in `Modules/Utilities/Postprocessing_module.F90`.

The successful selected-root **target** uses dense storage. This control uses
contracted scattering, a separate matrix construction and a separate spectrum.
The upstream entrypoint, double-precision executable/library pairing and build
pins are described in [the image guide](../../docker/ukrmol-plus/README.md).

## Finite small-model experiment

[`sparse-scattering-contract.json`](sparse-scattering-contract.json) fixes
the CAS(10,10) equilibrium l=4 dense reference, 40 retained neutral states,
dimension **8350**, both Pi sectors and scattering eigenpair counts
**128 → 512 → 2048**. Eigensolver tolerance is `1e-12`, with 1000 iterations
and a six-hour per-stage timeout. The reference's full common energy grid and
its extraction inputs are retained. These counts are controls, not recommended
production settings or an asserted convergence sequence.

[`sparse_scattering.py`](sparse_scattering.py) copies the immutable integral,
CONGEN, target-CI and target-property inputs from a qualified dense run into a
new diagnostic directory. The `--boundary-export mpi` route runs SCATCI with
its in-memory boundary exporter, followed by RSOLVE, EIGENP and RESON. The
retained `swinterf` route reproduces the original downstream rejection.
Every MPI rank receives the named `scatci.inp` argument; supplying
only MPI stdin leaves other ranks without the namelists. The original native
abort can return exit zero, so solver completion and emitted data are checked
in addition to the subprocess status.

The runner matches selected total energies to the dense spectrum, checks
absolute residuals and sign-aligned continuum coefficients for isolated roots,
and compares the complete phase grids modulo π. Near-degenerate boundary
subspaces need a separate comparison before release. Native fitted candidates
are retained, but **the independent fixed-window position/width verdict remains
unset** until that analysis is performed. Empty candidates do not mean zero width.

### Measured first solve and the partitioned-interface blocker

The original B1-only 128-root control measures **292.93 seconds / 1.228 GiB**
through its downstream rejection. PETSc reports `mpisbaij`, 11746514 used
structural entries, 13243850 allocated entries, a 256-vector Krylov basis and
largest-magnitude selection. All 128 emitted states match dense roots 1–128:
maximum energy error **5.97e-13 Hartree**, absolute residual **1.34e-10 Hartree**,
and sign-aligned continuum-coefficient error **7.88e-8**. These are successful
inner-region numerical checks, not a phase/position/width verdict.

SWINTERF then prints `TRUEAMP: NOT IMPLEMENTED FOR PARTITIONED R-MATRIX AND
REDUCED SET OF EIGENVECTORS SAVED ON THE CI FILE` and stops with exit zero.
RSOLVE rejects the missing boundary data. The original control and its exits
remain preserved; the 512/2048 successors of that owner were never launched.

The pinned outer source confirms a deeper issue: `source/libouter/swinterf.f`,
`CIVIO`, reads the Hamiltonian diagonal into its local `dg` array but does not
return it through its argument list. `READCIP` subsequently adds the core energy
to its separate, uninitialized `dgem` array; `ORDERD` and `TRUEAMP` consume that
array for the partitioned correction. Merely writing full coefficients would
not qualify this path. The pinned `doc/swinterf` also labels the partitioned
approach unimplemented in UKRmol+.

The successor uses native MPI `outer_interface` with `vecstore=3`, the unchanged
40-state map, radius and target property file. Its `write_boundary_data` emits
only the selected pole energies/amplitudes, with `ibut=0` and no partitioned
tail correction. This is an explicitly **uncorrected truncated pole sum**:
the omitted spectral contribution is set to zero and must be bounded through
dense-reference agreement and refinement. It uses the native surface-amplitude
normalization and outer-region solver. Different export paths must also agree
on channels, thresholds and sign/degeneracy-aware boundary data before release.

The first B1 MPI-export differential passes: boundary energies differ by at most
**5.97e-13 Hartree**, amplitudes by **3.46e-8**, and all static channel/threshold/
multipole records match byte-exactly. Its 128-pole phases differ by **1.56 rad**
modulo π, well above the 0.05-rad gate. The native pipeline completes, but its
driver initially rejects a zero-error stderr summary embedded in a fitted-value
block. The parser repair removes only the explicitly empty diagnostic summary;
raw original bytes and the failing analysis exit are retained.

### Independent omission diagnostic

The spectral boundary sum is evaluated with the selected pole energies and
channel amplitudes. Setting the complementary poles to zero changes the entire
energy-dependent R-matrix background, even if the selected eigenstates and their
boundary amplitudes are accurate. The outer-region transformation is nonlinear,
so compare complete propagated phases rather than treating a small pole error
as an observable bound.

`truncate_boundary` in the experimental replay module slices an energy-ordered
**dense-oracle** boundary file without changing channels, thresholds, multipoles,
radius, surface normalization or adding any correction. A full-count copy must
be byte-identical. The independent finite diagnostic uses **128/512/2048/4096/
6144/8192/8350** poles, both Pi sectors, the unchanged 99-point grid, and identical
1.8–3.3-requested-eV independent fits at background orders 1–4. The 8350-pole
endpoint must also reproduce the dense phases. It determines the omission error
and needed spectral coverage without paying for repeated eigensolutions.
These outer-only costs are recorded separately from native iterative controls.

## Read-only solver telemetry

[`scattering_telemetry.c`](scattering_telemetry.c) interposes the C `EPSSolve`
call without changing matrix or eigensolver options. After the native solve it
records the actual PETSc matrix type, Krylov basis dimensions, spectral selection,
convergence reason/count, stored/allocated entries and absolute eigenpair residuals.
Its headers and shared libraries must match the runtime image. For SBAIJ,
stored-entry counts describe the storage representation, including allocated
structural entries; they are not the full symmetric matrix's numerical density.
PETSc's reported memory is retained as telemetry; aggregate kernel/container
profiles supply the actual peak-memory measurement.

Build the sidecar in the SLEPc-enabled `build` image, which contains development
headers and binutils. Set `OMPI_CC=/opt/compiler/bin/gcc` for the relocated MPI
wrapper. A runtime image alone may lack system C headers or the assembler.
For example, with `BUILD_IMAGE` naming that compatible build image:

```bash
mkdir -p "$ARTIFACT_ROOT/telemetry"
docker run --rm --entrypoint /opt/ukrmolp/entrypoint.sh \
  -e OMPI_CC=/opt/compiler/bin/gcc \
  -v "$PWD:/source:ro" -v "$ARTIFACT_ROOT/telemetry:/telemetry" \
  "$BUILD_IMAGE" bash -c '
    mpicc -shared -fPIC -O2 -Wall -Wextra \
      -I/opt/ukrmolp/include -I/opt/ukrmolp/include/petsc \
      -I/opt/ukrmolp/include/slepc \
      /source/projects/ukrmol_co/scattering_telemetry.c \
      -o /telemetry/scattering_telemetry.so \
      -L/opt/ukrmolp/lib -lslepc -lpetsc -ldl
  '
```

Record the builder/runtime image IDs, engine/library/source hashes, compiler and
sidecar bytes. The local experiment builds in a temporary development container
and retains its downloaded header/binutils packages and checksums. The native
controls still run in the original pinned scientific image.

## Run and interpret a control

Set `ARTIFACT_ROOT` to a writable persistent root containing the completed dense
reference under `runs/`; select a new `CONTROL_NAME`, the source-built SLEPc
`IMAGE`, and a verified physical-core `CPU_GROUP`. The local contract uses four
ranks, a 16-GiB container, at least 32 GiB available RAM and 20 GiB free disk.
Check CPU ownership and SMT partners before each launch.

```bash
docker run --rm --entrypoint /opt/ukrmolp/entrypoint.sh \
  --cpuset-cpus="$CPU_GROUP" --memory=16g --memory-swap=16g --shm-size=1g \
  -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
  -e PYTHONPATH=/opt/qmodeling \
  -v "$PWD:/opt/qmodeling:ro" -v "$ARTIFACT_ROOT:/work" \
  -w /opt/qmodeling "$IMAGE" \
  python3 -m projects.ukrmol_co.sparse_scattering \
    --reference /work/runs/co-eq-ccdz-sa10-tight-cc40-matched \
    --workdir "/work/diagnostics/$CONTROL_NAME" \
    --telemetry /work/telemetry/scattering_telemetry.so \
    --eigenpairs 128 --ranks 4 --boundary-export mpi
```

Retain `config.json`, `comparison.json`, all per-sector solver/export inputs and
logs, CIDATA, boundary files, phases, telemetry, `resources.csv/json` and the
frozen source/build records. A nonzero driver exit stops the finite queue;
passing eigensolver checks with failing phases permits the next predeclared
spectral-refinement control, not a larger-model release.

Replay costs exclude generating the reused electronic/integral/configuration
inputs. Report those earlier costs separately. Missing residuals, silent dense
fallback, incomplete spectra/exports or inconsistent target lineage reject.
Observable agreement and stability on further refinement, including fixed-window
fits with unchanged backgrounds, precede a CAS(10,11) qualification. CAS(10,12)
then needs its own spectral refinement and a measured resource envelope with
headroom. None of these controls establishes electronic-model convergence or
replaces the geometry/extraction/continuity gates.
