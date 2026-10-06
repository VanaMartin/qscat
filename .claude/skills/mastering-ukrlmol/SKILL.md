---
name: mastering-ukrlmol
description: Use when building, deploying, profiling or troubleshooting UKRmol+ calculations in Docker on a local or SSH-accessible remote CPU host, or planning MPI and concurrent-geometry execution. Covers source pins, engine validation, persistent scratch, precision-library matching and reproducible artifacts.
---

# Mastering UKRmol+

Run UKRmol+ as an external scattering engine with reproducible inputs and
measured costs. The skill identifier preserves the name `mastering-ukrlmol`;
the software is UKRmol+.

Read `docker/ukrmol-plus/README.md` for build/run commands and engine pins.
Read the chosen experiment's README for its molecule, model, units, manifests
and qualification evidence. Keep host addresses and molecule-specific settings
in the experiment, supplied by the user or its continuation notes.

## Establish the execution contract

1. Determine the SSH host, remote repository checkout, persistent artifact root,
   image target, available physical CPUs, RAM and disk. Use explicit variables
   such as `REMOTE_HOST`, `REMOTE_REPO`, `RUN_ROOT` and `IMAGE` in shared examples.
2. On the host, inspect `lscpu -e=CPU,CORE,SOCKET,NODE`, `free -h`,
   `df -h "$RUN_ROOT"`, Docker availability and existing workloads. CPU IDs
   need not be consecutive physical cores; inspect topology before pinning.
3. Keep scratch and outputs on a persistent host filesystem, conventionally
   under `/home`, bind-mounted at `/work`. Check the image user's UID/GID and
   write access before launching. Put `TMPDIR` and `PSI_SCRATCH` inside the run.
   Integral and `fort.*` files also count toward the storage budget.
4. Choose new run and batch names. Completed and failed directories are evidence;
   preserve them rather than reusing an old directory for a retry.

## Build and validate the engine

Use the repository's UKRmol-specific Dockerfile. Its pinned upstream layer
supplies a compatible MPI/compiler/BLAS/ScaLAPACK/Psi4 environment. The source
build downloads SHA256-checked UKRmol-in/out archives; GBTOlib comes from the
inner-region archive so its interface matches the engine.

| Target | Role |
|---|---|
| `pilot` | Generic environment with upstream reference binaries |
| `build` | Compile/install UKRmol-in/out and bundled GBTOlib |
| `test` | Upstream serial and two-rank reference comparisons |
| `runtime` | Generic source-built double engine, gated through `test` |
| `co-pilot` | CO experiment on the reference engine |
| `co-source` | CO experiment on the tested source engine; default target |

Build on the intended CPU host. The tested toolchain is Linux x86-64; establish
new evidence before claiming support for another architecture. Record the
Dockerfile/build-helper hashes, source checksums, toolchain digest, image ID
and test logs. Mutable tags are convenient names, not provenance identifiers.
Pinned engine sources do not make apt/pip resolution or image bytes identical
across rebuilds.

Treat the upstream test stage as an engine gate, then compare a representative
molecule calculation with the reference engine. Upstream target tests do not
certify a scattering model's accuracy. Preserve the source-built provenance
and `/opt/ukrmolp/source-tests/` logs.

Keep precision variants coherent: select both `bin.double`/`lib.double` or
`bin.quad`/`lib.quad`. The source build replaces **double** only; quad remains
the upstream reference implementation. The entrypoint initializes the upstream
environment and defaults to double libraries.

## Execute finite, measurable jobs

Use one fresh container per calculation. Pin disjoint CPU groups, bound memory
and shared memory, and set `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1` and
`MKL_NUM_THREADS=1` when MPI provides the parallelism. Pass explicit ranks in
the experiment's job inputs. Avoid inferring an optimal layout from core count.

For the included CO experiment, `projects/ukrmol_co/batch.py` snapshots the
source read-only and records image/source hashes, commands, exit codes and
batch wall time. Use its explicit `--root`, `--image`, `--cpu-groups` and
`--memory` arguments. For a different molecule, reuse this execution pattern
with its own input deck and analyzer.

Measure a representative geometry at one, two and four ranks, then measure
independent jobs on the same physical-core allocation. Compare batch throughput,
aggregate CPU time, per-job kernel memory peak and total disk use. Compact
models may favour more independent workers; larger close-coupling models need
their own measurement. Extrapolate inner-region and per-energy outer-region
cost separately, label forecasts conditional, and remeasure after changing
active space, basis or channel count.

## Preserve a calculation's contract

Retain the exact input deck, script archive checksum, copied/instrumented
templates, configuration, stage timings, raw logs, resource records, eigenphases,
cross sections, target properties, K-matrices and native resonance output.
Instrument upstream command execution to stop on an engine failure: some
scripts print a failure and continue. Preserve the original batch exit code
when later analyzer repairs recover a successful engine calculation.

Label units at each boundary. In this pipeline, geometry is bohr, external
scattering grids are eV, cross sections are bohr², phases are radians, native
RESON positions/full widths are Rydberg, and saved target/resonance records
are Hartree. Verify other interfaces rather than assuming the same convention.

Check complete grids, finite/nonnegative cross sections, final-state sums,
expected symmetry degeneracy modulo π and electronic target-energy consistency.
Compare like electronic methods and integral approximations: density-fitted
Psi4 SCF is not an exact-integral UKRmol RHF reference.

## Diagnose failures before expanding the calculation

| Symptom | Action |
|---|---|
| Missing Perl modules or SciPy | Use the packaged dependency layer; the upstream image alone is incomplete for this workflow |
| Cold-cache HTTPS certificate failure | Check the installed CA bundle and `SSL_CERT_FILE`; the relocated Python has an upstream build-machine default path |
| MPI cannot find its compiler | Preserve the relocated `OMPI_FC`/`OMPI_CC` settings in the build stage |
| Negative kinetic-matrix eigenvalue | Inspect continuum containment, sphere, orthogonalization deletion and precision; preserve the failed input and rerun in a fresh directory |
| `NDIMX TOO SMALL`, `CDIMX TOO SMALL` | Increase all relevant CONGEN input workspaces, including `NODIMX`; available host RAM does not enlarge these input limits |
| Quad executable fails against double libraries | Select matching executable and shared-library variants |
| CC cross-section columns appear incomplete | Check for repeated blocks of final-state columns and excited-initial-state sections |
| No automatic resonance fit | Inspect phases and the search window; near threshold use pole/bound-state diagnostics, not a fabricated zero width |

Treat radius, deletion threshold, basis and precision fixes as model-dependent
choices requiring checks. A value that rescued one molecule's deck is not a
universal numerical prescription.

## Separate numerical stability from electronic qualification

Vary energy spacing/window, fit background, continuum angular cutoff, sphere
radius, deletion threshold and propagation. Separately vary target orbitals,
active/frozen spaces, atomic basis and explicit channel count. Inspect compressed,
equilibrium and stretched geometries and target-state/orbital continuity.
Compare published results critically, including empirical adjustments and
disagreement between papers.

Ground-state CASSCF orbitals are not a multi-spin/multi-irrep state average.
A numerically stable SEP result is not evidence that its target and scattering
correlation are balanced. An empty fit is not a bound-anion diagnosis. Symmetry
contributions need not be a full polar-molecule elastic cross section, and
fixed-nuclei data do not uniquely determine a local potential.

Report: what model ran; which checks passed; measured cost and batch layout;
where inputs/artifacts are preserved; convergence evidence; and the next
qualification gate. Use the experiment's chosen tolerances, with their physical
meaning and limitations stated. SSH/Docker batches are the implemented
interface here; HTTP and AWS execution require additional implementation.
