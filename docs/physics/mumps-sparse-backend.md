# The MUMPS complex-symmetric sparse backend

**Location:** `qscat.linalg` (`SparseLU`'s `backend=` dispatch, plus
`set_default_backend` / `get_default_backend` / `default_backend`); the
implementation is `qscat.linalg._mumps_backend` (`_MumpsBackend`,
`mumps_available`). **Provisioning:** `docker/base.Dockerfile` (system MUMPS +
synthesized pkg-config files); the `qscat[mumps]` optional extra
(`libs/qscat/pyproject.toml`). **Benchmark:** `benchmarks/mumps_vs_superlu.py`.
**Origin:** sub-project #8; design plan at
`docs/superpowers/specs/2026-07-2x-mumps-sparse-backend-*` and the task
briefs/reports under `.superpowers/sdd/`. **Units:** atomic units throughout,
as everywhere in the library.

## Why the sparse LU is the hot path

`docs/physics/nd-tensor-hamiltonian.md` already established the shape of the
cost: for the N₂ 2-D electron-scattering problem, assembling the Hamiltonian is
cheap (~0.1 s) and factorizing it is not. On the production N₂ 2-D deck --
electronic grid `q=8`, 33+15 elements, `n=335`; nuclear grid `q=14`, 23+10
elements, `n=428`; `N = 335·428 = 143,380`, `nnz = 3,276,450`, ~22.9
nonzeros/row -- one sparse LU factorization is the entire runtime, and the
back-substitution behind it is a rounding error by comparison.

The structure of the physics makes this worse, not better, in the way that
matters for planning. In a driven-equation (Lippmann-Schwinger / T-matrix)
scattering calculation, **every** final channel at a single collision energy
shares the same matrix `A = (E_tot·I − H)`: one factorization serves all
right-hand sides at that energy, so total cost scales with the number of
*energies* sampled, not the number of channels. `SparseLU` exists precisely to
factor once and solve many. That makes the per-energy factorization the one
number that decides whether a production energy sweep -- and by extension the
NO and F₂ decks behind it -- finishes in an afternoon or a week.

With scipy's SuperLU (the only backend before this sub-project), that number was
259.5 s and a 7.4 GB peak on the production deck (measured; see the table). A
multi-hundred-point energy sweep at ~260 s/point is on the order of a day per
molecule -- uncomfortably close to, and for the larger decks over, the "under an
hour for all of N₂/NO/F₂" bar this sub-project was chartered to clear.

## Why every matrix is complex-symmetric -- and why SuperLU can't use it

Under exterior complex scaling (ECS; see `docs/physics/femdvr-ecs.md` and the
"ECS consequences" section of `docs/physics/nd-tensor-hamiltonian.md`), the
assembled Hamiltonian `H` is **complex symmetric but never Hermitian**:
`H = Hᵀ ≠ H†`. Both matrices this backend is asked to factor inherit that
property directly, because both are `H` plus a scalar multiple of the identity:

- the **driven / resolvent matrix** `A = (E_tot·I − H)` of sub-project #6
  (`n2_2d_cross_section`): `Aᵀ = E_tot·I − Hᵀ = E_tot·I − H = A`.
- the **Crank-Nicolson propagation matrix** `A = (I + iH·dt/2)` of sub-project
  #7 (`n2_2d_td_cross_section`, via `make_sparse_cn_stepper`):
  `Aᵀ = I + iHᵀ·dt/2 = A`.

Adding a scalar diagonal preserves symmetry, so both are complex-symmetric,
`A = Aᵀ`. They are *not* Hermitian (the ECS tail and, for CN, the `i` see to
that), so none of the usual Hermitian/positive-definite fast paths apply -- only
the plain complex-symmetric one.

SuperLU is a general **unsymmetric** LU solver. It has no notion of `A = Aᵀ`: it
stores and factors both the `L` and the `U` triangle independently and computes
a fill-reducing ordering (COLAMD by default) that targets a general sparsity
pattern. Handing it a complex-symmetric matrix throws away, silently, exactly
the structure that would halve the work and the storage. There is no SuperLU
option that recovers it; exploiting `A = Aᵀ` requires a solver written for the
symmetric case.

## The MUMPS `SYM=2` backend

MUMPS is that solver. Its `SYM=2` matrix type is *general symmetric* (`A = Aᵀ`,
no definiteness assumed) -- exactly the complex-symmetric case here, as opposed
to `SYM=1` (symmetric positive definite, which these matrices are not) or
`SYM=0` (general unsymmetric). The backend drives it through python-mumps'
factor-once/solve-many `Context`, mirroring `SparseLU`'s own contract:

```python
ctx = mumps.Context()
ctx.set_matrix(sp.triu(A).tocsc(), symmetric=True)  # SYM=2: UPPER TRIANGLE ONLY
ctx.analyze()                                        # ordering='auto' -> SCOTCH
ctx.factor()                                         # factor once
x = ctx.solve(b)                                     # solve many (reuse ctx)
```

**The correctness trap: `SYM=2` reads only the upper triangle.** MUMPS assumes
the lower triangle from symmetry and never looks at it. So the backend passes
`sp.triu(A)`, not the full `A`: supplying the full matrix with `symmetric=True`
would double-count every off-diagonal entry and silently produce the wrong
factorization -- a plausible-looking wrong number, not a crash. `SparseLU`'s
`symmetric` flag (auto-detected as `A == A.T` via the cheap O(nnz)
`(abs(A - A.T)).max() == 0`, or set explicitly) selects `SYM=2` + upper-triangle
when true and falls back to `SYM=0` on the full matrix when false. The
differential test in `libs/qscat/tests/test_mumps_backend.py` actively guards
the trap: it confirms `sp.triu(A)` under `SYM=2` reproduces the full-matrix
SuperLU solve to machine precision, which it would not if the triangle handling
were wrong.

### Dispatch, and SuperLU as fallback + oracle

`SparseLU(A, backend=...)` selects the engine:

- **`"auto"`** (the default) -- use MUMPS if `mumps_available()`, else SuperLU.
  On a MUMPS-less box (a bare Mac) `"auto"` and `"scipy"` are identical in every
  observable way, so existing call sites are bit-for-bit unchanged there.
- **`"scipy"`** -- force SuperLU, always available (numpy/scipy only).
- **`"mumps"`** -- force MUMPS; raises a clear `RuntimeError` naming the
  `qscat[mumps]` extra if MUMPS is absent, rather than silently falling back.
- **process-wide override** -- `set_default_backend(name)` and the
  `default_backend(name)` context manager change what `"auto"` resolves to,
  process-wide and (for the context manager) scoped + exception-safe. An
  explicit `"scipy"`/`"mumps"` at a call site always wins over the override.
  This is the seam that forces a whole computation which builds `SparseLU`
  internally (e.g. `ve_cross_section_2d`, which exposes no `backend=` kwarg)
  onto one engine for a backend-equivalence check, with zero changes to the
  solver's signatures.

**SuperLU is deliberately kept as both the fallback and the differential
oracle.** It is the reference the MUMPS path is validated against, not dead
weight: `test_mumps_backend.py` matches the two backends' solves to ~7e-16
(machine precision), and the physics-level check
(`validation/n2/test_backend_equivalence.py`) recomputes a full #6 exact VE
cross section twice -- every internal `SparseLU` forced through SuperLU, then
through MUMPS -- and asserts `rtol = 1e-9`. Identical physics; only the cost
differs. Keeping the pure-Python SuperLU path means the core library stays
numpy/scipy-only and any MUMPS result can always be re-derived without MUMPS
present.

## The benchmark

`benchmarks/mumps_vs_superlu.py` builds the **real** N₂ 2-D driven matrix
`A = (E_tot·I − H_2D)` (`H_2D` from
`projects.n2_2d_cross_section.hamiltonian2d.build_h2d`, `E = 0.2 Ha`,
complex-symmetric) at three grids and factors + solves it under each backend.
The script **asserts nothing** about the speedup -- it measures. Each (grid,
backend) row runs in a **fresh subprocess** so peak RSS
(`resource.getrusage`'s `ru_maxrss`) is attributable to that one factorization;
it deliberately does **not** call `SparseLU.memory_bytes()`, which would
materialize SuperLU's `L`/`U` factors and add ~6 GB at production scale.
Environment: Docker `qmodeling-base:latest`, Linux aarch64, Python 3.12.12,
mumps_seq 5.5.1 / python-mumps 0.0.6.

| grid | backend | N | nnz | factor (s) | solve (s) | peak RSS (MB) | fill_factor | ordering | residual |
|---|---|---|---|---|---|---|---|---|---|
| working | scipy | 26,857 | 476,377 | 5.255 | 0.0266 | 479 | 44.62 | COLAMD | 1.81e-13 |
| working | mumps | 26,857 | 476,377 | 0.440 | 0.0070 | 171 | 6.89 | scotch | 1.30e-13 |
| td | scipy | 47,188 | 886,664 | 18.959 | 0.0622 | 1,410 | 58.33 | COLAMD | 5.09e-13 |
| td | mumps | 47,188 | 886,664 | 0.799 | 0.0121 | 260 | 7.26 | scotch | 3.36e-13 |
| production | scipy | 143,380 | 3,276,450 | 259.531 | 0.3871 | 7,420 | 93.35 | COLAMD | 3.16e-13 |
| production | mumps | 143,380 | 3,276,450 | 3.574 | 0.0507 | 809 | 7.95 | scotch | 2.91e-13 |

Per grid, MUMPS vs SuperLU (measured, not assumed):

- **working (27k):** factor **11.9×** faster, peak RSS **2.8×** smaller.
- **td (47k):** factor **23.7×** faster, peak RSS **5.4×** smaller.
- **production (143k):** factor **72.6×** faster (3.6 s vs 260 s), peak RSS
  **9.2×** smaller (0.8 GB vs 7.4 GB).

The residuals agree to ~1e-13 across backends at every grid: the two engines
compute the same solution, only the cost differs. The production factorization
drops from 260 s / 7.4 GB to 3.6 s / 0.8 GB -- which puts a multi-hundred-point
energy sweep for N₂/NO/F₂ comfortably inside the "under an hour" bar, with
room to spare, rather than skirting it.

## The mechanism -- two compounding causes, not one

It is tempting to read the table and credit the whole win to MUMPS's ordering,
because the `fill_factor` column jumps out. That reading is wrong. The MUMPS
advantage has **two** distinct, compounding causes:

1. **Symmetric single-triangle storage.** `SYM=2` factors only the upper
   triangle of `A = Aᵀ`; SuperLU stores and factors both `L` and `U`
   independently. That is roughly a **2× factor** in both arithmetic and
   storage before any ordering effect, and SuperLU cannot recover it -- it has
   no symmetric mode.

2. **A better fill-reducing ordering.** MUMPS's `analyze()` picks SCOTCH
   nested-dissection ordering, which holds `fill_factor` at ~7-8 at every scale
   (6.89 → 7.26 → 7.95). SuperLU's default COLAMD, ordering a general
   unsymmetric pattern, lets fill grow with N: 44.62 → 58.33 → 93.35. Lower
   fill means both **less arithmetic** (the factorization is superlinear in the
   number of factor entries) *and* **less memory** -- which is why both the
   time gap and the RSS gap widen together as N grows, rather than one or the
   other.

**Do not attribute the whole gap to ordering.** The two `fill_factor` numbers
are not even measured the same way and cannot be compared as a clean ratio:
SuperLU's is `(L.nnz + U.nnz) / A.nnz`, a two-triangle count, while MUMPS's is
its INFOG entries-in-factors over `A.nnz`, effectively a single-triangle count.
The `fill_factor` column therefore already **conflates** the storage halving
(cause 1) with the ordering improvement (cause 2); reading the ~44→7 and ~93→8
drops as pure ordering wins double-counts the triangle effect. The honest
statement is: MUMPS wins because it exploits symmetry to store one triangle
*and* orders that triangle better, and both effects grow with N.

## Honest caveats on the absolute numbers

- **SuperLU's absolute time is inflated by the container.** The 259.5 s
  production factorization here is roughly 2× the repo's historical
  Mac-native 128 s figure (recorded in `docs/physics/nd-tensor-hamiltonian.md`
  and `sparse_lu.py`'s docstring). That gap is Docker-on-Mac virtualization
  overhead, not a regression. The **speedup ratio is unaffected**: both
  backends run in the same container under the same overhead, so the 72.6× /
  9.2× figures are valid comparisons even though the SuperLU absolute is
  container-slowed.

- **SuperLU's production peak (7.4 GB) is lower than the docstring's 13.6 GB.**
  Those measure different things. The 13.6 GB figure in
  `nd-tensor-hamiltonian.md` / the historical spike is a factorization peak
  that *includes* materializing the `L`/`U` factors via `memory_bytes()` (a
  further ~+6 GB cache; see `sparse_lu.py`'s module docstring). The benchmark
  deliberately never calls `memory_bytes()`, so its 7.4 GB is the bare
  factorization high-water mark. Both are honest; they are not the same
  measurement. (The `sparse_lu.py` docstring's 13.6 GB figure is left as-is
  because it correctly describes the `memory_bytes()`-inclusive peak it is
  attached to; this note documents the distinction so the two numbers are not
  read as a contradiction.)

## Provisioning MUMPS

The clean path is the Docker base image, per Task 1:

- `docker/base.Dockerfile` installs Debian's `libmumps-seq-dev` (the
  sequential, non-MPI MUMPS) plus `libscotch-dev` (Debian's seq MUMPS links
  SCOTCH, not METIS, for its ordering).
- Debian ships MUMPS headers + libs but **no pkg-config `.pc` files**, whereas
  python-mumps discovers MUMPS via pkg-config under the conda-forge names
  `{d,z,c,s}mumps_seq`. The Dockerfile **synthesizes** those `.pc` files in the
  multiarch pkgconfig dir, pointing at the Debian libs; ELF `NEEDED` carries the
  transitive dependencies (`mumps_common_seq`, SCOTCH, LAPACK, gfortran), so
  `-l<name>_seq` links transitively. A build-time `pkg-config --exists` gate
  fails the base image if `dmumps_seq`/`zmumps_seq` don't resolve.
- `qscat[mumps]` (`python-mumps>=0.0.6`) is an **optional extra**; the core
  `qscat` dependencies stay numpy/scipy-only. The Docker `test` stage adds
  `--extra mumps` so the MUMPS backend is exercised, not skipped; the `runtime`
  stage deliberately omits it, keeping python-mumps out of the production image.

**The Mac dev box has no MUMPS**, and that is by design. All MUMPS work happens
in the container. Every MUMPS-touching test is `@skipif(not mumps_available())`,
so the Mac suite stays green (the MUMPS tests skip) while the same tests
**run and pass** in the Docker `test` image. The absence-path tests
(`backend="mumps"` errors clearly; `backend="auto"` falls back to SuperLU
bit-identically) are the mirror image -- they `@skipif(mumps_available())`, so
they run on the Mac and skip in the container.

## Lifecycle position and deferred levers

This is a **stage-4 optimization** in the qModeling lifecycle: an already-
validated capability (`SparseLU`) made faster on a proven hot path, with the
pure-Python SuperLU path kept as the differential oracle it is validated
against. It is not a Rust kernel -- the hot path here is a call into a mature
Fortran direct solver, where the right move is to *dispatch* to it, not to
reimplement LU.

Deliberately **not** done here, and left for later if a measurement justifies
them:

- **complex64 + iterative refinement** -- halve the factor storage again by
  factoring in single precision and refining back to double; needs a residual
  study to confirm the ECS-conditioned matrices tolerate it.
- **symbolic / numeric reuse across the TI energy sweep** -- the sparsity
  pattern of `(E_tot·I − H)` is *identical* at every energy; only the diagonal
  shifts. Reusing MUMPS's `analyze()` (symbolic) phase across energies, and
  possibly the numeric factorization structure, would cut the per-energy cost
  further. `SparseLU` currently re-analyzes per matrix.
- **MKL PARDISO as a second dispatch backend** -- a natural third `backend=`
  option on x86-64, sharing the same complex-symmetric contract; fits the
  existing dispatch seam with no API change.
- **the publishing / packaging pipeline** -- shipping `qscat[mumps]` as an
  installable wheel against a bundled or system MUMPS, beyond the Docker image.
- **a Rust non-LU (iterative / matrix-free) kernel** -- only if a problem
  outgrows what a direct sparse solver fits in RAM; at the sizes measured here,
  direct LU is demonstrably sufficient.

## Reproduce

In the Docker `test` image (or any container with system MUMPS + `qscat[mumps]`
installed):

```bash
uv run python -m benchmarks.mumps_vs_superlu --grids working td              # 27k + 47k
uv run python -m benchmarks.mumps_vs_superlu --grids working td production   # + 143k deck (~4.5 min, ~7.4 GB for SuperLU)
```

The driver measures each (grid, backend) in a fresh subprocess, writes the
Markdown table to `.superpowers/sdd/task-4-benchmark-table.md`, and prints it.

## See also

- `docs/physics/nd-tensor-hamiltonian.md` -- `SparseLU`, the Kronecker-sum
  Hamiltonian, the ECS complex-symmetry that makes `A = Aᵀ`, and the historical
  SuperLU-only production cost this backend improves on.
- `docs/physics/n2-2d-cross-section.md` -- sub-project #6, the driven
  `(E_tot·I − H)` solver whose factorization is the hot path.
- `docs/physics/n2-2d-td-cross-section.md` -- sub-project #7, the sparse
  Crank-Nicolson `(I + iH·dt/2)` factorization, the other complex-symmetric
  matrix this backend serves.
- `docs/physics/femdvr-ecs.md` -- the 1-D FEM-DVR-ECS grid and the exterior
  complex scaling that makes the Hamiltonian complex symmetric.
- `.superpowers/sdd/task-1-report.md` -- the verified MUMPS provisioning +
  `SYM=2` recipe; `.superpowers/sdd/task-4-report.md` -- the benchmark
  methodology and full results.
