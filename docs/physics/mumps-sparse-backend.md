# The MUMPS complex-symmetric sparse backend

**Location:** `qscat.linalg` (`SparseLU`'s `backend=` dispatch, plus
`set_default_backend` / `get_default_backend` / `default_backend`); the
implementation is `qscat.linalg._mumps_backend` (`_MumpsBackend`,
`mumps_available`). **Provisioning:** `docker/base.Dockerfile` (system MUMPS +
synthesized pkg-config files); the `qscat[mumps]` optional extra
(`libs/qscat/pyproject.toml`). **Benchmark:** `benchmarks/mumps_vs_superlu.py`.
**Origin:** the design rationale is recorded in
`docs/superpowers/specs/2026-07-26-mumps-sparse-backend-design.md`.
**Units:** atomic units throughout,
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
258 s and a 7.4 GB peak on the production deck (measured; see the table). A
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
`symmetric` flag (auto-detected, or set explicitly) selects `SYM=2` +
upper-triangle when true and falls back to `SYM=0` on the full matrix when
false. The differential test in `libs/qscat/tests/test_mumps_backend.py`
actively guards the trap: it confirms `sp.triu(A)` under `SYM=2` reproduces the
full-matrix SuperLU solve to machine precision, which it would not if the
triangle handling were wrong.

**The auto-detect must use a scaled tolerance, not exact equality -- or `SYM=2`
never engages on real physics.** This is a subtle and important point. The real
N₂ matrices are `A = Aᵀ` *mathematically*, but they are assembled by
Kronecker-sum reordering of float arrays, so `A - Aᵀ` is not bit-zero: on the
working deck, `max|A − Aᵀ| = 4.5e-13` against `max|A| = 1.3e4`, a **relative
asymmetry of ~3.6e-17** -- one ULP, pure round-off. An *exact*-equality detect
(`(abs(A − Aᵀ)).max() == 0`) therefore returns **False on every real N₂
matrix**, silently routing the MUMPS backend onto `SYM=0` (general unsymmetric,
full matrix) -- forfeiting the single-triangle storage that is the entire reason
the backend exists. So the auto-detect compares against a **scaled tolerance**:
`(abs(A − Aᵀ)).max() <= _SYM_RTOL · abs(A).max()` with `_SYM_RTOL = 1e-12`
(`qscat.linalg.sparse_lu`). That threshold sits ~5 orders of magnitude *above*
the real matrices' ~3.6e-17 relative asymmetry (a decisive accept) yet ~12
orders *below* the O(1) relative asymmetry of a genuinely non-symmetric matrix
(a decisive reject). The tightness is a correctness guard, not a nicety: `SYM=2`
takes the upper triangle as truth and reconstructs the lower from it, so a
truly-asymmetric matrix wrongly accepted would produce a *wrong* factorization.
A round-off-symmetric matrix is only perturbed at the ~1e-13 level, which a
backward-stable factorization absorbs; a zero matrix (`abs(A).max() == 0`) is
treated as trivially symmetric. An explicit `symmetric=True`/`False` always
overrides the auto-detect.

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
mumps_seq 5.5.1 / python-mumps 0.0.6. The **SYM** column records the MUMPS
matrix type actually driven -- `SYM=2` (complex-symmetric, single upper
triangle) when the auto-detect flagged the matrix symmetric, `SYM=0` (general
unsymmetric, full matrix) otherwise; SuperLU has no symmetric mode, so its cell
is `-`. Because the real N₂ matrices are round-off-symmetric and the auto-detect
now uses the scaled `_SYM_RTOL` tolerance (above), **every MUMPS row here runs
`SYM=2`**.

| grid | backend | SYM | N | nnz | factor (s) | solve (s) | peak RSS (MB) | fill_factor | ordering | residual |
|---|---|---|---|---|---|---|---|---|---|---|
| working | scipy | - | 26,857 | 476,377 | 4.936 | 0.0254 | 479 | 44.62 | COLAMD | 1.81e-13 |
| working | mumps | SYM=2 | 26,857 | 476,377 | 0.411 | 0.0052 | 147 | 3.54 | scotch | 1.63e-13 |
| td | scipy | - | 47,188 | 886,664 | 18.477 | 0.0620 | 1,409 | 58.33 | COLAMD | 5.09e-13 |
| td | mumps | SYM=2 | 47,188 | 886,664 | 0.789 | 0.0095 | 213 | 3.68 | scotch | 4.93e-13 |
| production | scipy | - | 143,380 | 3,276,450 | 258.116 | 0.3709 | 7,418 | 93.35 | COLAMD | 3.16e-13 |
| production | mumps | SYM=2 | 143,380 | 3,276,450 | 3.175 | 0.0388 | 625 | 4.00 | scotch | 4.37e-12 |

Per grid, MUMPS (`SYM=2`) vs SuperLU (measured, not assumed):

- **working (27k):** factor **12.0×** faster, peak RSS **3.3×** smaller.
- **td (47k):** factor **23.4×** faster, peak RSS **6.6×** smaller.
- **production (143k):** factor **81.3×** faster (3.2 s vs 258 s), peak RSS
  **11.9×** smaller (0.6 GB vs 7.4 GB).

The residuals agree to ~1e-12 or better across backends at every grid: the two
engines compute the same solution, only the cost differs. (The production MUMPS
residual, 4.4e-12, is a touch larger than SuperLU's 3.2e-13 because `SYM=2`
reconstructs the lower triangle from the upper of a matrix that was only
round-off-symmetric to begin with -- still ~4 orders of magnitude inside any
tolerance that matters.) The production factorization drops from 258 s / 7.4 GB
to 3.2 s / 0.6 GB -- which puts a multi-hundred-point energy sweep for N₂/NO/F₂
comfortably inside the "under an hour" bar, with room to spare, rather than
skirting it.

**Historical note -- these are the `SYM=2` numbers; an earlier revision of this
table reported `SYM=0`.** The first benchmark ran before the auto-detect was
fixed: the exact-equality symmetry check rejected the round-off-symmetric N₂
matrices, so the MUMPS rows silently ran `SYM=0` (general unsymmetric) and
reported *full*-factor `fill_factor` values (working 6.89, td 7.26, production
7.95) at 11.9× / 23.7× / 72.6× speedup and 2.8× / 5.4× / 9.2× memory. Switching
to the scaled-tolerance detect engaged the intended `SYM=2` single-triangle
mode, which roughly **halved** the `fill_factor` at every grid (6.89→3.54,
7.26→3.68, 7.95→4.00) and improved both the speedup and the memory ratio. The
speedup was always real; `SYM=2` makes it larger *and* makes the mechanism the
one the backend was chosen for.

## The mechanism -- two compounding causes, now separately measured

The MUMPS advantage has **two** distinct, compounding causes, and the two
`SYM=2`-vs-`SYM=0` benchmark runs let us *isolate* each one instead of
hand-waving. The key that makes the decomposition clean: a MUMPS `SYM=0` run
still uses MUMPS's SCOTCH ordering but stores **both** triangles, so its
`fill_factor` is a full-factor count directly comparable to SuperLU's, while the
`SYM=2` run adds single-triangle storage on top of the same ordering.

1. **A better fill-reducing ordering (cause 2, isolated by SYM=0 vs SuperLU).**
   Both SuperLU's `fill_factor` and MUMPS's `SYM=0` `fill_factor` are
   full-factor `(L+U)/A.nnz`-style counts, so they compare directly. SuperLU's
   default COLAMD, ordering a general unsymmetric pattern, lets fill grow with
   N: 44.62 → 58.33 → 93.35. MUMPS's `analyze()` picks SCOTCH nested-dissection
   ordering, which under `SYM=0` held the full-factor fill at ~7-8 at every
   scale (6.89 → 7.26 → 7.95, from the historical `SYM=0` run). That ~6× → ~12×
   fill reduction is **pure ordering** -- same two-triangle storage, better
   permutation.

2. **Symmetric single-triangle storage (cause 1, isolated by SYM=2 vs SYM=0).**
   Switching the *same* SCOTCH-ordered factorization from `SYM=0` to `SYM=2`
   drops the `fill_factor` by **roughly half** at every grid -- 6.89 → 3.54,
   7.26 → 3.68, 7.95 → 4.00 -- because `SYM=2` factors only the upper triangle
   of `A = Aᵀ`. This is the ~**2× factor** in arithmetic and storage that
   SuperLU structurally cannot recover (it has no symmetric mode), and it is now
   measured directly rather than asserted, by comparing the two MUMPS runs on
   identical matrices with identical ordering.

Lower fill means both **less arithmetic** (the factorization is superlinear in
the number of factor entries) *and* **less memory**, which is why both the time
gap and the RSS gap widen together as N grows. The honest, now-quantified
statement: SuperLU→MUMPS-`SYM=0` is the ordering win (~6-12× fill), and
`SYM=0`→`SYM=2` is the single-triangle win (~2× fill on top), and the shipped
default (`SYM=2`, ~3.5-4 fill vs SuperLU's ~44-93) delivers both.

Note the two shipped `fill_factor` columns are *not* the same measurement and
should not be divided as a raw ratio: SuperLU's is a two-triangle
`(L.nnz + U.nnz) / A.nnz`, while MUMPS `SYM=2`'s is its INFOG entries-in-factors
over `A.nnz`, a single-triangle count. The `SYM=0` intermediate above is exactly
what makes the decomposition rigorous: it supplies the apples-to-apples
full-factor count that pins the ordering effect before the triangle effect is
layered on.

## Honest caveats on the absolute numbers

- **SuperLU's absolute time is inflated by the container.** The 258 s
  production factorization here is roughly 2× the repo's historical
  Mac-native 128 s figure (recorded in `docs/physics/nd-tensor-hamiltonian.md`
  and `sparse_lu.py`'s docstring). That gap is Docker-on-Mac virtualization
  overhead, not a regression. The **speedup ratio is unaffected**: both
  backends run in the same container under the same overhead, so the 81.3× /
  11.9× figures are valid comparisons even though the SuperLU absolute is
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

The driver measures each (grid, backend) in a fresh subprocess and prints the
Markdown table; `--out PATH` also writes it to a file.

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
- `docker/base.Dockerfile` -- the verified MUMPS provisioning and the
  pkg-config shim the `SYM=2` build needs; the benchmark methodology and full
  results are in the "Benchmark" section above.
