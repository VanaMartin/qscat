# The N-dimensional tensor-product Hamiltonian

**Location:** `qscat.linalg` (`kron_sum`, `SparseLU`, `c_product`), `qscat.dvr`
(`kinetic_sparse`, `TensorGrid`, `kinetic_nd`, `potential_nd`, `hamiltonian_nd`).
**Origin:** port-scout archaeology of eMoScat's 2-D e-N2 machinery
(`.superpowers/sdd/n2-2d-exact-extraction.md`); design spec at
`docs/superpowers/specs/2026-07-22-nd-sparse-hamiltonian-design.md`
(sub-project #5). Consumed by sub-project #6, the exact 2-D electron-N2
scattering solver.
**Units:** atomic units throughout (energy in Hartree, length in Bohr, mass in
electron masses).

## Physical picture

Sub-projects #1-#4 built and validated a 1-D FEM-DVR-ECS stack and used it,
under the Local Complex Potential (LCP) approximation, to reduce the 2-D
electron-N2 vibrational-excitation problem to a 1-D nuclear-coordinate
calculation. The next step -- sub-project #6 -- drops that approximation and
solves the genuine 2-D problem, the same one eMoScat and Karel Houfek solved:
an electron at radial coordinate `r` scattering off N2 at internuclear
distance `R`, with the two coordinates coupled only through a potential
surface `V(r, R)`, no coupling in the kinetic energy. This document describes
the library layer that makes that 2-D (and, as it turns out, D-dimensional in
general) problem tractable, built and validated on its own terms *before* any
N2 physics touches it -- so that when a sub-project #6 result disagrees with
reference data, the disagreement can be attributed to physics, not to unproven
linear algebra.

## The construction

A Hamiltonian on a tensor product of D coordinate grids is separable in its
kinetic energy and diagonal in its potential exactly when it has the form

```
H(x_0, ..., x_{D-1}) = sum_d T_d(x_d)  +  V(x_0, ..., x_{D-1})
```

i.e. each kinetic term acts on one coordinate only, and the potential couples
the coordinates but never differentiates between them. Discretized on a
tensor-product basis (the outer product of D 1-D FEM-DVR-ECS bases, one per
coordinate), this becomes a **Kronecker sum plus a diagonal**:

```
H = sum_d  I x ... x T_d x ... x I  +  diag(V(x_0, ..., x_{D-1}))
```

where `T_d` is the 1-D kinetic-energy matrix on grid `d` (`qscat.dvr.kinetic`
or its sparse sibling `kinetic_sparse`) sandwiched between identity matrices
of every other grid's dimension, and `V` is evaluated pointwise at the tensor
grid's D-dimensional coordinate points and placed on the diagonal.

Two conditions make this valid, and both must hold or the construction is
simply the wrong matrix for the physics:

1. **Separable kinetic energy.** The kinetic operator must have no
   cross-derivative terms (`d^2/dx_i dx_j`, `i != j`); each `T_d` depends on
   coordinate `d` alone. This holds for the electron-N2 problem (independent
   radial coordinates `r` and `R`, each with its own reduced mass) but would
   fail for, e.g., a genuinely coupled bending/stretching normal-mode
   Hamiltonian.
2. **Diagonal (DVR) potential.** `V` must act as a pointwise multiplication in
   the chosen basis -- exactly the diagonal-potential approximation already
   used by `qscat.dvr.operators.hamiltonian` in 1-D (see
   `docs/physics/femdvr-ecs.md`), extended here to a potential of all D
   coordinates jointly. A DVR basis makes this an approximation whose error is
   controlled by how well the grid resolves `V`'s structure (element
   boundaries at potential discontinuities, sufficient quadrature order in
   regions of curvature), not an ad hoc truncation.

Given both, `kron_sum` plus a diagonal is not one useful representation among
several -- it is *the* matrix that this Hamiltonian discretizes to on this
basis, and every result in the "Validation" section below is a check that the
implementation actually produces it.

## Why it is dimension-general

`qscat.linalg.kron_sum` (the `qscat.linalg.kron` module) takes a sequence of
**arbitrary square sparse matrices**, not grids:

```python
def kron_sum(ops: Sequence[sp.spmatrix]) -> sp.csr_matrix:
    """sum_d  I x ... x ops[d] x ... x I,  for arbitrary len(ops)."""
```

It knows nothing about FEM-DVR, ECS, or physics -- only that each operand is
square and that the Kronecker-sum algebra applies. `D == 1` returns the single
operator unchanged (no special case is needed: with only one factor, the
identity blocks on either side both have dimension 1, so `I_1 x T_0 x I_1`
collapses to `T_0` by construction of the same loop that handles `D > 1`).

This has a concrete consequence beyond tidiness: a future angular-DVR,
finite-difference, or B-spline discretization for some coordinate composes
with FEM-DVR-ECS coordinates for the others at zero extra implementation
cost, because `kron_sum` never inspects what kind of operator it was handed.
`qscat.dvr.tensor.kinetic_nd` is the physics-flavored wrapper that supplies
FEM-DVR-ECS `T_d`'s specifically, but the underlying algebra layer does not
require it.

`qscat.dvr.tensor.TensorGrid` is the grid-aware layer one level up: it holds
D `FemDvrEcsGrid` instances, exposes `shape` (`(n_0, ..., n_{D-1})`), `size`
(`prod(shape)`), broadcastable coordinate arrays via `points()` so a potential
can be written as `V(r, R)` (or `V(x_0, ..., x_{D-1})` at any D) without
materializing a full meshgrid, a separable-state constructor `outer()`, and
`real_mask()` (below). `kinetic_nd`, `potential_nd`, and `hamiltonian_nd`
assemble the full construction on top of it:

```python
H = kinetic_nd(tgrid, masses) + diag(potential_nd(tgrid, V))
  = kron_sum([kinetic_sparse(g, m) for g, m in zip(tgrid.grids, masses)])
    + diag(V(*tgrid.points()))
```

## Index convention: C order, last axis fastest

A state on a D-dimensional tensor grid is a NumPy array of shape
`tgrid.shape = (n_0, ..., n_{D-1})`; `kron_sum`'s result acts on
`psi.ravel()` under NumPy's native C order, where the **last** axis is
fastest (`idx = i_{D-1} + n_{D-1}*(i_{D-2} + n_{D-2}*(...))`). This is the
convention that makes `reshape(tgrid.shape)` and `ravel()` round-trip without
any transposition, so it is the one this library uses throughout.

This is a **deliberate divergence** from eMoScat, which orders its 2-D basis
with the *first* coordinate (electronic `r`) fastest:
`idx = i_r + i_R * N_r` (`FemDvrEcsGrid2d.cpp:169`). The two conventions are
physically identical -- they differ only in how the same set of basis
functions is numbered -- but they are not interchangeable at the bit level:
a raw index dump compared directly against eMoScat's would look completely
scrambled despite both being correct. Anyone porting a reference vector,
matrix row, or debug printout from eMoScat must explicitly convert between
the two orderings; this is not a bug to chase if the numbers look
transposed relative to the C++ output.

## ECS consequences

Whenever any grid factor carries an exterior-complex-scaled (ECS) tail
(see `docs/physics/femdvr-ecs.md` for the 1-D theory), the assembled `H` is
**complex symmetric but never Hermitian**: `H = H^T != H^dagger`. Every
routine in this layer is written for that case unconditionally -- there is no
Hermitian fast path to fall into by accident, and no eigensolver or linear
solve in this library assumes conjugate symmetry.

Two things follow directly from `H` being merely symmetric, not Hermitian:

- **The c-product.** The natural inner product paired with a complex
  symmetric operator is the bilinear `c_product(a, b) = sum_i a_i b_i`, with
  **no** complex conjugation -- not NumPy's `vdot`, which conjugates its
  first argument and is the correct pairing only for a Hermitian operator.
  Using `vdot` under ECS produces a complex value with a plausible-looking
  magnitude and the wrong phase -- a quiet failure mode, not a crash. This
  has already bitten this repo once (sub-project #3's S-matrix, where the
  Hermitian convention produced a negative cross section), and even
  eMoScat's own reference implementation gets it formally wrong
  (`cblas_zdotc`, sesquilinear) and is saved only because every channel
  function it pairs against happens to be zeroed on the ECS tail. Naming the
  operation as `qscat.linalg.c_product` makes the correct choice explicit at
  every call site rather than relying on that kind of accidental
  cancellation.
- **`real_mask()`.** Under ECS, a driving term, an incident wavefunction, or a
  channel-projection function is only physically meaningful on the unscaled
  (real) region of the grid -- beyond the ECS pivot, the coordinate itself is
  a complex number, and evaluating a physical asymptotic form there is
  meaningless, not just imprecise. `TensorGrid.real_mask()` returns a flat
  boolean array, `True` exactly where **every** coordinate of a tensor-grid
  point lies in its own grid's unscaled region (`real_points <= R0`, ANDed
  across axes). Making this a property of `TensorGrid` itself, rather than
  something each physics routine has to remember to construct, is a
  deliberate defensive design choice: without it, forgetting to mask a driving
  term or a channel-projection integral silently contaminates the result with
  contributions from complex-valued coordinates the underlying asymptotic
  formula was never derived for -- again a plausible-looking wrong number,
  not a crash.

## What was actually validated

The full test suite lives in `libs/qscat/tests/test_kron_sum.py`,
`test_kinetic_sparse.py`, `test_tensor_grid.py`, `test_hamiltonian_nd.py`,
`test_sparse_lu.py`, and `test_nd_scale.py`, corresponding to V1-V6 of the
design spec.

**V1 -- `kron_sum` against dense `np.kron`.** Random small complex matrices at
D = 1, 2, 3, 4, with **unequal** per-axis dimensions (a square case would let
a transposed-index bug pass unnoticed). Exact to round-off
(`atol=1e-12`-`1e-14`).

**V2 -- `kinetic_sparse` against the dense `kinetic`.** Both share a private
`_element_block` helper for the per-element local `T_local` computation
(the `wze`/`dBF`/normalization/einsum block that the analytic
particle-in-a-box benchmark already pins), but the two functions
**deliberately keep independent bridge-accumulation code paths**: `kinetic`
scatter-adds each retained sub-block into a dense array (`T[...] += block`),
while `kinetic_sparse` emits COO triplets and lets `coo_matrix -> csr_matrix`
conversion sum the duplicate `(row, col)` entries at shared bridge indices.
Sharing the element-local math costs no coverage, since that part is already
exercised by the 1-D benchmarks; keeping the *accumulation* independent means
the dense-vs-sparse differential test still exercises two genuinely
different implementations of the bug-prone part -- bridge bookkeeping -- not
one implementation checked against itself. The two agree to round-off, and
`kinetic_sparse`'s nonzero count matches the analytic formula
`nnz = q^2 * tnel - 4q + 3 - tnel` (eMoScat `KineticEnergy.cpp:95`,
independently re-derived here) for every grid spec tested.

**V3 -- analytic separable benchmarks at D = 1, 2, 3.** Generality is
*exercised*, not merely asserted: every benchmark grid uses **unequal**
per-axis extents, element counts, quadrature orders, and (for the
oscillator) frequencies, so a transposed-axis bug cannot hide behind two
axes that happen to look alike.

- **D-dimensional particle in a box**, eigenvalues
  `sum_d n_d^2 pi^2 / (2 m_d L_d^2)`, checked via `kinetic_nd` alone (no
  potential). Measured relative error against the analytic sum, at the basis
  sizes actually used: **D=1: 3.4e-9**, **D=2: 3.8e-12**, **D=3: 5.4e-5**.
  The D=3 figure looks like an outlier next to D=1/D=2, but it is not a
  method problem: D=3 deliberately uses a coarse 2-element, 6-point-per-axis
  grid to keep a *dense* eigensolve of the full `N = n_0*n_1*n_2` matrix
  tractable in a routine test run. A convergence study (quadrature order
  q = 6 -> 8 -> 10, same element count) confirms the error is
  **grid-limited, not method-limited**: it falls from 5.4e-5 to 4.1e-8 to
  1.2e-11 -- the textbook exponential convergence of a correct spectral
  discretization, exactly the same signature `docs/physics/femdvr-ecs.md`'s
  Benchmark 1 uses to certify the underlying 1-D kinetic assembly.
- **D-dimensional harmonic oscillator**, eigenvalues
  `sum_d omega_d (n_d + 1/2)`, checked via `hamiltonian_nd` (kinetic plus a
  quadratic diagonal potential). Measured relative error: **D=1: 6.7e-9**,
  **D=2: 6.2e-4**, **D=3: 8.3e-4**. Unlike the box, these do not fall toward
  round-off as the basis grows arbitrarily, because the oscillator's
  Gaussian-tailed eigenfunctions are never exactly compactly supported: any
  finite hard-wall box leaves a residual truncation error, so D=2/D=3 sit at
  a few times `1e-4`-`1e-3`, consistent with the modest box half-widths and
  basis sizes used to keep the D=3 eigensolve (sparse shift-invert) fast.

**V4 -- D = 1 reproduces the existing 1-D stack bit-for-bit.**
`hamiltonian_nd(TensorGrid([g]), [mass], V)` and `kinetic_nd(TensorGrid([g]),
[mass])` are checked against the pre-existing dense
`qscat.dvr.operators.hamiltonian(g, V, mass)` and `qscat.dvr.kinetic.kinetic
(g, mass)` respectively, on an ECS grid, and match to **exactly 0.0** --
not merely within round-off tolerance, but bit-identical, because the D=1
tensor-product code path reduces algebraically to the same computation as
the pre-existing dense one. The practical consequence: this makes every
cross-section result already validated in sub-projects #1-#4 (which all sit
on top of the 1-D dense stack) a standing regression test on the new
N-dimensional code -- if a future change to `kron_sum`, `TensorGrid`, or
`kinetic_sparse` ever perturbs the D=1 case even slightly, this test catches
it immediately, long before it could reach any N2 physics.

**V5 -- `SparseLU` correctness and reuse.** Residual `norm(A@x - b) /
norm(b) < 1e-12` on a complex-symmetric ECS-derived matrix; a multi-RHS solve
matches looped single solves; one factorization reused across many
right-hand sides gives identical results to re-factorizing for each (as it
must, since nothing about the matrix changes between solves). A fixed-seed
300x300 complex-symmetric test matrix was also used to measure how
`fill_factor` (`(L.nnz + U.nnz) / A.nnz`) varies with the `permc_spec`
ordering scipy's `splu` is given: **NATURAL 18.20, COLAMD 17.60,
MMD_AT_PLUS_A 8.79** -- MMD_AT_PLUS_A roughly halves the fill relative to
the other two on this matrix. This is genuinely promising evidence that
`MMD_AT_PLUS_A` -- appropriate in principle for the structurally symmetric
sparsity pattern a Kronecker-sum Hamiltonian has -- is worth trying at
production scale, where the default COLAMD ordering measured a much larger
x93 fill-in (below). It is **not** a proven production result: it has only
ever been measured on one small, randomly generated 300x300 matrix, never on
the real N2 Hamiltonian, whose sparsity pattern (a 2-D FEM-DVR-ECS tensor
product with a ~23-nonzero-per-row band structure) is structurally quite
different from a dense random complex-symmetric matrix. Sub-project #6, which
actually needs to factorize at production scale, should measure `fill_factor`
on the real matrix before choosing an ordering, not assume this ratio
transfers.

**V6 -- production-scale smoke test** (`test_nd_scale.py`, marked `slow` and
excluded from the default `pytest` run). Assembles a `TensorGrid` at the
exact grid dimensions of the real eMoScat N2 2-D deck
(`reference/eMoScat/input/experimental/N2-model.json`) -- electronic grid
`q=8`, 33 real + 15 ECS elements, `n=335`; nuclear grid `q=14`, 23 real + 10
ECS elements, `n=428`; `N = 335*428 = 143,380` -- but with a generic
analytic potential (`1/(1+r^2) + 1/(1+R^2)`) rather than the real N2
potential surface, since `libs/qscat` must not import from `validation/` or
`projects/` and the N2-specific assembly is sub-project #6's job. It checks
matrix dimension, `nnz == 3,276,450` (matching eMoScat's own nonzero-count
formula, independently re-derived), and `max|H - H^T| < 1e-10`
(measured: `1.7e-13`) -- all **without factorizing** (see below for why).
Assembly alone costs about 0.1 s.

## Measured cost at production scale

Assembly is cheap; factorization is not. On the production N2 2-D grid
(`N = 143,380`, `nnz = 3,276,450`, 22.9 nonzeros/row):

| quantity | value |
|---|---|
| electronic grid | `q=8`, 33 real + 15 ECS elements -> `n = 335` |
| nuclear grid | `q=14`, 23 real + 10 ECS elements -> `n = 428` |
| N | 143,380 |
| nnz | 3,276,450 (22.9/row) |
| `max\|H - H^T\|` | 1.7e-13 |
| assembly time | ~0.1 s |
| `splu` factorization (COLAMD) | **128 s**, fill-in **x93**, L+U ~ 3.05e8 nnz |
| peak RSS during factorization | **13.6 GB** |
| back-substitution (per solve) | **440 ms** |

This is why V6 assembles but never factorizes: 128 s and 13.6 GB have no
place in a routine test run, but the assembly itself -- the part this
library is actually responsible for -- is fast and needs to be checked at
real scale, not just on toy grids.

The consequence that matters for sub-project #6's design: **the
factorization, not the back-substitution, is the dominant cost, and one
factorization serves every right-hand side at a fixed total energy.** In a
driven-equation (T-matrix) scattering calculation, all outgoing/final
channels at a single collision energy share the same `(H - E)` matrix, so
`SparseLU` needs to factor once per energy point and then back-substitute
once (440 ms) per channel -- cost scales with the number of *energies*
sampled, not with the number of vibrational channels computed at each
energy. A production run covering many energies is therefore dominated by
`number_of_energies * 128 s`, which is the actual planning constraint
sub-project #6 has to design around (e.g. subsampling the energy grid, or
tuning `SparseLU`'s `ordering` per the V5 measurement above), not the
per-channel back-substitution cost.

## Out of scope

Everything N2-specific -- the potential surface, vibrational channel
functions, the driven-equation T-matrix, and the cross section itself --
belongs to sub-project #6, not here. Also out of scope for this layer:
non-diagonal potentials, coupled partial waves, iterative/matrix-free
solvers (direct sparse LU is demonstrably sufficient at the sizes measured
above), sparse time propagation (`qscat.evolution.make_cn_stepper` remains
dense-only), and GPU/CUDA (deferred repo-wide; the port-scout confirmed
eMoScat's own 2-D path was never GPU either).

## See also

- `docs/superpowers/specs/2026-07-22-nd-sparse-hamiltonian-design.md` -- the
  design spec this sub-project implements.
- `.superpowers/sdd/n2-2d-exact-extraction.md` -- the port-scout archaeology
  of eMoScat's 2-D Hamiltonian assembly and production-scale measurements
  this document's cost table is drawn from.
- `docs/physics/femdvr-ecs.md` -- the 1-D FEM-DVR-ECS grid, kinetic
  assembly, and diagonal-potential approximation that each `T_d` and the
  D-dimensional `V` build on.
