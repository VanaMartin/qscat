# Sparse shift-invert eigensolver

`qscat.linalg.ShiftInvertEigs` — the `k` eigenvalues and eigenvectors of a large
sparse complex-symmetric matrix nearest a complex shift `sigma`.

Validated in **1-D only**. Nothing here is a claim about 2-D grids.

## Why shift-invert

A resonance is an **interior** eigenvalue: it sits in the lower half of the
complex energy plane, surrounded by the discretized rotated continuum, at
neither end of the spectrum. Krylov methods (Arnoldi/Lanczos) converge to
**extremal** eigenvalues, so applied to `H` directly they return the edges of the
spectrum and never the pole.

The standard remedy is the shift-invert spectral transform. Run Arnoldi on

```
OP = (A - sigma*I)^-1
```

whose extremal eigenvalues `1/(E - sigma)` correspond to `A`'s eigenvalues `E`
**nearest** `sigma`. ARPACK returns Ritz values of `OP`; SciPy converts them
back, so what comes out are eigenvalues of `A`.

## Why it is cheap here

The transform costs one sparse solve per matrix-vector product — which is
exactly `qscat.linalg.SparseLU`'s job, so this inherits the complex-symmetric
MUMPS `SYM=2` backend and its SuperLU fallback for free
(`docs/physics/mumps-sparse-backend.md`).

More usefully: `sigma*I` touches only the diagonal, so **`A - sigma*I` has the
same sparsity pattern for every shift**. `SparseLU.refactor` therefore applies —
on the MUMPS backend it reuses the symbolic analysis and skips the SCOTCH
ordering; on the scipy backend it re-runs `splu` (correct, but with no reuse).
A resonance hunt is a sweep of shifts through the complex plane, so it gets the
same discount the time-independent energy sweep gets from the identical trick
(`docs/physics/ti-energy-sweep-reuse.md`).

`ShiftInvertEigs` is a class rather than a function precisely so that the
factorization survives between shifts: the first `near(sigma)` factors, every
later one refactors.

## Conventions

Three choices that are easy to get wrong, and are pinned by tests.

**The shift sign is `A - sigma*I`.** SciPy's `OPinv` must solve
`(A - sigma*I) x = b` — *not* the driven solver's `E*I - H`. The wrong sign does
not raise; it returns eigenvalues **reflected about the shift**. Measured on a
diagonal spectrum `[0, 1, 2, 3, 10, 11, 12]` with `sigma = 9`:

| `OPinv` built from | returned |
|---|---|
| `A - sigma*I` (correct) | `10, 11` |
| `sigma*I - A` (wrong) | `7, 8` |

Both look entirely plausible. `test_shift_sign_convention_is_A_minus_sigma_I`
exists for this reason, and was itself verified to fail under the wrong sign.

**Eigenvalues come back sorted by `|E - sigma|`, nearest first** — not by
`Re E`. A shift-invert result is a local window around the shift, so distance
from the shift is its meaningful order. `qscat.dvr.eigen` returns a *whole*
spectrum and keeps its ascending-`Re E` order; neither convention should be
changed to match the other.

**Eigenvectors come back Euclidean-normalized** (`v†v = 1`), exactly as
`qscat.dvr.eigen` returns them. ECS observables need the bilinear
`qscat.linalg.c_product` normalization instead, and the region to normalize over
is the caller's decision (the LCP code, for instance, normalizes over the real
region only), so this class does not presume it.

## Validation

The oracle is dense `np.linalg.eig` on the same matrix.

- **Synthetic**: sparse complex-symmetric matrices (`n = 150`–`300`).
  Eigenvalues match to `rtol = 1e-9`; eigenvectors match to
  `|vᵀw| = 1` within `1e-6` after normalizing both to unit c-norm `sqrt(vᵀv)`
  (the right notion of "equal up to scale" for a complex-symmetric operator).
- **Physical**: the N₂ electronic FEM-DVR-ECS Hamiltonian at fixed `R = 2.02`
  (`T + diag(V_surface)`, `n = 113`, the same build `qscat.core.lcp` uses).
  Every eigenvalue the sparse solver returns at a given ECS angle is one the
  dense solver returns, to the same `rtol = 1e-9`; and the two-angle pole built
  from sparse spectra equals the one built from dense spectra to that tolerance.

`rtol = 1e-9` rather than something tighter is deliberate: the comparison runs
through a sparse factorization, and pinning sparse-solve agreement tighter has
failed CI on a different BLAS.

**The pole found this way is the physical one.** At `R = 2.02` the sparse
two-angle pole is `E = -0.661315 - 0.008333j` Ha, i.e.

| quantity | this solver | `docs/physics/n2-resonance.md` |
|---|---|---|
| `E_res` (relative to `v0`) | 2.441 eV | 2.445 eV |
| `Gamma` | 0.4535 eV | 0.455 eV |

Note the subtraction: an ECS eigenvalue is **absolute**, carrying
`v0(R) = -0.751 Ha`, while the literature quotes `E_res` measured from the
neutral curve. A shift seeded at the literature value rather than at
`v0 + E_res` finds nothing.

## Measured working range

On that N₂ electronic Hamiltonian (`n = 113`), seeding at
`sigma = E_pole + offset·(1+i)` and asking whether the pole appears among the
`k` returned:

| offset (Ha) | k=2 | k=4 | k=6 | k=8 | k=16 |
|---|---|---|---|---|---|
| 0.001 | found | found | found | found | found |
| 0.01 | found | found | found | found | found |
| 0.05 | found | found | found | found | found |
| 0.1 | found | found | found | found | found |
| 0.2 | found | found | found | found | found |
| 0.5 | **raises** | found | found | found | found |
| 1.0 | **raises** | **raises** | **raises** | found | found |
| 2.0 | absent | absent | absent | found | found |

Three things this says:

1. The seed is forgiving. A shift `0.2 Ha` away — an order of magnitude larger
   than `Gamma = 0.0167 Ha` — still finds the pole at any `k >= 2`. A BO/LCP
   level is a far better guess than that, which is what makes the 2-D plan's
   "seed from `resonance_levels`" step credible.
2. When it fails, it **raises**. The failure mode at a distant shift with a small
   Krylov space is `ConvergenceError` ("ARPACK did not converge … 1 of 2
   eigenvalues converged"), not a plausible wrong answer. The one silent outcome
   — `offset = 2.0`, small `k` — is not an error at all: the pole genuinely is
   not among the `k` nearest eigenvalues to a shift that far away.
3. `k = 8` was reliable everywhere tested. `k = 2` is not worth using.

**Cost note.** On the scipy backend a second `near()` is not faster than the
first (measured 1.9 ms vs 2.3 ms at `n = 113` — noise at this size), because
`SparseLU.refactor` re-runs `splu` there. The reuse is a MUMPS-backend property
and it is the large-matrix case that will show it; no speedup should be claimed
from these 1-D numbers.

## Limits

- **1-D only.** The primitive has not been run on a 2-D tensor Hamiltonian. The
  N₂ electronic matrix is `n = 113`; the N₂ 2-D working deck is ~143k and H₂⁺ is
  ~1.15M, where the factorization, not ARPACK, is expected to dominate.
- **Selection is not solved here.** Within a single shift-invert window several
  rotated-continuum eigenvalues are *narrower* in `|Im E|` than the pole itself,
  so "the narrowest state" is not a selector. Only the two-angle criterion
  (`qscat.ecs.find_resonance_pole` / `match_angle_stable`) separates them, and
  the tests here use it exactly that way. A shift parked in the rotated continuum
  at `sigma = -0.30 - 0.40j` yields **no** angle-stable state, which is the
  correct and useful negative.
- **Near-singularity was not observed to be a problem.** `A - sigma*I` is
  ill-conditioned by construction — that is the amplification mechanism — and at
  every offset tested down to `0.001 Ha` the results matched dense to `1e-9`.
  MUMPS `SYM=2` behaviour at a shift pathologically close to an eigenvalue has
  not been probed; these measurements ran on the scipy/SuperLU backend.
