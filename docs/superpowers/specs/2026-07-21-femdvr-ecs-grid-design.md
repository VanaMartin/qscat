# FEM-DVR-ECS Radial Grid + Kinetic Operator — Design Spec (sub-project #1)

**Date:** 2026-07-21
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — implementation pending
**Lifecycle:** `qm-method-lifecycle` stage 1 (Design). Sub-project #1 of the electron–N₂
resonance port (sub-project #2 = the resonance θ-stabilization solver, built on this).

## Context

This is qModeling's first foundational numerical capability: a **finite-element
Discrete-Variable-Representation radial grid with Exterior Complex Scaling**
(Rescigno–McCurdy FEM-DVR-ECS) and its kinetic-energy operator. Every downstream method
reuses it — the N₂ resonance pole solver next, then the time-independent and
time-dependent cross-section solvers. It is ported from `reference/eMoScat`
(`FemDvrEcsGrid`, `KineticEnergy`, `DiscreteStates`); the full extraction is in
`.superpowers/sdd/femdvr-ecs-extraction.md` (to be promoted to `docs/physics/`).

It is validated **entirely on analytic benchmarks independent of N₂**, so that when the
N₂ resonance (sub-project #2) later looks right, the grid underneath is already known
trustworthy.

## Method

- **FEM-DVR:** the radial domain is split into elements; within each, Gauss-Lobatto
  quadrature nodes/weights define a DVR basis. Adjacent elements share their boundary node
  via a **bridge function** whose global weight is the sum of the two elements' local
  Lobatto weights at that node (C⁰ continuity, L² orthonormal). The two outermost points
  are dropped (Dirichlet ψ=0).
- **Kinetic operator:** `T = −(1/2m) d²/dz²`, assembled per element as
  `T_ij = (1/2m) Σ_l wze[l]·dBF(i,l)·dBF(j,l)` with `dBF = dLp/hz` and each basis function
  normalized by `1/√(bridge-summed global weight)`. Bridge-corner diagonal terms accumulate
  across the shared node.
- **Diagonal-potential DVR:** `V_ij ≈ V(x_i)·δ_ij` (exact for the orthonormal Lobatto basis;
  documented as an assumption, not a bug).
- **Exterior Complex Scaling:** `z(x) = x` for `x ≤ R0`, `z(x) = R0 + (x−R0)·e^{iθ}` beyond.
  The pivot R0 is an element boundary (= sum of the real-region element lengths); elements
  past it have length scaled by `e^{iθ}` (θ in degrees). The Jacobian enters through
  `hz = ½·(complex element length)`.
- **Eigenproblem:** `H = T + diag(V)` is complex-**symmetric but non-Hermitian**; solved as
  a **standard** eigenproblem (DVR basis is orthonormal — no metric) via `np.linalg.eig`,
  eigenvalues sorted by ascending real part. (numpy normalizes eigenvectors to `v†v=1`; the
  ECS c-product convention is `vᵀv=1` — documented; re-normalization only needed if
  bit-matching the Fortran, not for eigenvalues.)

Units: atomic (Hartree, bohr). Electron mass = 1.

## Interface

```
projects/femdvr_ecs/
├── spec.py        ElementSpec(length: float, angle_deg: float = 0.0)
│                  GridSpec(quadrature: int, elements: list[ElementSpec], x_min: float = 0.0)
│                    - complex elements (angle_deg != 0) must be contiguous at the END (ECS tail);
│                      validated on construction. Pivot R0 = sum of real-region lengths.
├── grid.py        FemDvrEcsGrid(spec) ->
│                    .points   (complex ndarray, shape (nb,))  — ECS-scaled DVR points
│                    .weights  (complex ndarray, (nb,))        — bridge-summed weights
│                    .real_points (real ndarray, (nb,))        — unscaled coords (bookkeeping)
│                    .R0 (float), .n (int = nb)
├── kinetic.py     kinetic(grid, mass: float) -> complex ndarray (nb, nb)
├── operators.py   hamiltonian(grid, V, mass) -> ndarray   # V: callable V(z) or ndarray at grid.points
│                  eigen(H) -> (E: complex ndarray sorted by Re, vecs: ndarray)
└── tests/         analytic benchmarks (below)
```

Gauss-Lobatto nodes/weights come from `numpy.polynomial.legendre` (no scipy dependency;
not the reference's hand-rolled QL solver).

## Validation (the `numerical-validation` gate — all N₂-independent)

Do not promote to `qscat` until all pass.

- **B1 — Particle in a box (θ=0).** Real elements spanning `[0, L]`, `V=0`, `m=1`. Exact
  `E_n = n²π²/(2 L²)`. Check n=1…5 to `rtol ≤ 1e-6` at sufficient Lobatto order, and confirm
  **spectral convergence** (error falls fast as `quadrature` rises). Isolates kinetic
  assembly, bridge coupling, and Dirichlet trimming (no potential, no ECS).
- **B2 — Harmonic oscillator (θ=0).** `V = ½ m ω² (x−x_c)²`, box wide enough that the low
  levels are unaffected by the walls, `x_c` centered. Exact `E_n = ω(n+½)`. Check n=0…4 to
  `rtol ≤ 1e-6`. Exercises the diagonal-potential DVR and the `1/√w` normalization (a
  bridge-weight bug shows as a boundary-localized energy shift).
- **B3 — ECS continuum rotation (θ>0, V=0).** Free particle on a grid with a real region
  `[0, R0]` plus a complex tail at angle θ. The discretized continuum eigenvalues lie along
  the ray `arg(E) ≈ −2θ`. Check that mid-spectrum eigenvalues satisfy `|arg(E) + 2θ|` within
  a few degrees. A clean, oracle-free validation of the ECS Jacobian / `e^{−2iθ}` kinetic
  scaling.
- **B4 — Bound-state θ-independence.** Square well `V = −V0` on `[0, a]` (0 beyond), deep
  enough for a bound level `E_b < 0`. Diagonalize on grids with two different θ; `E_b` must
  agree to `rtol ≤ 1e-6` (bound states are ECS-invariant). Confirms ECS does not perturb the
  discrete spectrum.

(The **resonance-pole** validation — a model barrier shape resonance vs an independent
transcendental solve, plus θ-stabilization — is deliberately deferred to **sub-project #2**,
which builds the stabilization scan that `eMoScat` does not automate.)

## Promotion target

On validation, promote:
- `FemDvrEcsGrid` + `kinetic` + `hamiltonian`/`eigen` → **`libs/qscat/qscat/dvr/`**.
- ECS contour/transform helpers (the `z(x)` map, angle handling) → **`libs/qscat/qscat/ecs/`**.
Bring the benchmark tests along; update `CLAUDE.md`/`docs` to describe the new capability.
The pure-Python implementation is the permanent differential-test oracle (no Rust unless
profiling later shows a real hot path — out of scope here).

## Out of scope (this sub-project)

- The N₂ resonance θ-stabilization scan and V_d(R)/Γ(R) generation (sub-project #2).
- The eMoScat JSON grid-spec adapter (`parse_points` 10%-merge, `uniform_increment` growth) —
  added in sub-project #2 where matching N₂'s exact element boundaries matters.
- Sparse/RowCompressed kinetic; left-side ECS (`complex_negative`); the graded-angle
  "special" path.
- Any Rust kernel (dense numpy is the correct end state for this layer).

## Verification

- `uv run pytest projects/femdvr_ecs` green (B1–B4).
- B1/B2 converge to `rtol ≤ 1e-6` and show spectral convergence with Lobatto order.
- B3 eigenvalues align with `arg(E) ≈ −2θ`; B4 bound level θ-invariant to `rtol ≤ 1e-6`.
- Runs on CPU; no N₂ data or physics involved.
