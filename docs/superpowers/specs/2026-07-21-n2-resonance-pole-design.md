# N₂ Electronic Resonance Pole Solver — Design Spec (sub-project #2)

**Date:** 2026-07-21
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — implementation pending
**Lifecycle:** `qm-method-lifecycle` stage 1. Sub-project #2 of the electron–N₂ resonance
port; builds on sub-project #1 (`qscat.dvr` FEM-DVR-ECS grid, merged).

## Context

Compute the electron–N₂ ²Π_g resonance curves **E_res(R), Γ(R), V_d(R)** by finding the
θ-stabilized ECS complex-eigenvalue pole of the fixed-R electronic Hamiltonian, on the
`qscat.dvr` grid. This produces the inputs the cross-section solvers (sub-project #3+) need,
and flips the N₂ benchmark harness **B1** check (`validation/n2/`) from PENDING to a real
pass/fail. Algorithm mapped by port-scout (`.superpowers/sdd/n2-lcp-model-extraction.md`,
`make_vres`); note eMoScat does **not** automate the pole identification — we design it,
using **two-angle matching**.

## Method

- **Electronic Hamiltonian at fixed R** (electron mass = 1, atomic units):
  `H_el(R) = kinetic(grid, 1) + diag(V_eff_el(grid.points, R))`, where
  `V_eff_el(r,R) = l(l+1)/(2 r²) + V_int(r,R)`, `l = 2`, `V_int(r,R) = −λ(R)·exp(−α_c r²)`,
  `λ(R)` and the N₂ parameters exactly as validated in `validation/n2/model.py`
  (single source of the N₂ potential; params in `validation/n2/config.json`).
- **Pole finder — two-angle matching:** diagonalize `H_el(R)` on two ECS grids that differ
  ONLY in scaling angle (θ_a=35°, θ_b=44°). Continuum eigenvalues rotate with θ; the
  resonance pole does not. The pole is the eigenvalue lying at ~the same complex energy in
  both spectra — i.e. `argmin |E_a − E_b|` over pairs within a resonance-energy search
  `window` (a box in the complex plane around the expected pole). Return the matched value
  (mean of the two, or θ_a's) plus the matching residual `|E_a − E_b|` as a stability
  diagnostic.
- **R-scan with continuation:** for each nuclear R in a grid, seed the search `window` from
  the previous R's pole (start R₀ seeded from the known resonance region ≈ 0.08–0.10 Ha).
  `E_res(R) = Re(E_pole)`, `Γ(R) = max(0, −2·Im(E_pole))`, `V_d(R) = V₀(R) + E_res(R)`.
- **Simplification (documented):** skip eMoScat's electron-affinity bootstrap — seed the R₀
  window directly from the expected resonance region; rely on continuation for the scan.

## Interface

```
projects/n2_resonance/
├── potential.py   v0(R), lam(R), v_int(r,R), v_eff_el(r,R)  — N₂ potentials, params from
│                  validation/n2/config.json (matches validation/n2/model.py exactly)
├── grid_n2.py     n2_electronic_grid(angle_deg, *, r_pivot, r_max, ...) -> FemDvrEcsGrid
│                    hand-built GridSpec: dvr_order 8, real region to r_pivot + complex tail
│                    at angle_deg (approximating N2.json's electronic grid).
├── pole.py        electronic_hamiltonian(R, grid) -> ndarray
│                  find_pole(R, grid_a, grid_b, window) -> (E_pole: complex, residual: float)
│                  resonance_curve(R_grid, grid_a, grid_b) -> (E_res[R], Gamma[R], V_d[R])
└── tests/         validation (below)
```

The **general** two-angle pole finder (angle-independent of N₂) is promoted to
`qscat.ecs.find_resonance_pole(eigs_a, eigs_b, window)`; the N₂ potential/grid stay in the
project.

## Validation (`numerical-validation` gate)

- **V1 (key) — resonance at equilibrium.** `E_res(R₀=2.01943)` ∈ **[2.3, 2.5] eV** and
  `Γ(R₀)` ∈ **[0.35, 0.55] eV** (port-scout prototype ≈2.44/0.46 eV; literature 2.3–2.4 eV /
  0.4–0.5 eV). Converted to Hartree in the test.
- **V2 — pole is genuinely θ-stable & converged.** The matching residual `|E_a − E_b|` at R₀
  is small (≪ Γ, e.g. < 1e-3 Ha), and `E_pole(R₀)` is stable (to a few %) under a change in
  grid resolution (Lobatto order / element count) — a real resonance, not a discretization
  artifact.
- **V3 — smooth curves.** `E_res(R)`, `Γ(R)` are smooth over R ∈ [≈1.6, 3.0] bohr (no
  mode-hopping discontinuities); `Γ(R) ≥ 0`.
- **Harness wiring.** Update `validation/n2/experiment.py` B1 to compute `E_res(R₀)` via this
  solver and assert the literature window — flipping B1 from PENDING to PASS (or FAIL if
  wrong). The N₂ harness must still run in the CPU docker.

## Promotion target

On validation: `qscat.ecs.find_resonance_pole` (general two-angle matcher). The N₂-specific
`potential.py`/`grid_n2.py`/`resonance_curve` stay in `projects/n2_resonance/` for now (they
become part of the N₂ model package when the cross-section solvers land). Pure Python; no
Rust (out of scope). Update `CLAUDE.md`/`docs` for the promoted `qscat.ecs` addition.

## Out of scope (this sub-project)

- The time-independent and time-dependent **cross-section** solvers (sub-project #3+) — this
  produces V_d(R)/Γ(R) only.
- The electron-affinity bootstrap; the eMoScat JSON grid-spec adapter; the graded-angle path.
- Reproducing eMoScat's exact production element boundaries (a converged hand-built grid
  gives the correct physical pole regardless).

## Verification

- `uv run pytest projects/n2_resonance` green (V1–V3).
- `E_res(R₀)` in [2.3, 2.5] eV, `Γ(R₀)` in [0.35, 0.55] eV; residual small; curves smooth.
- The N₂ harness B1 flips PENDING → PASS: `uv run python validation/n2/experiment.py` shows
  B1 PASS; `docker run --rm qmodeling:runtime python validation/n2/experiment.py` still exits 0.
- Runs on CPU; uses `qscat.dvr`/`qscat.ecs`.
