# N₂ Time-Independent VE Cross-Section Solver — Design Spec (sub-project #3)

**Date:** 2026-07-21
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — implementation pending
**Lifecycle:** `qm-method-lifecycle` stage 1. Sub-project #3; builds on #1 (`qscat.dvr`) and
#2 (`qscat.ecs` resonance pole → V_d(R)/Γ(R)).

## Context

Compute the electron–N₂ vibrational-excitation (VE) cross sections σ_{0→v'}(E) via the
**time-independent Local Complex Potential (LCP)** method — a resolvent/driven-equation solve
on the nuclear FEM-DVR-ECS grid — and compare them, at a handful of anchor points, to
**Karel Houfek's `CSVE.V00.J00` golden data**. This flips the 6 PENDING **C5** anchors in the
N₂ benchmark harness (`validation/n2/`) to real checks. Full extraction (formulas, signs,
normalization) in `.superpowers/sdd/ti-cross-section-extraction.md`.

**Model-difference caveat (important).** eMoScat's LCP has no time-independent path; the
formulas below are port-scout's exact Laplace-transform equivalent of the *implemented*
time-dependent LCP. Houfek's data is from a time-independent **2-D** model (explicit electron
+ nuclear coordinates), a *different approximation* than our 1-D LCP. So agreement is expected
to be **"quite close" but not exact** — the C5 tolerance is loose and documented accordingly
(see Validation).

## Method (port-scout-verified formulas, atomic units, bohr²)

- **Neutral vibrational states:** diagonalize `T_nuc(μ) + diag(V₀(R))` → (ε_v, χ_v). μ=12766.36.
- **Doorway:** `d_v(R) = √(Γ(R)/2π)·χ_v(R)`.
- **Driven equation** (one `np.linalg.solve` per collision energy E, initial channel v):
  `[E_tot − T_nuc(μ) − diag(V_d(R)) + i·diag(Γ(R))/2]·ξ = d_v`, with `E_tot = E + ε_v`,
  `V_d(R) = V₀(R) + E_res(R)` from sub-project #2.
- **VE S-matrix & cross section:** `S_{v'←v}(E) = ⟨d_{v'}|ξ⟩` (dot product in the
  1/√w-normalized DVR basis); `σ_{v→v'}(E) = 4π³·|S|²/(2E)`. Open-channel guard: `E_tot − ε_{v'} > 0`, else 0.
- (DA channel formula is available in the extraction but **out of scope** this iteration —
  anchors are VE only.)
- **V_d(R)/Γ(R) source: recompute per nuclear-R** (user choice) — call the sub-project #2
  two-angle pole finder at each nuclear grid R. For the **complex ECS-tail** R points, the
  electronic potential is analytic in R (and Γ→0 beyond the ≈2.4 bohr crossing), so V_d/Γ
  continue cleanly; a smooth analytic-continuation fallback is used if direct complex-R
  pole-finding is numerically delicate.

## Interface

```
projects/n2_ti_cross_section/
├── nuclear_grid.py   n2_nuclear_grid() -> FemDvrEcsGrid  (dvr_order 14, real 0→12 bohr + 35° tail)
├── vibrational.py    vibrational_states(grid, mu, n) -> (eps[n], chi[n, :])   # T_nuc(mu)+V0
├── vres.py           vres_on_grid(grid) -> (Vd[nb], Gamma[nb])  # recompute pole per nuclear-R
├── cross_section.py  ve_cross_section(grid, mu, Vd, Gamma, eps, chi, v_init, vprimes, E) -> sigma[...]
│                       (builds H_res, solves driven eqn, forms S, returns sigma in bohr²)
└── tests/            validation (below)
```

Reuses `qscat.dvr` (grid, kinetic → T_nuc(μ)), `qscat.ecs`, and the sub-project #2 pole finder
(`projects/n2_resonance/`). Promotion of any reusable piece (e.g. a generic
`resolvent_cross_section`) to `qscat` decided at the end.

## Validation

**Part 1 — internal correctness (Houfek-independent; the real bug gate):**
- Neutral vibrational spacing `ε₁ − ε₀ ≈ N₂ ω_e` (~0.29 eV / ~0.0107 Ha) within a few %.
  - **Note (superseded):** this real-N₂ target was superseded by a maintainer
    decision to accept eMoScat's model potential as-is (see
    `.superpowers/sdd/task1fix-report.md`). The shipped Part-1 check instead
    validates the FEM-DVR solver against the **analytic Morse spectrum of
    eMoScat's own potential** (matched to ~1e-14), not against real N₂. Under
    eMoScat's model, `ε₁−ε₀ ≈ 0.0124` Ha — ~16% above real N₂'s `ω_e` — which
    is a documented property of the model potential, not a solver error; the
    criterion as written above is not what the solver is actually checked
    against.
- σ ≥ 0 for all computed anchors; channel-opening thresholds correct (`σ=0` when `E_tot−ε_{v'}<0`).
- Resonance enhancement: the VE cross section is largest in the ~2–3 eV region (the ²Π_g resonance).

**Part 2 — cross-model comparison (the C5 flip):**
- Compute σ_{0→v'}(E) at the 6 C5 anchor (E, v') coordinates.
- Compare to Houfek's `CSVE.V00.J00` values. **Report the actual per-anchor agreement**, then
  set the C5 pass tolerance to a **documented** bound reflecting the LCP-vs-2D model difference
  (expected: same order of magnitude + correct trend; the exact factor finalized after the first
  run, e.g. within a factor of ~2, and clearly labeled a cross-*model* check — NOT exact
  agreement, NOT 5%).
- **Harness wiring:** update `validation/n2/experiment.py` C5 from PENDING to a real PASS/FAIL
  at the documented tolerance, with the detail string noting it's an LCP-vs-2D comparison.
  Must still run in the CPU docker; degrade to a labeled FAIL (not a crash) on solver error.

## Out of scope

- The DA channel; full σ(E) curves over all channels/energies (anchor scope only).
- The time-dependent solver (a separate future sub-project; would cross-check LCP internally).
- Reproducing Houfek's 2-D model itself.
- Rust (dense `np.linalg.solve` is the correct end state for a toy).

## Verification

- `uv run pytest projects/n2_ti_cross_section` green (Part-1 internal checks + the anchor
  computation).
- Vibrational spacing ≈ N₂ ω_e; σ≥0; thresholds correct; resonance in ~2–3 eV.
- The 6 C5 anchors compared to Houfek's data, agreement reported, tolerance documented; N₂
  harness C5 flips PENDING → PASS at that tolerance; `validation/n2/experiment.py` exits 0;
  docker run still exits 0.
- Uses `qscat.dvr`/`qscat.ecs`; CPU-only.
