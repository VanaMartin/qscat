# N₂ Time-Dependent VE Cross-Section Solver — Design Spec (sub-project #4)

**Date:** 2026-07-22
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — implementation pending
**Lifecycle:** `qm-method-lifecycle` stage 1. Sub-project #4; builds on #1 (`qscat.dvr`),
#2 (`qscat.ecs` pole → V_d/Γ), #3 (the TI cross section, which is the differential oracle here).

## Context

Compute the electron–N₂ VE cross sections σ_{0→v'}(E) by **time-dependent wavepacket
propagation** (Crank-Nicolson) of the Local Complex Potential model, and cross-check them
against the already-validated **time-independent** result (sub-project #3). This flips the last
PENDING check — Group **D** — in the N₂ benchmark harness. Time-dependent extraction (doorway,
propagation, energy transform) is in `.superpowers/sdd/n2-lcp-model-extraction.md`.

**Why this is lower-risk than #3:** the TD and TI methods are Laplace/Fourier transforms of the
same physics. With the `(1/i)` normalization,
`S_TD(E) = (1/i)∫₀^∞ e^{i E_tot t}⟨d_{v'}|e^{−iH_res t}|d_v⟩ dt = ⟨d_{v'}|(E_tot−H_res)⁻¹|d_v⟩ = S_TI(E)`.
So **TD must converge to TI** (up to time-step/finite-time discretization) — an *exact*
differential oracle, not just a loose cross-model comparison.

## Method (atomic units, bohr²)

Reuses #3's setup verbatim: the nuclear FEM-DVR-ECS grid, vibrational states (ε_v, χ_v),
V_d(R)/Γ(R) (`vres_on_grid`), and the doorways `d_v = √(Γ/2π)·χ_v`. Then:

- **Initial wavepacket:** `ψ(0) = d_{v_init}`.
- **Crank-Nicolson propagation** under the time-independent, non-Hermitian
  `H_res = T_nuc(μ) + diag(V_d − iΓ/2)`:
  `(I + iH_res·dt/2)·ψ_{n+1} = (I − iH_res·dt/2)·ψ_n`.
  H_res is constant ⇒ LU-factor `A = I + iH_res·dt/2` ONCE, reuse each step (one solve/step).
  The `−iΓ/2` term makes ‖ψ‖ decay (the resonance depletes); propagate until the norm is
  negligible.
- **Correlation & energy transform:** accumulate `c_{v'}(t_n) = ⟨d_{v'}|ψ(t_n)⟩` (c-product, no
  conjugate — same convention as #3); then
  `S_{v'}(E) = (1/i)·Σ_n w_n·e^{i(E+ε_{v_init})t_n}·c_{v'}(t_n)·dt` (Simpson/Filon weights `w_n`),
  and `σ_{v→v'}(E) = 4π³·|S|²/(2E)`, `0` if `E_tot − ε_{v'} ≤ 0`.

## Interface

```
projects/n2_td_cross_section/
├── propagator.py       make_cn_stepper(H, dt) -> stepper(psi) -> psi_next   (LU precomputed once)
├── td_cross_section.py td_ve_cross_section(grid, mu, Vd, Gamma, eps, chi, v_init, vprimes,
│                          E, *, dt, n_steps) -> sigma[...]  (propagate + correlate + transform)
└── tests/              validation (below)
```

Reuses `projects.n2_ti_cross_section.{nuclear_grid, vibrational, vres, cross_section}` (the last
for the TI oracle) and `qscat.dvr` (T_nuc). The generic **Crank-Nicolson propagator** promotes to
`qscat.evolution` (the currently-empty `evolution` subpackage — its natural home).

## Validation

**V1 (key) — TD matches TI.** At the 4 gated VE anchors (and a few extra (E,v')), TD σ agrees
with the TI σ from sub-project #3 to a tight tolerance (e.g. `rtol ≤ 10%`, limited only by
dt/finite-T discretization). This is the exact differential oracle.

**V2 — propagation convergence.** Halving `dt` and/or increasing `n_steps` (longer T) changes σ
negligibly; ‖ψ(T)‖ has decayed to a small fraction of ‖ψ(0)‖ (resonance depleted).

**V3 — Houfek anchors.** The 4 gated anchors vs Houfek's `CSVE.V00.J00` land within the same
factor-3 cross-model bound as TI (they should, since TD≈TI).

**Harness wiring.** Update `validation/n2/experiment.py` Group D from PENDING to a real check:
compute TD σ at the gated anchors, assert (a) agreement with TI within the V1 tolerance and
(b) the Houfek factor-3 bound; degrade to a labeled FAIL (not a crash) on error; still runs in
the CPU docker.

## Out of scope

- The DA channel; full σ(E) curves (anchor scope only, like #3).
- Chebyshev propagation (CN chosen); adaptive time-stepping.
- Reproducing eMoScat's exact Filon-quadrature constants (the TI oracle, not eMoScat's TD
  numbers, is what TD is checked against).
- Rust (dense CN with a cached LU is the correct toy end state).

## Verification

- `uv run pytest projects/n2_td_cross_section` green (V1 TD≈TI, V2 convergence, V3 Houfek).
- TD σ matches TI σ at the gated anchors within the stated tolerance; convergence in dt/n_steps.
- N₂ harness Group D flips PENDING → PASS: `python -m validation.n2.experiment` shows D PASS,
  exit 0; docker run still exits 0.
- CN propagator promoted to `qscat.evolution`, mypy-clean; `CLAUDE.md`/docs updated.
- Uses `qscat.dvr`/`qscat.evolution`; CPU-only.
