# Resonance-aware nuclear mesh + iterative 2-D refinement — Design Spec

**Date:** 2026-07-28
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — spec for review
**Lifecycle:** `qm-method-lifecycle` — extends the `qscat.tuning` discretisation tuner (merged
2026-07-28) by closing its documented "Genuine finding #3": the a-priori nuclear mesh, built from
`v0(R)` alone, is blind to the interaction-region structure and is ~5× under-converged for a 2-D
dissociation observable until manually refined.

## Context

`qscat.tuning.propose_grid` builds the nuclear mesh from `k(R) = √(2μ(E_max − v0(R)))` — the
neutral/ion channel curve **alone**. The interaction `V_int(r,R)` is not used for the nuclear mesh
at all. But the DA/DR dynamics happen on the **anion/resonance curve `V_d(R)`**, whose crossing +
the resonance width `Γ(R)` peak live in `V_int`/λ(R) around R≈2.5–2.7 bohr — a region `v0(R)` is
smooth through. So the a-priori nuclear grid under-resolves exactly where the dissociation flux is
generated: F₂'s σ_DA comes out 0.31 vs the converged 1.66 until the nuclear grid is refined. The
1-D probes miss it (they test the *smooth* bound vibrational states and the *asymptotic* outgoing
wave); only the skill's 2-D spot-check catches it, manually.

The fix (agreed): the heavy nuclei move **adiabatically on the electronic curves** — `v0(R)` (VE)
AND `V_d(R)` (dissociation). Drive the nuclear mesh from the **worst-case over both curves plus the
Γ(R) coupling region**, using the machinery we already have (`qscat.core.lcp.local_complex_potential`
/ `find_resonance_pole`-per-R from sub-project B). Add an **iterative 2-D refinement loop** as the
general fallback for any structure the a-priori guess misses.

## Goal

Make `propose_grid`'s a-priori nuclear grid **2-D-converged for the dissociation observable on the
first pass** (or after a couple of cheap iterations), by feeding the interaction into the guess —
and, because the hand decks were deliberately conservative, land the re-tuned F₂/H₂⁺ grids
**~10–20% SMALLER than those decks while being correct** (adaptive redistribution: fine at the
resonance + the fast outgoing wave, coarse in the smooth well/asymptotic regions).

## Key efficiency constraint (first-class)

The per-R resonance scan (~2 electronic diagonalizations per R) is the tuner's real cost. It MUST
be **sparse in the far channel region and dense only where `V_int` is non-negligible**:
- **Interaction region** `[R_lo, R_hi]` — where `max_r |V_int(r,R)|` (equivalently `|λ(R)|`) exceeds
  a small fraction of its peak. Sample `V_d(R)/Γ(R)` at a *reasonable* density here (enough to
  resolve the crossing — the continuation walk needs consecutive R anyway).
- **Far channel space** (`R > R_hi`, `R < R_lo`) — `V_int → 0`, so `V_d → ` the asymptotic anion
  energy and `Γ → 0`. Compute this **once** (a single asymptotic solve); the nuclear mesh there is
  driven by `v0(R)` + the outgoing-channel wavenumber alone. This reuses
  `local_complex_potential`'s existing continuation-freeze: stop the pole scan once past the
  interaction region with `Γ` below threshold, freeze at the asymptote.

## Approach

1. **Interaction-region detector** (`qscat.tuning`): `interaction_region(model, ...) -> (R_lo, R_hi)`
   from `|V_int|`/`λ(R)` vs a fraction of its peak. Pure-ish (uses `model.v_int`/`model.lam`).
2. **Efficient resonance-curve sampler** (`qscat.tuning`): `resonance_curve(model, elec_grid_a,
   elec_grid_b, R_samples) -> (V_d, Γ)` — reuse the sub-project-B pole continuation, but sample
   `R_samples` DENSE only inside `[R_lo, R_hi]` and take a SINGLE far/asymptotic value (freeze
   beyond). Returns `V_d(R)`, `Γ(R)` on the sample set (interpolatable).
3. **Multi-curve nuclear mesh** (extend `qscat.tuning.mesh`/`propose`): the nuclear `k(R)` profile
   = worst-case (max) over `k_{v0}(R)` and `k_{V_d}(R) = √(2μ(E_max − V_d(R)))`, with extra local
   refinement where `Γ(R)` is large (halve elements there). Only in `[R_lo, R_hi]` does `V_d` differ
   from the asymptote, so the extra density lands exactly at the resonance. `propose_grid` gains a
   `resonant=True`/`channel="dissociation"` path (VE-only stays `v0`-driven — its nuclear
   wavefunction is the smooth bound χ_v, already resolved).
4. **Iterative 2-D refinement loop** (`qscat.tuning` + the skill): the general fallback —
   `refine_to_2d_convergence(model, g_r, g_R, observable, energy, rtol, max_iter)`: run the 2-D
   observable on `(g_r, g_R)` vs a once-refined variant; where they disagree, LOCALLY refine the
   coordinate/region driving the disagreement (the marginal 2-D `k`-field, or a resonance-region
   halving) and repeat until converged or `max_iter`. Model-agnostic; catches structure the
   adiabatic-curve heuristic misses. The skill's step 6 calls it.
5. **Re-tune F₂ and H₂⁺** (validation): show the resonance-aware a-priori grids are 2-D-converged
   for σ_DA/σ_DR (the spot-check passes on/near the first pass) AND ~10–20% smaller than the hand
   decks.

## Deliverables

- **D1** `interaction_region` + `resonance_curve` (efficient sparse-far/dense-interaction sampler,
  reusing `local_complex_potential`).
- **D2** the multi-curve (`v0` ⊔ `V_d` ⊔ Γ-refined) nuclear mesh path in `propose_grid`.
- **D3** `refine_to_2d_convergence` iterative loop + the skill wiring.
- **D4** re-tuned F₂/H₂⁺ grids: 2-D-converged for σ_DA/σ_DR, ~10–20% fewer points than the decks;
  gated (the F₂ 2-D spot-check now passes on the resonance-aware a-priori grid, not only after a
  manual refine).

## Validation

- `interaction_region(F2)` brackets ~[1.5, 4] bohr (where λ(R) is significant); `resonance_curve`
  does the dense scan only there + one far value (assert the far region used ≤ a few solves).
- The resonance-aware nuclear mesh places its finest elements at the F₂ crossing (R≈2.5–2.7), not
  spread by `v0`; the `@slow` **2-D DA spot-check passes on the FIRST resonance-aware
  `propose_grid`** (σ_DA within rtol of the converged value — the fix for finding #3).
- Re-tuned F₂/H₂⁺: `grid_cost` ≤ ~0.8–0.9× the committed deck (the 10–20% reduction) AND the 2-D
  observable converged — both hold together (correct AND smaller).
- `refine_to_2d_convergence` converges the F₂ σ_DA from the old `v0`-only guess in ≤ a couple of
  iterations (the general fallback still works even without the adiabatic heuristic).

## Sub-project decomposition (tasks for the plan)

1. `interaction_region` (from `|V_int|`/`λ(R)`) + tests.
2. `resonance_curve` efficient sampler (dense-interaction / single-far, reusing
   `local_complex_potential`) + tests (asserts the far region is computed sparsely).
3. Multi-curve nuclear mesh in `mesh`/`propose` (worst-case `v0`⊔`V_d`, Γ-refinement) + tests.
4. `refine_to_2d_convergence` iterative loop + skill wiring.
5. Re-tune F₂/H₂⁺ validation: 2-D-converged a-priori grids, 10–20% smaller, gated.

## Out of scope

- **Multi-channel / non-adiabatic (coupled) nuclear curves** — the worst-case-over-adiabatic-curves
  model suffices for the current single-resonance diatomics/H₂⁺.
- **Re-gridding production runs onto the new grids** — the tuner EMITS configs; adoption stays opt-in.
- **The electronic-mesh worst-case-over-R** improvement (a separate, smaller defensive fix) — note
  it, but the finding-#3 gap is nuclear, so this sub-project targets the nuclear mesh.

## Verification

- `uv run pytest -q -m "not slow"` pass; the `@slow` F₂ 2-D spot-check passes on the resonance-aware
  a-priori grid. `uv run mypy libs/qscat/qscat` 0; `uv run ruff check .` clean.
- The resonance scan is demonstrably sparse-far/dense-interaction (measured solve count).
- F₂/H₂⁺ re-tuned grids are correct (2-D-converged) AND ~10–20% smaller than the hand decks —
  gated. `docs/physics/discretisation-tuning.md` "finding #3" updated to "resolved"; `CLAUDE.md` +
  the skill updated.
