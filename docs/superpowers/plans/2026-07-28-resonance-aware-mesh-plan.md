# Resonance-aware nuclear mesh + iterative 2-D refinement Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the discretisation tuner's finding #3 — make `propose_grid`'s a-priori NUCLEAR mesh resolve the resonance/interaction region by driving it from the adiabatic resonance curve `V_d(R)/Γ(R)` (sampled DENSELY only where `V_int` is non-negligible, ONCE in the far channel space), add an iterative 2-D refinement loop as the general fallback, and re-tune F₂/H₂⁺ to be 2-D-converged AND ~10–20% smaller than the hand decks.

**Architecture:** The heavy nuclei move adiabatically on the electronic curves — `v0(R)` (VE) and the anion/resonance curve `V_d(R)` (dissociation). The nuclear mesh's `k(R)` becomes the worst-case over both curves, refined where `Γ(R)` peaks. `V_d/Γ` come from sub-project B's `local_complex_potential` pole continuation, refactored to run over an arbitrary R-set so a new `qscat.tuning.resonance` sampler can call it DENSELY inside the interaction region and take a single far/asymptotic value. `refine_to_2d_convergence` is the model-agnostic fallback the skill's step 6 calls. See `docs/superpowers/specs/2026-07-28-resonance-aware-mesh-design.md`.

**Tech Stack:** Python ≥3.12, NumPy/SciPy, `qscat.core.lcp` (`local_complex_potential` + the extracted walk), `qscat.core.dissociation` (`anion_electronic_states`, `da_cross_section`), `qscat.ecs.find_resonance_pole`, `qscat.tuning` (`analyze`/`mesh`/`propose`/`probes`), `qscat.dvr`, `qscat.model`. pytest, mypy --strict over `libs/qscat/qscat`, ruff.

## Global Constraints

- **Atomic units.**
- **Efficiency is FIRST-CLASS:** the per-R resonance scan (~2 electronic diagonalizations per R) runs DENSELY only inside the interaction region `[R_lo, R_hi]` (where `V_int` is non-negligible); the far channel space is a SINGLE asymptotic value (`V_d → ε_asym`, `Γ → 0`). A test measures the solve count and asserts the far region is not densely scanned.
- **Interaction strength is protocol-universal:** use `max_r |model.v_int(r, R)|` as the R-profile of interaction strength (works for neutral AND ionic models — NOT `model.lam`, which the narrowed `ResonanceModel` protocol dropped and `IonicResonanceModel` lacks).
- **The lcp refactor is BEHAVIOR-PRESERVING:** extracting the pole-continuation walk from `local_complex_potential` must leave it bit-identical — gated by the existing `libs/qscat/tests/test_lcp.py` (incl. the N₂ `vres` differential-oracle test).
- **VE stays `v0`-driven:** the resonance-aware path is for the DISSOCIATION channel (DA/DR). A VE-only nuclear mesh keeps `k(R)=√(2μ(E_max−v0(R)))` (its bound `χ_v` is already resolved). `propose_grid` gains an opt-in `resonant`/`channel` flag; default behavior unchanged.
- **`qscat.tuning` may import `qscat.core`/`qscat.model`** (tooling). mypy --strict over `libs/qscat/qscat` clean; ruff clean; the E741/E731 per-file-ignore convention applies.

## File Structure

- `libs/qscat/qscat/core/lcp.py` (modify) — extract `resonance_pole_walk`; `local_complex_potential` calls it.
- `libs/qscat/qscat/tuning/resonance.py` (create) — `interaction_region`, `resonance_curve`.
- `libs/qscat/qscat/tuning/mesh.py` (modify) — a `combined_profile`/worst-case-`k` helper.
- `libs/qscat/qscat/tuning/propose.py` (modify) — the resonant nuclear path.
- `libs/qscat/qscat/tuning/refine2d.py` (create) — `refine_to_2d_convergence`.
- `libs/qscat/qscat/tuning/__init__.py` (modify) — exports.
- Tests: `libs/qscat/tests/test_tuning_resonance.py`, `test_tuning_refine2d.py`; extend `test_tuning_propose.py`.
- `validation/tuning/test_resonance_aware.py` (create) — the F₂/H₂⁺ re-tune gate.
- `.claude/skills/discretisation-tuner/SKILL.md`, `docs/physics/discretisation-tuning.md`, `CLAUDE.md` (modify).

---

### Task 1: `interaction_region` (`qscat.tuning.resonance`)

**Files:**
- Create: `libs/qscat/qscat/tuning/resonance.py`
- Modify: `libs/qscat/qscat/tuning/__init__.py`
- Test: `libs/qscat/tests/test_tuning_resonance.py`

**Interfaces:**
- `interaction_region(model, *, r_probe=..., R_max=..., frac=0.02, n=400) -> tuple[float, float]` — the `[R_lo, R_hi]` where the interaction strength `s(R) = max over r of |Re(model.v_int(r, R))|` exceeds `frac × s.max()`. Returns the outermost such bracket (the R-window where `V_int` is non-negligible).

**Design notes:** `s(R) = max_r |Re(v_int(r_probe, R))|` over a fixed electronic probe set `r_probe` (a modest grid, e.g. `linspace(0.1, 15, 200)`) — universal (uses only `model.v_int`, in the protocol). Scan `R` on `linspace(1e-3, R_max, n)`; `R_lo`/`R_hi` = first/last R with `s(R) ≥ frac·s.max()`. This brackets where the coupling lives (F₂: ~[1.5, 4] bohr).

- [ ] **Step 1: Write the failing test**

```python
# libs/qscat/tests/test_tuning_resonance.py
from __future__ import annotations

import numpy as np
from qscat.model import F2, N2
from qscat.tuning import interaction_region


def test_f2_interaction_region_brackets_the_coupling():
    R_lo, R_hi = interaction_region(F2)
    assert 0.5 < R_lo < 2.5 and 3.0 < R_hi < 8.0     # ~[1.5,4]-ish, where lambda(R) is significant
    assert R_lo < R_hi


def test_region_is_where_vint_is_nonnegligible():
    # outside [R_lo,R_hi], max_r|v_int| is < a few % of its peak
    R_lo, R_hi = interaction_region(N2, frac=0.05)
    r = np.linspace(0.1, 15, 200)
    s = lambda R: np.max(np.abs(np.real(N2.v_int(r[:, None], np.array([[R]])))))
    peak = max(s(R) for R in np.linspace(0.5, 6, 60))
    assert s(R_hi + 1.0) < 0.05 * peak * 1.5          # well outside -> small
```

- [ ] **Step 2–5:** run→fail; implement `interaction_region` in `resonance.py` + export; run→pass; mypy+ruff; commit `feat(tuning): interaction_region (where V_int is non-negligible)`.

---

### Task 2: Extract the pole walk + `resonance_curve` (efficient sampler)

**Files:**
- Modify: `libs/qscat/qscat/core/lcp.py` (extract `resonance_pole_walk`)
- Modify: `libs/qscat/qscat/tuning/resonance.py` (`resonance_curve`)
- Modify: `libs/qscat/qscat/tuning/__init__.py`
- Test: `libs/qscat/tests/test_lcp.py` (must still pass), `libs/qscat/tests/test_tuning_resonance.py`

**Interfaces:**
- In `lcp.py`: `resonance_pole_walk(model, R_descending, elec_grid_a, elec_grid_b, seed_window, *, re_half_width, im_half_width, resid_tol) -> tuple[NDArray[float64], NDArray[float64]]` — the continuation walk over a descending-R array, returning `(E_res, Gamma)` aligned with `R_descending` (the SAME accept/recenter/freeze logic currently inlined in `local_complex_potential`). `local_complex_potential` is refactored to build its descending real-R array + seed window and call this.
- In `qscat.tuning.resonance`: `resonance_curve(model, elec_grid_a, elec_grid_b, *, R_max, n_dense=25, region=None) -> tuple[NDArray, NDArray, NDArray]` — returns `(R_samples, V_d, Gamma)`. `region = region or interaction_region(model)`; build `R_samples` DENSE (`n_dense` points) inside `[R_lo, R_hi]` plus a FEW points down to a small inner R, plus a SINGLE far point at `R_max` (the asymptote); run `resonance_pole_walk` over them descending. `V_d = Re(E_pole)`; the far point supplies the asymptote (`Γ→0` there).

**Design notes:** the refactor must be behavior-preserving — `local_complex_potential`'s result is unchanged (gated by `test_lcp.py`, incl. `test_matches_n2_vres_oracle`). The efficiency win: `resonance_curve`'s `R_samples` are dense only in `[R_lo,R_hi]` + one far value, so the scan does ~`n_dense + few` solves, NOT hundreds. `V_d(R)` outside `[R_lo,R_hi]` is the frozen asymptote (the walk's freeze already does this once `Γ` is small). Seed window from `anion_electronic_states(elec_grid_a, model, R_max, 1)`.

- [ ] **Step 1: Write the failing test** (append to `test_tuning_resonance.py`)

```python
def _elec_grids():
    from qscat.core.grids import electronic_grid
    return (electronic_grid(r_max=16.0, order=7, n_complex=6, angle_deg=35.0),
            electronic_grid(r_max=16.0, order=7, n_complex=6, angle_deg=44.0))


def test_resonance_curve_dense_interaction_sparse_far():
    from qscat.model import F2
    from qscat.tuning import resonance_curve, interaction_region
    ga, gb = _elec_grids()
    R_lo, R_hi = interaction_region(F2)
    R, Vd, G = resonance_curve(F2, ga, gb, R_max=22.0, n_dense=20)
    # most samples land inside the interaction region; the far region is sparse (~1 pt near R_max)
    inside = (R >= R_lo) & (R <= R_hi)
    assert inside.sum() >= 15                     # dense inside
    assert (R > R_hi + 1.0).sum() <= 3            # sparse far
    assert np.all(np.isfinite(Vd)) and np.all(G >= 0.0)
    # Gamma peaks inside the interaction region (the resonance), ~0 far
    assert G[inside].max() > 10 * (G[R > R_hi + 1.0].max() if (R > R_hi + 1.0).any() else 0.0) + 1e-12
```

- [ ] **Step 2–6:** run→fail; extract `resonance_pole_walk` from `lcp.py` (refactor `local_complex_potential` to call it — VERIFY `test_lcp.py` still passes, incl. the N₂ oracle); implement `resonance_curve`; run→pass; mypy+ruff; commit `feat(tuning): resonance_curve (dense-interaction/sparse-far V_d(R)/Gamma sampler) + lcp walk extraction`.

---

### Task 3: Multi-curve resonance-aware nuclear mesh (`mesh`, `propose`)

**Files:**
- Modify: `libs/qscat/qscat/tuning/mesh.py` (a worst-case-`k` combiner)
- Modify: `libs/qscat/qscat/tuning/propose.py` (the resonant nuclear path)
- Modify: `libs/qscat/qscat/tuning/__init__.py`
- Test: `libs/qscat/tests/test_tuning_mesh.py`, `test_tuning_propose.py`

**Interfaces:**
- `mesh.combined_profile(profiles: list[PotentialProfile]) -> PotentialProfile` — a `PotentialProfile` whose `k(x)` is the elementwise MAX over the inputs (interpolated onto a common `x`), `kappa` the elementwise MIN, turning points/singularities unioned. (Lets the mesh resolve the worst-case of several curves.)
- `propose_grid(model, coordinate, energy_range, *, rtol=1e-3, incident=None, channel="ve") -> FemDvrEcsGrid` — `channel="dissociation"` (or `resonant=True`) makes the NUCLEAR path build `k(R)` from the worst-case of `analyze_potential(model.v0, …)` AND `analyze_potential(V_d_interp, …)` (from `resonance_curve`, mass μ), with extra local refinement where `Γ(R)` is large. `channel="ve"` (default) is unchanged (`v0`-only).

**Design notes:** in `propose._nuclear_adapter` (or a resonant sibling): call `resonance_curve(model, elec_a, elec_b, R_max=x_max)` → `(R, V_d, Γ)`; build `V_d_of_R` by interpolation (constant-extrapolate the far asymptote); `profile_v0 = analyze_potential(model.v0, 0, x_max, μ, e_max)`; `profile_vd = analyze_potential(V_d_of_R, 0, x_max, μ, e_max)`; `combined = combined_profile([profile_v0, profile_vd])`; feed to `optimal_real_mesh`; then HALVE elements overlapping the `Γ`-peak region (reuse `mesh`'s refinement post-pass keyed on the Γ>threshold R-interval). Because `V_d` differs from its asymptote only inside `[R_lo,R_hi]`, the extra density lands exactly at the crossing. Keep the electronic path unchanged.

- [ ] **Step 1: Write the failing test**

```python
def test_resonant_nuclear_mesh_refines_the_crossing():
    from qscat.model import F2
    from qscat.tuning import propose_grid, interaction_region
    g_ve = propose_grid(F2, "nuclear", (0.01, 0.05))                      # v0-only
    g_res = propose_grid(F2, "nuclear", (0.01, 0.05), channel="dissociation")
    R_lo, R_hi = interaction_region(F2)
    # the resonance-aware grid packs MORE real points into [R_lo, R_hi] than the v0-only one
    def frac_in(g):
        rp = g.real_points[g.real_points < g.R0]
        return ((rp >= R_lo) & (rp <= R_hi)).sum() / max(rp.size, 1)
    assert frac_in(g_res) > frac_in(g_ve)
```

- [ ] **Step 2–6:** run→fail; implement `combined_profile` + the resonant nuclear path (+ the Γ-refinement); run→pass; also confirm `channel="ve"` default is UNCHANGED (a test that `propose_grid(F2,"nuclear",…)` == the pre-change grid); mypy+ruff; commit `feat(tuning): resonance-aware nuclear mesh (worst-case v0 ⊔ V_d + Gamma refinement)`.

---

### Task 4: `refine_to_2d_convergence` iterative loop (`refine2d`) + skill

**Files:**
- Create: `libs/qscat/qscat/tuning/refine2d.py`
- Modify: `libs/qscat/qscat/tuning/__init__.py`, `.claude/skills/discretisation-tuner/SKILL.md`
- Test: `libs/qscat/tests/test_tuning_refine2d.py`

**Interfaces:**
- `refine_to_2d_convergence(observable, g_r, g_R, *, rtol=1e-2, max_iter=4) -> tuple[FemDvrEcsGrid, FemDvrEcsGrid, dict]` — `observable(g_r, g_R) -> float` (a scalar the caller closes over, e.g. σ_DA at one energy). Compare `observable(g_r, g_R)` vs a once-`refine`d NUCLEAR variant (then electronic); whichever refinement moves the observable more is the under-resolved coordinate → adopt that refinement; repeat until `|Δ|/|value| < rtol` or `max_iter`. Return the converged `(g_r, g_R)` + a `detail` dict (iterations, per-iter values, which coordinate refined). Model-agnostic — the general fallback for structure the adiabatic heuristic misses.

**Design notes:** reuse `probes.refine`. The observable is a closure (keeps `refine2d` model-free). Cap `max_iter`; if it doesn't converge, return the best + `detail["converged"]=False` (a real signal, don't loop forever). Wire into the skill's step 6: the 2-D spot-check becomes "call `refine_to_2d_convergence`; if it needed refinement, the a-priori grid was insufficient — report the converged grid + the cost delta."

**Harness for the `@slow` test:** copy the F2 DA harness from `validation/tuning/test_emoscat_decks.py::test_f2_2d_da_cross_section_spot_check` (it already builds the F2 `(g_r, g_R)`, the anion `eps`/`chi` via `anion_electronic_states` + `vibrational`, and calls `da_cross_section` — reuse that assembly verbatim; only the grids and the loop differ). Read that test first.

- [ ] **Step 1: Write the failing test** (`@slow` — one 2-D observable; assemble the F2 harness by copying the spot-check test above)

```python
@pytest.mark.slow
def test_refine_converges_f2_da_from_coarse_guess():
    # Assemble F2 (g_r, coarse g_R) + eps/chi exactly as the spot-check test does.
    # observable(gr, gR) = da_cross_section(TensorGrid([gr, gR]), F2, eps, chi, 0, 0.03)[0]
    # A deliberately coarse nuclear guess -> refine_to_2d_convergence lifts sigma_DA
    # to the converged value in a few iterations.
    obs = lambda gr, gR: float(da_cross_section(TensorGrid([gr, gR]), F2, eps, chi, 0, 0.03)[0])
    g_r2, g_R2, detail = refine_to_2d_convergence(obs, g_r, g_R_coarse, rtol=0.05, max_iter=4)
    assert detail["converged"] and obs(g_r2, g_R2) > 1.0    # lifted from ~0.3 toward 1.66
```

- [ ] **Step 2–6:** run→fail; implement `refine2d` + the skill step-6 wiring; run the `@slow` test FOREGROUND (a few 2-D solves, minutes — be patient, do NOT background/Monitor); mypy+ruff; commit `feat(tuning): refine_to_2d_convergence iterative loop + skill step-6 wiring`.

---

### Task 5: Re-tune F₂/H₂⁺ — 2-D-converged AND ~10–20% smaller (validation)

**Files:**
- Create: `validation/tuning/test_resonance_aware.py`
- Modify: `docs/physics/discretisation-tuning.md`, `CLAUDE.md`, `.claude/skills/discretisation-tuner/SKILL.md`
- Test: the above.

**Design notes:** the payoff gate. Reuse the F2 harness from `validation/tuning/test_emoscat_decks.py::test_f2_2d_da_cross_section_spot_check` (grids, `eps`/`chi`, `da_cross_section`) — this test is the resonance-aware counterpart: it swaps the deck/`v0`-only nuclear grid for `propose_grid(..., channel="dissociation")` and asserts convergence on the first pass. For F₂: `g_R = propose_grid(F2, "nuclear", (0.01,0.05), channel="dissociation")`; build the 2-D grid with a proposed electronic grid; assert (a) the `@slow` 2-D DA spot-check now CONVERGES on the FIRST resonance-aware a-priori grid (σ_DA within ~rtol of the once-refined value — the finding-#3 fix, contrast with the old `v0`-only grid which was ~5× off, σ_DA≈0.31 vs converged 1.64); (b) `grid_cost(g_R)["n_points"] ≤ 0.9 × deck_nuclear_n` (the 10–20% reduction vs the committed F₂ nuclear deck; read `deck_nuclear_n` from the same test's deck grid). For H₂⁺: the same on the proxy nuclear grid (2-D-representative + smaller than the proxy deck). Update `docs/physics/discretisation-tuning.md` "Genuine finding #3" from "known limitation" to "RESOLVED (resonance-aware mesh)"; note the iterative loop as the general fallback; update `CLAUDE.md` + the skill.

- [ ] **Steps:** write the `@slow` F₂ (and H₂⁺-proxy) gate (2-D-converged + smaller-than-deck); run FOREGROUND (minutes); if the 2-D spot-check does NOT converge on the first resonance-aware pass (a real finding — the adiabatic heuristic wasn't enough), report with numbers and fall back to documenting the iterative loop as the mechanism; update docs/CLAUDE/skill; mypy+ruff; commit `feat(tuning): resonance-aware F2/H2+ re-tune (2-D-converged + ~10-20% smaller) + docs`.

---

## Verification (whole sub-project)

- `uv run pytest -q -m "not slow"` pass; the `@slow` refine-loop + re-tune gates pass; `test_lcp.py` still passes (the lcp refactor is behavior-preserving).
- `uv run mypy libs/qscat/qscat` 0; `uv run ruff check .` clean.
- The resonance scan is demonstrably sparse-far/dense-interaction (the Task-2 solve-count test).
- The resonance-aware F₂ nuclear grid packs its density at the R≈2.5–2.7 crossing; the `@slow` 2-D DA spot-check CONVERGES on the first resonance-aware `propose_grid` (contrast the old `v0`-only ~5× gap); the re-tuned F₂/H₂⁺ grids are ~10–20% smaller than the hand decks AND 2-D-converged.
- `channel="ve"` default behavior is unchanged. `docs/physics/discretisation-tuning.md` finding #3 → resolved; `CLAUDE.md` + skill updated.

## Out of scope (this plan)

- **Electronic-mesh worst-case-over-R** (a separate defensive fix — finding #3 is nuclear).
- **Multi-channel / non-adiabatic (coupled) nuclear curves.**
- **Re-gridding production runs** onto the new grids (the tuner emits configs; adoption stays opt-in).
