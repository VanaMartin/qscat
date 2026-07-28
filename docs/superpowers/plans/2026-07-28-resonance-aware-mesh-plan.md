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
- `interaction_region(model, *, r_probe=None, R_max=8.0, frac=0.02, n=400) -> tuple[float, float]` — the `[R_lo, R_hi]` window over which the coupling `s(R) = max over r of |Re(model.v_int(r, R))|` is still TRANSITIONING between its low-R and high-R regimes. `R_lo`/`R_hi` bracket the middle `(1−2·frac)` of `s`'s own `[min,max]` rise on the scan.

**Design notes (CORRECTED — see the spec's "Premise correction"):** `s(R) = max_r |Re(v_int(r_probe, R))|` over `r_probe = linspace(0.1, 15, 200)` — universal (uses only `model.v_int`). Scan `R` on `linspace(1e-3, R_max, n)`. `λ(R)` (hence `s(R)`) is a SIGMOID saturating to a nonzero plateau on both ends, so a raw `s ≥ frac·s.max()` threshold degenerates to the scan edges — do NOT use it. Instead: `span = s.max()−s.min()`; `R_lo` = first R with `(s−s.min) ≥ frac·span`, `R_hi` = last R with `(s−s.min) ≤ (1−frac)·span`; if `span ≤ 0` return the full scanned domain. This brackets the transition zone where `V_d(R)` is non-saturated (F₂: ≈[0.66, 3.03] bohr, containing the R≈2.5–2.7 crossing).

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


def test_region_ends_where_coupling_has_saturated():
    # outside [R_lo,R_hi], s(R) has SATURATED (stopped changing) — NOT that it is
    # small in absolute terms (the sigmoid plateau is large). This is the property
    # the downstream resonance_curve sampler relies on (freeze at one far value).
    R_lo, R_hi = interaction_region(N2, frac=0.05)
    r = np.linspace(0.1, 15, 200)
    s = lambda R: np.max(np.abs(np.real(N2.v_int(r[:, None], np.array([[R]])))))
    peak = max(s(R) for R in np.linspace(0.5, 6, 60))
    assert abs(s(R_hi + 1.0) - s(R_hi)) < 0.05 * peak * 1.5   # saturated past R_hi
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

### Task 3 (REVISED): Resonance-aware nuclear mesh — crossing super-refinement + exit-wave order

> **This supersedes the original k-merge design.** A validated diagnostic (2026-07-28) showed the
> worst-case-`k` merge is INERT: `min_len=0.15` floors the exit-region elements, so the resonant
> grid (614 pts) came out ≈ the v0 grid (609) and would stay ~5× off on σ_DA. The real levers
> (matching how the eMoScat deck converges with FEWER points than brute refinement) are: **(1)
> resolve the fast exit wave with a high ORDER, (2) super-refine the narrow crossing, (3) trim the
> real extent.** The prior commits `c4f2404`/`248d451` (the k-merge + `combined_profile`) are
> SUPERSEDED — this task builds on the still-good scaffolding they added (`channel` param, the
> `elec_grids`/`resonance_n_dense` kwargs, the ValueError guards) but REPLACES the mesh logic and
> REMOVES `combined_profile` if it ends up unused (YAGNI).

**Validated numbers (F2, use as the design anchors — a subagent should reproduce, not trust blindly):**
- Crossing `R* = 2.598` = the `Re(V_d(R)) − v0(R)` sign-change (lands on `F2.R_c=2.595`). `argmax Γ`
  is WRONG (Γ has a frozen-plateau artifact at the walk's inner edge → `R≈0.06`).
- `V_d` asymptote ≈ −0.127 Ha; `K_exit = √(2μ(E_max − V_d_asym)) ≈ 78` (μ≈17316). λ_exit ≈ 0.080 bohr.
- Points/wavelength = `order·λ_exit / L`. Order-6 at L=0.15 → 3.2 ppw (UNDER-resolved). Order-14 at
  L=0.15 → 7.5 ppw (resolved). Deck uses q=14 + 0.024-bohr crossing elements + ~10.7-bohr extent.

**Files:**
- Modify: `libs/qscat/qscat/tuning/mesh.py` (a local-refinement + resolution-order helper; remove `combined_profile` if unused)
- Modify: `libs/qscat/qscat/tuning/propose.py` (the revised resonant nuclear path)
- Modify: `libs/qscat/qscat/tuning/__init__.py`
- Test: `libs/qscat/tests/test_tuning_mesh.py`, `test_tuning_propose.py`, and a `@slow` 2-D convergence test in `test_tuning_propose.py`

**Interfaces:**
- `mesh.refine_elements_in_window(real_lengths, x_min, R_lo, R_hi, target_len) -> list[float]` — subdivide (span-preservingly) every real element overlapping `[R_lo, R_hi]` until each is ≤ `target_len`; elements outside are untouched. (Reuse the `_clamp_lengths_span_preserving` subdivide logic — the local, override-`min_len` refinement for the crossing.)
- `mesh.order_for_wavenumber(k, element_len, *, target_ppw=6.0, orders=(6,8,10,14)) -> int` — the smallest `order` in `orders` with `order·(2π/k)/element_len ≥ target_ppw` (i.e. resolves wavenumber `k` at that element length); falls back to `max(orders)` if none qualifies.
- `propose_grid(model, coordinate, energy_range, *, rtol=1e-3, incident=None, phase_coeff=None, channel="ve", elec_grids=None, resonance_n_dense=25) -> FemDvrEcsGrid` — `channel="dissociation"` (nuclear only) takes the revised resonant path below. `channel="ve"` (default) UNCHANGED.

**Design notes (the revised resonant nuclear path):**
1. `R, Vd, Gamma = resonance_curve(model, ga, gb, R_max=x_max, n_dense=resonance_n_dense)` (build/accept `elec_grids` as before).
2. **Exit wave:** `Vd_asym = Re(Vd)` at the largest sampled R; `K_exit = √(2·μ·max(e_max − Vd_asym, e_max))` (the fast outgoing dissociation wavenumber). Set the ECS-tail `channel_k = K_exit` (physical — the wave the tail must absorb). Build the base real mesh from the **v0 profile** (as the ve path) — but choose the DVR order via `order_for_wavenumber(K_exit, min_len)` instead of the point-count-min sweep, so the exit wave is resolved (F2 → order 14).
3. **Crossing:** `R* = ` the outermost `Re(Vd) − v0` sign-change (interpolated); if none, fall back to `R[argmin(Re(Vd) − v0)]`. `δ` from the Γ-closing width — the R-range around `R*` where `Γ > 0.1·Γ_significant` (exclude the frozen inner plateau: restrict to `R ≥ R_lo` from `interaction_region`), clamped to a sane `[0.15, 0.6]` bohr half-width. Super-refine `[R*−δ, R*+δ]` to `target_len ≈ 0.03` bohr via `refine_elements_in_window` (this OVERRIDES `min_len` locally — that is the point).
4. **Extent:** derive `x_max` for the resonant path from where both `v0` and `V_d` have flattened to their asymptote (within a small tol) plus a fixed outgoing-wave margin (a few bohr) before the ECS tail — NOT the blunt `_NUCLEAR_X_MAX_DEFAULT=18`. Target ≈ the deck's ~10–11 bohr for F2. If a principled derivation is hard, expose it as a parameter defaulting to a reduced value (e.g. 12.0) and note it — the ECS tail does the absorbing, so the real exit region need only host the outgoing wave a few λ past the interaction region.
5. Assemble the `GridSpec` (real refined lengths + ECS tail at the chosen order/angle) as the ve path does.

- [ ] **Step 1 — unit tests (FAST, no models):**
```python
def test_order_for_wavenumber_resolves_fast_wave():
    from qscat.tuning import order_for_wavenumber
    import numpy as np
    # K_exit=78, lambda≈0.08; at 0.15-bohr elements order 6 gives 3.2 ppw (too few), 14 gives 7.5
    assert order_for_wavenumber(78.0, 0.15, target_ppw=6.0) == 14
    assert order_for_wavenumber(5.0, 0.15, target_ppw=6.0) <= 8       # slow wave, low order ok

def test_refine_elements_in_window_only_local_and_span_preserving():
    from qscat.tuning import refine_elements_in_window
    import numpy as np
    lengths = [0.5]*10                                   # span 5.0, from x_min=0
    out = refine_elements_in_window(lengths, 0.0, 2.0, 3.0, 0.1)
    assert abs(sum(out) - 5.0) < 1e-12                   # span preserved
    assert max(out) <= 0.5 + 1e-12                       # nothing coarsened
    # elements straddling [2,3] are now <=0.1; those outside stay 0.5
    assert min(out) <= 0.1 + 1e-12
```
- [ ] **Step 2 — VE-unchanged regression (FAST):** `test_ve_channel_default_unchanged` — `propose_grid(F2,"nuclear",(0.01,0.05))` `.n` and `.points` identical with/without `channel="ve"` (the default path is byte-unchanged). Keep the ValueError-guard tests (`test_propose_grid_rejects_unknown_channel`, `..._rejects_dissociation_electronic`).
- [ ] **Step 3 — the crossing/order test (small elec grids for speed):** `propose_grid(F2,"nuclear",(0.01,0.05),channel="dissociation",elec_grids=<small>,resonance_n_dense=10)` returns a grid whose (a) order == `order_for_wavenumber(K_exit,min_len)` (≥10), and (b) the smallest real element is `≤ 0.05` bohr AND sits within `[R*−δ, R*+δ]` around `R*≈2.6` (the super-refined crossing) — NOT at the inner wall. Mark `@pytest.mark.slow` if the resonance scan pushes it over ~30s.
- [ ] **Step 4 — the LOAD-BEARING `@slow` 2-D convergence gate (this proves the mechanism):**
```python
@pytest.mark.slow
def test_resonant_nuclear_grid_converges_f2_da():
    # Copy the F2 DA harness from
    # validation/tuning/test_emoscat_decks.py::test_f2_2d_da_cross_section_spot_check
    # (electronic grid, anion eps/chi, da_cross_section). Build g_R via the resonant path.
    # The resonant a-priori grid must give sigma_DA within ~15% of the once-refined value
    # (converged ≈1.6), i.e. CONVERGED on the first pass — contrast the v0 grid's ~0.31 (5x off).
    ...
    sig_base = da_cross_section(TensorGrid([g_r, g_R]), F2, eps, chi, 0, 0.03)[0]
    sig_ref  = da_cross_section(TensorGrid([g_r, refine(g_R)]), F2, eps, chi, 0, 0.03)[0]
    assert abs(sig_base - sig_ref) / abs(sig_ref) < 0.15      # converged on the first resonant pass
    assert sig_base > 1.0                                     # lifted off the v0 grid's ~0.31
```
- [ ] **Step 5:** run all; `uv run mypy libs/qscat/qscat`; `uv run ruff check .`; commit `feat(tuning): resonance-aware nuclear mesh — crossing super-refine + exit-wave order (supersedes k-merge)`. If `combined_profile` is now unused, remove it + its test in the same commit (note it in the message).

---

### Task 4: `refine_to_2d_convergence` iterative loop (`refine2d`) + skill

**Files:**
- Create: `libs/qscat/qscat/tuning/refine2d.py`
- Modify: `libs/qscat/qscat/tuning/__init__.py`, `.claude/skills/discretisation-tuner/SKILL.md`
- Test: `libs/qscat/tests/test_tuning_refine2d.py`

**Interfaces:**
- `refine_to_2d_convergence(observable, g_r, g_R, *, rtol=1e-2, max_iter=4) -> tuple[FemDvrEcsGrid, FemDvrEcsGrid, dict]` — `observable(g_r, g_R) -> float` (a scalar the caller closes over, e.g. σ_DA at one energy). Each iteration: evaluate the observable on a once-`refine`d NUCLEAR variant and a once-`refine`d ELECTRONIC variant; whichever moves the observable MORE (larger `|Δ|/|value|`) is the under-resolved coordinate → adopt that refinement and record the step. STOP when the larger of the two relative moves is `< rtol` (converged — record NO step for that final check) or when `max_iter` adopted steps is hit (not converged). Return `(g_r, g_R)` (refined) + `detail` with EXACTLY these keys: `detail["converged"]: bool`, `detail["iterations"]: list[dict]` where each adopted step is `{"coordinate": "nuclear"|"electronic", "value": float, "rel_move": float}`, and `detail["final_value"]: float`. So an already-converged input returns `iterations == []` and `converged is True`; a never-converging one returns `len(iterations) == max_iter` and `converged is False`. Model-agnostic — the general fallback for structure the adiabatic heuristic misses.

**Design notes:** reuse `probes.refine`. The observable is a closure (keeps `refine2d` model-free). Cap `max_iter`; if it doesn't converge, return the best + `detail["converged"]=False` (a real signal, don't loop forever). Wire into the skill's step 6: the 2-D spot-check becomes "call `refine_to_2d_convergence`; if it needed refinement, the a-priori grid was insufficient — report the converged grid + the cost delta."

**Test strategy — the loop LOGIC is tested cheaply (a real F2-DA loop is 20–40 min of solves, uncheckable in-harness).** The observable is a closure, so a SYNTHETIC one exercises the whole loop deterministically and fast. Test: (1) convergence + right-coordinate adoption, (2) the `max_iter` cap, (3) already-converged (0 iterations). The real F2-DA integration is a `@slow` marker only (its per-solve cost is documented; the Task-3 gate already proved the resonant grid converges F2 DA).

- [ ] **Step 1: Write the failing tests (FAST synthetic observable — deterministic, no 2-D solves):**
```python
# The observable "converges" as the NUCLEAR grid gains points; the electronic grid
# is already fine (refining it does nothing). A closure over g.n mimics a real
# observable's approach to its converged limit, so the loop must (a) pick the nuclear
# coordinate, (b) stop when |Δ|/|value| < rtol.
def _obs_factory(exact=1.66, scale=200.0):
    return lambda g_r, g_R: exact - scale / g_R.n     # rises toward `exact` as g_R refines

def test_refine_adopts_nuclear_and_converges():
    from qscat.tuning import propose_grid, refine_to_2d_convergence
    from qscat.model import F2
    g_r = propose_grid(F2, "electronic", (0.01, 0.05))
    g_R = propose_grid(F2, "nuclear", (0.01, 0.05))
    obs = _obs_factory()
    g_r2, g_R2, detail = refine_to_2d_convergence(obs, g_r, g_R, rtol=1e-2, max_iter=6)
    assert detail["converged"]
    assert g_R2.n > g_R.n and g_r2.n == g_r.n          # refined nuclear, left electronic alone
    assert all(step["coordinate"] == "nuclear" for step in detail["iterations"])

def test_refine_caps_at_max_iter_when_never_converging():
    from qscat.tuning import propose_grid, refine_to_2d_convergence
    from qscat.model import F2
    g_r = propose_grid(F2, "electronic", (0.01, 0.05))
    g_R = propose_grid(F2, "nuclear", (0.01, 0.05))
    # An observable that keeps changing by a fixed relative amount never converges.
    flip = {"v": 1.0}
    def obs(g_r, g_R):
        flip["v"] *= -2.0
        return flip["v"]
    _, _, detail = refine_to_2d_convergence(obs, g_r, g_R, rtol=1e-3, max_iter=3)
    assert detail["converged"] is False and len(detail["iterations"]) == 3

def test_refine_already_converged_is_zero_iterations():
    from qscat.tuning import propose_grid, refine_to_2d_convergence
    from qscat.model import F2
    g_r = propose_grid(F2, "electronic", (0.01, 0.05))
    g_R = propose_grid(F2, "nuclear", (0.01, 0.05))
    _, _, detail = refine_to_2d_convergence(lambda a, b: 3.14, g_r, g_R, rtol=1e-3, max_iter=4)
    assert detail["converged"] and len(detail["iterations"]) == 0
```
- [ ] **Step 1b (optional `@slow`, documented — do NOT block on it running in-harness):** a `test_refine_converges_f2_da_from_coarse_guess` marked `@pytest.mark.slow` that closes `observable` over the real `da_cross_section` (harness from `validation/tuning/test_emoscat_decks.py::test_f2_2d_da_cross_section_spot_check`). Note in a comment that a full run is ~20–40 min (multiple 2-D solves) and is not run in the fast suite — the fast synthetic tests above gate the loop logic.

- [ ] **Step 2–6:** run→fail; implement `refine2d` + the skill step-6 wiring; run the `@slow` test FOREGROUND (a few 2-D solves, minutes — be patient, do NOT background/Monitor); mypy+ruff; commit `feat(tuning): refine_to_2d_convergence iterative loop + skill step-6 wiring`.

---

### Task 5: Re-tune F₂/H₂⁺ — 2-D-converged AND ~10–20% smaller (validation)

**Files:**
- Create: `validation/tuning/test_resonance_aware.py`
- Modify: `docs/physics/discretisation-tuning.md`, `CLAUDE.md`, `.claude/skills/discretisation-tuner/SKILL.md`
- Test: the above.

**Design notes:** the payoff gate. **VERIFIED NUMBERS (controller, 2026-07-28 — use these, don't re-derive the expensive convergence):**
- **F₂**: `propose_grid(F2,"nuclear",(0.01,0.05),channel="dissociation")` → **n=1000** (order 14) vs the **974-pt deck** = **1.027× (deck-PARITY)**; σ_DA(E=0.03)=**1.6562**, CONVERGED (matches deck 1.66, finding-#3 refine² 1.658). The old v0-only grid = 609 pts but σ_DA≈0.31 (5× off — the finding-#3 bug).
- **H₂⁺**: resonant → **n=489** (order 8) vs the **510-pt proxy deck** = **0.959× (~4% SMALLER)**. Convergence NOT laptop-verifiable (2-D DR ~1.15M unknowns full — Docker/MUMPS); size + build-success is the laptop-verifiable part.

**The honest finding (state it plainly in the docs, do NOT hide it):** the "10–20% smaller than the deck" expectation does NOT hold for F₂ — the eMoScat deck is a near-optimal expert hand-tuning, and the old propose_grid WAS smaller (609) only because it was *under-converged* (finding #3). To CONVERGE you need ~deck size; the resonance-aware tuner reaches it **automatically** (deck-parity for F₂, ~4% under for H₂⁺). The deliverable is **convergence + automation at deck-competitive size**, not a point reduction.

**Gates (honest):** (a) SIZE (FAST, cheap — no solve): `F2_resonant.n <= 1.05 * F2_deck.n` (deck-parity) and `H2P_resonant.n <= H2P_proxy_deck.n` (no larger; in fact smaller). Read `F2_deck` from `validation.tuning.test_emoscat_decks.CONFIGS["F2"].da_grid().grids[1]` and `H2P_proxy_deck` from `validation.h2plus.config.proxy_grid().grids[1]`. (b) CONVERGENCE (F₂, `@slow`, documented — the controller verified σ_base=1.6562; write the test but it need not run in the fast suite): copy the `test_f2_2d_da_cross_section_spot_check` harness, build `g_R` via `channel="dissociation"`, assert `|σ_base − σ_refined|/σ_refined < 0.15` and `σ_base > 1.0`. H₂⁺ convergence = a documented Docker/MUMPS deferral (like the existing H₂⁺ 2-D handling). Update `docs/physics/discretisation-tuning.md` "Genuine finding #3" → "RESOLVED (resonance-aware mesh)" with the mechanism + these numbers + the honest size finding; note `refine_to_2d_convergence` as the general fallback; update `CLAUDE.md`'s `qscat.tuning` entry.

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
