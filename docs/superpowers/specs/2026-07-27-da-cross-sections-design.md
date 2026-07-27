# Dissociative-attachment (DA) cross sections — LCP + exact-2D — Design Spec

**Date:** 2026-07-27
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — implementation pending
**Lifecycle:** `qm-method-lifecycle` — stage 2–4 (toy → validate) for a NEW capability (the DA
exit channel), on top of the existing 1-D LCP TD (`projects/n2_td_cross_section`) and the
promoted exact-2-D engine (`qscat.core`). Follows the DA-framework diagnosis
(`docs/physics/diatomic-ve-cross-sections.md`).

## Context

So far the library computes only the **vibrational-excitation (VE)** exit channel
(`e⁻ + AB(v=0) → AB⁻* → e⁻ + AB(v')`). The transient anion `AB⁻*` has a second fate —
**dissociative attachment (DA)**: `e⁻ + AB(v=0) → AB⁻* → A + B⁻`, the anion breaking apart with
the electron captured, outgoing flux in the **nuclear** coordinate. DA is a major channel for
NO and (dominantly, exothermically) F₂; it is closed for N₂ in the measured range. We need
σ_DA(E) for **both** the approximate **LCP** (local-complex-potential) model and the **exact
2-D** model — the LCP is the approximation under test, the exact 2-D is the oracle
([[research-program-is-method-validation]]).

**eMoScat implements DA for both models correctly**, and is the source for both formulas:
- LCP DA: `source/module_LCP.cpp` + `source/ModelLCP/SMatrix.cpp` (`Make_Vres`, the S-matrix).
- Exact-2-D DA: `source/Model2d/{MultiTestFunction2d,TestFunction2d,Potentials2d}.cpp` (the
  nuclear-coordinate test function).

The DA-framework diagnosis already established: the anion has exactly one bound electronic
state per molecule; thresholds `ε_e(R→∞) − eps[0]` are **F₂ −0.069 Ha (exothermic/open), N₂
+0.502 Ha (closed), NO +0.172 Ha** — validated against the physics and the eMoScat decks (N₂
DA disabled). The naive exact-2-D TI driven-equation DA fails structurally (fragment state
invalid at small R; ECS absorbs the surface flux) — so the exact-2-D DA is done the eMoScat
way (**time-dependent**, nuclear test function), NOT via a TI volume T-matrix.

## The physics (both from eMoScat)

### σ normalization (shared, two equivalent conventions)
The LCP writes σ per channel as **4π³|S_T|²/(2E)** (`SMatrix.cpp`, `S_T` = the T-matrix-like
transition amplitude), the exact-2-D as **π|S − δ|²/(2E)** (`MultiTestFunction2d.cpp`, `S` = the
S-matrix). These are the SAME for an inelastic channel (`S = δ − 2πi·S_T`, so π|S−δ|² = 4π³|S_T|²)
— exactly the identity already used in `qscat.core.driven` (4π³|T|²) vs `qscat.core.time_dependent`
(π|S−δ|²). DA is off-diagonal (δ=0). Incident mass = 1 (electron). Threshold shift per channel;
a channel contributes only where `E + shift > 0`.

### LCP DA (the approximate model)
The nuclear wavepacket lives on the **complex resonance curve**:
`H_res = T_nuc(μ) + diag(V_res(R) − i·Γ(R)/2)` — already in `projects/n2_td_cross_section`
(`H_res = kinetic(grid, mu) + diag(Vd − 0.5j·Gamma)`), doorway `d_v(R) = √(Γ(R)/2π)·χ_v(R)`.
- **Propagation:** `ψ(t) = e^{−iH_res t} d_0` (the incident v=0 doorway), `‖ψ‖` decays as the
  resonance autodetaches (the `−iΓ/2` sink).
- **VE (existing):** `S_VE,v'(E) = (1/i)∫ e^{iE_tot t} ⟨d_v'|ψ(t)⟩ dt`.
- **DA (new):** the flux that **survives to the dissociation boundary** `R = X` (the outer
  grid edge). Per `SMatrix.cpp::close_multistep`:
  `S_DA(E) = √(K/2πμ)·e^{−iK·X} · ∫ e^{iE_tot t} ψ(X, t) dt`, `K = √(2μ(E + shift_DA))`,
  `shift_DA = eps[0] − ε_e(R→∞)` (so `E + shift_DA` = the outgoing nuclear kinetic energy).
  I.e. DA reads the propagating wavepacket's amplitude at the last real grid point (`buffer_[DA][i]
  = psi.f(grid.nr())`) and energy-transforms it with the outgoing-nuclear normalization.
- **V_res(R), Γ(R):** the R-dependent electronic **resonance pole** — at each R, the fixed-R
  electronic ECS Hamiltonian's complex resonance eigenvalue: `V_res(R) = Re, Γ(R) = −2·Im`
  (`Make_Vres`). Reuse `qscat.ecs.find_resonance_pole` per R (the repo already finds the N₂
  resonance pole). `ε_e(R→∞)` (the DA threshold) is the anion bound electronic state at the
  dissociation limit — the machinery already prototyped in the DA diagnosis.

### Exact-2-D DA (the oracle)
The exact version — the full electron–nuclear coupling, no LCP reduction — via the **existing
`qscat.core` TD propagation** (order-3 Padé under `H_2D`) with a **nuclear-coordinate outgoing
test function**:
- **DA test function:** `Φ_DA(r,R) = φ_e(r)·g_out(R)`, `φ_e` = the anion bound electronic state
  at the dissociation limit (fragment B⁻), `g_out(R)` = an outgoing Gaussian in R placed at
  **large R** (eMoScat NO/F₂ decks: R₀ ≈ 7.5–9.7 bohr) — so it samples Ψ only where the
  fragment state is valid.
- **Correlation:** `c_DA(t) = ⟨Φ_DA | Ψ(t)⟩` (c-product), recorded alongside the VE
  correlations in the same propagation.
- **Transform:** the same Tannor-Weeks form as VE, but with `η_out_DA(E) = ⟨g_out(R) |
  F^out_nuc(k_R, R)⟩`, `F^out_nuc = √(μ/2πk_R)·e^{ik_R R}` the outgoing nuclear function (a
  plane wave, not a Riccati-Hankel), and the SAME `η_in` (electronic incident). `σ_DA =
  π|S_DA|²/2E`.
- **Cost:** the nuclei are heavy (μ ~ 1.4–1.7×10⁴), so dissociation is slow — the propagation
  must be long enough for the nuclear wavepacket to reach `g_out`'s position; a T- and
  box-convergence study is part of the deliverable.

## Deliverables

**D1 — LCP DA** (small, in `projects/n2_td_cross_section` / a diatomic LCP module): add the
DA boundary-flux channel to the existing doorway propagation + energy transform. Produces
`σ_DA(E)` alongside `σ_VE(E)` for N₂/NO/F₂ from the LCP model.

**D2 — V_res(R)/Γ(R) construction** (if not already reusable): the R-dependent complex
resonance curve via `find_resonance_pole` per R + the anion asymptotic state, so the LCP has
its inputs for NO/F₂ (N₂'s may already exist).

**D3 — exact-2-D DA** (`qscat.core.time_dependent` + a nuclear test function): the DA test
function `φ_e ⊗ g_out(R)`, `η_out_DA`, and σ_DA via the existing TD propagation; T/box
convergence for the slow nuclear dissociation. The exact oracle.

**D4 — validation + figures:** LCP-vs-exact σ_DA (the LCP-under-test-vs-2-D-oracle comparison,
the scientific output), unitarity/flux checks, and per-molecule VE+DA figures. Cross-check
against eMoScat's own DA output where feasible.

## Interface (sketch)

```python
# LCP DA (extend the 1-D doorway transform)
def lcp_da_sigma(grid_R, mu, Vd, Gamma, chi, eps, eps_e_inf, E, *, dt, n_steps): ...
    # propagate d_0 on H_res; S_DA = sqrt(K/2πμ) e^{-iKX} ∫ e^{iE_tot t} ψ(X,t) dt; σ=4π³|S|²/2E

# exact-2-D DA (qscat.core.time_dependent)
def da_test_function(tgrid, phi_e, *, r0_out, p0_out, sigma_out): ...   # φ_e(r) ⊗ g_out(R)
def eta_outgoing_nuclear(grid_R, k_R, mu, *, r0_out, p0_out, sigma_out): ...  # <g_out|√(μ/2πk_R)e^{ik_R R}>
# td_ve_cross_section gains DA channels (nuclear test functions) alongside VE
```

## Validation

- **Thresholds:** already validated (F₂ open/exothermic, N₂ closed, NO 0.172) — the DA σ must
  respect them (zero below threshold, onset at the right E).
- **F₂ signature:** strong exothermic DA dominant at low E (the famous F₂ case); N₂: σ_DA ≈ 0
  in [0, 0.2] (channel closed).
- **Unitarity / flux conservation:** `Σ_v'|S_VE,v'|² + |S_DA|² ≤ 1` (deficit = ECS absorption);
  the total captured flux is conserved between VE + DA.
- **LCP vs exact-2-D:** the primary scientific comparison — where the LCP DA holds vs where the
  exact 2-D departs (the same "approximation under test" framing as the VE benchmark).
- **eMoScat cross-check:** compare σ_DA against eMoScat's `CSDA` output where a run is available
  (documented tolerance, not exact — a cross-model check).

## Sub-project decomposition

- **A — LCP DA** (D1 + D2): the boundary-flux channel on the existing doorway TD; the
  R-dependent V_res/Γ inputs for NO/F₂. Cheap, self-contained, validated by unitarity +
  thresholds. Do first.
- **B — exact-2-D DA** (D3): the nuclear test function in `qscat.core` TD + the slow-nuclear
  T/box convergence. The oracle. Larger.
- **C — comparison + figures** (D4): LCP-vs-exact σ_DA, VE+DA per-molecule figures, docs.

Each its own spec → plan → execute → merge.

## Out of scope

- **H₂⁺ dissociative recombination (DR)** — the same machinery with the neutral's MANY bound
  electronic states (a Rydberg series → cutoff) + the Coulomb tail; a later ionic sub-project.
- **The TI driven-equation DA** — diagnosed as structurally wrong for the exact 2-D (fragment
  state validity + ECS absorption); the exact-2-D DA is time-dependent, per eMoScat.
- **Rotational (J>0) / multiple partial waves / the angular coupled-channel DA.**

## Verification

- `uv run pytest -q -m "not slow"` pass; `uv run mypy libs/qscat` 0; `uv run ruff check .` clean.
- LCP σ_DA respects thresholds + unitarity; F₂ strong, N₂ ≈ 0 — gated in tests.
- exact-2-D σ_DA converged in T/box, respecting thresholds/unitarity, compared to the LCP.
- Committed VE+DA figures per molecule; `docs/physics/` DA note extended with results;
  `CLAUDE.md` updated. eMoScat `CSDA` cross-check recorded where run.
```
