# Dissociative-attachment (DA) + dissociative-recombination (DR) cross sections — Design Spec

**Date:** 2026-07-27
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — implementation pending
**Lifecycle:** `qm-method-lifecycle` — stage 2–4 (toy → validate) for a NEW capability (the
dissociation exit channels), on top of the promoted exact-2-D engine (`qscat.core`) and the
existing 1-D LCP TD (`projects/n2_td_cross_section`).

## Context

The library computes only the **vibrational-excitation (VE)** exit channel so far. The transient
anion `AB⁻*` has a second fate — **dissociative attachment (DA)**: `e⁻ + AB(v=0) → AB⁻* → A + B⁻`,
outgoing flux in the **nuclear** coordinate. The ionic analogue is **dissociative recombination
(DR)**: `e⁻ + AB⁺ → A + B` (H₂⁺). We need these for both the exact-2-D model (the oracle) and the
approximate LCP model (the approximation under test — [[research-program-is-method-validation]]).

**eMoScat computes DA and DR exactly and TIME-INDEPENDENTLY** — the source for the method:
- **exact-2-D TI DA/DR:** `source/time_independent_model.cpp` (`compute_2d_noncoupled_model`,
  H₂⁺ DR with multiple Rydberg channels; the same T-matrix does DA for one channel).
- **LCP DA:** `source/module_LCP.cpp` + `source/ModelLCP/SMatrix.cpp` (the doorway wavepacket).

**Correction to an earlier note:** a prior prototype of mine called the exact-2-D TI DA
"structurally broken." That was WRONG — the prototype used the *VE* interaction `V_int` in the
DA T-matrix, giving a ~10⁶ unitarity violation. eMoScat uses the **rearrangement interaction**
`V_DR = H − H_final` (below); with it, the TI DA is correct (verified: σ_DA is O(1) bohr²,
mostly within the unitarity cap for F₂; N₂/NO closed). There is NO structural obstacle to a TI
DA — it is the natural, exact route, and it generalizes directly to DR.

Thresholds already validated: `ε_e(R→∞) − eps[0]` = **F₂ −0.069 Ha (exothermic/open), N₂ +0.502
(closed), NO +0.172** — one bound anion electronic state each; H₂⁺ DR has a Rydberg *series* of
neutral electronic states.

## The physics (from eMoScat `time_independent_model.cpp`)

### Exact-2-D TI DA/DR (the oracle)
The **same driven equation as VE**, projected onto the dissociation channel with the correct
rearrangement interaction:
1. **Driven solve:** `Ψ₊ = Ψ_i + (E_tot·I − H)⁻¹ (V_int · Ψ_i)` — identical to
   `qscat.core.driven.ve_cross_section` (reuse it via `return_wavefunction`). Incident
   `Ψ_i = channel_vector(k, χ_0, l)` (neutral) or a **Coulomb** incident `coulomb::sF_en` (H₂⁺).
2. **Rearrangement interaction:** `V_DR(r,R) = V_int(r,R) + v0(R) − V_int(r, R→∞)` = `H − H_final`
   (H minus the final-channel asymptotic Hamiltonian: electron bound in `V_int(r,∞)`, free nuclei
   on `v0`). This — not `V_int` — is the operator in the DA/DR T-matrix.
3. **Exit channel(s):** for each bound electronic state `φ_e^(n)(r)` at the dissociation limit
   (DA: n=0, the single anion state; DR: n = the Rydberg series, cut off), with `E_n = ε_e^(n)`,
   `K_n = √(2μ(E_tot − ε_e^(n)))`, the exit channel is `Φ_n(r,R) = φ_e^(n)(r) · F^nuc_{K_n,0}(R)`,
   `F^nuc` = the energy-normalized regular nuclear Bessel (l=0, mass μ; eMoScat `bessel::s_jEn`,
   `= √(2μ/πK)·sin(K R)`), masked to the real region.
4. **T-matrix + σ:** `T_n = ⟨Φ_n | V_DR | Ψ₊⟩` (c-product), `σ_n = 4π³|T_n|²/(2E)` (= π|S|²/2E, the
   same T-vs-S convention as VE; open only where `E_tot − ε_e^(n) > 0`).

Reuses everything: the sparse driven solve (`SparseLU.refactor` sweep), `channel_vector`, the
anion-electronic-state eigensolve (already prototyped), the σ prefactor. The additions are `V_DR`,
the nuclear Bessel, and the per-channel exit projection.

**Numerical caveat (a real part of the work):** the nuclei are heavy (μ ~ 10³–10⁴), so `K` is
large (F₂: K~50, wavelength ~0.13 bohr) — the nuclear grid must resolve the fast outgoing Bessel.
A finer nuclear real region (or a mass-adapted outgoing representation) and a convergence study
are needed; the prototype's near-threshold unitarity overshoot is this resolution effect.

### LCP DA (the approximation)
The 1-D reduction — the doorway `d_0(R) = √(Γ(R)/2π)·χ_0` propagated on `H_res = T_nuc + diag(V_res
− iΓ/2)` (already in `projects/n2_td_cross_section`). DA = the flux surviving to the dissociation
boundary: `S_DA(E) = √(K/2πμ)·e^{−iKX}·∫ e^{iE_tot t} ψ(R_max,t) dt` (`SMatrix.cpp`), `σ = 4π³|S|²/2E`.
Needs the R-dependent `V_res(R)/Γ(R)` (the fixed-R electronic resonance pole, via
`qscat.ecs.find_resonance_pole`).

## Deliverables

**D1 — exact-2-D TI DA** (`qscat.core.driven` extension): `V_DR`, the nuclear Bessel exit channel,
the anion-state helper, and `da_cross_section(tgrid, model, ...)` returning σ_DA(E) via the reused
driven solve. Includes the nuclear-grid convergence for the fast outgoing wave. **The oracle.**
**D2 — LCP DA**: the boundary-flux channel on the existing 1-D doorway + the R-dependent V_res/Γ
inputs for NO/F₂. **The approximation.**
**D3 — exact-2-D DR (H₂⁺)**: the same D1 T-matrix looped over the Rydberg electronic channels +
the Coulomb incident/exit (`qscat.special`/`coulomb::sF_en`, the ECS Coulomb functions). Deferred
after DA if the Coulomb tail needs its own groundwork, but it is the *same method*, not a new one.
**D4 — validation + figures**: unitarity, thresholds, F₂-strong/N₂-closed, LCP-vs-exact σ_DA, and
per-molecule VE+DA figures; eMoScat cross-check where runnable.

## Interface (sketch)

```python
# qscat.core: exact-2-D TI DA (reuses the driven solve)
def anion_electronic_states(model, electronic_grid, R_inf, n_states=1): ...   # (eps_e[], phi_e[])
def v_dr_diag(tgrid, model): ...                       # V_int + v0(R) - V_int(r, R_inf), flat
def da_cross_section(tgrid, model, eps, chi, v_init, E, *, ordering="COLAMD"): ...
    # Ψ+ from the driven solve; per anion channel: T = <phi_e (x) F_nuc(K) | V_DR | Ψ+>; σ=4π³|T|²/2E
# qscat.special: energy-normalized regular nuclear Bessel (l=0, mass μ)
def riccati_bessel_en_mass(R, K, l, mu): ...           # sqrt(2 mu K/pi) R j_l(K R)
# LCP DA: add the boundary-flux channel to the 1-D doorway transform
def lcp_da_sigma(grid_R, mu, Vd, Gamma, chi, eps, eps_e_inf, E, *, dt, n_steps): ...
```

## Validation

- **Thresholds** (validated): σ_DA = 0 below threshold; onset at the right E.
- **Unitarity:** `Σ_v'|S_VE,v'|² + Σ_n|S_DA,n|² ≤ 1` (deficit = ECS absorption) — the primary
  numerical-correctness gate (the prototype's near-threshold overshoot must resolve away with a
  finer nuclear grid).
- **Signatures:** F₂ strong exothermic DA; N₂ σ_DA ≈ 0 (closed) in [0, 0.2].
- **LCP vs exact-2-D:** the scientific output — where the LCP DA holds vs departs.
- **eMoScat cross-check:** vs `output/H2+/sigma.txt` (DR) / the 2-D `CSDA` where a run is available.

## Sub-project decomposition

- **A — exact-2-D TI DA** (D1): `V_DR` + nuclear Bessel + the driven-solve reuse + the nuclear-grid
  convergence; unitarity/threshold gates; F₂/N₂/NO σ_DA. The oracle; do first (it is the corrected
  method and the highest-value piece).
- **B — LCP DA** (D2): the boundary-flux channel + R-dependent V_res/Γ. The approximation.
- **C — comparison + figures** (D4). **D — H₂⁺ DR** (D3): the Rydberg-channel loop + Coulomb tail.

Each its own spec → plan → execute → merge.

## Out of scope (for the first sub-projects)

- **The Coulomb tail machinery** for H₂⁺ (ECS Coulomb functions, `coulomb::sF_en`) — sub-project D.
- **Rotational (J>0) / multiple electron partial waves / angular coupled-channel DA.**
- A TD exact-2-D DA — unnecessary; the TI route (D1) is exact and cheaper (no long propagation).

## Verification

- `uv run pytest -q -m "not slow"` pass; `uv run mypy libs/qscat` 0; `uv run ruff check .` clean.
- exact-2-D σ_DA converged (nuclear grid), respecting thresholds + unitarity; F₂ strong, N₂ ≈ 0 —
  gated in tests.
- LCP σ_DA respects thresholds/unitarity; LCP-vs-exact comparison + committed VE+DA figures;
  `docs/physics/` DA note extended with results; `CLAUDE.md` updated; eMoScat cross-check recorded.
