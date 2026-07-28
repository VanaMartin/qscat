# H₂⁺ model port + dissociative recombination (DR) — Design Spec

**Date:** 2026-07-28
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — spec for review
**Lifecycle:** `qm-method-lifecycle` — a NEW model (the first **ionic** target) + a new exit channel
(**dissociative recombination**), built on the promoted `qscat.core`/`qscat.model` engine and the
exact-DA machinery (`da_cross_section`, `anion_electronic_states`, `riccati_bessel_en_mass`).

## Context

Everything so far is **neutral** electron–diatomic scattering: the scattering electron sees a
neutral molecule, so the asymptotic channel functions are free Riccati-Bessel/Hankel. **H₂⁺ is an
ion.** The electron sees the **+1 core charge**, so it moves in a **−1/r Coulomb tail**; the
asymptotic functions become **Coulomb** functions, the incident wave is a long-range Coulomb wave
(the electronic grid must reach **~1300 bohr**), and the process is **dissociative recombination**
(DR): `e⁻ + H₂⁺(v) → H + H`, captured into a Rydberg resonance that dissociates on the neutral
curve. There is a **Rydberg series** of neutral electronic exit channels (cut off at a measurable
number, ~3). This model is deliberately **not laptop-runnable at full size** (the electronic grid
alone is huge); it targets the Docker/MUMPS path.

Extracted from eMoScat (`source/coulomb.cpp`, `Model2d/CoupledModel2d.cpp`,
`time_independent_model.cpp`'s `compute_2d_noncoupled_model`, `input/experimental/H2p.json`) via
port-scout. This is sub-project **D** of the DA/DR design spec
(`2026-07-27-da-cross-sections-design.md`), scoped here to **model infrastructure + the full DR
cross section**, validated analytically + on small proxies + Docker-ready (no required full run).

## The physics (atomic units; extracted formulas)

### Coulomb special functions (the ionic generalization of Riccati-Bessel)
Energy-normalized Coulomb functions, Sommerfeld parameter `η = m·z/k`, argument `ρ = k·x`:
- Regular: `F_en(x,k,z,m,l) = √(2m/(πk))·F_l(η, ρ)`; outgoing: `H⁺_en = √(2m/(πk))·(G_l + i F_l)`.
- **z→0 reduction:** `F_l(0,ρ) = ρ·j_l(ρ)`, so `F_en(x,k,0,m,l) = √(2mk/π)·x·j_l(kx) =`
  the existing `riccati_bessel_en` (mass m) — the direct differential-test hook.
- `F`/`G`/`H±` for **complex** ρ (ECS-rotated x) via the Barnett COULCC algorithm; in the DR run
  `η = m·z/k` is **real** (z=−1, k real), only ρ is complex. **Backend: `mpmath.coulombf`/
  `coulombg`** (arbitrary precision, accept complex args — matches COULCC functionally).
- Two eMoScat gotchas to NOT replicate: (i) its `sH1` wrapper has a copy-paste bug (returns F, not
  G+iF) — define `H⁺ = G + iF` correctly; (ii) COULCC's `ifail` is ignored — our port checks
  mpmath convergence. `H⁺ = G + iF` outgoing, `H⁻ = G − iF` incoming (confirmed).

### The H₂⁺ potential (`Model2d/CoupledModel2d.cpp`, `time_independent_model.cpp`)
- **σ capture interaction** (the H₂⁺ analogue of `v_int`), `x=r` electronic, `y=R` nuclear:
  `Q(R) = (a₂ − R − a₃R⁴)/7`, `S(R) = tanh(R/a₄)⁴`, `E(r) = e^{−r²/3}/r`,
  `V_int(r,R) = −a₁·(1 − tanh Q(R))·S(R)·E(r)`, with `a₁=1.6435, a₂=6.2, a₃=0.0125, a₄=1.15`.
  (Singular `1/r` at r=0; the grid starts at r>0.)
- **Ion nuclear curve** (a pure Morse — the `1/R` proton repulsion is folded in, NOT explicit):
  `v0(R) = h2p_v_zero(R) = V₀(e^{−2α(R−R₀)} − 2e^{−α(R−R₀)})`, `V₀=0.1027, R₀=2.0, α=0.69`
  (min −0.1027 Ha at R=2, → 0 as R→∞). The initial vibrational state lives in this well.
- **Electronic surface:** `surface(r,R) = v0(R) + V_int(r,R) + ℓ(ℓ+1)/(2r²) − 1/r`, **ℓ=1**
  (p-wave), the **−1/r** electron–core Coulomb attraction (the ionic term). `μ=918.25`, electron
  mass 1, charge `z=−1`.

### The DR cross section (`compute_2d_noncoupled_model`)
`H_2D = −½∂²_r − (1/2μ)∂²_R + surface(r,R)` (evaluated on the **complex** ECS coordinates — we fix
eMoScat's real-part `// FIXME`). Then:
1. **Initial vibrational state:** diagonalize `T_nuc(μ) + diag(v0(R))` → `E_vib = eps[v_init]`,
   `χ_{v_init}` (tail zeroed) — exactly `vibrational_states`.
2. **Rydberg exit states:** diagonalize the electronic Hamiltonian at `R_inf = nuclear R0`,
   `−½∂²_r + surface(r, R_inf)` → its **bound** eigenstates (electron bound to the +1 core via
   −1/r) = the Rydberg channels `φ_e^(n)`, `E_ryd(n) = ε_e^(n)` (n = ri_min..ri_max, ~3) — exactly
   `anion_electronic_states` with `n_states` = the channel count (below the continuum threshold).
3. **Coulomb incident:** `Ψ_i = F_en(r, k, z=−1, m=1, ℓ) · χ_{v_init}` (masked) — the Coulomb
   generalization of `channel_vector` (charge-aware).
4. **Driven solve:** `Ψ₊ = Ψ_i − (E_tot·I − H_2D)⁻¹ V_int Ψ_i`, `E_tot = E + E_vib`, `V_int` =
   the σ capture interaction — reuses the `ve_cross_section` driven sweep (`SparseLU.refactor`).
5. **Rearrangement interaction:** `V_DR(r,R) = V_int(r,R) + v0(R) − V_int(r, R_inf)` (→0 as
   R→∞) — identical in form to the DA `v_dr_diag`.
6. **Per Rydberg channel n** (open iff `E_tot − E_ryd(n) > 0`): `E_DR = E_tot − E_ryd(n)`,
   `K = √(2μ E_DR)`, exit `Φ_n = φ_e^(n)(r)·F^nuc_{K,0}(R)`, `F^nuc` = `riccati_bessel_en_mass`
   (l=0, mass μ, truncated on the far tail); `T_n = ⟨Φ_n | V_DR | Ψ₊⟩` (**c-product** — eMoScat
   uses a conjugated dot `zdotc`, but the ECS-correct choice is the c-product, as everywhere in
   qscat; a validation check confirms the rotated-tail difference is negligible),
   `σ_n = 4π³|T_n|²/(2E)`.

This is `da_cross_section` **generalized**: (a) a Coulomb incident instead of Bessel, (b) a LOOP
over `n_channels` Rydberg exit states instead of one anion state. The T-matrix, `V_DR`, driven
solve, and nuclear-Bessel exit are the SAME.

## Deliverables

**D1 — Coulomb special functions** (`qscat.special.coulomb`): `coulomb_f_en`, `coulomb_g_en`,
`coulomb_h1_en` (energy-normalized, complex-arg, mpmath), validated by z→0→Bessel, the Wronskian
`F'G − FG' = k`, and known real values.
**D2 — the ionic model** (`qscat.model`): a `charge` attribute added to the `ResonanceModel`
protocol (0 for the neutral diatomics, −1 for H₂⁺); a new ionic model form (H₂⁺'s Morse `v0` +
σ-capture `v_int` + `−1/r` Coulomb surface) with an `H2P` registry entry; the neutral
`DiatomicResonanceModel` gains `charge=0` (unchanged behavior).
**D3 — Coulomb channel functions** (`qscat.core.channels`): `channel_vector` becomes charge-aware
— `z=0` keeps the fast `riccati_bessel_en` path (neutrals unaffected), `z≠0` uses `coulomb_f_en`.
**D4 — the DR cross section** (`qscat.core`): `dr_cross_section(tgrid, model, eps, chi, v_init,
n_channels, E, ...)` — the Rydberg-loop + Coulomb-incident generalization of `da_cross_section`,
reusing `v_dr_diag`, `anion_electronic_states` (n_states=n_channels), the driven sweep, and the
nuclear Bessel.
**D5 — discretization + Docker-readiness**: the H₂⁺ grids (electronic real→1300 + exp-ECS tail 5°;
nuclear real→14 + exp-ECS tail 22°; order 8) via `segmented_grid`/the exp-tail helper, in a
per-model config; a **reduced-grid smoke run** that exercises the whole DR path on a laptop, and
the full deck wired for Docker/MUMPS (not run here).

## Interface (sketch)

```python
# qscat.special.coulomb
def coulomb_f_en(x, k, z, m, l): ...   # sqrt(2m/pi k) F_l(m z/k, k x); z=0 -> riccati_bessel_en
def coulomb_h1_en(x, k, z, m, l): ...  # outgoing G_l + i F_l, energy-normalized
# qscat.model: ResonanceModel gains `charge: int` (property); new ionic model + H2P registry entry
# qscat.core.channels
def channel_vector(tgrid, k, chi_v, l, *, charge=0): ...   # Coulomb when charge != 0
# qscat.core (dissociation/recombination)
def dr_cross_section(tgrid, model, eps, chi, v_init, E, *, n_channels=3, ordering="COLAMD"): ...
    # Psi+ from the driven solve (Coulomb incident); per Rydberg channel n:
    #   T = <phi_e^(n) (x) F_nuc(K_n) | V_DR | Psi+>;  sigma_n = 4 pi^3 |T|^2 / 2E
```

## Validation (given non-laptop full size)

- **Coulomb functions (analytic):** z→0 reproduces `riccati_bessel_en`/`riccati_hankel_en` to
  ~1e-10; the Wronskian `F_l' G_l − F_l G_l' = k` (mass-1) / `= 1` in ρ; a few tabulated
  Abramowitz-Stegun values; behavior on a complex (ECS-rotated) argument is finite/smooth.
- **The model:** `v0`/`v_int`/`surface` match the extracted formulas elementwise; `charge=−1`; the
  Rydberg electronic states at R_inf are bound (real, below the −1/r continuum) and number ≥
  n_channels.
- **DR well-posedness (small proxy):** on a REDUCED grid (electronic r_max ~60, nuclear ~14 — big
  enough to hold a couple of Rydberg states and the incident wave, small enough for a laptop
  SuperLU solve), `dr_cross_section` is finite, ≥0, respects each channel's threshold (σ_n=0 below
  `E_ryd(n)−eps[v_init]`), and the c-product-vs-conjugated-dot difference is negligible. This is a
  well-posedness/threshold gate, NOT a converged cross section (the real grid is 1300 bohr).
- **Docker/MUMPS:** the full deck builds and a smoke run starts under the container (a couple of
  energies), demonstrating the path; a converged σ_DR(E) curve is an explicit follow-on.
- **No independent golden data** ships (eMoScat's `output/H2+/sigma.txt` is absent from the
  snapshot); the exact solver is the oracle, as for the neutral DA.

## Sub-project decomposition (tasks for the plan)

1. **Coulomb special functions** (D1) — `qscat.special.coulomb` + analytic tests.
2. **The `charge` protocol + ionic model** (D2) — `ResonanceModel.charge`, `DiatomicResonanceModel.charge=0`, the ionic model form + `H2P` registry.
3. **Charge-aware `channel_vector`** (D3) — Coulomb dispatch; neutrals unchanged (regression-gated).
4. **`dr_cross_section`** (D4) — Rydberg-loop + Coulomb-incident generalization of `da_cross_section`.
5. **H₂⁺ discretization + config** (D5a) — the grids (exp-ECS tail) + per-model config + the reduced-proxy grid.
6. **Small-proxy DR validation + Docker-ready + docs** (D5b) — the well-posedness/threshold gate on the reduced grid, the full deck Docker-wired, a smoke run, and the physics-note section.

## Out of scope (this sub-project)

- **A converged full-size σ_DR(E) curve / figure** — needs the 1300-bohr Docker/MUMPS run; the
  small-proxy gate + Docker-ready deck are delivered, the converged run is a follow-on.
- **The π channel** (`p_pi_potential`) — the committed DR is σ-only.
- **Optimizing the Coulomb functions** (Rust / a direct COULCC port) — mpmath first; optimize if
  it becomes the bottleneck.
- **Rotational (J>0), coupled-channel (non-adiabatic) DR.**

## Verification

- `uv run pytest -q -m "not slow"` pass; `uv run mypy libs/qscat/qscat` 0; `uv run ruff check .` clean.
- Coulomb functions match the analytic limits; the ionic model matches the extracted formulas and
  keeps the neutral models' behavior (charge=0); `channel_vector` charge=0 is bit-unchanged.
- `dr_cross_section` is well-posed on the reduced proxy, respects per-channel thresholds, and the
  c-product choice is justified; the full H₂⁺ deck builds + smoke-runs under Docker/MUMPS.
- `docs/physics/` gains an H₂⁺-DR section; `CLAUDE.md` updated; the DA/DR spec's sub-project D
  marked delivered.
