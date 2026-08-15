# H₂⁺ dissociative recombination (the first ionic model)

**Location:** the ionic pieces live in `qscat` (`qscat.special.coulomb`, `qscat.model.H2P` /
`IonicResonanceModel`, `qscat.core.dr_cross_section`); the deck + validation in
`validation/h2plus/` (`config.py` — the grids; `dr.py` — the driver; `test_dr.py` — the gate).
**Origin:** sub-project D of the DA/DR design spec
(`docs/superpowers/specs/2026-07-27-da-cross-sections-design.md`), designed in
`docs/superpowers/specs/2026-07-28-h2plus-dr-design.md` and ported from eMoScat via port-scout.
**Units:** atomic units.

## What this is — and why it is different

Everything before H₂⁺ was **neutral** electron–diatomic scattering. **H₂⁺ is an ion**: the
scattering electron sees the **+1 core charge**, i.e. a long-range **−1/r Coulomb tail**. That one
fact changes three things:

1. **The channel functions become Coulomb**, not free Riccati-Bessel. The incident wave is a
   long-range Coulomb wave — so the electronic grid must reach **~1300 bohr**.
2. **The process is dissociative recombination** (DR): `e⁻ + H₂⁺(v) → H + H`. The electron is
   captured into a Rydberg resonance of the neutral H₂, which dissociates on the neutral curve.
3. **There is a Rydberg *series* of exit channels** (the neutral's bound electronic states,
   accumulating at the −1/r continuum edge), cut off at a measurable number (`N_CHANNELS = 3`).

This is deliberately the **first non-laptop-scale** model — the full 2-D deck is ~1.15 million
unknowns (see Discretization). It targets the Docker/MUMPS path; the test suite validates it
analytically and on a reduced proxy.

## Coulomb special functions (`qscat.special.coulomb`)

The charge-z generalization of `riccati_bessel_en`/`riccati_hankel_en`. Energy-normalized regular /
irregular / outgoing Coulomb functions, Sommerfeld parameter `η = m·z/k`, argument `ρ = k·x`:

  `F_en(x,k,z,m,l) = √(2m/(πk))·F_l(η,ρ)`,  `G_en` likewise,  `H1_en = √(2m/(πk))·(G_l + i F_l)`.

Backed by **mpmath** (`coulombf`/`coulombg`), which accepts **complex** ρ — required for
ECS-rotated arguments. At **z=0**, `F_l(0,ρ) = ρ·j_l(ρ)`, so `coulomb_f_en(·,·,0,m,l)` reduces
EXACTLY to `riccati_bessel_en` (mass m) — the direct differential-test hook, verified to ~1e-16.
(Note: `coulomb_h1_en(z=0) = i·riccati_hankel_en`, since `G_l(0,ρ)=−ρ y_l(ρ)` while `F_l=+ρ j_l(ρ)`
give `G+iF = i(j_l+i y_l)ρ` — the code keeps the physical `H⁺=G+iF`; only the *test* carries the
`i`.) Two eMoScat gotchas were NOT replicated: its `sH1` wrapper had a copy-paste bug (returned F,
not G+iF), and its wrappers ignored `coulcc`'s `ifail`.

## The ionic model (`qscat.model.H2P`)

The model layer was generalized for ions: the `ResonanceModel` protocol gained a `charge`
attribute (0 for the neutral diatomics — unchanged — and −1 for H₂⁺) and shed the engine-unused
`lam`. `H2P` is an `IonicResonanceModel` with the extracted eMoScat form (μ=918.25, ℓ=1,
charge=−1):

- **ion Morse** `v0(R) = V₀(e^{−2α(R−R₀)} − 2e^{−α(R−R₀)})`, `V₀=0.1027, R₀=2.0, α=0.69` (the
  initial vibrational state lives here; the `1/R` proton repulsion is folded into this single
  Morse, not explicit);
- **σ-capture** `v_int(r,R) = −a₁(1−tanh Q(R))·S(R)·(e^{−r²/3}/r)`, `Q=(a₂−R−a₃R⁴)/7`, `S=tanh(R/a₄)⁴`,
  `a₁=1.6435, a₂=6.2, a₃=0.0125, a₄=1.15`;
- **surface** `= v0(R) + v_int(r,R) + ℓ(ℓ+1)/(2r²) + charge/r` — the `charge/r = −1/r` is the ionic
  electron–core Coulomb attraction.

Adding an ion is data + validation, no engine changes (the same lesson as the neutral molecules).

## The DR cross section (`qscat.core.dr_cross_section`)

`dr_cross_section` is `da_cross_section` **generalized**, reusing `v_dr_diag`,
`anion_electronic_states`, `riccati_bessel_en_mass`, and the driven Lippmann-Schwinger sweep. The
only new physics is (a) a **Coulomb incident** (`channel_vector(..., charge=−1)`) and (b) a **loop**
over the Rydberg exit channels instead of one anion state:

1. `Ψ₊ = Ψ_i − (E_tot·I − H_2D)⁻¹ V_int Ψ_i`, `E_tot = E + eps[v_init]`, Coulomb `Ψ_i` (potentials
   on the **complex** ECS coordinate — eMoScat's real-part `// FIXME` is fixed);
2. Rydberg states `φ_e^(n)`, `E_ryd(n) = eps_e^(n)` from `anion_electronic_states(…, n_states=N)`
   (they are bound below the −1/r continuum — the same bound-state solver);
3. `V_DR = V_int + v0(R) − V_int(r,R_inf)` (the rearrangement interaction);
4. per open channel n (`E_tot > E_ryd(n)`): `T_n = ⟨φ_e^(n)·F^nuc_{K_n,0} | V_DR | Ψ₊⟩` (c-product),
   `σ_n = 4π³|T_n|²/2E`.

**The c-product is the ECS-correct choice** where eMoScat used a conjugated dot (`zdotc`); the
port validates they agree — the relative difference is **≈3.4×10⁻¹²** on the proxy (the
rotated-nuclear-tail contribution is negligible there), so the convention question is settled.

## Discretization (`fem_grid_exp_tail`, `validation/h2plus/config.py`)

The Coulomb tail forces a huge electronic grid. A new builder `fem_grid_exp_tail` (like
`segmented_grid` but with an exponential-growth ECS tail — reusing the `_ecs_tail` helper) builds
the eMoScat deck:

| | real region | ECS tail | 2-D size |
|---|---|---|---|
| electronic | → **1300 bohr** (n=1406) | 5°, exp-growth ×25 | |
| nuclear | → 14 bohr (n=818) | 22°, exp-growth ×25 | ~**1.15 M unknowns** |

That is firmly **Docker/MUMPS** territory — not laptop-runnable. `config.proxy_grid()` (electronic
→60 bohr) and an even smaller test grid give a laptop-feasible **well-posedness/threshold** gate,
not a converged cross section.

## Validation (analytic + small-proxy + Docker-ready)

- **Coulomb functions**: z→0 reproduces `riccati_bessel_en`/`riccati_hankel_en` (~1e-16); a known
  mpmath value; finite on complex ECS arguments.
- **The model**: `v0`/`v_int`/`surface` match the extracted formulas; `charge=−1`; the neutrals keep
  `charge=0` and `channel_vector(charge=0)` is byte-identical to before.
- **DR (small proxy, `@slow`)**: σ_DR finite, ≥0, correctly shaped; a genuinely **closed** Rydberg
  channel returns exactly 0 (the 3rd channel, threshold ≈0.0426 Ha, above the probe energies); and
  the **c-product vs conjugated-dot** agreement (≈3.4e-12) justifies the convention. This is a
  well-posedness gate, NOT a converged σ_DR — the real grid is 1300 bohr.
- **Docker/MUMPS**: `validation/h2plus/dr.py`'s `main()` is the full-deck smoke path (a couple of
  energies); the **converged σ_DR(E) curve** (needs the container + MUMPS) is now delivered — see
  the next section. No independent golden data ships (eMoScat's `output/H2+/sigma.txt` is absent
  from the snapshot), so — as for the neutral DA — the exact solver is the oracle.

## The converged full-size σ_DR(E) curve (delivered)

Run the full 1300-bohr deck (~1.15 M unknowns) under MUMPS via the reproducible generator
`validation/h2plus/dr_curves.py` (Docker/MUMPS-only, NOT in the test suite — same convention as
`dr.py`'s `main()`):

```
uv run python -m validation.h2plus.dr_curves
```

It writes two figures + their sidecar data (`.npz`/`.csv`) into `docs/physics/figures/`:

- `h2plus-dr-cross-section.png` — the coarse `config.energy_grid()` sweep (0.001..0.050 Ha), first
  `N_CHANNELS = 3` Rydberg exit channels, linear axes.
- `h2plus-dr-cross-section-shortrange.png` — 200 log-spaced energies across the DR1 resonance in
  [0.005, 0.007] Ha, log–log (the accuracy figure).

![H2+ DR cross section, short range (log–log)](figures/h2plus-dr-cross-section-shortrange.png)

The DR1 (n=0) channel peaks at **E ≈ 6.31×10⁻³ Ha, σ ≈ 1.54×10⁻³ bohr²** above a ~10⁻¹⁰
background; DR2 (n=1) is ~10⁻⁶; DR3 (n=2) is closed in this window (threshold ≈ 0.0426 Ha). The
solve runs in ~8 s/energy on the `sadaharu` host with the OpenMP MUMPS backend.

## Follow-ons

The π channel (`p_pi_potential`); optimizing the mpmath Coulomb functions (a Rust/COULCC port) if
they become the bottleneck; rotational / coupled-channel (non-adiabatic) DR.
