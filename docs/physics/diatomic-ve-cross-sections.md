# NO and F₂ exact-2D VE cross sections (the model port)

**Location:** `validation/diatomic/` (`config.py` — per-molecule grid/energy config;
`curves.py` — the exact-2D TI oracle driver + committed figures; `test_diatomic.py`).
**Computed entirely through the promoted library** `qscat.core` + `qscat.model` — the first
consumers of the generalization beyond N₂ (sub-project A). **Origin:** sub-projects B (NO) and
C (F₂) of the diatomic VE-scattering spec. **Units:** atomic units.

## What this is

The exact 2-D time-independent vibrational-excitation cross section σ_{0→v'}(E) for **NO** and
**F₂** — the same model and method as N₂ (`H = −½∂²_r − (1/2μ)∂²_R + v0(R) + l(l+1)/2r² −
λ(R)e^{−α_c r²}`), differing only in parameters (`qscat.model.{NO,F2}`). Adding each molecule
was **data + validation, no new solver code**: a `qscat.model` registry entry + a
`validation/diatomic/config.MoleculeConfig` grid/energy entry + `compute_ti_curve` calling
`qscat.core.driven.ve_cross_section`.

## The oracle (no independent golden data)

Unlike N₂ (gated against Houfek's independent `CSVE.V00.J00`), **no independent cross-section
data ships for NO/F₂** — so the exact-2D TI solver *is* the reference (the research program's
stance: the exact solution is the oracle, the LCP/TD approximations are under test). The
committed σ(E) curves are that oracle:

![NO exact-2D VE cross section](figures/no-2d-ti-cross-section.png)
![F₂ exact-2D VE cross section](figures/f2-2d-ti-cross-section.png)

Both show clear **boomerang** oscillation structure from a low-lying shape resonance, decaying
smoothly at higher E:

| | partial wave l | α_c | neutral vib spacing eps₁−eps₀ | resonance window |
|---|---|---|---|---|
| N₂ | 2 (d-wave) | 0.40 | 0.0124 Ha | ~0.07–0.10 Ha (broad) |
| NO | 1 (P-wave) | 1.00 | 0.0091 Ha | ~0.02–0.05 Ha (sharp) |
| F₂ | 1 (P-wave) | 3.00 | 0.0039 Ha | ~0.01–0.04 Ha (very sharp, near threshold) |

NO and F₂ have **lower, sharper** resonances than N₂: NO's ²Π shape resonance sits low, and F₂
— weakly bound (D₀ = 0.06 Ha, ~1.6 eV), a strong dissociative-attachment system — has an
extremely sharp near-threshold resonance with boomerang features only ~0.004 Ha wide.

**Convergence:** the N₂-style FEM-DVR-ECS grid (electronic r_max = 16, nuclear r_max = 22) is
converged for NO and F₂ as well — the exact-2D σ(E) is unchanged to <1 % at electronic
r_max = 16/24/32 for NO, so the sharp low-E swings are *genuine* resonance structure, not grid
noise.

## The three molecules side by side

![Exact-2D VE cross sections: N2 vs NO vs F2](figures/diatomic-ve-comparison.png)

N₂'s broad d-wave resonance vs NO/F₂'s sharp low-lying P-wave ones, all from the same model
and the same `qscat.core` solver — only the parameters differ.

## Dissociative attachment (DA) — the second exit channel

Beyond vibrational excitation (VE — the electron re-emitted, `e⁻ + AB(v=0) → AB⁻* → e⁻ +
AB(v')`), the transient anion `AB⁻*` can **dissociate**: `e⁻ + AB(v=0) → AB⁻* → A + B⁻`. The
outgoing flux is then in the **nuclear** coordinate (R→∞) rather than the electronic one. eMoScat
measures it with a second test function on the nuclear coordinate (`Model2d/MultiTestFunction2d.cpp`,
`TestFunction2d.cpp`, `Potentials2d.cpp`): the DA exit channel is

  Φ_DA(r,R) = φ_e(r) · F^out_R(R),

where **φ_e(r) is the anion's bound electronic state at the dissociation limit** — the bound
eigenstate of the electronic Hamiltonian `−½∂²_r + v0(R_∞) + l(l+1)/2r² − λ(R_∞)e^{−α_c r²}`
(`Neutral2dPotential` at the outer nuclear edge; λ(R_∞) → λ_inf) — and `F^out_R(R)` is the
outgoing nuclear wave `√(μ/2πk_R)·e^{ik_R R}`, `k_R = √(2μ(E_tot − ε_e))`. The cross section is
`σ_DA = π|S_DA|²/2E` with S_DA from the same Tannor-Weeks transform / driven-equation projection
as VE. This is the eMoScat convention (deck: the first, nuclear-coordinate test function; "1
transversal eigenstate" = the single bound anion electronic state).

**Thresholds (computed here, validated against the physics):** `ε_e(R_∞) − eps[0]` in collision
energy —

| | anion ε_e (Ha) | DA threshold (E_coll) | status |
|---|---|---|---|
| N₂ | −0.243 | +0.502 Ha | **closed** in the measurement range — matching eMoScat *disabling* N₂'s DA |
| NO | −0.060 | +0.172 Ha | opens above the resonance |
| F₂ | −0.127 | **−0.069 Ha (exothermic)** | **open at all E>0** — matching F₂'s famously *strong* dissociative attachment |

There is exactly one bound anion electronic state per molecule (a real eigenvalue below a
complex ECS continuum), matching "1 transversal eigenstate." **F₂'s exothermic DA and N₂'s
closed channel are the key physical validations of the setup.**

**Status of the DA cross-section numerics (open — a focused research sub-project).** The
threshold and branching physics above are validated, but a correct σ_DA(E) is not a
normalization tweak — two natural TI extractions were tried and *diagnosed* as insufficient:

1. **Volume T-matrix `⟨φ_e(R_∞)·F_R | V_int | Ψ₊⟩`** — wrong: the fragment electronic state
   `φ_e(R_∞)` is valid only *asymptotically* (large R), but the integrand concentrates at small
   R (where `V_int` and `Ψ₊` are large and the fragment state is the wrong electronic state),
   giving a `|S_DA|² ~ 10⁶` unitarity violation.
2. **Surface flux extraction** (project `Ψ₊` onto `φ_e`, read the outgoing nuclear current at a
   surface) — wrong: the ECS contour *absorbs* the outgoing flux by design, so there is no
   clean asymptotic outgoing wave at a real-region surface (the current there is ~0).

The correct treatment needs the **R-dependent anion electronic state `φ_e(r;R)`** (the adiabatic
anion curve — the electronic Hamiltonian's bound state at *each* R, which evolves with R), with
the DA nuclear motion governed by the **anion's own potential curve** `ε_e(R)`, not a free plane
wave. That is the full 2-D DA theory (the local-complex-potential / boomerang regime). It is a
well-scoped follow-on — the framework, the anion-state machinery, and the validated thresholds
here are what it builds on; the eMoScat TD DA output (`Model2d/MultiTestFunction2d.cpp`, the
nuclear test function) is the cross-check.

## Dissociative recombination (DR) — H₂⁺, many channels

The same machinery generalizes to **dissociative recombination** for a molecular *ion*
(`e⁻ + AB⁺ → A + B`): the exit channels are the **multiple** bound electronic states of the
neutral AB at large R (a Rydberg series → in principle infinitely many, cut off at a
measurable number). This is the H₂⁺ case — deferred with the Coulomb tail (the ionic
`sphHankel1En` Coulomb branch exists in eMoScat), but structurally it is "DA with N electronic
channels instead of 1."

## Time-dependent route (next)

The N₂ time-dependent route (order-3 Padé + Tannor-Weeks, `qscat.core.time_dependent`) matches
the exact TI to ~1–2 % across N₂'s broad resonance. For NO and especially F₂ the resonances are
**much sharper and lower** — their ~0.004–0.01 Ha boomerang features are at or below the
finite-time-propagation resolution `2π/T`, so a TD-vs-TI reproduction to the N₂-level tolerance
needs a dedicated long-propagation convergence study (the near-threshold limit already documented
for N₂'s own low-E edge — see `docs/physics/n2-2d-td-cross-section.md` and issue #1). That study
is the natural follow-on; the exact-2D TI oracle delivered here is what it would be gated against.
