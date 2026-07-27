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

## Time-dependent route (next)

The N₂ time-dependent route (order-3 Padé + Tannor-Weeks, `qscat.core.time_dependent`) matches
the exact TI to ~1–2 % across N₂'s broad resonance. For NO and especially F₂ the resonances are
**much sharper and lower** — their ~0.004–0.01 Ha boomerang features are at or below the
finite-time-propagation resolution `2π/T`, so a TD-vs-TI reproduction to the N₂-level tolerance
needs a dedicated long-propagation convergence study (the near-threshold limit already documented
for N₂'s own low-E edge — see `docs/physics/n2-2d-td-cross-section.md` and issue #1). That study
is the natural follow-on; the exact-2D TI oracle delivered here is what it would be gated against.
