# Nonlocal resonance model — vibrational excitation

Date: 2026-08-18

## Purpose

Extend `qscat.core.nrm` from dissociative attachment to **vibrational excitation**,
completing the time-independent nonlocal resonance model of Houfek, Rescigno &
McCurdy, Phys. Rev. A **77**, 012710 (2008). This is "spec 2", scoped out of
`2026-08-17-nrm-ti-da-design.md` at the time.

VE is not merely another observable. It is the process for which **PRA 77 publishes
curves for the nonlocal model itself**, and it is the only in-repo way to test the
nonlocal potential `F(E,R,R')` against a published oracle.

That matters because of how spec 1 ended. On F₂ the model reproduced the exact 2-D
oracle to 0.06–0.33 %; on NO it collapsed by 5–8 orders, and the collapse could not be
resolved — an equation-by-equation audit found the code correct, and eight hypotheses
and two mechanisms were killed by measurement (`docs/physics/nonlocal-resonance-model.md`
§8). The investigation ran out of oracles precisely because the paper publishes **no**
NO or N₂ dissociative-attachment cross section at all. VE reverses that: it is O(1)
rather than exponentially sensitive near threshold, and the paper plots it for every
molecule in this study.

## Reference literature

Primary: `reference/literature/houfek-2008-pra77-012710.md` (PRA 77, 012710). **Every
equation number below refers to that paper.** Eq. (55)–(61) and (69) were transcribed
into that note during spec 1; Eq. (28)–(38) — the T-matrix decomposition this spec
implements — are **not yet anchored there** and must be added as part of this work,
under the `mastering-references` skill.

Supporting: `reference/literature/domcke-1991-physrep208-97.md` — PRA 77 disagrees with
its Eq. (4.14) on the coupling in the resonant T-matrix (`V_dk` vs `V_{d,−k̄}`,
p. 012710-4). We implement PRA 77's form; that disagreement bears directly on Eq. (31)
below, so the note must be re-read before the conjugation convention is fixed in code.

## The method

The VE T-matrix splits into background and resonant parts:

```
T^VE_{vi→vf} = T^bg + T^res                                          Eq. (28)

T^res_{vi→vf} = ⟨χ_vf | V*_{dk_f} | Ψ_d⁺⟩                             Eq. (31)

T^bg_{vi→vf} = ⟨χ_vf φ⁻_{k_f}| V_int |χ_vi 𝒥_{k_i}⟩
               − ⟨χ_vf | V*_{dk_f} 𝒥_{dk_i} | χ_vi⟩                   Eq. (37)

𝒥_{dk}(R) = ∫ dr φ_d(r;R) 𝒥_k(r)                                     Eq. (38)

σ_{vi→vf}(E) = 4π³ |T^VE|² / k_i²
```

`Ψ_d⁺` is the solution of the nuclear equation Eq. (52) — **the same object
`nrm_da_cross_section` already computes**. `V_dk⁺` is the existing `coupling.v_dk_plus`,
evaluated at the final momentum `k_f` where `k_f²/2 = E_tot − ε_vf`.

So `T^res` costs one coupling evaluation per final channel on top of machinery that
already exists and is gated. `T^bg` is the new work.

### The conjugation in Eq. (31) and (37)

`V*_{dk_f}` is written with a conjugate. Spec 1 established (and the equation audit
confirmed against p. 012710-6) that under ECS the pairing is the **bilinear c-product,
without conjugation**, and that PRA 77's Eq. (59)/(60) fix that convention. The same
reasoning applies here, and PRA 77's own footnote at p. 012710-4 — its disagreement with
Domcke's Eq. (4.14) — is precisely about this element. The implementation uses the
c-product form; the physics note must record the choice and its justification, because it
is the single easiest place in this spec to be silently wrong.

### `φ⁻` — the one open construction

`T^bg` needs the P-space continuum with **incoming** boundary conditions,
`φ⁻_{k_f}(r;R)`, at every nuclear `R`. Eq. (34)–(36) give `φ_k⁻ = (φ_k⁺)*` for a real
discrete state in the radial case, and `φ_k⁻ = (φ_{−k̄}⁺)*` in general.

Under exterior complex scaling conjugation is not innocent: it maps an outgoing wave to
an incoming one but also conjugates the contour. The production path is therefore
**`φ⁻ = (φ⁺)*` restricted to the real region**, where the contour is real and where the
integrand of Eq. (37) lives — `V_int` is a Gaussian in `r`, and spec 1 measured couplings
decaying to ~1e-12 by the real-region edge.

That restriction is an assumption, so it is **checked, not asserted**: a second genuine
driven solve with an incoming reference wave serves as the differential oracle, and the
two must agree wherever `V_int` has support. If they disagree, this spec stops and
reports rather than choosing the convenient one.

## Scope

**In scope**

- `T^res` (Eq. 31), `T^bg` (Eq. 37), `𝒥_dk` (Eq. 38), `φ⁻`, and `σ_VE`.
- Transitions **0→0 and 0→1**.
- Discrete-state choices **A** (`PhysicalDiscreteState`) and **B**
  (`AsymptoticDiscreteState`).
- Molecules **N₂ and F₂**.
- A committed N₂ figure, laid out like the paper's Fig. 4.
- Eq. (28)–(38) added to the tracked reference note.

**Out of scope**

- 0→8 and other high overtones (needs more vibrational states and a wider window).
- NO (a follow-on; see "What this does not settle").
- The time-dependent NRM (spec 3).
- Discrete-state choice C.
- Any attempt to resolve the NO **DA** collapse. This spec may inform it and must not
  pretend to close it.

## Which comparisons have a published anchor

Not every combination is plotted in the paper, and the spec must not imply otherwise:

| Molecule / transition | Choice A | Choice B |
|---|---|---|
| N₂ 0→0 | Fig. 4 (top) | Fig. 8 (top-left) |
| N₂ 0→1 | Fig. 4 (middle) | **not plotted** — Fig. 8 omits it, stating results of all calculations are "practically the same in this particular case" (p. 012710-10) |
| F₂ 0→1 | Fig. 6 (top) | Fig. 8 (bottom-left) |
| F₂ 0→0 | **not plotted** | **not plotted** |

Three of four pairs have a published anchor. F₂ 0→0 is gated against the exact solver like the
others but carries no *external* corroboration, and the note must say so. N₂ 0→1's *absence* from Fig. 8 is itself a
published statement worth testing: all four routes should agree there.

## Architecture

New module `libs/qscat/qscat/core/nrm/vibrational_excitation.py`:

| Function | Responsibility |
|---|---|
| `j_dk(elec_grid, phi_d, R_values, energy, ell)` | Eq. (38): `∫dr φ_d(r;R) 𝒥_k(r)`, one entry per `R` |
| `t_resonant(nuclear_grid, chi_f, v_dk_f, psi_d)` | Eq. (31), a contraction over the nuclear grid |
| `t_background(...)` | Eq. (37), both terms |
| `nrm_ve_cross_section(...)` | assembles `σ_{v→v'}(E)`; `include_background: bool = True` |

One addition to `scattering.py`: `scattering_state_minus`, the incoming-boundary
solution, used as the oracle for the conjugation shortcut and available for the general
case.

**`include_background` is a first-class parameter, not a debug flag** — the paper's whole
argument is the difference between the two curves, and the figure plots both.

Everything else is reused unchanged: `nrm_ingredients`, `nonlocal_operator`,
`solve_nuclear`, `v_dk_plus`, `discrete_state`. The cross-section normalisation is
`driven.py`'s existing `4π³|T|²/k_i²`, so the exact and NRM curves are compared on
identical footing rather than through two conventions that happen to agree.

`qscat.core` must not import `qscat.model` at runtime; the guard added in spec 1's final
fix wave (`pkgutil.walk_packages`) covers the new module automatically.

## Validation

**No figure digitisation.** The exact 2-D solver is already anchored to Houfek's
published `CSVE.V00.J00` data at rtol 1e-3 (`validation/n2`), and PRA 77's claim is that
nonlocal+background reproduces the exact result. Testing NRM+bg against our own exact
solver therefore inherits the published anchor instead of re-deriving it from figure
pixels.

Four checks, in increasing strength:

1. **`Ψ_d⁺` consistency.** The VE and DA paths must produce the identical `Ψ_d⁺` for the
   same molecule, choice and energy — they call the same solver, and a divergence means
   one of them is mis-wiring the nuclear equation.
2. **`φ⁻` agreement.** The conjugation shortcut against the genuine incoming-boundary
   solve, over the region where `V_int` has support. This gates the shortcut; a failure
   blocks the layer.
3. **Background is load-bearing, and its importance decreases with inelasticity.**
   PRA 77 states both (p. 012710-8, and Fig. 4's own panels). So `σ` with and without
   `T^bg` must differ materially for 0→0, and differ *less* for 0→1. This is a
   falsifiable published prediction, not a self-consistency check.
4. **The result: NRM+bg against the exact 2-D oracle** for all four molecule/transition
   pairs (N₂ 0→0, N₂ 0→1, F₂ 0→0, F₂ 0→1), both discrete-state choices. The oracle is
   our own exact solver, which exists for every transition — what F₂ 0→0 lacks is the
   *paper's* curve, not an oracle, so it is gated like the rest and only its external
   corroboration is missing. Bands recorded from the first converged run and
   justified in the physics note — **never preset**, and a sentinel left in the gate is
   a task failure. This follows spec 1's precedent exactly.

Additionally, N₂ 0→1 should show all four routes agreeing, which is what Fig. 8's
omission asserts.

## Cost

Measured deck sizes: N₂'s VE deck is 107 electronic × 251 nuclear (26,857 unknowns);
F₂'s is larger. For N₂, `F(E)` is ~106 inversions of a 251×251 complex matrix — roughly
0.2 s per energy, against 7.5 s for F₂'s 974-node DA deck. `T^bg` adds ~251 solves of
107×107 per final channel per energy, which is negligible beside `F`.

A dense N₂ sweep for the figure is therefore affordable, which is why the figure is
specified as a curve rather than a handful of anchors. F₂ VE runs on its own `ti` deck
and is sized during implementation; if it proves expensive, the F₂ gate drops to anchors
and the note says so.

## Deliverables

- `libs/qscat/qscat/core/nrm/vibrational_excitation.py` + `scattering_state_minus`.
- Tests for each layer, and a `validation/` comparison for N₂ and F₂.
- `docs/physics/figures/n2-ve-nrm-vs-exact.png` — N₂ 0→0 and 0→1, four curves each
  (exact / LCP / nonlocal / nonlocal+bg), laid out like Fig. 4.
- Eq. (28)–(38) added to `reference/literature/houfek-2008-pra77-012710.md`.
- `docs/physics/nonlocal-resonance-model.md` extended with a VE section: the measured
  comparison, the conjugation convention and its justification, the `φ⁻` construction and
  its checked assumption, and which pairs have published anchors.
- `qscat-run`: the `nrm` method extended to the `ve` observable.

## What this does not settle

The NO **dissociative-attachment** collapse remains open. If N₂ and F₂ VE both reproduce
the exact oracle, that is evidence the nonlocal potential and the T-matrix machinery are
sound, which narrows the NO DA question — but it does not close it, and the physics note
must not claim otherwise. Running **NO VE**, where the paper does publish curves
(Fig. 5, Fig. 8), is the natural follow-on and is deliberately out of scope here.

## Risks

- **`φ⁻` may not reduce to a conjugation** on this discretisation. Check 2 catches it;
  the fallback is the genuine second solve, at roughly double the `T^bg` cost — which is
  small in absolute terms.
- **The conjugation convention in Eq. (31)/(37)** is the easiest thing to get silently
  wrong, and σ depends on `|T|²`, which hides phase errors — spec 1 measured σ_DA as
  blind to the coupling phase. Mitigation: `T^res` and `T^bg` are compared *before*
  squaring, in the layer tests.
- **F₂ VE cost** is unmeasured. If it dominates, the F₂ gate degrades to anchors rather
  than the spec expanding.
