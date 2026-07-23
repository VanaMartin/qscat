# N₂ Exact 2-D Time-Independent VE Cross-Section — Design Spec (sub-project #6)

**Date:** 2026-07-23
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — implementation pending
**Lifecycle:** `qm-method-lifecycle` stages 1–3. Builds on #1 (`qscat.dvr`/`qscat.ecs`),
#3 (vibrational states, the 1-D LCP TI solver it is compared against) and #5
(`qscat.linalg`, the N-dimensional sparse Hamiltonian layer).

## Context

Sub-projects #1–#4 solved the electron–N₂ vibrational-excitation problem twice, both times
under the **Local Complex Potential** approximation: the 2-D problem was reduced to a 1-D
nuclear calculation via a resonance pole, giving `V_d(R)` and a local width `Γ(R)`. The
results agree with Karel Houfek's golden data within a documented factor-3 cross-model
bound, with **two honestly-documented structural limitations**:

1. **The elastic channel is resonance-only.** The LCP T-matrix omits non-resonant background
   (potential) scattering, so σ_elastic agrees only near the resonance.
2. **The local width has no electron-energy dependence**, giving the model the wrong
   (non-Wigner) threshold law, so σ diverges as ~1/E toward every channel opening.

This sub-project drops the approximation entirely and solves the genuine 2-D problem — the
same one Houfek solved, and the model the LCP was derived *from*. Port-scout archaeology
(`.superpowers/sdd/n2-2d-exact-extraction.md`) established that eMoScat's 2-D potential
surface is *exactly* `v0(R) + v_int(r,R) + l(l+1)/2r²`, i.e. the functions already ported
and verified in `validation/n2/model.py`. **No new potential physics.**

**Why this is worth doing:** the 2-D calculation is a full scattering calculation, so its
elastic T-matrix contains the background scattering the LCP omits, and its width emerges
with the correct energy dependence rather than being imposed. Both documented NOTEs should
close. Whether they do is a falsifiable prediction of this sub-project.

### What is actually being measured (read this before designing any check)

**The model is given, not studied.** This potential surface is a *testbed*, chosen because
it can be solved exactly — it is deliberately not an attempt to describe real N₂. (This is
why the model's `D_0 ≈ 2×` real N₂'s dissociation energy was accepted as-is in sub-project
#3 rather than corrected.) The research goal is to find **where standardly used techniques
actually fail**.

That inverts the usual direction of validation, and it determines what this sub-project's
deliverable is:

- **Houfek's data validates our solver.** It is an independent implementation of the same
  model and method, so agreement means our 2-D code is right. It is *not* the physical
  truth being sought.
- **Once verified, the exact 2-D result becomes the oracle**, and the LCP approximation is
  the thing **under test**.
- **The LCP-vs-exact difference is therefore the scientific result**, not an error to be
  minimized. A large, well-characterized discrepancy is a successful outcome — it is the
  answer to "where does this standard technique fail?"

Consequently no model parameter may ever be tuned to improve agreement with real-molecule
data, and V5/V6/V7 below (approximation vs exact) rank *above* V4 (solver vs Houfek) in
scientific importance, even though V4 must pass first for the rest to mean anything.

## Method (atomic units, bohr²)

**The Hamiltonian**, on a tensor product of two FEM-DVR-ECS grids (electronic `r`,
nuclear `R`), assembled by `qscat.dvr.hamiltonian_nd`:

```
H_2D = -(1/2)d²/dr²  -  (1/2μ)d²/dR²  +  v0(R) + l(l+1)/(2r²) - λ(R)e^{-α_c r²}
```

with `μ = 12766.36`, `l = 2`, a single fixed partial wave, no rotation. Both coordinates
carry an ECS tail: scaling `r` imposes the outgoing-electron boundary condition, scaling
`R` absorbs dissociative-attachment flux. `H_2D` is complex **symmetric**, never Hermitian.

**The driven (two-potential Lippmann–Schwinger) equation.** The entrance-channel
Hamiltonian is `H_0 = T_r + l(l+1)/2r² + T_R + v0(R)`; the *only* perturbation is the
electron–molecule interaction

```
V_int(r,R) = -λ(R)·exp(-α_c r²)
```

Note carefully that `v0(R)` and the centrifugal term belong to the channel, **not** to the
perturbation — putting them in `V_int` would be a physics error, not a convention choice.
Then, for collision energy `E` and initial vibrational state `v` (`E_tot = E + ε_v`):

```
Ψ_i   = F_{E,l}(r) ⊗ χ_v(R)                     masked to the unscaled region
Ψ_sc  = (E_tot - H_2D)⁻¹ · [V_int · Ψ_i]        one sparse LU per energy
Ψ⁽⁺⁾  = Ψ_i + Ψ_sc
T_{v→v'} = c_product( χ_{v'} ⊗ F_{E',l},  V_int · Ψ⁽⁺⁾ )     masked
σ_{v→v'}(E) = 4π³|T|²/k²,   k² = 2E,  k'² = 2(E_tot - ε_{v'})
```

`σ = 0` for closed channels (`E_tot - ε_{v'} ≤ 0`).

**The channel functions.** `F_{E,l}` is the **energy-normalized** regular free radial
solution — for a neutral target the Coulomb function reduces to a Riccati–Bessel function:

```
F_{E,l}(r) = sqrt(2/(πk))·(kr)·j_l(kr) = sqrt(2k/π)·r·j_l(kr)      (m = 1)
```

matching eMoScat's `sphBesselJEn` (`source/bessel.cpp:50`) and `sF_en`
(`source/coulomb.cpp:75`). `scipy.special.spherical_jn` supplies `j_l`.

**Two conventions that must not be gotten wrong** (both are documented traps this repo has
already been bitten by):

- **DVR coefficients.** The FEM-DVR basis is pre-normalized by `1/sqrt(w)`, so a *function*
  becomes coefficients as `c_j = f(r_j)·sqrt(w_j)` using the **bridge-summed complex**
  weight — `qscat.dvr.TensorGrid.sqrt_weights()`, added in #5 for exactly this. With both
  sides in coefficient form, `c_product` *is* the quadrature integral, with no extra weights.
- **ECS masking.** `Ψ_i` and every channel projection must be zeroed outside the unscaled
  region (`TensorGrid.real_mask()`). A projection that extends onto the scaled tail is
  meaningless. eMoScat does this explicitly and only gets away with a Hermitian dot product
  (`cblas_zdotc`) because of it.

**Elastic and inelastic share one formula.** Houfek's elastic column uses `π|S−1|²/k²`
while ours uses `4π³|T|²/k²`. With `S = 1 − 2πiT`, `|S−1|² = 4π²|T|²`, so the two are
identical. Unlike the LCP case, the 2-D elastic T-matrix includes background scattering.

## Interface

```
projects/n2_2d_cross_section/
├── electronic_grid.py   n2_electronic_grid(*, r_max, angle_deg, order, n_complex, ...)
│                        -> FemDvrEcsGrid   (parametrized: the convergence study varies it)
├── channels.py          riccati_bessel_en(r, k, l) -> F_{E,l} values
│                        channel_vector(tgrid, F_vals, chi) -> masked DVR coefficients
├── driven.py            build_h2d(tgrid, ...) -> csr;  solve_driven(H, E_tot, rhs) -> Psi_sc
└── cross_section_2d.py  ve_cross_section_2d(..., E, vprimes) -> sigma[...]
```

Reuses `qscat.dvr` (`TensorGrid`, `hamiltonian_nd`, `real_mask`, `sqrt_weights`),
`qscat.linalg` (`SparseLU`, `c_product`), `projects.n2_ti_cross_section`
(`n2_nuclear_grid`, `vibrational_states`) and `validation.n2.model` (`v0`, `lam`, `v_int`).
Nothing new is promoted to `qscat` in this sub-project — the reusable layer was #5.

**χ_v renormalization.** `qscat.dvr.eigen` returns eigenvectors with numpy's Hermitian
`v†v = 1`. Under ECS the correct normalization is the c-product `vᵀv = 1`. For strictly
bound vibrational states these differ negligibly, but the solver must renormalize explicitly
rather than rely on that.

## Validation

**V1 — free-particle limit.** With `λ ≡ 0` the perturbation vanishes, so `T` must be
identically zero and σ = 0 to round-off. Catches sign, masking and normalization errors
independently of any reference data.

**V2 — energy normalization of `F_{E,l}`.** Verify `F` satisfies the free radial equation on
the grid, and that its normalization matches the analytic energy-normalized convention.

**V3 — convergence study (a required deliverable, not a side note).** σ at a fixed anchor
must be stable against: electronic box size `r_max`, **ECS angle θ**, DVR order, and ECS tail
length; and against the nuclear grid parameters. **θ-independence is the sharpest available
check** — a converged ECS result must not move when the contour rotates. Report a
convergence table; choose the smallest grid where σ is stable to ~1% and use *that* for the
benchmark. This deliberately redoes the box/angle study eMoScat performed but never
documented.

**V4 — the six benchmark anchors vs Houfek. This is the gate, not the goal.** The same
`reference.ANCHOR_COORDS` the LCP solver is measured on. Because this is the same model and
method Houfek used, agreement should be far better than the LCP's factor-3 bound. V4 passing
is what earns the exact solver the right to be used as an oracle in V5–V7. **A tolerance
will be set from the converged result, not assumed in advance** — see Open Questions.

**V5 (primary deliverable) — LCP vs exact, head to head.** At all six anchors, report
`σ_LCP / σ_exact` with the exact result as the reference. Current LCP-vs-Houfek ratios are
0.44, 0.77, 0.83, 1.01 at the GATED anchors; those become secondary once a verified exact
solver exists. **The size and structure of this discrepancy is the sub-project's scientific
output.** A large, well-characterized failure of the LCP is a successful result.

**V6 — do the LCP's two documented structural failures show up as predicted?** The elastic
anchor (0.2 Ha, v'=0) and the near-threshold anchor (0.02 Ha, v'=1) are where the LCP is
*expected* to fail, for reasons already derived: no background scattering, and an
energy-independent width giving the wrong threshold law. Quantify the failure against the
exact result at both. **This is a measurement, not a pass/fail requirement** — making it a
gate would create pressure to fudge it.

**V7 — nuclear dynamics, not just integrated cross sections.** The interesting differences
between an approximate and an exact treatment need not show up in σ, which integrates a lot
away. Compare the **nuclear-coordinate density** of the driven solution: project the exact
`Ψ⁽⁺⁾(r,R)` onto `R` (over the unscaled region) and compare its shape against the LCP
solver's 1-D driven solution `ξ(R)` at the same energy. Report the comparison as a figure
and a few summary numbers (e.g. centroid and width in `R`). This is exploratory — there is
no pass/fail — but it is where a local-width approximation is most likely to visibly break,
and it is the natural bridge to the further model systems and higher-dimensional models
planned for later sub-projects.

**Harness wiring.** A new group **E** in `validation/n2/experiment.py`, guarded like B1/C5/D1
so a solver error becomes a labeled FAIL rather than a crash. Whether group E runs at full
size in Docker depends on the converged grid cost — see Open Questions.

## Out of scope

- The dissociative-attachment (DA) channel.
- Full σ(E) curves and the boomerang structure. The API takes an array of energies from the
  start, so a curve run is later a matter of compute time only.
- The time-dependent 2-D route (wavepacket + Tannor–Weeks). A possible sub-project #7; it
  would need a sparse propagator, since `qscat.evolution.make_cn_stepper` is dense-only.
- Coupled partial waves; rotation; anything beyond the single fixed `l = 2`.
- Rust. The cost is entirely inside SuperLU.
- Promoting anything new into `qscat`.
- **Other model systems** (F₂, NO, H₂⁺ — decks exist in `reference/eMoScat/input/`) and
  **models above two dimensions**, both of which are planned follow-on work. They are out of
  scope here, but they are the reason the #5 library layer is dimension-general and why
  anything N₂-specific in this sub-project belongs under `projects/`, never in `libs/qscat`.

## Open questions (to be resolved by the work, not guessed now)

1. **The V4 tolerance.** Set it from the converged result and the residual model differences,
   and document the reasoning — do not pick a number in advance and tune the grid until it
   passes. If agreement is worse than the LCP's, that is a finding to investigate and report,
   not to paper over.
2. **Docker/harness cost.** If the converged grid is cheap enough, group E runs at full size.
   If not, the harness runs a reduced grid with a documented looser tolerance while the
   converged result lives in the test suite or a one-off study. Decide from the measurement.
3. **Ordering.** `MMD_AT_PLUS_A` roughly halved fill versus COLAMD on a small random matrix
   (8.79 vs 17.60). Measure on the real Hamiltonian and pick accordingly — this is a
   speed/memory choice only and cannot affect correctness.

## Verification

- `uv run pytest projects/n2_2d_cross_section libs/qscat validation/n2 -q` → all pass.
- `uv run mypy libs/qscat` → 0 errors. `uv run ruff check .` → clean.
- `uv run python -m validation.n2.experiment` → group E reports real PASS/FAIL/NOTE rows;
  the pre-existing **19 PASS / 0 PENDING / 2 NOTE / 0 FAIL** must not regress (the 2 NOTEs
  may become PASS if V6 succeeds — that would be the headline result).
- `docker/build.sh test` → passes.
- Convergence table and the V4/V5/V6 comparison recorded in `docs/physics/n2-2d-cross-section.md`.
