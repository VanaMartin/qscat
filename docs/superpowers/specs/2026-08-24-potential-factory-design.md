# Potential factory — model surfaces fitted to real molecules

Date: 2026-08-24

## Purpose

Every `ResonanceModel` in `qscat.model` is a hand-tuned testbed: Houfek, Rescigno
& McCurdy chose the N₂-like, NO-like and F₂-like constants so that the fixed-`R`
resonance *resembles* the real molecule's (`E_res(R_0) ≈ 2.4 eV`, `Γ(R_0) ≈
0.46 eV` for N₂), while the neutral curve is a Morse with `D_0` ≈ 2× real N₂'s and a
vibrational ladder that is off (`docs/physics/n2-resonance.md`).

This spec defines a **factory**: given a real diatomic (or diatomic ion) and its
published resonance data, produce a `ResonanceModel` whose neutral curve `V_0(R)`
and resonance curve `V_res(R) = E_res(R) − iΓ(R)/2` — and, beyond the pole, the
energy-dependent width `Γ̃(ε, R)` — match that molecule as closely as the
two-dimensional model form permits, and **say how closely**. The program's purpose
is unchanged (`docs/physics/potential-factory-options.md`): the exact 2-D solution
stays the oracle and the approximations stay the thing under test. A closer model
makes the verdicts about those approximations transferable to the real molecule;
no parameter is ever tuned to an experimental observable.

The first target is **O₂** (Alt & Houfek, Phys. Rev. A **103**, 032829 (2021)),
because its nonlocal model is published *complete* — `V_0`, the anion curve,
`A(R)`, `B(R)`, `α`, `c` — in exactly the format the factory consumes.

The survey of methods and data sources that led to this design, including the
feasibility spike, is `docs/physics/potential-factory-options.md`; this spec does
not repeat it.

## Decisions taken in brainstorming (2026-08-24)

| Question | Decision |
|---|---|
| Target tier for v1 | **T3** — the published nonlocal-model functions. T2 (raw R-matrix eigenphase sums) is kept as a first-class *slot* in the target and as a recorded future direction (an own or vendored R-matrix capability), not built now. |
| First molecule | **O₂** (²Π_g, `l = 2`). |
| Neutral-curve fidelity | **Richer analytic form** (EMO), not Morse, not a spline. |
| The "3-D" the contract must extend to | **Both** meanings — electron angle `θ_e` (coupled partial waves) and a second nuclear coordinate — so the contract is generic in the coordinate tuple. |
| Tolerances | **Calibrated on N₂ and NO, not chosen a priori** — see "Tolerance budget" below. |

## Reference literature

- `reference/literature/houfek-2006-pra73-032721.md`, `houfek-2008-pra77-012710.md`
  — the model form and the N₂/NO/F₂ constants the factory must reproduce in the
  round-trip test.
- Alt & Houfek, Phys. Rev. A **103**, 032829 (2021) — **not yet a
  `reference/literature` note**; writing it (with page-anchored locators for §II–III,
  Eqs. 20–31, Tables I–II) under the `mastering-references` skill is the first task
  of the plan. Every equation number below refers to that paper unless stated.
- Čížek & Houfek, *Nonlocal theory of resonance electron–molecule scattering*
  (CRC 2011, ch. 4) — the general parametrisation `V_dk(R) = Σ f_i(ε) g_i(R)` and the
  threshold law `Γ_l(ε) ∝ ε^{l+1/2}` (§4.3.2, Eq. 4.113–4.116); background for the
  target format, cited but not required as a note.
- `reference/literature/houfek-2008-pra77-012710.md` Eq. (21), (68) — `V_dk⁺` and
  `Γ = 2π|V_dk⁺|²`, the repo's forward model for the T3 loss
  (`qscat.core.nrm.coupling`).

## The target format (what a molecule *is* to the factory)

A `Target` is a tiered, provenance-carrying description. Tiers are independent
optional slots; the factory fits whichever are present and reports on each.

```
Target
  molecule      name, mu, ell, charge, symmetry label of the resonance
  coordinates   tuple of named nuclear axes -- ("R",) for a diatomic; a second
                axis or an electronic angle appears here for 3-D targets
  neutral  (T0) constants {R_e, D_e, omega_e, omega_e x_e, ...} and/or a
                tabulated curve V_0(R) with its uncertainty
  resonance(T1) V_ion(R) = V_0 + E_res (table or callable), Gamma(R) local width,
                EA of the fragment (the asymptote V_ion(inf) - V_0(inf) = -EA)
  coupling (T3) the energy-dependent width:  Gamma~(eps, R) = 2 pi eps^alpha A(R)
                exp(-B(R) eps)  with A(R), B(R), alpha, and background c;
                an energy window [eps_min, eps_max]  (Alt & Houfek Eq. 24-27)
  eigenphase(T2) OPTIONAL tables delta(eps; R) -- reserved; a loader only, no loss
  provenance    a locator per slot (paper, table/figure, page) in the
                mastering-references form
```

Each curve is a callable on the coordinate tuple. For a 3-D target the same object
carries two nuclear axes (a resonance *surface*) or a partial-wave-resolved coupling
(the `Λ`-block K-matrix of `docs/physics/angular-coupled-channels.md`); the fitter
does not change.

For an **ion** (`charge = −1`) the `resonance`/`coupling` slots are replaced by a
`rydberg` slot — quantum defects `μ(R)` or the bound Rydberg curves — with
`qscat.core.bo.electronic_curves` as the forward model. Out of scope for v1; the
slot exists so the contract does not have to change.

## The ansatz (what the factory can emit)

`FlexibleDiatomicModel`, implementing `ResonanceModel`, a superset of
`DiatomicResonanceModel`:

```
v0(R)       = EMO:  D_e * [ (1 - exp(-beta(R) (R - R_e)))^2 - 1 ]
              beta(R) = sum_i beta_i y_p(R)^i,  y_p = (R^p - R_e^p)/(R^p + R_e^p)
              (Le Roy; N_beta = 0 recovers Morse exactly, zero at infinity as now)

V_int(r,R)  = -lambda(R) exp(-alpha(R) r^2)                        (well)
              + beta_b(R) exp(-alpha_b (r - r_b)^2)                 (barrier shell, optional)

lambda(R), alpha(R), beta_b(R):  the existing sigmoid form
              f_inf + f_0 / (1 + exp(f_1 (R - R_f)))  times a low-order
              polynomial in y_p(R) -- analytic, bounded, entire in R
```

Constraints the ansatz enforces by construction:

- **ECS analyticity.** Every term is entire in `r`; in `R` the only poles are those
  of `y_p` at `|R| = R_e` in the complex plane, which the nuclear ECS tail (pivot
  `R_0 > R_e`) never approaches. The factory computes and publishes the maximal safe
  nuclear/electronic ECS angles for the emitted parameters (the same attribute
  `IonicResonanceModel.max_nuclear_ecs_angle_deg` already carries).
- **Threshold law.** `ell` comes from the resonance's symmetry and is not fitted;
  `Γ̃ ∝ ε^{l+1/2}` is then exact, and the factory reports the exponent as a check
  rather than a residual. A dipole exponent `α = √(d + ¼)` (polar targets) is **out
  of scope** — the ansatz cannot produce it, and the report must say so when a
  target carries one.
- **Special case.** With `N_beta = 0`, `alpha(R) ≡ alpha_c`, `beta_b ≡ 0`, the
  ansatz *is* `DiatomicResonanceModel`; `qscat.model.N2/NO/F2` are exact points of
  the parameter space. That is what makes the round-trip test (Validation (a)) an
  exact oracle.

Why two free functions plus a shell: the spike showed the bare well `(λ, α)` is a
genuine two-parameter family — at fixed `E_res = 2.3 eV`, `l = 2`, `α` tunes `Γ`
from 1.67 eV (`α = 0.2`) to 0.16 eV (`α = 0.8`) — so `λ(R)` **and** `α(R)` are
needed to follow both `E_res(R)` and `Γ(R)`; Houfek's constant `α_c` can follow one.
The shell is what shapes the off-resonance phase and the `B(R)` fall-off (T3), and
gives an `l = 0` system a barrier at all.

## The forward models (all existing code)

| Tier | Model quantity | Computed by |
|---|---|---|
| T0 | `v0(R)` | the ansatz itself (closed form) |
| T1 | `E_res(R)`, `Γ(R)` | `qscat.core.lcp.resonance_pole_walk` → `qscat.ecs.find_resonance_pole` on an `ecs_angle_family` pair |
| T1 | `EA` asymptote | `qscat.core.dissociation.anion_electronic_states` at `R_max` |
| T3 | `Γ̃(ε, R) = 2π|V_dk⁺(ε,R)|²` | `qscat.core.nrm.coupling.v_dk_plus` + `gamma_from_coupling`, discrete state `AsymptoticDiscreteState` (choice B — R-independent, formally exact for this model, `docs/physics/nonlocal-resonance-model.md`) |
| T3 | `V_d(R)` (Eq. 29 of Alt & Houfek: `V_0 + E_res − Δ̃`) | `qscat.core.nrm.ingredients.nrm_ingredients` |
| T2 | `δ(ε; R)` of the model (for the future loss) | `qscat.core.nrm.scattering.scattering_state` → asymptotic phase; **not built in v1** |

The T3 comparison is therefore like with like: the paper's `Γ̃` *is* `2π|V_dε|²`
of a nonlocal model, and the repo computes the same object for the 2-D model.

The gradient, when a stage needs one, is the Hellmann–Feynman formula under the
c-product: `∂E_pole/∂V(r_i) = ψ_i² / (ψᵀψ)` (`qscat.linalg.c_product`), so one
eigenvector gives the exact sensitivity of the pole to every potential value, and
the chain rule gives it per ansatz parameter. No autodiff dependency.

## Staged fitting

Each stage seeds the next. A stage that cannot meet its tolerance **stops and
reports**; the factory never returns a model from a later stage that violated an
earlier one.

1. **T0 — neutral curve.** Closed-form least squares of the EMO to the tabulated
   curve, or to `(R_e, D_e, ω_e, ω_e x_e)` when only constants exist. `mu` from the
   isotopic masses. Output: `v0`, the vibrational ladder (via
   `qscat.core.vibrational`), and the ladder's error against the published one.
2. **T1 — pole curves.** Per nuclear node `R_j` (descending, as the pole walk
   does), solve the 2×2 problem `(λ, α) ↦ (E_res, Γ)` for the target
   `(V_ion(R_j) − V_0(R_j), Γ(R_j))` by Newton on the c-product gradient, seeded by
   continuation from the neighbouring node. This needs a **continuation tracker in
   parameter space** — the spike showed a bracket scan on `λ` loses the branch for
   `α ≳ 1` — with the same discipline as the `R`-walk (accept only angle-stable
   poles with residual `≪ Γ`; freeze-and-flag on breakdown). Then fit the analytic
   `λ(R)`, `α(R)` forms to the per-node solutions (global: differential evolution or
   Nelder–Mead as the paper used; local: least squares), and **re-verify** the pole
   curve of the smoothed model — the smoothing is where the residual actually
   lives. The asymptote `EA` is a hard constraint on `λ_inf`.
3. **T3 — energy-dependent width.** Global least squares of the model's
   `Γ̃(ε,R)` against `2π ε^α A(R) e^{−B(R)ε}` over the window, on a grid of
   `(ε, R)`; `beta_b(R)` (the shell) enters here, `λ, α` may re-polish within the T1
   tolerance. The loss is on `log Γ̃` so the threshold region and the peak count
   alike.
4. **Report.** `FitReport`: the parameter table, per-tier residuals (RMS and
   maximum, in the units the tolerance is quoted in), the attainability verdict per
   tier (`met` / `not met: <which feature>` / `not attempted`), the ECS angle bounds,
   the crossing `R_c` and the DA threshold sign (must equal the molecule's), and the
   provenance copied through from the `Target`.

The output is (a) a `FlexibleDiatomicModel` instance, (b) a serialised parameter
file (YAML/JSON) with the report, and (c) a `MoleculePreset` for `apps/qscat-run`
— the model **plus grids built by `qscat.tuning.propose_grid` and the
`discretisation-tuner` loop** (grids are per-potential; a copied deck is never
acceptable). `qscat-run` selects models by registry name (`molecule: O2`), so a
fitted molecule becomes a registry entry exactly as a hand-tuned one does.

## Tolerance budget — calibrated on N₂ and NO (the decision from review)

Tolerances on the *curves* are not chosen a priori. They are derived from what they
do to the *observables*, measured on the two molecules where the exact 2-D solver
already gives resonance levels and cross sections and where the hand-tuned model is
an exact point of the ansatz: N₂-like and NO-like.

**Procedure (a sensitivity study, itself a deliverable):**

1. Take `qscat.model.N2` (and `NO`). Compute its exact observables on the
   converged deck: `driven.ve_cross_section` (elastic + first two excitations over
   the published window), `dissociation.da_cross_section` where open (NO), the BO
   levels `lcp.resonance_levels`, and the exact poles `exact_resonance_states`.
2. Perturb one curve feature at a time by a controlled amount and recompute:
   `E_res(R)` shifted by `δE` (uniform, and a tilt); `Γ(R)` scaled by `(1 + δΓ)`
   (uniform, and only near the crossing); `R_c` moved by `δR`; `v0`'s `ω_e` by
   `δω`; `B(R)` (via the shell) by `δB`. The perturbations are realised as ansatz
   parameter changes, so every perturbed model is a legitimate `ResonanceModel` and
   the perturbed curves are *measured*, not assumed.
3. Score each perturbed observable against the unperturbed one with the metrics
   the repo already has: peak positions in units of a resonance width
   (`qscat.core.assignment.peak_positions` / `peak_alignment`), peak heights
   (relative), integrated cross section (relative), resonance-level positions and
   widths.
4. **The budget** is the largest curve deviation that keeps every observable metric
   within a stated observable-level tolerance. That tolerance is the one number
   chosen by hand, and it is chosen at the observable level because that is where
   "reproduces the behaviour correctly" is defined: *peak positions within
   0.25 of a width, peak heights and integrated σ within 10 %, level positions
   within 0.1 of their width, level widths within 10 %* — proposed, to be
   confirmed when the sensitivity numbers exist.

The result is a per-feature budget table — it will not be a single percentage,
because the observables are far more sensitive to `R_c`, to `Γ` near the crossing
and to the ladder than to `Γ` at short `R`. That table becomes the O₂ acceptance
criteria, and the same table applies to every later molecule until re-measured.

**Push-back recorded.** "Reproduce the resonance and cross sections correctly"
must mean *correctly relative to the exact 2-D solution of the fitted model*, not
relative to experiment. The sensitivity study is oracle-based and never touches
measured cross sections; if it did, tolerances would be tuned to the very
observables the program uses as verdicts, and the LCP/NRM comparisons on the fitted
model would no longer be independent. Two consequences the study must state
explicitly: (i) a tolerance calibrated on N₂/NO's `l = 2`/`l = 1` boomerang
structure is a *transfer* to O₂ (also `l = 2`, similar mass) and will need
re-measuring for a light, broad-resonance system such as H₂; (ii) NO's exact DA is
exponentially sensitive near threshold (`docs/physics/nonlocal-resonance-model.md`
§8), so NO's DA enters the budget through peak positions and the threshold, not
through relative magnitude, or it would forbid every fit.

## Validation

- **(a) Round-trip differential oracle.** Extract `qscat.model.N2`'s own T0/T1/T3
  curves with the forward models above, package them as a `Target`, run the
  factory, and recover Table I (`houfek-2006`, p. 032721-5) — an exact known
  answer, and the test that the fitter, the tracker and the report are right
  before any real data is touched. The gate is on the **`v0` constants**
  (`D_0`, `α_0`, `R_0`, well-determined by the curve) and on the **curves**
  `λ(R)`, `α(R)` and the pole curve over the fitted range — not on the sigmoid's
  individual constants, which are degenerate there (N₂'s `R_λ = −27.98`,
  `λ_1 = 1.06` make the sigmoid nearly a single exponential over `R ∈ [1.6, 3]`,
  so many constant sets give the same curve). Repeat for NO and F₂ (different
  `l`, different crossing geometry).
- **(b) O₂ acceptance** against the budget table above, on the published T0/T1/T3
  data. The report is the deliverable even when a tier is `not met` — a documented
  limit of the 2-D form is a result.
- **(c) Grid and 2-D convergence.** The emitted `MoleculePreset`'s grids come
  from the tuner, and the 2-D spot-check (`refine_to_2d_convergence`) must pass
  on the fitted model before its observables are quoted.
- **(d) `physics-reviewer`** before promotion: units, the `Γ`-support condition
  (`Γ(R) ≠ 0` only where `V_0 < V_ion`), the DA threshold sign, the ECS angle bounds.
- **(e) Spurious-pole guard.** The tracker's acceptance test (angle-stable *and*
  residual `≪ Γ` *and* `Re E` above a threshold floor) is unit-tested against the
  fake near-threshold match the spike produced at `(0.11 eV, 0.26 eV)` for `l ≤ 1`.

## Data acquisition for O₂ (explicit, because it is the risk)

Alt & Houfek publish Table II (the seven nonlocal constants) and Table I (EA(O),
EA(O₂), D₀); the potential curves `V_0(R)`, `V_ion(R)` and `Γ(R)` are **figures**
(Fig. 2, Fig. 4). O₂ therefore runs in two phases (decision 2026-08-24):

1. **Image match first.** Digitise Fig. 2 (and Fig. 4 for the T3 shape) with a
   stated precision that enters the T0/T1 residuals as an uncertainty floor, and
   fit to that. The purpose of this phase is to learn *what the 2-D form can do*
   on a real target, not to produce the final O₂ model.
2. **Benchmark on the real data.** The tabulated curves are requested from the
   authors (the second author is the user's own group); the acceptance verdict in
   "Validation (b)" is taken against those tables, and the phase-1 fit is
   re-run on them.

T0 constants from spectroscopy (NIST / Huber–Herzberg; `R_e`, `ω_e`, `ω_e x_e`,
`D_0 = 5.165 eV` from Table I) are the fallback for `v0`. All locators go through
`mastering-references`.

**N₂ and NO are inputs, not benchmarks against nature.** Their *calculated* curves
from the existing 2-D models are the factory's input in the round-trip test and the
sensitivity study; reproducing them well is the requirement. Nothing in v1 is
matched to an experiment — this is theoretical work. Getting closer to
experimental data is a stated later direction ("Kept open"), not part of this
spec.

## Placement and lifecycle

- Toy: `projects/potential_factory/` (`target.py`, `ansatz.py`, `tracker.py`,
  `fit.py`, `report.py`, `sensitivity.py`), following the qm-method-lifecycle.
- Promotion: `libs/qscat/qscat/factory/` — a **sibling of `qscat.model`** that may
  import `qscat.model` and `qscat.core`; `qscat.core` stays import-clean
  (`test_core_no_model_import.py` unchanged). `FlexibleDiatomicModel` itself joins
  `qscat.model`.
- Validation home: `validation/factory/` (round-trip gates, the sensitivity study,
  the O₂ report); `validation/diatomic/config.py` gains nothing — O₂'s grids are
  tuner-built.
- Rust: none foreseen; the cost is the eigensolves and sparse factorizations
  already optimised.

## Kept open (recorded, not built)

- **T2 — raw R-matrix eigenphase sums.** `Target.eigenphase` and a loader exist
  from v1; the *capability* to compute fixed-nuclei eigenphase sums for a real
  molecule (UKRmol+ as a dependency, or an own one-electron/model R-matrix) is a
  separate future direction noted in `docs/physics/potential-factory-options.md`.
  The model-side phase `δ(ε;R)` is cheap to add from `scattering_state` when that
  day comes.
- **SUSY/Bargmann oracle.** An exact per-`R` construction of a local potential with
  the target pole, to bound what *any* single-`l` local potential can carry of a T3
  target. Not in v1.
- **Flexible surface + gradient (route D).** Only if the ansatz fails a target that
  matters; the c-product gradient above is the enabling piece.
- **Ions.** The `rydberg` slot; `IonicResonanceModel` as the ansatz.
- **Closer to experiment.** At some later point the fitted models may be
  compared with measured cross sections (as Alt & Houfek do for O₂ VE). That is a
  *comparison*, made after the fit, never a fitting target — and it is not part
  of v1.
- **3-D.** Coordinate-tuple contract in place from v1; Legendre components
  `v_λ(r,R)` as extra ansatz terms, or a second nuclear axis, when
  `angular-coupled-channels` or a triatomic model is built.

## Out of scope for v1

Polar molecules (dipole threshold law), ions, any 3-D ansatz, the R-matrix
capability, the SUSY oracle, and any fit to an experimental observable.
