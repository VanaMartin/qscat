# NO and F₂ exact-2D VE cross sections (the model port)

**Location:** `validation/diatomic/` (`config.py` — the per-molecule eMoScat nuclear
decks; `test_diatomic.py` — the cross-section gate). The per-molecule *curve and figure*
drivers this note's committed figures came from were retired into `apps/qscat-run`, so
the curves are now produced from a config (e.g.
`apps/qscat-run/examples/f2-da-lcp-vs-exact.yaml`) rather than from solver code here.
**Computed entirely through the promoted library** `qscat.core` + `qscat.model` — the first
consumers of the generalization beyond N₂ (sub-project A). **Origin:** sub-projects B (NO) and
C (F₂) of the diatomic VE-scattering spec. **Units:** atomic units.

## What this is

The exact 2-D time-independent vibrational-excitation cross section
$\sigma_{0\to v'}(E)$ for **NO** and **F₂** — the same model and method as N₂,

$$
H_\mathrm{2D} = -\frac{1}{2}\frac{\mathrm{d}^2}{\mathrm{d}r^2}
  - \frac{1}{2\mu}\frac{\mathrm{d}^2}{\mathrm{d}R^2}
  + v_0(R) + \frac{l(l+1)}{2r^2} - \lambda(R)\,e^{-\alpha_c r^2}
$$

differing only in parameters (`qscat.model.{NO,F2}`). Adding each molecule
was **data + validation, no new solver code**: a `qscat.model` registry entry + a
`validation/diatomic/config.MoleculeConfig` grid/energy entry + a driver calling
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

**Which grid these figures use.** The curves above are computed on each molecule's own
**eMoScat production deck** (`grid: {preset: emoscat}`, what `apps/qscat-run` uses):
NO `(132, 597)` and F₂ `(132, 974)` electronic × nuclear points. They are *not* on the
N₂-style shared grid.

**Convergence.** An earlier r_max study on the **shared N₂-style grid** (electronic
r_max = 16, nuclear r_max = 22) found σ(E) unchanged to <1 % at electronic r_max = 16/24/32
for NO — evidence that the sharp low-E swings are genuine resonance structure rather than grid
noise. That measurement stands, but it was made on the shared grid and is *not* a convergence
statement about the per-molecule decks these figures now use. No equivalent r_max sweep has
been run on the production decks; the case for them rests instead on their being eMoScat's own
convergence-tested decks and on the discretisation tuner's reproduce-and-beat comparison
against them (`docs/physics/discretisation-tuning.md`).

The structure is, independently, not a sampling artefact: these curves are sampled at
0.001 Ha — 97 points for F₂, 117 for NO — against boomerang features this note measures at
~0.004 Ha wide, so each feature carries several samples.

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
outgoing nuclear wave $\sqrt{\mu/2\pi k_R}\;e^{i k_R R}$, with
$k_R = \sqrt{2\mu(E_\mathrm{tot} - \varepsilon_e)}$. The cross section is
$\sigma_\mathrm{DA} = \pi|S_\mathrm{DA}|^2/2E$ with $S_\mathrm{DA}$ from the same Tannor-Weeks transform / driven-equation projection
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

**The DA cross section — a TIME-INDEPENDENT driven-equation T-matrix (eMoScat
`time_independent_model.cpp`).** eMoScat computes DA (and H₂⁺ DR) exactly and time-*in*dependently,
via the SAME driven equation as VE ($\Psi_+ = \Psi_i + (E-H)^{-1} V_\mathrm{int} \Psi_i$) but projected onto the DA
exit channel with the **rearrangement interaction**

  $$V_\mathrm{DR}(r,R) = V_\mathrm{int}(r,R) + v_0(R) - V_\mathrm{int}(r, R\to\infty)$$

i.e. $H - H_\mathrm{final}$, NOT `V_int`,

and the DA T-matrix
$T_\mathrm{DA} = \langle \phi_e(r)\,F^\mathrm{nuc}_{K_R,0}(R) \,\vert\, V_\mathrm{DR} \,\vert\, \Psi_+ \rangle$, $\sigma_\mathrm{DA} = 4\pi^3|T_\mathrm{DA}|^2/2E$ —
where `φ_e` is the anion bound electronic state at the dissociation limit and `F^nuc` is the
energy-normalized regular nuclear Bessel ($l=0$, mass $\mu$), $K_R = \sqrt{2\mu(E_\mathrm{tot} - \varepsilon_e)}$. (An earlier
prototype of mine used `V_int` instead of `V_DR` and got a ~10⁶ unitarity violation — that was
the bug, NOT a structural obstacle to a TI DA.) With `V_DR`, σ_DA is O(1) bohr² and within the
unitarity cap `π/2E` for F₂; N₂/NO closed (correct). Implemented in
`qscat.core.dissociation` (`anion_electronic_states`, `v_dr_diag`, `da_cross_section`).

**The discretisation must be per-molecule.** DA's outgoing flux is in the NUCLEAR coordinate, and
the heavy nuclei make the exit wave $F^\mathrm{nuc} = \sqrt{2\mu/\pi K}\,\sin(K_R R)$ oscillate fast (F₂:
`K_R ≈ 58`, wavelength ~0.107 bohr). On the single N₂-style nuclear grid (1.0-bohr outer
elements) σ_DA did NOT converge — it swung `0.16 → 26 → 0.54 → 2.3 → 4.0` bohr² as the nuclear
*quadrature* was raised, because increasing points-per-element (p-refinement) cannot resolve an
oscillatory integrand; element density (h-refinement) must. The fix is eMoScat's **already-tested
per-molecule grids** (`reference/eMoScat/input/{NO,F₂}/grids.txt`), whose nuclear grids are far
finer over the dissociation region (NO: 37×0.2 bohr over [1.6, 9.0]; F₂: 40×0.2 bohr over
[2.7, 10.7] plus 0.024 bohr around the 2.5–2.7 resonance). On the F₂ eMoScat grid (nuclear
`n = 974`, `R0 = 10.7`) σ_DA(E=0.03) = **1.66 bohr²**, **stable to < 0.002 %** under a
quadrature-order refinement (nuc_quad 14→16) of that already-fine mesh — the coarse-grid `25.99`
was a pure quadrature artifact. (A full element-doubling *h*-refinement at this ~10⁵-unknown deck
size exceeds a laptop's memory; confirming it needs the Docker+MUMPS backend — a follow-on.) `qscat.core.grids.segmented_grid`
builds any such deck; `validation.diatomic.config.MoleculeConfig.da_grid()` carries the NO/F₂
decks. This distinction is why VE converged on the coarse nuclear grid but DA did not: **VE's
outgoing flux is electronic** (needs fine *electronic* resolution, already validated at
r_max = 16), **DA's is nuclear** — they stress different coordinates.

The exact-2D TI σ_DA(E) oracle curves on those grids (F₂ exothermic; NO opening ~0.17 Ha):

![F₂ exact-2D TI dissociative attachment](figures/f2-2d-ti-da-cross-section.png)
![NO exact-2D TI dissociative attachment](figures/no-2d-ti-da-cross-section.png)

**Automatic discretisation.** The eMoScat decks above are hand-tuned; `qscat.tuning` (the
`discretisation-tuner` skill) now computes a grid from the potential curve alone — an
equidistribution mesh bounding the per-element de Broglie phase, an h/p quadrature sweep, and a
double-ECS-safe tail — calibrated and gated against exactly these N₂/NO/F₂ decks, including F₂'s
K≈58-78 DA wave. See docs/physics/discretisation-tuning.md.

**H₂⁺ DR** is the same T-matrix looped over the neutral's MANY bound electronic (Rydberg) states
+ a Coulomb incident (`coulomb::sF_en`). See the DA design spec.

## Dissociative recombination (DR) — H₂⁺, many channels

The same machinery generalizes to **dissociative recombination** for a molecular *ion*
(`e⁻ + AB⁺ → A + B`): the exit channels are the **multiple** bound electronic states of the
neutral AB at large R (a Rydberg series → in principle infinitely many, cut off at a
measurable number). This is the H₂⁺ case — deferred with the Coulomb tail (the ionic
`sphHankel1En` Coulomb branch exists in eMoScat), but structurally it is "DA with N electronic
channels instead of 1."

## Local-complex-potential (LCP) DA — the approximation under test

The exact-2D DA above is the ORACLE. The **local-complex-potential (LCP)** model is the
*approximation*: it reduces the full electron–nuclear problem to a **1-D nuclear** problem on a
complex potential `V_d(R) − iΓ(R)/2`, where the fixed-R electronic resonance pole gives
$V_d(R) = \operatorname{Re} E_\mathrm{pole}(R)$, $\Gamma(R) = \max(0, -2\operatorname{Im} E_\mathrm{pole}(R))$. Implemented model-independently in
`qscat.core.lcp`: `local_complex_potential` finds `E_pole(R)` by two-angle ECS matching
(`qscat.ecs.find_resonance_pole`) of `−½∂²_r + model.surface(r,R)`, **seeded from the bound anion
state at R_inf** (`anion_electronic_states`) and continued inward (validated against the N₂
`vres_on_grid` oracle to ~1e-5). `lcp_da_cross_section` is the **time-independent resolvent** form
(the T→∞ limit of eMoScat's `ModelLCP/SMatrix.cpp` doorway propagation, cheaper and
confound-free): from the doorway $d = \sqrt{\Gamma/2\pi}\;\chi_{v_0}$, solve
$\psi_\mathrm{sc} = (E_\mathrm{tot}\mathbb{1} - H_\mathrm{res})^{-1} d$, and the
DA amplitude is the outgoing flux at the boundary $X$,
$S_\mathrm{DA} = \sqrt{K/2\pi\mu}\;\psi_\mathrm{sc}(X)$,
$\sigma_\mathrm{DA} = 4\pi^3|S_\mathrm{DA}|^2/2E$.

**Two subtleties are decisive** (each collapsed σ_DA by many orders when wrong): (1) the nuclear
grid must be the **fine per-molecule eMoScat deck** — the K≈58 outgoing dissociation wave is
unresolved on the coarse shared grid (σ_DA drops ~36 orders); (2) the boundary observable is the
wavefunction **value** $\psi_\mathrm{sc}(X) = \psi_\mathrm{sc}^\mathrm{coeff}[b]/\sqrt{w_b}$, NOT the raw DVR coefficient (a √w
boundary-weight factor, ~27× in σ). These are the same lessons as the exact DA (per-molecule grid)
plus a DVR-normalization one.

**Result (F₂).** The LCP's error is **systematic and energy-dependent — there is no plateau of
agreement.** Recomputed 2026-08-17 at 41 energies (0.010–0.050 Ha, step 0.001):

![F₂ σ_DA: LCP vs exact-2D oracle](figures/f2-2d-da-lcp-vs-exact.png)

| E (Ha) | σ_DA LCP | σ_DA exact-2D | LCP/exact |
|---|---|---|---|
| 0.010 (threshold) | 1.410 | **5.366** | 0.263 |
| 0.020 | 1.56 | 3.36 | 0.47 |
| 0.030 | 1.471 | 1.656 | **0.888** |
| 0.040 | 1.02 | 0.72 | 1.43 |
| 0.050 | 0.490 | 0.282 | 1.736 |

The exact σ_DA falls monotonically by a factor of 19 across this range while the LCP stays
nearly flat, so the ratio sweeps **0.263 → 1.736**, crossing unity near E ≈ 0.032. Earlier
versions of this note described the LCP here as "a good approximation away from threshold",
and `CLAUDE.md` quoted "~11% of exact" — both were reading the single crossing point as if it
were a characteristic agreement. It is not: the LCP under-predicts below ~0.03 and
over-predicts above it, and only *passes through* good agreement on the way.

**That sweep is a real property of F₂'s deck, and it was checked rather than assumed.**
The LCP's `V_d`/`Γ` come from an ECS resonance-pole walk, which is not guaranteed to be
independent of the electronic box, so the walk was rerun at
`r_max` = 16/32/48/64/80/96 on both molecules. On F₂ **five of the six agree to ~1.8 %**
— the shipped deck's curve is determined, and only `r_max = 96` breaks, to a ratio range
of [0.027, 3.51]. On NO the spread is **3.98 × 10⁴**, non-monotone across the whole
range, with the larger boxes producing deep oscillatory minima spanning 1e-6 to 1e-1
where the shipped deck is smooth and flat.

The distinction matters, because the two failures are not the same kind: **on F₂ the
LCP's `V_d`/`Γ` are determined and wrong** by the documented energy-dependent factor;
**on NO they are undetermined** — a second defect, in NO's pole walk, independent of the
volume-T-matrix one and still unfixed.

**Those extra walks are not plotted, deliberately.** They are not alternative estimates
of the same quantity: on NO the walk does not converge, so they are failed computations,
and drawing them beside a real curve would invite a reader to average them or to read
the spread as an uncertainty band. What the ladder establishes is a yes/no — whether the
method determines `V_d`/`Γ` on that molecule at all — which is a sentence, not a curve.
Each figure's legend carries its own verdict, and
`validation/diatomic/da_figure.py --ladder` recomputes the spread on demand.

(The older 13-point curve was not wrong — its endpoints match this dense run exactly — it was
simply too sparse to show that the ratio never settles.)

The **physically-sensible departures** the comparison exists to expose (the exact solver is
the oracle, the LCP is under test): (a) it **under-predicts the near-threshold spike** — the
exact σ_DA rises sharply as E→threshold while the LCP stays smooth (0.263 at E=0.010);
(b) for the sibling VE channel, the LCP **elastic omits the non-resonant background** that the
exact elastic T-matrix contains (`driven.py` documents this) — a known *qualitative* LCP
limitation. The LCP VE cross section itself is not computed in *this* sub-project, so no
quantitative LCP-VE band is recorded in this note. It is computed elsewhere:
`validation/diatomic/ve_nrm.py` runs it (as `qscat.core.lcp.
lcp_ve_cross_section`, the graduated 1-D solver) alongside the exact and nonlocal routes, and
`docs/physics/nonlocal-resonance-model.md` §8.4 publishes the measured
LCP-over-exact bands for N₂ and F₂. N₂'s DA channel is closed (threshold +0.5 Ha), so
LCP and exact both give ≈0 there — a consistency sanity, no figure.

**NO — where the LCP fails outright.** 151 energies (0.150–0.300 Ha, step 0.001),
recomputed 2026-08-24 with the flux extraction. The exact σ_DA peaks at
**1.349 × 10⁻⁹ bohr² at E = 0.1720** and decays **smoothly, by ten orders of
magnitude**, to 8.47 × 10⁻²⁰ at E = 0.300. The LCP does not decay: it stays near
10⁻⁴ across the whole range. The ratio therefore runs from **9.70 × 10⁴** near the
peak to **4.84 × 10¹⁴** at the top of the range, **never crossing unity** — unlike
F₂, where it sweeps through 1 near E ≈ 0.032.

Three earlier claims about this curve are withdrawn, all of them artefacts of the
volume-form extraction (`nonlocal-resonance-model.md` §7.2.1):

- the "sharp spike at threshold, peak 0.0925 bohr²" — there is no spike; the true
  peak is 1.349e-9 and the curve is smooth;
- "thirteen orders of magnitude" and the ratio "1.8 × 10⁹" — both were reading the
  cancellation residue, not the cross section;
- **the oscillations**, which this note previously read as structure. They were
  structure in the residue. The converged curve has none.

The far tail is real, not a floor: flux extraction at `r_max` 16 against 48 agrees
to 5–6 significant figures even at 8.4568e-20 (E = 0.300), where the volume route's
two boxes differ by 5460×.

The LCP curve on the figure below carries its own, separate caveat — NO's ECS pole
walk does not converge in `r_max`, so that curve is one member of a family the method
does not pin down. Its legend says so.

![NO σ_DA: LCP vs exact-2D oracle](figures/no-2d-da-lcp-vs-exact.png)

This is a far stronger statement than the "shown for completeness, not a quantitative agreement
claim" hedge this note previously carried, and it is the more useful result: the local-complex-
potential approximation does not reproduce the exponential suppression of dissociative
attachment away from threshold at all. On a 12-point curve at a 0.0136 Ha step that behaviour
was not visible; on the dense curve it is the dominant feature.

**The LCP's verdict is observable-dependent — do not generalise either way.** The results above
are about **σ_DA cross sections** for F₂ and NO, and they are unflattering. Measured on a
different observable the same approximation does very well: for N₂'s **resonance levels**, the
BO + local approximation agrees with the exact (non-BO) 2-D poles to sub-meV in both position
and width (`docs/physics/exact-2d-resonances.md`). Neither result licenses the other. An
approximation that reproduces where a resonance sits need not reproduce how much flux leaves
through a particular exit channel, and here it does not.

## Time-dependent route — validated on N₂, not attempted for NO or F₂

The N₂ time-dependent route (order-3 Padé + Tannor-Weeks, `qscat.core.time_dependent`) matches
the exact TI oracle to a median of 0.2 % across 161 energies spanning 0.060–0.220 Ha, with 93 %
of in-window points inside 5 % — see `docs/physics/n2-2d-td-cross-section.md`.

NO and F₂ have **sharper, lower-lying** resonances than N₂, so a TD route for them has not been
attempted and remains the natural follow-on, gated against the exact-2D TI oracle delivered here.

**What this note used to say, and why it was withdrawn.** It previously argued that NO's and F₂'s
"~0.004–0.01 Ha boomerang features are at or below the finite-time-propagation resolution `2π/T`",
so TD-vs-TI agreement at the N₂ level should not be expected without a long-propagation study.
That reasoning was inherited from a claimed finite-T resolution limit in the N₂ note which has
since been **disproved**: those measurements came from an order-1 Crank-Nicolson propagator, and
at order-3 Padé the same energies track the oracle (E = 0.06/0.08/0.12 moved from ratios of
0.229/0.575/0.348 to 0.997/0.958/1.015). There is no established `2π/T` barrier to inherit, so
this note no longer offers one as a reason to expect poor agreement — or as a reason not to try.

A genuine caution does remain, and it is a different one: **NO's cross sections really do carry
sharp structure.** Its v'=1 excitation is measured to swing 34.76 → 6.13 bohr² across a single
0.002 Ha step (`validation/diatomic/test_diatomic.py`). Resolving *that* is a sampling question
for whoever runs the TD study, not a predicted failure of the method.
