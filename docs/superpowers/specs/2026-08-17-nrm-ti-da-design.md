# Nonlocal resonance model — TI core and dissociative attachment

Date: 2026-08-17

## Purpose

Implement the **nonlocal resonance model (NRM)** — the Feshbach projection-operator
approximation to the 2-D electron-diatomic model — as a third method alongside the
exact 2-D solver (`qscat.core.driven`/`qscat.core.dissociation`) and the local
complex potential (`qscat.core.lcp`).

This is the research program's standard shape. The exact 2-D solver is the oracle;
the approximation is under test. The repo already measures the LCP against that
oracle and has documented where it fails: it misses the non-resonant background in
VE-elastic, and its energy-*independent* width gets DA threshold behavior wrong
(`docs/physics/diatomic-ve-cross-sections.md`). The NRM sits strictly between LCP
and exact — it keeps the energy dependence and the nonlocality the LCP throws away
— so the question this capability answers is: **how much of the LCP's documented
error is bought back by nonlocality alone, and where does even the NRM break down?**

PRA 77's own answer is that it depends on how the discrete state `φ_d(r;R)` is
chosen, and that for DA the "intuitive" physical choice can fail through a
Born-Oppenheimer breakdown while an R-independent choice stays near-exact. This spec
reproduces that comparison inside qscat, on the repo's own exact oracle rather than
on the paper's figures.

Scope of *this* spec is the NRM core plus **dissociative attachment**. Vibrational
excitation needs the background T-matrix and is spec 2; the time-dependent NRM is
spec 3.

## Reference literature and notation

Primary, already tracked:

- K. Houfek, T. N. Rescigno, C. W. McCurdy, *Probing the nonlocal approximation to
  resonant collisions of electrons with diatomic molecules*, Phys. Rev. A **77**,
  012710 (2008) — `reference/literature/houfek-2008-pra77-012710.md`. **Every
  equation number in this spec refers to this paper unless stated otherwise.**
- K. Houfek, T. N. Rescigno, C. W. McCurdy, Phys. Rev. A **73**, 032721 (2006) —
  `reference/literature/houfek-2006-pra73-032721.md`, the 2-D model itself.

To be ingested as part of this work, under the `mastering-references` skill:

- W. Domcke, *Theory of resonance and threshold effects in electron-molecule
  collisions: the projection-operator approach*, Phys. Rep. **208**, 97 (1991) —
  the canonical NRM formulation. Needed because PRA 77 explicitly disagrees with
  its Eq. (4.14): the coupling in the resonant T-matrix must be `V_d,-k̄` (or
  `V*_dk` for a real discrete state), not `V_dk`, and PRA 77 states Domcke's form
  "was, in our opinion, used incorrectly" (p. 012710-4, following Eq. 31). We
  implement PRA 77's form; the note records the disagreement so the choice is
  traceable at the line of code that makes it.
- P. L. Gertitschke, W. Domcke, *Time-dependent wave-packet description of
  dissociative electron attachment*, Phys. Rev. A **47**, 1031 (1993) — the
  time-dependent nonlocal treatment of DA. Ingested now because spec 3 needs it and
  the PDF is in hand; not cited by spec 1's code.

**Notation.** This spec uses PRA 77's symbols throughout, which collide with
qscat's in one place. PRA 77's `V_d(R)` (Eq. 20) is the *discrete-state potential*
`V_0(R) + ⟨φ_d|H_el|φ_d⟩`. qscat's existing `Vd` in `qscat.core.lcp` is the real
part of the LCP curve, `V_0(R) + Re E_res(R)`. The two only *almost* coincide, and
only for the physical discrete-state choice (p. 012710-8). The new code therefore
names PRA 77's quantity `V_d_discrete` in any context where both are in scope, and
the module docstring states the collision explicitly.

## Scope

**In scope**

- The NRM ingredients: discrete state `φ_d(r;R)`, projected electronic spectrum
  `E_n(R)`, discrete-continuum couplings `V_dn(R)` and `V_dk⁺(R)`.
- The nonlocal potential `F(E,R,R')` (Eq. 53), assembled per Eq. (60)–(61).
- The nuclear wave equation (Eq. 52) and `σ_DA` (Eq. 54).
- Two discrete-state choices: **A** the "intuitive"/physical state (Sec. VI A) and
  **B** the R-independent state `φ_b` (Sec. VI B).
- Molecules **F₂ and NO**, on the per-molecule nuclear decks already in the repo.
- A `qscat-run` method `nrm` producing σ_DA(E) artifacts.

**Out of scope**

- Vibrational excitation: the resonant T-matrix (Eq. 31), the background T-matrix
  (Eq. 37), and the final-state projections. Spec 2. *Note the seam is not where it
  first appears* — see "The P-space continuum is not deferrable" below.
- Discrete-state choice **C**, the "compact" `λ_spec(R)` state (Sec. VI C,
  eMoScat's `MakePhiD_spec`). The paper itself calls it model-specific tuning
  rather than a general recipe.
- N₂. Added in spec 2 with VE, where the published Houfek data anchors it.
- The time-dependent NRM. Spec 3.
- Any matrix-free / iterative variant of the `F` application (see "Cost").

## Reference implementation: untrusted

`reference/eMoScat` contains an NRM module — `source/module_NRM.cpp` (723 lines),
`include/module_NRM.h`, and per-molecule `input/{N2,NO,F2,coupled}/NRM.txt`.

**eMoScat never delivered the NRM as a working capability.** Nothing in that module
is treated as correct, complete, or validated by this spec. It is read the way one
reads someone's abandoned notebook: it may suggest *which* steps exist, never that a
step is right. This is a stronger position than "port carefully" — no line of the
new implementation may be justified by "eMoScat does it this way", and no test may
assert agreement with it. **PRA 77 is the sole specification.** Where the C++ and
the paper differ, the paper wins without further argument.

The internal evidence agrees with that judgement:

- `TimeIndependentSolution` is a half-finished stub: the `V_dn` multiplications into
  the Green operator are commented out (lines 669–670), every term inside `M(n)`
  except `E_t` is commented out (line 655), and the right-hand side carries a
  `FIXME` (line 703). It computes nothing meaningful.
- `build_potentials` and `build_coupling` *look* complete but carry `FIXME` markers
  and commented-out alternative formulations at exactly the delicate points (the
  `PHP` triple product, lines 409–421 and 571–575).
- `build_coupling` contradicts the paper on the ECS treatment of `φ_k⁺` — see "Open
  question" below. Under this section's rule that is not a puzzle to adjudicate but
  a settled matter: the paper's ECS treatment is the specification, and eMoScat's
  real-region variant is at most a hypothesis to test if the specified form fails
  its own check.
- `python/nrm_test.py` is exploratory scratch work on eigenvalue-curve reordering,
  not a test.

Consequently the validation plan below rests entirely on the paper and on qscat's
own validated code. Per `CLAUDE.md`, nothing in `reference/` is imported, built, or
edited.

## Method

For each nuclear DVR point `R_j`, with `H_el(r;R) = T_r + l(l+1)/2r² + V_int(r,R)`
(Eq. 17):

1. **Discrete state.** `φ_d(r;R_j)`, normalized, square-integrable, satisfying
   `φ_d(r;R) → φ_b(r)` as `R → ∞` (Eq. 67). Both choices below satisfy this by
   construction.
2. **Projectors.** `Q(R_j)_{ik} = √w_i φ_d(r_i;R_j) φ_d(r_k;R_j) √w_k` (Eq. 58),
   `P = 1 − Q` (Eq. 57). Note the **bilinear** outer product — no conjugation, since
   under ECS the natural inner product is the c-product
   (`qscat.linalg.c_product`), consistent with the paper's own statement that the
   complex-scaled `P H_el P` is symmetric and "we have to use for the wave functions
   the scalar product defined without complex conjugation" (p. 012710-6).
3. **Discrete-state potential.** `V_d(R_j) = V_0(R_j) + ⟨φ_d|H_el|φ_d⟩` (Eq. 20).
4. **Projected spectrum.** Solve `(P H_el P) φ_n = E_n(R_j) φ_n` (Eq. 56). The
   operator is complex **symmetric**, not Hermitian, so eigenvectors are normalized
   under the c-product. Eigenvalues lie in the fourth quadrant.
5. **Couplings.** `V_dn(R_j) = Σ_{ik} √w_i φ_d(r_i) H_el(r_i,r_k) φ_n(r_k) √w_k`
   (Eq. 59).

Then per energy `E`:

6. **Nonlocal potential.** `F(E,R_i,R_j) = Σ_n √W_i V_dn(R_i) M(n)⁻¹_{ij} V_dn(R_j)
   √W_j` (Eq. 60), with `M(n)_{ij} = [E − T_R − V_0(R) − E_n(R)]_{ij}` (Eq. 61).
   `T_R` is the nuclear kinetic matrix; `V_0` and `E_n` are diagonal in the nuclear
   DVR. No singularity arises because `E_n` is complex while `E` is real (p.
   012710-6).
7. **Nuclear equation.** `[E − T_R − V_d(R)] Ψ_d⁺(R) − ∫ F(E,R,R') Ψ_d⁺(R') dR' =
   V_dk_i⁺(R) χ_vi(R)` (Eq. 52).
8. **Cross section.** `σ_DA_vi(E) = (2π²/k_i²)(K_DA/μ) lim_{R→∞} |Ψ_d⁺(R)|²`
   (Eq. 54).

### Discrete-state choices

**A — "intuitive"/physical (Sec. VI A).** At each `R`, take the resonance energy
`E_res(R)` from the ECS pole, then solve the fixed-nuclei electron **scattering**
problem at the real energy `Re E_res(R)` and truncate the resulting wave function
smoothly with

```
f(r) = 1 − 1/(1 + e^{−(r − r_d)}),    r_d = 10 a_0      Eq. (69), p. 012710-8
```

multiplying by `e^{−iδ(R)}` to make it real and normalizing to unity. Where
`E_res(R)` is real and negative the electron is bound and `φ_d` is that bound state
directly. This is **not** the ECS pole eigenvector — `qscat.core.lcp.
resonance_eigenstate` supplies `E_res(R)` but the state itself comes from a
real-energy scattering solve, the same primitive the coupling layer needs.

**B — R-independent (Sec. VI B).** `φ_d(r;R) = φ_b(r)` for all `R`: the electronic
bound state of the isolated atom, i.e. the `R → ∞` limit of choice A. One electronic
eigenproblem for the whole calculation. Implemented directly as that limit of choice
A, which makes Eq. (67) hold by construction rather than by a separate potential.

### The P-space continuum is not deferrable

The right-hand side of Eq. (52) carries `V_dk_i⁺(R)`, the coupling at the **real**
incident electron energy — not one of the discretized `V_dn`. PRA 77 is explicit
(p. 012710-6): these "are at specific real electron energies and hence must be
evaluated directly using Eq. (21) where the background continuum function
`φ_k⁺(r;R)` is obtained by solving Eq. (18) in the electronic DVR basis under
exterior complex scaling."

So the P-space continuum solve belongs to spec 1 even though the *background
T-matrix* does not. Eq. (18) is a scattering problem for `P H_el P` at real energy
`k²/2` with incident wave `𝒥_k(r)`, which in driven (scattered-wave) form is

```
φ_k⁺ = 𝒥_k + (P H_el P − E)⁻¹ [−(P H_el P − E) 𝒥_k]
V_dk⁺(R) = ⟨φ_d| H_el |φ_k⁺⟩                      Eq. (21)
```

with `𝒥_k` the Riccati-Bessel function of order `l` (`qscat.special.radial`) and all
products bilinear. This is the standard driven-equation form the repo already uses
in `qscat.core.driven`, applied here in the electronic coordinate with the `P`
projection.

### The `φ_k⁺` boundary treatment is specified, and checked

PRA 77 states that `φ_k⁺(r;R)` is obtained by solving Eq. (18) in the electronic DVR
basis **under exterior complex scaling** (p. 012710-6). That is what we implement.

The reason this gets its own section is that it is easy to get wrong and its error
propagates everywhere downstream — `V_dk⁺` is both the right-hand side of Eq. (52)
and the outgoing projection. So it is checked independently rather than assumed,
using Eq. (68):

```
Γ(E,R) = 2π |V_dk⁺(R)|²,   E = k²/2       Eq. (68), p. 012710-8
```

Evaluated at `E = Re E_res(R)`, this must reproduce the `Γ(R)` that
`qscat.core.lcp.local_complex_potential` already computes from the ECS pole — two
independent routes to the same physical width, one of them already validated in the
repo. **This is a gate, not a diagnostic**: if the specified ECS form fails it, the
coupling layer is blocked until the disagreement is understood, and no downstream
layer is built on top of an unverified `V_dk⁺`. The outcome is recorded in the
physics note either way.

## Architecture

New package `libs/qscat/qscat/core/nrm/`. A package rather than a flat module:
`qscat/core/lcp.py` is already 40 KB and this capability is comparable, with six
genuinely separable layers.

| File | Responsibility | Depends on |
|---|---|---|
| `discrete_state.py` | `DiscreteState` protocol — `phi_d(R) -> NDArray` on the electronic grid, plus a `phi_b` asymptote. `PhysicalDiscreteState` (choice A), `AsymptoticDiscreteState` (choice B) | `scattering.py`, `qscat.core.lcp.resonance_pole_walk` |
| `scattering.py` | `scattering_state(h_el, E, ell, grid)` — the real-energy fixed-nuclei scattering solve, shared by choice A and by the coupling | `qscat.dvr`, `qscat.special.radial`, `qscat.linalg.SparseLU` |
| `ingredients.py` | `NrmIngredients` (`V_d_discrete`, `E_n`, `V_dn`) + its builder: the per-R loop, `Q`/`P`, the `PHP` eigenproblem, state tracking, tail freezing | `discrete_state.py`, `qscat.linalg.c_product` |
| `coupling.py` | `V_dk⁺(R;E)` at real energies (Eq. 21) | `scattering.py`, `discrete_state.py` |
| `nonlocal_potential.py` | `nonlocal_operator(ingredients, nuclear_grid, model, E) -> NDArray` (Eq. 60–61) | `ingredients.py` |
| `dissociation.py` | `nrm_da_cross_section(...)` — the Eq. (52) solve and Eq. (54) | all of the above |

The layering matters for testability: `ingredients` is a pure function of (grids,
model, `φ_d`) and is computed **once per molecule per choice**, then reused across
the whole energy sweep — the same analyze-once/solve-per-energy shape as
`SparseLU.refactor` in `qscat.core.driven`. Every layer above has a check that does
not require the layer below (see Validation).

`qscat.core` must not import `qscat.model` at runtime; the existing
`test_core_no_model_import.py` guard covers the new package automatically since it
scans `qscat.core`.

### Two numerical details that will silently produce wrong answers

- **The `PHP` null mode.** `P` annihilates `φ_d`, so `P H_el P` has a spurious zero
  eigenvalue whose eigenvector is `φ_d` itself. It must be excluded from the sum
  over `n` in Eq. (60), which the paper's Eq. (56) basis implies but does not spell
  out. We identify it by eigenvalue magnitude *and* by asserting its eigenvector
  overlaps `φ_d` above a threshold — never by index, which is not stable across `R`.
- **Adiabatic state tracking across `R`.** The `n`-th eigenvalue of `PHP` at
  adjacent `R` must be matched by nearest-energy continuation, not by the
  eigensolver's ordering, or `E_n(R)` and `V_dn(R)` acquire discontinuities that
  corrupt `F`. We walk `R` inward from the outermost real point (where the states
  are cleanly ordered and `φ_d → φ_b`) matching by minimum
  `|E_n(R_{i+1}) − E(R_i)|`, with a continuity assertion on `E_n(R)` as a guard.
  Avoided-crossing regions are where this is expected to be hardest; the guard is
  what turns a silent corruption into a visible failure.
- **ECS tail.** The ingredients are electronic quantities evaluated at a nuclear
  coordinate; the paper's `φ_d` is "real and localized in the inner region where the
  electronic coordinate is not complex scaled" (p. 012710-6, following Eq. 59), and
  Eq. (67) makes `V_dk(R) → 0` as `R → ∞`. So `E_n`, `V_dn` and `V_d` are computed
  on the real nuclear region and continued into the complex tail at their
  outermost-real values. The continuation rule is checked, not assumed: the tail
  must not carry appreciable `V_dn`, and the `σ_DA` result must be insensitive to
  where the real region ends.

## Cost

Measured against the actual decks (`validation/diatomic/config.py`): the electronic
grid is 132 points for both molecules; the nuclear deck is 974 points (819 real) for
F₂ and 597 (507 real) for NO.

- **Ingredients**, once per molecule per choice: 819 dense 132×132 complex
  eigenproblems ≈ 5×10¹⁰ flops — seconds.
- **`F(E)`**, per energy: 131 inversions of a 974×974 complex matrix ≈ 3×10¹¹ flops
  ≈ 10 s, 15 MB peak. A 41-energy sweep ≈ 7 minutes.
- **Nuclear solve**, per energy: one dense 974×974 solve — negligible beside `F`.

This stays inside the repo's CPU-runnable-on-a-laptop invariant, so `F` is formed
**explicitly**, exactly as Eq. (60) writes it. That is the whole point: an explicit
`F` can be plotted, its local limit inspected numerically, and its assembly checked
term by term against the paper. A matrix-free variant (`F` applied as
`Σ_n D_n M(n)⁻¹ D_n ψ` inside a preconditioned iterative solve) is 10–50× cheaper
and may be needed for spec 2 or larger grids; `nonlocal_operator`'s signature is
kept narrow enough that it can be swapped behind the same call site, but it is not
built here.

## Validation

Three checks, ordered by what they isolate. The first two test the new code against
*existing validated code*, so a failure localizes to the new layer.

There is deliberately no check against eMoScat's own NRM output: `build_potentials`
writes `Vdn.bin`/`Vd.bin`/`En.bin`, but no run output is committed anywhere in
`reference/eMoScat` (verified — the tree contains the four `NRM.txt` inputs and no
`.bin` files at all). The C++ is a porting reference for the algorithm, not a
numerical oracle.

1. **`Γ` consistency (coupling layer).** `2π|V_dk⁺(R)|²` at `E = Re E_res(R)`
   against `lcp.local_complex_potential`'s `Γ(R)`, over the real nuclear region.
   Tests `scattering.py` + `coupling.py` with no `F` and no nuclear solve involved.
   Also the decider for the ECS-vs-real-region open question above.
2. **Local-limit bridge (nuclear-solve layer).** Replacing `F` by its local limit
   `−(i/2)Γ(E,R)δ(R−R')` must make `nrm_da_cross_section` reproduce
   `lcp.lcp_da_cross_section` to solver tolerance on the same grid. This is a true
   differential oracle: it catches sign, weight (`√W`), normalization and
   boundary-condition errors that an order-of-magnitude agreement with the exact
   solver would hide.
3. **Exact-2D gate (the result).** `σ_DA^NRM(E)` against
   `dissociation.da_cross_section` for F₂ and NO across the energy range the
   existing LCP comparison uses, for both discrete-state choices, reported alongside
   the existing LCP ratios.

**What counts as success for check 3 is a measurement, not a threshold.** The repo's
documented LCP behavior is that σ_LCP/σ_exact sweeps 0.263 → 1.736 across
0.010–0.050 Ha on F₂ and fails outright on NO (ratio 1.8×10⁹ away from threshold).
The NRM's job is to do better, and PRA 77 predicts *how much* better differs by
choice: B near-exact, A degraded by Born-Oppenheimer breakdown, worst for the wider
resonances. The gate is therefore written the way the repo's other cross-model gates
are — a recorded ratio curve with an asserted band set from the first converged run,
plus an assertion that choice B beats the LCP by a stated margin at the anchors.
Any band is chosen after the first run and justified in the physics note; none is
invented here.

### Test placement

- `libs/qscat/tests/test_nrm_ingredients.py`, `test_nrm_coupling.py`,
  `test_nrm_dissociation.py` — checks 1–3 and unit-level behavior, fast paths only.
- The full-deck F₂/NO sweeps are `@pytest.mark.slow`, following
  `validation/tuning/test_emoscat_decks.py`'s precedent: a ~7-minute sweep is over
  any per-group harness budget.
- `validation/diatomic/` gains the NRM-vs-LCP-vs-exact comparison as the
  reported-artifact layer.

## Deliverables

- `libs/qscat/qscat/core/nrm/` (six modules above) with tests.
- `docs/physics/nonlocal-resonance-model.md` — the method, the two discrete-state
  choices, the ECS-vs-real-region finding, the measured σ_DA comparison, and the
  `V_d` naming collision.
- Reference notes `reference/literature/domcke-1991-physrep208-97.md` and
  `gertitschke-1993-pra47-1031.md`.
- `qscat-run` support for `methods: [ti, lcp, nrm]` on the DA observable, plus an
  example config producing the F₂ three-way σ_DA overlay figure.
- A committed figure `docs/physics/figures/f2-da-nrm-vs-lcp-vs-exact.png`.

## Risks

- **The physical discrete state (choice A) may be delicate to construct.** It needs
  a scattering solve at `Re E_res(R)` at every `R`, a phase rotation, and a cutoff
  at `r_d = 10 a₀` on a grid whose real region ends at `r = 16 a₀`. If it proves
  unstable, choice B alone still delivers a working NRM and the A-vs-B contrast
  moves to spec 2 — the spec degrades gracefully rather than blocking.
- **The ECS-vs-real-region question may not resolve cleanly** against `Γ(R)`. That
  outcome blocks the coupling layer and must be understood before proceeding, since
  everything downstream multiplies through `V_dk⁺`.
- **F₂'s deck is large enough that the sweep is minutes, not seconds.** Iteration on
  the full deck will be slow; unit tests run on a reduced grid and the full sweep is
  `@slow`.
