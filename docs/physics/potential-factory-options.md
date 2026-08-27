# Potential factory: building model surfaces that match real molecules — options survey

**Status:** the survey behind the factory — methods, data sources, one
feasibility spike, and the decisions taken on it (last section). The first
implementation stage, the fitter round-tripped on N₂/NO/F₂, is
[`potential-factory.md`](potential-factory.md); the sensitivity budget and the O₂ fit are not built.
**Relates to:** `qscat.model` (`DiatomicResonanceModel`, `IonicResonanceModel`),
`qscat.core.lcp` (the fixed-`R` pole walk that would be the forward model),
`qscat.core.nrm` (the nonlocal model whose *inputs* are the natural target
format), `docs/physics/angular-coupled-channels.md` (the 3-D direction this
must stay compatible with).
**Units:** atomic units unless a value is explicitly given in eV.

## The question

Every model in `qscat.model` is a **given** testbed — a Morse neutral curve, a
sigmoid interaction strength $\lambda(R)$ and a Gaussian electron–molecule well
$V_{int}(r,R) = -\lambda(R)\,e^{-\alpha_c r^2}$ — whose constants were *hand-chosen* so that
the fixed-`R` resonance resembles a real molecule's (Houfek, Rescigno & McCurdy
2006 for N₂/NO, 2008 for F₂; see `reference/literature/`). The resemblance is
qualitative: the N₂-like model's `D_0` is ≈2× real N₂'s, its vibrational
spacing is off, and only `E_res(R_0) ≈ 2.4 eV`, $\Gamma(R_0) \approx 0.46\text{ eV}$ were matched
(`docs/physics/n2-resonance.md`).

The direction is a **factory**: given a real diatomic (or diatomic ion) and
published data about it, produce a `ResonanceModel` whose neutral curve
`V_0(R)` and resonance curve $V_{res}(R) = E_{res}(R) - i\Gamma(R)/2$ match that
molecule as closely as the two-dimensional model form permits — first in the
existing 2-D `(r, R)` form, later in a 3-D form. The purpose is unchanged from
the rest of the program: the exact 2-D solution stays the oracle and the
approximations (LCP, NRM, BO levels) stay the thing under test; a closer model
only makes the *verdicts* about those approximations transferable to the real
molecule. No parameter is ever tuned to improve agreement with experiment on
an observable.

## What "matching a real molecule" can mean — the target hierarchy

The 2-D model is a one-electron, single-partial-wave, *local* potential. A real
molecule's fixed-`R` electron scattering has many partial waves mixed by the
non-spherical field, exchange, a polarisation tail $-\alpha_d/2r^4$, and (for polar
molecules) a dipole tail. "As close as possible" therefore means: match the
**single-channel projection** the nonlocal resonance theory itself works
with. In increasing richness:

| Tier | Target | Source of data | What in the model it constrains |
|---|---|---|---|
| T0 | Neutral curve `V_0(R)`: `R_e`, `D_e`, $\omega_e$, $\omega_e x_e$, or a full RKR / *ab initio* curve | Spectroscopy (NIST, Huber–Herzberg); Le Roy-style EMO/MLR fits; MRCI curves | `v0(R)` — needs a form richer than Morse (Morse fits 3 constants) |
| T1 | Anion / resonance curve: `V_ion(R) = V_0 + E_res(R)` and local width $\Gamma(R)$; the asymptote `V_ion(∞) = V_0(∞) − EA(atom)`; the crossing `R_c` where $\Gamma \to 0$ | *Ab initio* bound-state curves (MRCI/CCSD(T)) for the bound region; R-matrix / complex-Kohn / CAP resonance curves for the metastable region | `V_int(r,R)` at each `R`: at least two free parameters per `R` (see the spike) |
| T2 | The fixed-`R` **eigenphase sum** $\delta(\varepsilon; R)$ of the resonant symmetry over an energy window, i.e. the energy-dependent width $\tilde\Gamma(\varepsilon,R)$, its threshold exponent $\varepsilon^{l+1/2}$ (Wigner) and the background phase | R-matrix (UKRmol) eigenphase sums — the exact input Alt & Houfek (2021) used for O₂ | `V_int(r,R)`'s *shape* in `r`, not just its depth/range; the model's `l` |
| T3 | The full nonlocal-model functions `V_d(R)`, $V_{d\varepsilon}(R)$ (or `A(R)`, `B(R)`, $\alpha$) | Published NRM fits (N₂, H₂, HCl, HBr, HF, F₂, O₂, …) | Same as T2, but already reduced to smooth functions of `R` — the most convenient published format |
| T4 | Observables (VE/DA cross sections) | Experiment | **Not a target.** Matching observables is circular with the program's purpose |

T3 is the sweet spot: the published NRM literature has already done the hard
reduction from many-electron scattering to $(V_0, V_d, V_{d\varepsilon})$, and the 2-D
model's fixed-`R` electronic problem is *exactly* the single-channel picture
those functions describe. T1 is what the current models were tuned to by eye;
T2 is what T3 was fitted from.

### The published recipe for the target data (verified source)

Alt & Houfek, Phys. Rev. A **103**, 032829 (2021), §III (read for this
survey; not yet a `reference/literature` note) builds the O₂ nonlocal model as:

1. `V_0(R)` and the bound anion curve from CASSCF/MRCI (aug-cc-pVQZ), the anion
   curve **shifted** so the asymptote reproduces the experimental EA(O) =
   1.461 eV (Table I).
2. Fixed-nuclei R-matrix (UKRmol) eigenphase sums $\delta_{\text{sum}}(\varepsilon; R)$ at
   `R = 1.80 … 2.25 a₀` (Fig. 3), fitted with RESON to a Breit–Wigner form
   for `E_res(R)`, $\Gamma(R)$ (Eq. 20–21).
3. The nonlocal model fitted to the same eigenphase sums with
   $\delta = \delta_{bg} + \delta_{res}$, $\delta_{bg} = c\,\varepsilon^\alpha$, $\tilde\Gamma(\varepsilon,R) = 2\pi\,\varepsilon^\alpha A(R)\, e^{-B(R)\,\varepsilon}$,
   `A(R) = (a₀ + a₁R) e^{a₂R}`, `B(R) = b₀ + b₁R`, $\alpha = l + \tfrac12$ with `l = 2`
   (Eq. 22–27), by Nelder–Mead on the mean-squared eigenphase error;
   $V_d(R) = V_0 + E_{res} - \tilde\Delta(E_{res}, R)$ (Eq. 29). Table II lists the seven
   fitted constants.

That is a complete, published target-data format for a factory: $(V_0(R), V_{ion}(R), A(R), B(R), \alpha, c)$. Čížek & Houfek's review chapter (*Low-Energy
Electron Scattering from Molecules, Biomolecules and Surfaces*, ch. 4, §4.3.2,
also read) states the same parametrisation generally — $V_{dk}(R) = \sum_i f_i(\varepsilon) g_i(R)$, threshold law $\Gamma_l(\varepsilon) \propto \varepsilon^{l+1/2}$ (Eq. 4.114), dipole
exponent $\alpha = \sqrt{d + \tfrac14}$ for polar targets (Eq. 4.115–4.116) — and notes that
direct *ab initio* `V_d`, `V_{dk}` over the full $(\varepsilon, R)$ range exist only for
H₂ and HeH⁺.

## The inverse step: from target curves to `V_int(r; R)`

At fixed `R` the electronic problem is one-dimensional radial scattering at a
single `l` in `V_int(r;R) + l(l+1)/2r²`. Its complete scattering data is the
phase shift $\delta_l(k)$ for all `k > 0` plus the bound states, and the **fixed-`l`
inverse problem is well posed** (Gel'fand–Levitan / Marchenko). That is the
right frame; the fixed-*energy* inverse problem (Newton–Sabatier) is both the
wrong data type (many `l`, one `E`) and mathematically unsound (Ramm 2001, see
references). Six routes, in decreasing directness:

### A. Parametric ansatz + least squares (recommended first step)

Keep an analytic, ECS-safe form — Gaussians and exponentials, entire in `r` —
and fit its parameters per `R` (or as smooth functions of `R`) so the pole of
the fixed-`R` problem lands on the target $E_{res}(R) - i\Gamma(R)/2$, optionally
also fitting $\delta(\varepsilon;R)$ over a window. The forward model already exists
(`qscat.core.lcp.resonance_pole_walk` → `qscat.ecs.find_resonance_pole`);
the fit needs a *continuation* tracker in parameter space rather than a
bracket scan (the spike below shows why).

**Feasibility spike (2026-08-24, throwaway script, not repo code).** For the
pure Gaussian well $-\lambda\, e^{-\alpha r^2}$ on the `electronic_grid(r_max=16, order=7,
n_complex=6)` pair at 35°/44°, root-finding $\lambda$ at fixed $\alpha$ so that
`E_res` hits a target, and reading off $\Gamma$:

| `l` | `E_res` (eV) | $\alpha$ → $\Gamma$ (eV) |
|---|---|---|
| 2 | 2.3 | 0.20 → 1.67, 0.30 → 0.57, 0.40 → 0.40, 0.60 → 0.23, 0.80 → 0.16 |
| 2 | 4.0 | 0.20 → 3.82, 0.30 → 2.07, 0.40 → 1.41, 0.60 → 0.83, 0.80 → 0.57 |
| 2 | 1.0 | 0.20 → 0.24, 0.30 → 0.08, 0.80 → 0.02 |
| 1 | 0.3 / 1.0 / 2.3 / 4.0 | $\alpha = 2.0$: 0.06 / 0.38 / 1.39 / 3.38 |

So at fixed `E_res` the range $\alpha$ tunes $\Gamma$ over an order of magnitude —
$(\lambda, \alpha)$ is a genuine two-parameter family covering real N₂'s
`(2.3 eV, 0.4 eV)` at $(\lambda, \alpha) \approx (5.0, 0.4)$ (Houfek's values) and F₂/NO-like
`p`-wave pairs. Two consequences for the design:

- Houfek's one-free-function form ($\lambda(R)$ with constant $\alpha_c$) can follow
  `E_res(R)` **or** $\Gamma(R)$, not both. A factory needs at least two free
  functions of `R` — $\lambda(R)$ and $\alpha(R)$ — and probably a third term (a
  repulsive Gaussian shell or a polarisation-like tail) to shape the
  off-resonance phase (T2) and to give `s`-wave (`l = 0`) systems a barrier at
  all: a bare Gaussian well has **no** `l = 0` resonance, and the probe's
  apparent `l = 0`/`l = 1` "poles" at `(0.11 eV, 0.26 eV)` were a spurious
  angle-stable match near threshold, the known weakness of
  `find_resonance_pole` (`docs/physics/h2plus-resonance-states.md`).
- Empty cells in the spike were bracket failures of the crude $\lambda$ scan (for
  $\alpha \ge 1.2$ the `d`-wave well must be much deeper than the scanned $\lambda \le 14$),
  not physical limits. The production tool must track the pole by continuation
  in $(\lambda, \alpha)$ — the same discipline as the `R`-walk.

**Pros:** cheapest; reuses the pole walk, the `ResonanceModel` protocol, the
tuner; the result is a small dataclass a reader can print. **Cons:** the
attainable set is whatever the ansatz spans — the factory must *report* the
residual per tier rather than promise a match.

### B. Exact per-`R` construction: rational S-matrix (Bargmann) potentials via SUSY/Darboux

A potential whose Jost function is a prescribed rational function of `k` can be
written in closed form (Bargmann 1949; Marchenko with a degenerate kernel); the
modern constructive route is a chain of supersymmetric (Darboux)
transformations, each adding one `S`-matrix pole — bound, virtual, or a
resonance pair $\pm k_r - i\kappa$ — with the threshold behaviour of the chosen `l`
built in and singular `r^{−2}` cores handled explicitly (Sparenberg & Baye,
Phys. Rev. C **55**, 2175 (1997); Baye & Sparenberg, J. Phys. A **37**, 10223
(2004)). This gives, at each `R`, an **exact** local potential with the target
pole and a controlled background — a mathematically clean answer to T1 and a
constructive one to T2 (a rational fit of $\delta(\varepsilon;R)$ → poles → potential).

**Pros:** exact; analytic (rational-exponential) hence ECS-continuable;
a strong *oracle* for checking how much of T2 an ansatz in A actually captures.
**Cons:** the potential's shape is dictated, not chosen — long exponential
tails and singular cores appear, `R`-smoothness must be imposed by smoothing
the *pole trajectories*, and the form does not drop into `DiatomicResonanceModel`
without a new model class. Best used as a reference/seed for A, not as the
shipped surface.

### C. Numerical fixed-`l` inversion (Marchenko / Gel'fand–Levitan) from δ(ε;R)

Take the R-matrix eigenphase sum as data, extend it to all `k` with a
Jost-function / R-matrix-pole parametrisation (Fabrikant's single-pole
R-matrix is the "R-matrix with a radius constraint" reading of the user's
proposal — Čížek & Houfek §4.1.1 note it is equivalent to the projection-operator
model), and integrate the Marchenko equation numerically. **Pros:** uses the
richest data (T2) directly. **Cons:** the extension to all energies dominates
the answer; the reconstructed potential is a table, not analytic, so it must be
re-fitted to an ECS-safe form anyway — which collapses C into A with a better
loss function. Not recommended as a separate route.

### D. Flexible surface + gradient optimisation (the "PINN" reading)

Parametrise `V_int(r,R)` richly — a Gaussian/exponential basis in `r` with
spline or small-network coefficients in `R`, or a small network of analytic
activations — and minimise a loss over the tiers with gradients. The essential
observation is that **no autodiff framework is needed**: for the
complex-symmetric ECS Hamiltonian the pole's sensitivity is the
Hellmann–Feynman formula with the bilinear c-product,
$\partial E_{pole}/\partial V(r_i) = \psi_i^2 / (\psi^T\psi)$ (`qscat.linalg.c_product`), so one
eigenvector gives the exact gradient with respect to *every* potential value at
once, and the gradient with respect to any parameter follows by the chain rule.
Phase shifts at real energy have the analogous formula through the scattering
solution. **Pros:** generalises to more targets (T2's window, several
resonances, the bound anion region) and to 3-D with no change of principle.
**Cons:** the flexibility is the risk — a network is not entire in `r` in any
controlled way, so ECS continuation can blow up; the design must confine
flexible parts to `r < R_0^{ECS}` and keep the tail analytic. A literal PINN
(residual of the Schrödinger equation as loss) is the wrong tool here: the
forward solver is already exact and cheap; only the *outer* optimisation is
wanted. The 2023–2026 "PINN inverse scattering" papers (neutron–α `P₃/₂`
potentials from phase shifts; RBF-network fixed-energy inversion) are this
route with an unnecessary network — see references.

### E. Bayesian / Monte-Carlo parameter inference

MCMC over the ansatz parameters with the tier residuals as a likelihood gives
posterior widths — which parameters the data pins and which are free. Useful
as a **UQ layer on top of A or D**, not as the fitting engine: each likelihood
evaluation is a pole walk (seconds), so a chain of 10⁴–10⁵ samples is hours,
acceptable once per molecule. Simulated annealing / differential evolution is
the same cost and the right global optimiser for A's start.

### F. Rejected: fixed-energy inversion (Newton–Sabatier)

Requires phase shifts for *all* `l` at one energy, and the scheme is shown
inconsistent in general (Ramm, math-ph/0105021). The model fixes one `l`; the
data varies $\varepsilon$. Wrong on both counts.

## Constraints the repo imposes on any surface

- **ECS analyticity.** `r` and `R` are complex on the tails; every term must
  be entire (or at least analytic in the sector) and bounded under rotation.
  `IonicResonanceModel` already documents an angle bound forced by an `R⁴`
  term (`max_nuclear_ecs_angle_deg = 22.5`); a factory must compute and
  publish the same bound for whatever it emits, and prefer forms with none.
- **Threshold law.** The model's `l` sets $\Gamma \propto \varepsilon^{l+1/2}$; it must be chosen
  from the resonance's symmetry (N₂ ²Π_g → `l = 2`; NO ³Σ⁻, F₂ ²Σ_u⁺, H₂
  ²Σ_u⁺ → `l = 1`; HCl ²Σ⁺ → `l = 0`, a virtual state, needs a barrier term or
  a dipole-like tail to exist at all). A polar target's $\varepsilon^{\sqrt{d+\tfrac14}}$ exponent
  cannot be produced by a centrifugal barrier; matching it needs an explicit
  long-range term and is out of scope for the first version.
- **The $\Gamma$-support condition.** $\Gamma(R) \neq 0$ only where `V_0(R) < V_ion(R)`
  (Váňa & Houfek 2017); the crossing `R_c` is a *derived* quantity of the fit,
  and its position relative to the neutral turning points decides whether DA
  is open — a fitted model must reproduce the molecule's DA threshold sign.
- **Grids are per-potential.** A new surface needs its own discretisation; the
  factory's output must feed `qscat.tuning.propose_grid` and the
  `discretisation-tuner` loop, not a copied deck.
- **`qscat.core` never imports `qscat.model`.** The factory is a *producer* of
  `ResonanceModel` instances and belongs beside `qscat.model` (or under
  `projects/` while a toy), never inside `qscat.core`.

## Candidate molecules and where their target data is published

Values are order-of-magnitude reminders from memory of the literature, **not
verified locators** — each must go through the `mastering-references` skill
before it is used as a target.

| System | Resonance / symmetry | `l` | Tier available | Published source of target data |
|---|---|---|---|---|
| N₂ / N₂⁻ | ²Π_g shape, ≈2.3 eV, Γ ≈ 0.4 eV at `R_e` | 2 | T3 | Berman, Estrada, Cederbaum & Domcke, Phys. Rev. A **28**, 1363 (1983) — the first full NRM; R-matrix `E_res(R)`, $\Gamma(R)$ in later work (Laporta et al. 2014, arXiv:1402.3814) |
| O₂ / O₂⁻ | ²Π_g, bound at `R_e`, resonant at short `R` | 2 | **T3, complete** | Alt & Houfek, PRA **103**, 032829 (2021), Table II + Figs. 2–4 |
| CO / CO⁻ | ²Π shape ≈ 1.5–2 eV, broad | 1 | T1 | Laporta, Cassidy, Tennyson & Celiberto (2012), arXiv:1206.2268 — R-matrix `E_res(R)`, $\Gamma(R)$ |
| H₂ / H₂⁻ | ²Σ_u⁺, ≈3 eV, very broad; DA to H⁻ | 1 | T3 | Čížek, Horáček & Domcke, J. Phys. B **31**, 2571 (1998); Horáček et al., PRA **70**, 052712 (2004), PRA **73**, 022701 (2006) |
| HCl / HCl⁻ (also HBr, HF) | ²Σ⁺ virtual state + dipole; DA to Cl⁻ | 0 (dipole) | T3 | Čížek, Horáček & Domcke, PRA **60**, 2873 (1999); HF: J. Phys. B **36**, 2837 (2003) |
| F₂ / F₂⁻ | ²Σ_u⁺, bound at `R_e`, resonant at short `R`; exothermic DA | 1 | T3 | Brems, Beyer, Nestmann, Peyerimhoff & Domcke, J. Chem. Phys. **117**, 10635 (2002) |
| NO / NO⁻ | ³Σ⁻ (also ¹Δ, ¹Σ⁺), near-threshold | 1 | T1–T3 | Trevisan, Houfek, Zhang, Orel, McCurdy & Rescigno, PRA **71**, 052714 (2005) |
| H₂⁺ / H₂ (ion, DR) | Rydberg series, quantum defects $\mu(R)$ | 0/2 | T0 + defects | exact ion curve; Hvizdoš et al. 2018 (already in `reference/literature`) |
| HeH⁺ (ion, DR) | Rydberg series | 0 | T3 | Movre & Meyer 1997 (cited by Čížek & Houfek §4.3.2 as one of two systems with *ab initio* `V_{dk}` over the full range) |

For an **ion**, the target changes: `v0` is the ion core, and the
$\sigma$-capture well must reproduce the neutral's Rydberg quantum defects $\mu(R)$
(equivalently the bound Rydberg curves), not a resonance width. The
`IonicResonanceModel` form and `qscat.core.bo.electronic_curves` already give
the forward model; the same fitter applies with a different loss.

## Kept open: computing the T2 data in this repository

The T2 tier assumes fixed-nuclei eigenphase sums from an *ab initio* R-matrix
code. A recorded future direction (decision of 2026-08-24) is to be able to
produce them here — either by vendoring/depending on an R-matrix suite
(UKRmol+ is the one Alt & Houfek used) or by an own implementation, at least
for one-electron model targets. Until then the factory's `Target.eigenphase`
slot is a loader for externally computed tables. The model-side counterpart —
the 2-D model's own fixed-`R` phase $\delta(\varepsilon;R)$ — is cheap from
`qscat.core.nrm.scattering.scattering_state` and is what a T2 loss would
compare against.

## The 3-D question

"3-D" has two different meanings, and the factory's target format differs:

1. **Electron angle $\theta_e$** (`docs/physics/angular-coupled-channels.md`):
   $V(r, \theta_e, R) = \sum_\lambda v_\lambda(r,R) P_\lambda(\cos\theta_e)$. The single-`l` model becomes the
   approximation under test. The target then is the fixed-`R` **`K`-matrix in
   the $\Lambda$ block** (several `l`), or at least the eigenphase sum with its
   partial-wave composition — R-matrix codes provide exactly that. Route A/D
   generalise by fitting $v_\lambda(r,R)$ for $\lambda = 0, 2, \ldots$.
2. **Two nuclear coordinates** (a triatomic, or a diatomic with rotation): the
   target becomes a resonance *surface* `E_res(R₁,R₂)`, $\Gamma(R_1,R_2)$; CO₂'s
   ²Π_u is the classic case. Data are sparser and the exact solver's cost
   climbs steeply (dense eigensolves are already out; sparse shift-invert is
   the path).

Either way the factory's contract should be "a set of target curves/surfaces +
a loss per tier + a `ResonanceModel` out", so that the dimension is a property
of the ansatz and the target, not of the fitter.

## Recommendation (for the design discussion)

Build **A** with **B as its oracle** and **E as an optional UQ layer**; defer
**D** until A's ansatz demonstrably cannot reach a T2 target on a molecule
that matters; drop **C** and **F**. Start from the tier-T3 data that is
already reduced in the literature (O₂ is complete and published by the
supervisor's group; N₂ and H₂ are the classic NRMs), and make the factory's
*first* deliverable a feasibility map — for each candidate molecule, the
per-tier residual the 2-D form can reach — rather than a fitted model. That
map is itself a result of the program: it says which real-molecule features
the two-dimensional model can carry at all.

## Decisions taken on this survey (2026-08-24)

- **Target tier: T3**, the published nonlocal-model functions (`V_0`,
  `V_ion = V_0 + E_res`, $\tilde\Gamma(\varepsilon,R) = 2\pi\,\varepsilon^\alpha A(R)\, e^{-B(R)\varepsilon}$). They are the
  richest data already reduced to smooth functions of `R`, and the 2-D model's
  fixed-`R` electronic problem is exactly the single-channel picture they
  describe. T2 (raw R-matrix eigenphase sums) stays a reserved slot in the
  target, with the capability to compute them here a recorded later direction
  (see "Kept open" above).
- **First molecule: O₂**, because Alt & Houfek (2021) publish its nonlocal
  model complete (Table II), it shares N₂'s `l = 2`, and it is bound at `R_e`
  and resonant at short `R`, so it exercises the crossing logic. Phase 1 is an
  image match of the published figures; the authors' tabulated curves, once
  obtained, are the benchmark.
- **Neutral curve: EMO** (Le Roy), not Morse (cannot carry a real ladder) and
  not a spline (not analytic on the ECS tail).
- **Contract generic in the coordinate tuple**, so both later "3-D" readings
  (electron angle; a second nuclear coordinate) change the ansatz and the
  target, never the fitter.
- **Tolerances are measured, not chosen**: a per-feature budget from a
  sensitivity study on N₂/NO against the exact 2-D solver, with a single
  hand-chosen observable-level tolerance. The exact solution of the fitted
  model is the oracle; nothing is compared with experiment. Until that study
  exists, the fitter's `Tolerances` are placeholders and say so.

## References found in this survey (2026-08-24)

Read for this survey (PDFs, not yet `reference/literature` notes):

- V. Alt, K. Houfek, *Resonant collisions of electrons with O₂ via the
  lowest-lying ²Π_g state of O₂⁻*, Phys. Rev. A **103**, 032829 (2021).
  <https://doi.org/10.1103/PhysRevA.103.032829> — §II–III, Eqs. 20–31, Tables
  I–II, Figs. 2–4.
- M. Čížek, K. Houfek, *Nonlocal theory of resonance electron–molecule
  scattering*, ch. 4 in *Low-Energy Electron Scattering from Molecules,
  Biomolecules and Surfaces* (CRC, 2011), §4.1.1, §4.2, §4.3.1–4.3.2.

Seen only as search results (titles/abstracts; to be verified before citing):

- J.-M. Sparenberg, D. Baye, *Inverse scattering with singular potentials: a
  supersymmetric approach*, Phys. Rev. C **55**, 2175 (1997); D. Baye,
  J.-M. Sparenberg, *Inverse scattering with supersymmetric quantum mechanics*,
  J. Phys. A **37**, 10223 (2004). <https://iopscience.iop.org/article/10.1088/0305-4470/37/43/014>
- *Single- and coupled-channel radial inverse scattering with supersymmetric
  transformations*, arXiv:1401.0439; *Unified Wronskian formulation of inverse
  scattering with supersymmetric quantum mechanics*, arXiv:2508.19022.
- *The Jost function and Siegert pseudostates from R-matrix calculations at
  complex wavenumbers*, Eur. Phys. J. A (2024).
  <https://link.springer.com/article/10.1140/epja/s10050-024-01316-4>
- A. G. Ramm, *Analysis of the Newton–Sabatier scheme for inverting fixed-energy
  phase shifts*, arXiv:math-ph/0105021 — the inconsistency result behind F.
- *Constructing inverse potentials from scattering phase shifts using
  physics-informed neural networks: application to neutron–alpha scattering*,
  arXiv:2605.02264; *Fixed-energy inverse scattering with radial basis function
  neural networks*, PTEP 2023, 113A01. — route D with a network.
- *Neural network potentials facilitating accurate complex scaling for molecular
  resonances*, Phys. Chem. Chem. Phys. (2024). <https://pubs.rsc.org/en/content/articlelanding/2024/cp/d4cp02452d>
- R. J. Le Roy, *dPotFit* (JQSRT 2017) and *betaFIT* — EMO / MLR / DELR forms for
  T0. <https://uwaterloo.ca/atmospheric-chemistry-experiment/leroy-programs>;
  Araújo & Ballester, *A comparative review of 50 analytical representations of
  potential energy interaction for diatomic systems*, Int. J. Quantum Chem.
  (2021). <https://onlinelibrary.wiley.com/doi/full/10.1002/qua.26808>
- Laporta et al., *Electron-impact resonant vibration excitation cross sections
  and rate coefficients for carbon monoxide*, arXiv:1206.2268; Laporta et al.,
  N₂, arXiv:1402.3814.
- Čížek, Horáček, Domcke, HCl: Phys. Rev. A **60**, 2873 (1999); H₂: Phys. Rev.
  A **70**, 052712 (2004), **73**, 022701 (2006); HF: J. Phys. B **36**, 2837
  (2003). Brems et al., F₂: J. Chem. Phys. **117**, 10635 (2002).
