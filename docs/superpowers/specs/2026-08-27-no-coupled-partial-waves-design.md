# NO coupled partial waves — a two-centre anisotropic model and its resonance poles

Date: 2026-08-27

## Purpose

Every resonance model in this repository fixes a single electronic partial wave.
The interaction $-\lambda(R)e^{-\alpha_c r^2}$ is isotropic in the electron's
angular coordinate, so one $l$ is not an approximation inside the model — it is
the whole model. That makes one question unanswerable: **how much does the
decomposition into a fixed-$l$ resonant model change the answer?**

This spec builds the smallest model in which that question has an answer. The
electron–molecule well is moved off the origin onto the two nuclei, which makes
the interaction anisotropic and couples partial waves; the multi-wave solution
is then the oracle and the fixed-$l$ solution is the approximation under test —
the same relationship the repository already has between the exact 2-D solver
and the LCP.

**NO is the target because it is the only molecule in the registry where the
coupling can be strong.** In a linear molecule the Legendre component
$v_\lambda$ couples $l \leftrightarrow l\pm\lambda$, and a homonuclear molecule
admits only even $\lambda$. N₂, O₂ and F₂ therefore reach their first coupling
partner across $\Delta l = 2$, over a centrifugal barrier that suppresses it —
which is why single-wave models have worked for them for decades. NO is
heteronuclear: $\Delta l = 1$ is allowed, and the resonance can leak into a
neighbouring wave without paying a barrier. NO also has a validated exact-2-D VE
*and* DA baseline, so the follow-on cross-section work has two observables to
land on.

Nothing here is fitted to, or compared with, experiment. The anisotropy is
geometric — the wells sit on the nuclei — not calibrated to a measurement.

The captured direction this spec discharges is
`docs/physics/angular-coupled-channels.md`; that note's framing, its warning
against truncating on "which wave hosts the resonance", and its cost estimate
are all carried forward here (the cost estimate is corrected — see
"The block assembler").

## Decisions taken in brainstorming (2026-08-27)

| Question | Decision |
|---|---|
| Molecule | **NO** — the only registry molecule admitting $\Delta l = 1$ coupling. |
| Anisotropy form | **Two-centre Gaussian well**, wrapping `qscat.model.NO` — not a new parameter set. |
| $\Lambda$ block | **$\Lambda = 1$** ($\pi$), $l = 1,2,3,\dots$ — the physical block, and it has no barrier-free channel to host spurious states. $\Lambda$ stays a parameter of the model; $\Lambda = 0$ is not run. |
| How the anisotropy is fixed | **A continuation sweep is the result.** $s: 0 \to 1$ (one centre $\to$ wells on the nuclei) and then $\kappa: 0 \to 0.5$. No value has to be declared "the NO one" — both endpoints are defined by geometry. **Amended after measurement:** the walk stops early, where $\Gamma$ exceeds $\varepsilon$ — see "Where the continuation actually ends". |
| Deliverable | **Resonance poles**, not cross sections. Fixed-$R$ electronic poles as a screen, then the full 2-D coupled poles. |
| Build shape | **Screen-gated, toy-first.** Machinery in `projects/`, campaigns in `validation/`, nothing into `qscat` in this spec. |
| $\lambda$ truncation | **None needed** — the Gaunt coefficient makes the $\lambda$ sum finite and exact once $N_l$ is fixed. |

## The model

### Form

$$
V_{\rm int}(\vec r, R) \;=\; -\tfrac{1+\kappa}{2}\,\lambda(R)\,
  e^{-\alpha_c\,|\vec r - d\,\hat z|^2}
\;-\; \tfrac{1-\kappa}{2}\,\lambda(R)\,
  e^{-\alpha_c\,|\vec r + d\,\hat z|^2},
\qquad d = s\,R/2 .
$$

Write $\lambda_A = \tfrac{1+\kappa}{2}\lambda(R)$ and
$\lambda_B = \tfrac{1-\kappa}{2}\lambda(R)$ for the two well depths; both
wells sit at the same distance $d$ from the centre, on opposite sides.

The coupled model **wraps** a `ResonanceModel` rather than restating its
constants, so the embedding below is structural and cannot drift:
$\lambda(R)$, $\alpha_c$, $\mu$ and the neutral curve $v_0(R)$ come from
`qscat.model.NO` unchanged, and the channel ladder starts at $l = $ `NO.ell`
$= 1$, which coincides with $\Lambda = 1$'s own minimum.

Two parameters, both with a geometric meaning:

- **$s \in [0, 1]$** — how far the wells sit from the molecular centre. $s = 0$
  puts both at the origin, where the sum collapses to
  $-\lambda(R)e^{-\alpha_c r^2}$ identically, *for any $\kappa$*. $s = 1$ puts
  them on the nuclei.
- **$\kappa \in [0, 1)$** — the amplitude asymmetry. $\kappa = 0$ is a symmetric
  two-centre well: only even $\lambda$ survive, and within $\Lambda = 1$ the
  resonance at $l = 1$ can reach only $l = 3$. This is the *homonuclear
  control* — the N₂/O₂/F₂ situation — reached without changing models.
  $\kappa > 0$ switches the odd $\lambda$ terms on and opens
  $l = 1 \leftrightarrow l = 2$.

The molecular axis is a symmetry axis of this potential, so $\Lambda$ is exactly
conserved and each $\Lambda$ block is an independent problem. Restricting to
$\Lambda = 1$ costs nothing.

### The coupled radial potential

$$
V_{ll'}(r,R) \;=\; \sum_\lambda c^\lambda_{ll'}\, v_\lambda(r,R),
\qquad
c^\lambda_{ll'} \;=\; \langle Y_{l\Lambda} | P_\lambda | Y_{l'\Lambda} \rangle
$$

with

$$
c^\lambda_{ll'} = (-1)^\Lambda \sqrt{(2l+1)(2l'+1)}
\begin{pmatrix} l & \lambda & l' \\ 0 & 0 & 0\end{pmatrix}
\begin{pmatrix} l & \lambda & l' \\ -\Lambda & 0 & \Lambda\end{pmatrix}.
$$

The 3-$j$ symbols vanish unless $|l - l'| \le \lambda \le l + l'$, so once $N_l$
is fixed the $\lambda$ sum is **finite and exact**. The note's open question —
whether the $\lambda$ truncation converges independently of the $l$ truncation
— does not arise in this model: there is one truncation parameter, $N_l$.

**That identity explains the finiteness; it is not the production route.**
$V_{ll'}$ is computed as a single angular quadrature against orthonormal
associated Legendre factors,

$$
V_{ll'}(r,R) = \int_{-1}^{1} \Theta_{l\Lambda}(x)\, V_{\rm int}(r,x,R)\,
  \Theta_{l'\Lambda}(x)\,{\rm d}x , \qquad x = \cos\theta,
$$

with $\Theta_{lm}$ normalized so that $\int \Theta_{lm}\Theta_{l'm}\,{\rm d}x =
\delta_{ll'}$. This is the same number by the expansion above, and it is
preferred for three reasons: an isotropic potential returns the Kronecker delta
to round-off, so the $s=0$ embedding gate is an identity rather than a
tolerance; it needs no Wigner 3-$j$ implementation, and `sympy` is not a
dependency of this repository; and it survives a change of well shape. The
Legendre components $v_\lambda$ remain as a **diagnostic** — they are what
gates 3 and 5 below are stated on.

### How $v_\lambda$ is computed

**Primary: Gauss–Legendre quadrature in $\cos\theta$**, 64 nodes,

$$
v_\lambda(r,R) = \frac{2\lambda+1}{2}\int_{-1}^{1}
  V_{\rm int}(r,\cos\theta,R)\, P_\lambda(\cos\theta)\, {\rm d}\cos\theta .
$$

Three reasons this is the primary route rather than the closed form:

1. **Complex-safe for free.** The quadrature nodes are real $\cos\theta$ values;
   $r$ and $R$ carry the ECS phase and pass straight through. Gaussians are
   entire, so the analytic continuation onto the ECS tail is unproblematic.
2. **The $s = 0$ collapse is exact.** At $s = 0$ the integrand is
   $\theta$-independent, and Gauss–Legendre integrates $P_\lambda$ against a
   constant to exactly zero for $\lambda > 0$. The embedding gate below is
   therefore a round-off identity, not a tolerance.
3. It survives a change of well shape without new analysis.

**Numerical requirement.** The Gaussian must be evaluated at the shifted
argument, $\exp(-\alpha_c|\vec r \mp d\hat z|^2)$, and **never** as the product
$e^{-\alpha_c(r^2+d^2)}\,e^{\,2\alpha_c r d\cos\theta}$: with NO's
$\alpha_c = 1.0$, $r$ up to 16 bohr and $d$ up to ~3 bohr, the second factor
reaches $e^{96}$ against a first factor of $e^{-265}$. The product is
representable; the factors are not.

**Oracle: the closed form**, used in tests only, at real argument:

$$
v_\lambda(r,R) = -(2\lambda+1)\,\big[\lambda_A + (-1)^\lambda \lambda_B\big]\,
  e^{-\alpha_c(r-d)^2}\, \tilde i_\lambda(2\alpha_c r d),
$$

which makes the $\kappa = 0$ symmetry manifest: there $\lambda_A = \lambda_B$,
the bracket is $\tfrac{1}{2}\lambda(R)\,[1 + (-1)^\lambda]$, and every odd
$\lambda$ vanishes identically.

where $\tilde i_\lambda(z) = e^{-z} i_\lambda(z) = \sqrt{\pi/2z}\;$`ive(`$\lambda
+ 1/2, z$`)` is the exponentially scaled modified spherical Bessel function
(`scipy.special.ive`). The $e^{-z}$ scaling is what makes the identity
$e^{-\alpha_c(r^2+d^2)}i_\lambda(2\alpha_c r d) = e^{-\alpha_c(r-d)^2}\tilde
i_\lambda(2\alpha_c r d)$ evaluable.

### The three-way comparison

Moving the wells apart changes the **monopole** $v_0$ as well as creating the
higher $\lambda$ components. That monopole shift is not coupling and must not be
counted as coupling. Every point of the continuation is therefore computed three
ways:

| Label | What it is |
|---|---|
| **full** | all $V_{ll'}$, $l = 1\dots N_l$ — the multi-wave oracle |
| **fixed-$l$** | the $l = 1$ diagonal block alone, at the *same* $(s, \kappa)$ — the approximation under test |
| **anchor** | $s = 0$ — must reproduce `qscat.model.NO` exactly |

The reported result is **full minus fixed-$l$ at matched $(s,\kappa)$**, which is
precisely "how the decomposition into a fixed-$l$ resonant model affects the
output". The anchor is the regression gate that replaces Houfek's external
certification, which no longer applies once the anisotropy is on.

## The block assembler

State vectors are **channel-outermost**:
$\Psi = [\psi_{l=1}(r,R),\, \psi_{l=2}(r,R),\, \dots]$. Then

$$
H_{ll'} = \delta_{ll'}\Big[T_r \oplus T_R
  + \mathrm{diag}\big(v_0(R) + \tfrac{l(l+1)}{2r^2}\big)\Big]
  + \mathrm{diag}\, V_{ll'}(r,R),
$$

so every diagonal block is exactly what `qscat.dvr.hamiltonian_nd` already
builds, every off-diagonal block is a `scipy.sparse.diags`, and the assembly is
one `scipy.sparse.bmat`. One block filler serves both campaigns: the fixed-$R$
electronic case is the same code with `qscat.dvr.kinetic_sparse` in place of
`kinetic_nd` and no nuclear axis.

Two consequences worth stating because they are better than
`docs/physics/angular-coupled-channels.md` predicted:

- **Nonzeros do not grow as $N_l^2$.** Because the off-diagonal blocks are
  diagonal, $\mathrm{nnz} = N_l\cdot\mathrm{nnz}(H_{2D}) + (N_l^2 - N_l)\cdot N$.
  With FEM-DVR kinetic bands at ~20–30 nonzeros per row, $N_l = 4$ costs ~4× the
  dimension and only ~1.2× that again in nonzeros. LU fill-in remains the binding
  constraint and must be measured, not assumed.
- **At $s = 0$ the matrix is block-diagonal to round-off**, and its $l = 1$
  block is `NO.hamiltonian(tgrid)`.

## Phase 1 — the electronic screen (laptop)

### Setup

- **Electronic grid**: `electronic_grid(r_max=16.0, order=8, n_complex=8)` —
  deliberately NOT `validation/diatomic/config.py`'s NO deck. See "Where the
  continuation actually ends": the published 6-element tail cannot represent the
  broadened resonance, and using it would measure the deck rather than the
  model.
- **Pole location**: `qscat.ecs.find_resonance_pole` on two electronic ECS
  angles, **44° and 52°** — again not the published pair (35°/44°), for the same
  measured reason — but requesting $k > 1$ eigenpairs near
  the shift, because *a second angle-stable pole entering the window is the
  single-discrete-state assumption failing*, and that is the finding being
  hunted. `find_resonance_pole` assumes one pole; the screen must therefore
  inspect the matched set, not just take its first element.
- **$R$ sample**: the real nuclear nodes of NO's deck
  (`nuc_real = ((1, 1.0), (1, 1.6), (37, 9.0))`) restricted to
  $R \in [1.6, 6.0]$ bohr. Nodes where the anion state has dropped below
  threshold are "no target" and are reported as such, following the existing
  crossing-slice convention.
- **Continuation**: $s \in \{0.0, 0.1, \dots, 1.0\}$ at $\kappa = 0.3$, then
  $\kappa \in \{0.0, 0.1, 0.2, 0.3, 0.4, 0.5\}$, each walked in $s$ from 0
  until the stop condition fires. Each solve is seeded
  from the previous point's pole, so the walk is continuously connected to the
  known $s = 0$ answer.
- **Truncation**: $N_l \in \{1, 2, 3, 4\}$ at every point, plus $N_l = 5$ at
  whichever point moved most.

### Outputs

1. $E_{\rm res}(R)$ and $\Gamma(R)$ for **full** and **fixed-$l$** across the
   continuation, with the $N_l$ convergence study.
2. The pole *count* in the window at each point — the split test.
3. **The cheap observable.** Both curves fed into
   `qscat.core.lcp.lcp_ve_cross_section` turn the curve difference into a
   cross-section difference with no coupled scattering solver at all, on
   $E \in [0.02, 0.10]$ Ha in steps of 0.002 (41 energies) across NO's
   resonance, which sits at 0.02–0.05 Ha.

   **VE only, not DA.** NO's dissociative-attachment channel opens at
   $+0.172$ Ha, above the resonance, so $\sigma_{\rm DA}$ is a
   $10^{-19}\,a_0^2$ tail there and the LCP is documented to miss it by five
   to seven orders of magnitude. A quantity that wrong cannot discriminate a
   few-percent change in $\Gamma(R)$.
   This is an approximate route — it is the LCP, whose own error against the
   exact 2-D solver is documented and energy-dependent — but it is a *differential*
   comparison on one grid with one approximation applied to both sides, so it
   measures the effect of the coupling rather than the quality of the LCP.

### The gate

Phase 2 runs **only if at least one** of these holds:

- **(a)** a second angle-stable pole enters the window anywhere on the
  continuation;
- **(b)** $\max_R |\Gamma_{\rm full} - \Gamma_{\rm fixed}| / \Gamma_{\rm fixed}
  > 0.05$ at $\kappa = 0.5$ and the largest $s$ **both** branches reached — not
  at $s = 1$, which the walk need not reach, and not at each branch's own
  endpoint, which would silently compare different $s$;
- **(c)** the LCP vibrational-excitation cross sections built from the two
  curves differ by more than 5 % at any sampled energy, elastic or
  $0 \to 1$.

The 5 % threshold is chosen against two known scales: the resonance curves in
this repository converge to $10^{-9}$–$10^{-7}$ Ha, so 5 % is far above the
numerical floor; and the approximations already in production (LCP vs exact)
depart by factors, so an effect below 5 % could not change any conclusion those
approximations are used for.

**If the gate does not open, that is the deliverable.** Phase 2 is written up as
deliberately not run, and the result is recorded as a measured negative: the
fixed-$l$ decomposition is sound for a NO-like model over the full geometric
range of the anisotropy. This outcome must be reported as prominently as the
positive one.

## Phase 2 — the 2-D coupled poles (container / MUMPS)

Runs at the $(s, \kappa)$ the screen selected — the point of largest effect, and
the $\kappa = 0$ homonuclear control at the same $s$ for contrast.

**What needs no change.** `qscat.core.exact_resonance_states` touches its grids
only to validate `ndim` and to read `el_0, nu_0`; everything else goes through
`model.hamiltonian(tgrid)`. A coupled model satisfying the existing
`ResonanceModel` shape and returning the block matrix therefore slots in as-is,
with `qscat.core.grids.ecs_angle_family` supplying the three-grid family
unchanged. The two-angle stability criterion is unaffected by the channel
structure.

**What does need work.** The verification layer is single-channel:

- `qscat.core.bo`'s reference states $\phi_j(r;R)\chi_v(R)$ must be embedded in
  the $l = 1$ channel (zero elsewhere), or the overlap taken channel-wise and
  summed. Which of the two is used must be stated with the result, because they
  answer different questions — "does the coupled pole still look like the
  single-wave state" versus "is it a resonance at all".
- `real_weight`'s mask must be tiled over channels.

**Seeds** come from Phase 1's *coupled* curve through
`qscat.core.lcp.resonance_levels`. Seeds are passed in, so the exact solver
still never calls the approximation it exists to measure.

**Truncation**: $N_l = 2$ first, since that already doubles NO's exact-2-D deck;
$N_l = 3$ as the convergence check only if $N_l = 2$ shows an effect. Peak RSS
and LU fill are recorded, since they decide whether the follow-on cross-section
spec is affordable at all.

**Deliverable**: coupled versus fixed-$l$ pole positions and widths at the
selected parameters, with the angle-stability residuals and the overlap verdicts
both reported — angle stability alone is necessary and not sufficient, as the
H₂⁺ campaign established.

## Files

Mirroring the potential-factory split (machinery in `projects/`, campaigns in
`validation/`):

| File | Responsibility |
|---|---|
| `projects/no_coupled_channels/anisotropy.py` | $v_\lambda$ by quadrature, the Gaunt coefficient table, and the `ive`-based closed-form oracle |
| `projects/no_coupled_channels/blocks.py` | the channel-block filler, shared by the 1-D and 2-D assemblers |
| `projects/no_coupled_channels/model.py` | `TwoCentreCoupledModel` — wraps a `ResonanceModel`, exposes `hamiltonian(tgrid)` and the fixed-$R$ electronic assembler |
| `validation/coupled/screen.py` | Phase 1 campaign and its report |
| `validation/coupled/exact2d.py` | Phase 2 campaign and its report |
| `validation/coupled/figures.py` | the pole-trajectory figure and the curve overlay |

`projects/` must not import `validation/` (`tests/test_layering.py`). Nothing
enters `qscat` in this spec; promotion is a separate decision once the block
filler has had a second consumer.

## Validation gates

In order of how much their absence would invalidate the result:

1. **$s = 0$ reproduces `qscat.model.NO` exactly.** The coupled $H$ is
   block-diagonal to round-off and its $l = 1$ block equals
   `NO.hamiltonian(tgrid)`. This replaces Houfek's certification and it is an
   identity, not a tolerance.
2. **The screen's $s = 0$ curve equals the published NO curve** produced by
   `validation.factory.base_experiments --molecule NO --stage curves`, to that
   campaign's own convergence tolerance.
3. **$v_\lambda$ against the closed-form oracle** at real argument, over the
   $(r, R)$ range of the deck and $\lambda \le 10$ (the largest
   $l + l'$ reachable at $N_l = 5$).
4. **Quadrature convergence**: 64 nodes against 128, agreeing to $10^{-12}$
   relative to $\max|v_\lambda|$.
5. **$\kappa = 0$ gives exactly zero odd-$\lambda$ components** — checked on the
   coefficients, not on an observable.
6. **$N_l$ convergence** on the reported pole, never on which wave hosts the
   resonance. The note's warning stands: the non-resonant waves carry the
   background whose interference shapes the profile, and truncating on
   "resonance-hosting" would bake in an approximation of the same character as
   the one under test.
7. **Complex symmetry** ($H = H^T$, never Hermitian) and ECS angle stability, as
   everywhere else in the repository.

## Where the continuation actually ends

Measured before Phase 1 was written, on the coupled model at $\kappa = 0.3$,
$N_l = 2$, with no acceptance cut applied:

| $s$ | $R=2.0$: $\varepsilon$ / $\Gamma$ / ratio | $R=2.5$ | $R=3.0$ |
|---|---|---|---|
| 0.0 | 0.054 / 0.037 / 0.69 | bound | bound |
| 0.3 | 0.084 / 0.079 / 0.94 | 0.035 / 0.020 / 0.56 | 0.044 / 0.028 / 0.64 |
| 0.4 | 0.103 / 0.114 / 1.11 | 0.071 / 0.062 / 0.87 | 0.089 / 0.092 / 1.04 |
| 0.7 | 0.154 / 0.280 / 1.82 | 0.149 / 0.288 / 1.93 | 0.155 / 0.387 / 2.50 |
| 1.0 | 0.164 / 0.478 / 2.91 | — | — |

Three findings, all of which change this spec rather than merely informing it.

**The resonance energy stays accessible; the resonance does not.**
$\varepsilon = E_{\rm res} - v_0(R)$ never leaves 0.05–0.16 Ha, straddling the
0.02–0.10 Ha window the cross sections probe. But $\Gamma/\varepsilon$ crosses
1 at $s \approx 0.35$–$0.45$ at every $R$ and reaches ~1.8–2.9 by $s = 0.7$–1.
At $s = 1$ the width is 0.478 Ha — five times the whole probe window — so the
state cannot produce resolvable structure there whatever a pole finder reports.
**The walk therefore stops on $\Gamma/\varepsilon \ge 1$**, and where it stops
is a result, not a failure.

**Stopping on pole-finder failure instead would have measured the deck.** On the
published NO deck (6 tail elements, 35°/44°) the pole is "lost" at $s = 0.5$. On
this spec's deck (8 elements, 44°/52°) the two-angle residual at the *same*
physical point ($s = 0.4$, $R = 2.0$) is 3.0e-9 against 6e-4 — five orders — and
the pole is followed to $s = 1$. Nothing physical happens at $s = 0.5$. The
44°/52° pair exceeds the usual 45° ceiling, which is safe here because the
binding constraint is that the Gaussian stay bounded on the contour
($\mathrm{Re}(z^2) \ge 0$), and that is a joint condition on angle *and* tail
extent: this contour holds $\min \mathrm{Re}(z^2) = 0$, while 50° with a tail
run out to $|z| = 1111$ overflows the potential to 1.75e259.

**A spurious state sits in the search window.** A near-threshold state
($\varepsilon \approx 0.001$, $\Gamma \approx 0.006$) is picked by a
nearest-to-seed rule whenever the true pole becomes unrepresentable, silently
and with a plausible value. Its residual is 7e-4 against the genuine pole's
~1e-9, so the campaign guards on residual and records a point above the cut as
*no target* rather than as a pole.

**And the crossing runs the other way.** At $R \ge 2.5$ the $s = 0$ state is
*bound* ($\Gamma = 0$, residual $10^{-12}$); the anisotropy is what makes it
resonant. The continuation crosses the bound-to-resonant threshold from below,
the opposite direction to what this spec first assumed.

## What this does not establish

- **Nothing about real NO.** The anisotropy is geometric — the wells sit on the
  nuclei — not fitted to any measured or computed property of the molecule. The
  result is about the fixed-$l$ *reduction*, on a model chosen to make that
  reduction testable.
- **No coupled cross section.** The observable route here is the LCP applied to
  both curves, and only in vibrational excitation. The exact coupled scattering solution — a multichannel
  `channel_vector`, a T-matrix in $l \otimes v$ space, and the energy sweep —
  is deliberately out of scope and is the follow-on spec.
- **$\Lambda = 0$ is not run**, so the strongest possible non-resonant
  background (a barrier-free $l = 0$ channel) is absent from the comparison.
- **No rotation.** $\Lambda$ is exactly conserved here because the model is
  fixed-axis and rotationless. Coriolis mixing between $\Lambda$ blocks is a
  real effect in a real molecule and is absent by construction, consistently
  with every other model in the repository.
- **Nothing is compared with experiment**, here or in the follow-on.

## Open questions deliberately left to the follow-on

- Whether the single-wave model *re-fitted* to match the multi-wave one does
  better than the single-wave model *as used* — this spec measures the second.
  The first is the fairer test of "how good can the approximation be made", and
  it needs the cross-section route to be meaningful.
- Whether the time-dependent route extends unchanged. Crank–Nicolson and the
  Padé steppers do not care that $H$ is block-structured, so it probably does,
  which would give a second independent numerical route — but "probably" is not
  a validated claim and it is not tested here.
