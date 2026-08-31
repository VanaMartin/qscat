# The coupled cross section: does the fixed-$l$ reduction change the observable?

Date: 2026-08-28

## Purpose

The preceding phase measured what the fixed-partial-wave reduction costs in the
*resonance curve*: at $s=0.3$, $\kappa=0.5$ the coupled and fixed-$l$ models
differ by 58 % in $\Gamma(R)$ with the channel truncation converged to 0.2 %.
It could not measure what that costs in an *observable*. The route it used —
applying the curve difference to the shipped local-complex-potential baseline —
is valid only while that difference stays a perturbation of the width it
modifies, and past $s \approx 0.2$ it is twice that width, at which point
$\Gamma$ clamps to zero across the doorway and there is no cross section left
to compare.

So the observable question is answerable only by solving the coupled scattering
problem exactly. That is this spec: $\sigma_{\rm VE}(E)$ from the exact 2-D
driven equation with the partial waves coupled, against the same quantity with
one partial wave, on the same deck.

Nothing is fitted to, or compared with, experiment. The anisotropy is
geometric — the wells sit on the nuclei — not calibrated to any measured
property of NO.

## What the preceding phase already implies

Measured at $s = 0.3$, $\kappa = 0.5$ over the sampled $R$:

| | $\varepsilon$ range (Ha) | $\Gamma$ range (Ha) |
|---|---|---|
| fixed-$l$ ($N_l = 1$) | 0.040 – 0.150 | **0.0239** – 0.328 |
| coupled ($N_l = 4$) | 0.020 – 0.111 | **0.0082** – 0.135 |

**The coupling pulls the resonance down and makes it narrower** — three times
narrower at the narrow end, which is the opposite of what "more channels means
more decay routes" would suggest. It is consistent with the preceding phase's
other observation that raising $N_l$ reduced the number of $R$ at which the
pole could not be resolved at all (20, 17, 11, 8 of 41 for $N_l = 1,2,3,4$).

That sets up the prediction this spec exists to test. NO's neutral vibrational
spacing is 7.1–8.5 mHa (it falls with $v$; 8.5 mHa over the lowest levels). The
coupled model's narrowest width, 8.2 mHa, is **comparable to that spacing** —
the boomerang condition, where oscillations survive. The fixed-$l$ model's
narrowest width, 23.9 mHa, is about three times the spacing, which washes them
out.

> **Prediction:** the fixed-$l$ reduction should produce a substantially
> *smoother* cross section than the coupled model. If it holds, the reduction
> does not merely get $\Gamma$ wrong — it destroys resolvable structure.

This is a prediction, not a result. It is written down before the run so that
confirming it is a test rather than a description.

## Decisions taken in brainstorming (2026-08-28)

| Question | Decision |
|---|---|
| Which observable | **Vibrational excitation.** NO's dissociative-attachment threshold is $+0.172$ Ha, above the whole sweep, so DA is shut and the LCP misses that tail by five to seven orders — it cannot discriminate anything. |
| Which comparison | **Same entrance, both exits reported.** Entrance in $l = \Lambda = 1$ for both models; the coupled one reports the total (summed over exit partial waves) and the $l\to l$ restricted amplitude. Both come from one solve. |
| Deck | NO's production deck, $132 \times 597 = 78{,}804$ per channel. |
| Truncation | $N_l = 4$ for the oracle, $N_l = 3$ as the convergence check. **Not $N_l = 2$** — measured, $N_l = 2$ is ~21 % from converged. |
| Energy mesh | Threshold-aware: 0.25 mHa background over $[0.002, 0.150]$ Ha plus 21 points at 0.05 mHa around each of the 20 vibrational thresholds. 1008 energies. |
| Where it runs | **sadaharu, with MUMPS.** Not negotiable — see "Why it cannot run here". |
| Where the code lives | Toy stage: `projects/no_coupled_channels/` and `validation/coupled/`. Nothing enters `qscat`. |

## Why it cannot run here

Measured on this laptop, building the coupled 2-D Hamiltonian on NO's deck and
factorising $E_{\rm tot}\mathbb{1} - H$ with SuperLU:

| $N_l$ | unknowns | nnz | factor | solve |
|---|---|---|---|---|
| 1 | 78,804 | 1,795,536 | 36 s | 0.30 s |
| 2 | 157,608 | 3,736,782 | 208 s | 1.11 s |

Doubling the unknowns cost 5.8× in factorisation time — superlinear fill-in —
which extrapolates to roughly 20 minutes per factorisation at $N_l = 4$
(315,216 unknowns). And `SparseLU.refactor` gives **no reuse on the SuperLU
backend**: it re-runs `splu`, so every energy would pay that in full. A
1008-energy sweep is therefore not merely slow here, it is infeasible.

MUMPS changes both terms: it was measured 72× faster in factor time and 9×
smaller in peak RSS on exactly this class of ECS complex-symmetric matrices,
and `factor(reuse_analysis=True)` skips the ordering, so the sweep pays the
analysis once and a numeric factorisation per energy. sadaharu has 32 cores and
123 GB, against a working set of a few GB at this size.

**Cost estimate:** ~15 s per energy at $N_l = 4$, so 4.2 h for the coupled
model, ~1.2 h for fixed-$l$ (four times cheaper per solve), ~3 h for the
$N_l = 3$ convergence check. About nine hours in total — an overnight run.
The estimate is from the SuperLU scaling and the published MUMPS ratio; the
first thing the campaign does is measure one energy and report the real figure
before committing to the sweep.

## The energy mesh

$E \in [0.002, 0.150]$ Ha. That span covers the near-threshold region, both
models' full resonance range (coupled 0.020–0.111, fixed-$l$ 0.040–0.150), and
twenty vibrational thresholds. DA stays shut throughout.

Two scales have to be resolved and they differ by an order of magnitude:

- **Threshold cusps** are non-analytic *at* each channel opening, so they need
  points bracketing the threshold tightly rather than a fine mesh everywhere.
- **Overlapping resonances.** The vibrational levels are spaced 7–8.5 mHa while
  the widths run 8–135 mHa, so the resonances genuinely overlap and their
  interference structure is spread across the range, on the scale of a width.

Hence background plus clusters: 0.25 mHa background (30+ points across the
narrowest width) and 21 points at 0.05 mHa spanning $\pm 0.5$ mHa around each
threshold. 1008 energies, against 2961 for a uniform 0.05 mHa mesh that would
cost three times as much to resolve twenty places.

**Implementation note:** the mesh generator must deduplicate with a *tolerance*,
not by exact equality. A first pass rounding to 9 decimals left pairs of
energies 0.0004 mHa apart — harmless numerically, but they are wasted solves at
15 s each.

Channels reported: elastic and $0 \to 1 \ldots 4$.

## Architecture

Three pieces. The first two are new; the third is the campaign.

### A channel vector that knows its block

`qscat.core.channels.channel_vector` builds $F_{E,l}(r)\chi_v(R)$ as a vector of
length $N_r N_R$. A coupled state is $N_l$ such blocks, channel-outermost, so
the coupled channel vector embeds that in a chosen block and zeros the rest.

The kinematics stay simple: $k = \sqrt{2(E - \varepsilon_v)}$ depends on $v$
alone. The partial wave changes the Bessel order in $F_{E,l}$, not the momentum,
so every channel at a given $v$ shares one $k$.

### A coupled driven solve

Mirrors `qscat.core.driven.ve_cross_section`: build $H$ from
`CoupledModel.hamiltonian(tgrid)`, form the driving term from the entrance
vector in block $l_i$, solve $(E_{\rm tot}\mathbb{1} - H)\Psi^+ = V\psi_i$ once
per energy with `SparseLU` built at the first energy and `refactor`ed
thereafter, and project onto each $(v', l_f)$ with the non-conjugated
c-product. The post-form T-matrix and the $\sigma = 4\pi^3|T|^2/2E$
normalisation are unchanged from the single-channel case.

It returns two cross sections per $(v', E)$:

- **total**, $\sum_{l_f} |T_{v'l_f,\,v l_i}|^2$ — what an angle-integrated
  measurement sees, and the like-for-like partner of the fixed-$l$ model's
  single exit;
- **restricted**, $l_f = l_i$ only — which isolates how the coupling changes the
  entrance amplitude through virtual excursions into other waves.

Their difference is the flux redistributed into other partial waves, and it is
free.

This duplicates roughly forty lines of sweep boilerplate from `driven.py`. That
is deliberate at toy stage: generalising a shipped solver before the coupled
shape has been used twice is what the lifecycle exists to prevent. Promotion —
either a channel-aware `driven.ve_cross_section` or a shared sweep helper — is
the follow-on, and the duplication should be called out in the note rather than
left for a reader to notice.

### The campaign

Runs at $(s, \kappa) = (0.3, 0.5)$: the anisotropy where the preceding phase's
oracle is converged to 0.2 % and the width difference is 58 %, over all 41
comparable $R$. Three models over the full mesh, on one deck — $N_l = 1$ (fixed-$l$),
$N_l = 3$, $N_l = 4$.

The two identity gates below add single-energy runs only, at parameters the
sweep does not otherwise visit: $s = 0$ (any $N_l$) and $\kappa = 0$ at
$N_l = 2$ against $N_l = 1$. Neither is a fourth sweep; together they cost a
handful of solves.

## Validation gates

In order of how much their absence would invalidate the result:

1. **$s = 0$: the coupled cross section equals the fixed-$l$ one exactly.** At
   zero anisotropy the models are the same Hamiltonian, so this is an identity,
   not a tolerance — and it now runs end-to-end through the solve, the
   projection and the normalisation rather than stopping at the Hamiltonian as
   the preceding phase's embedding gate did. One energy.
2. **$\kappa = 0$: $N_l = 2$ equals $N_l = 1$.** Only even Legendre components
   survive a symmetric well, so $l = 1$ cannot reach $l = 2$ at *any*
   anisotropy. Also one energy, also an identity.
3. **$N_l = 3$ against $N_l = 4$** over the whole sweep — the truncation must be
   converged where the answer is quoted, and the preceding phase measured
   $N_l = 2$ as inadequate for this.
4. **The single-channel path reproduces the shipped solver.** With $N_l = 1$ and
   $s = 0$ the coupled route must agree with `qscat.core.driven.ve_cross_section`
   on the same deck and energies — which is what certifies the forty duplicated
   lines.
5. Complex symmetry and the usual ECS checks, as everywhere.

## What this does not establish

- **Nothing about real NO.** The anisotropy is geometric, not fitted, so the
  result is the cost of the reduction in a NO-like model at a chosen
  anisotropy — not the error in NO.
- **One entrance channel.** Both models enter in $l = 1$. A physical cross
  section averages over an entrance distribution the model does not define.
- **One $(s, \kappa)$.** The preceding phase's sweep showed the width difference
  is 55–65 % across $s = 0.2$–0.5, so the point is representative rather than
  special — but this spec measures the observable at one anisotropy only.
- **$\Lambda = 0$ is not run**, so the barrier-free $l = 0$ background never
  enters.
- **No rotation.** $\Lambda$ is exactly conserved because the model is
  fixed-axis; Coriolis mixing is absent by construction, as in every model in
  this repository.
- **No dissociative attachment**, which is energetically shut for NO over this
  range.

## Open questions deliberately left to the follow-on

- Whether the coupled solver should be promoted into `qscat.core` as a
  channel-aware `ve_cross_section`, and whether the single-channel path should
  then be routed through it. That decision wants a second consumer, which this
  spec does not provide.
- Whether the 2-D coupled *resonance states* — the quasi-bound levels of the
  coupled $(r, R)$ problem, which the preceding phase's spec had gated as its
  Phase 2 — say the same thing as the cross section about where the structure
  is. Cheaper per solve than a sweep, and the natural cross-check if the
  smoothing prediction holds.
