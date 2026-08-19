# Angular (θ) extension: coupled partial-wave channels — research direction

**Status:** *not implemented.* This is a captured research direction, recorded
2026-07-23 for later work. No code, spec, or plan exists yet; the brainstorming
around the physics and the model construction is still to be done.
**Relates to:** `docs/physics/n2-resonance.md`,
`docs/physics/n2-cross-section.md`, `docs/physics/nd-tensor-hamiltonian.md`.
**Units:** atomic units throughout.

## The question to answer

Every N₂ model in this repo so far fixes a single electronic partial wave: the
²Π_g shape resonance is represented by one `l` (the `d`-wave), and the
electron-molecule interaction `-lambda(R) * exp(-alpha_c * r^2)` is isotropic
in the electron's angular coordinate. The direction recorded here adds `theta`
(the electron's polar angle relative to the molecular axis) by expanding in
spherical harmonics and solving the resulting **coupled-channel** problem.

The purpose is not a better N₂ number. It is the same purpose as the rest of
the program: **validate an approximation by computing the thing it
approximates.** Concretely —

> Treat the single-partial-wave model as an *approximation* of the multi-wave
> model, and measure how the fixed-`l` decomposition into a resonant model
> affects the observable.

That is a well-posed validation question of exactly the same shape as
"LCP vs. the exact 2-D solution": one model is the oracle, the other is under
test, and the difference between them is the result. It is currently
unanswerable, because with one channel there is nothing to compare against.

## Why the decomposition is sound

For a linear molecule in the body-fixed frame, the projection `Lambda` of the
electronic angular momentum on the molecular axis is **exactly conserved**
(D_inf_h symmetry). Restricting to a single `Lambda` block therefore costs
nothing — it is a symmetry, not an approximation. Within that block, expanding
the interaction in Legendre components and projecting onto spherical harmonics
gives

```text
V(r, theta_e, R) = sum_lambda v_lambda(r, R) * P_lambda(cos theta_e)

V_{l l'}(r, R)   = sum_lambda <Y_{l Lambda}| P_lambda |Y_{l' Lambda}> * v_lambda(r, R)
```

which is **exact until the sum over `l` is truncated**. For the ²Π_g resonance,
`Lambda = 1` and `l = 2, 4, 6, ...`. The truncation order `N_l` then becomes an
ordinary convergence parameter, studied exactly like box size `R_max`, the ECS
radius `R_0`, and the ECS angle `theta` already are — which is the standard
this repo already holds itself to.

## The truncation criterion — the one thing to get right

The tempting rule is *"keep the partial waves that host a resonance, discard
the rest."* **Do not use that rule.** The non-resonant waves carry the
**background** scattering, and it is the interference between the resonance and
that background which shapes the profile (the Fano-like asymmetry). Discarding
them would reintroduce precisely the deficiency the exact-2-D work exists to
expose — a missing elastic background is one of the two documented NOTEs
against LCP.

The criterion must be **observable convergence**: keep adding `l` until
`sigma_{0->v'}(E)` stops moving by more than tolerance. The two criteria give
different answers, and the resonance-hosting one quietly bakes in an
approximation *of the same character as the one under test*, which would make
the result circular.

## What it would buy

A single-`l` model **cannot, by construction**, test whether the
single-discrete-state Feshbach assumption underlying both LCP and the nonlocal
resonance model survives when the resonance can leak across several partial
waves. Multi-channel leakage is a genuinely new failure mode — not a more
accurate number — and it is invisible to every model currently in the repo.

## Consequences to plan for

1. **It breaks a library assumption.** `qscat.dvr.potential_nd` /
   `hamiltonian_nd` assume a **diagonal** potential; a coupled `V_{l l'}` is
   exactly what is not diagonal. `qscat.linalg.kron_sum` is unaffected (the
   channel index enters as a dimension carrying a zero kinetic operator), so
   the shape is

   ```
   H = I_l (x) (T_r (+) T_R) + sum_{l,l'} E_{l l'} (x) diag(V_{l l'}(r, R))
   ```

   i.e. a **block/channel-coupled assembler alongside** the diagonal one —
   still sparse, still Kronecker-structured. An addition, not a rewrite, but
   one worth designing deliberately rather than bolting on.

2. **Cost scales as `N_l^2` in nonzeros.** Dimension grows by `N_l`, nnz by
   roughly `N_l^2`. At `N ~ 1e5` with four partial waves that is `4e5`
   unknowns and ~16x the nonzeros. Sparse-LU **fill-in** is where this bites,
   which makes the `MMD_AT_PLUS_A` ordering result load-bearing and may force
   iterative/preconditioned solvers from optional to necessary.

3. **It adds channels, not dimensions.** This remains a 2-D continuum problem
   `(r, R)` with `N_l` channels. It is *complementary to*, not the same as,
   the higher-dimensional direction (a second nuclear coordinate — a bend, or
   a triatomic stretch — giving a true `(r, R_1, R_2)` problem). Only the
   latter exercises the D-general library layer.

4. **It costs the external oracle.** The current interaction is isotropic and
   generates **no coupling at all**, so the anisotropy would have to be
   *designed* (at minimum a `lambda = 0` and a `lambda = 2` component). That
   is legitimate — the model is chosen, not derived — but Houfek's independent
   `CSVE.V00.J00` data would no longer certify the solver. The cheap
   replacement check: **switching the anisotropy off must reproduce the
   current single-channel model exactly.**

## Lead worth scouting first

`reference/eMoScat` contains `source/Model2d/CoupledModel2d.cpp`,
`include/Model2d/CoupledModel2d.h`, and `input/coupled/` decks — coupled-model
machinery that exists but was apparently **never exercised for N₂**. Run
`port-scout` over it before designing the coupling structure; it may already
encode the intended `v_lambda` form and channel bookkeeping.

## Open questions for the brainstorming session

- What functional form for `v_lambda(r, R)`? How is the anisotropy strength
  chosen, and against what is it calibrated?
- How many `lambda` terms, and does the `lambda` truncation converge
  independently of the `l` truncation?
- Does the resonance pole stay a single pole under coupling, or split? (The
  two-spectrum matcher `qscat.ecs.find_resonance_pole` assumes one.)
- Is the comparison run at fixed model parameters, or is the single-wave model
  re-fitted to best match the multi-wave one? (The second is the fairer test of
  "how good can the approximation be made," the first of "how good is it as
  used.")
- Does the TD route (`docs/physics/n2-td-cross-section.md`) extend
  unchanged? Crank-Nicolson does not care that `H` is block-structured, so it
  probably does — which would give a second independent numerical route for
  free.
