# N₂ vibrationally-elastic/inelastic cross section: time-dependent (Crank-Nicolson) method

**Location:** `qscat.evolution.make_cn_stepper` (the general propagator),
`projects/n2_td_cross_section/` (`propagator.py` -- thin re-export,
`td_cross_section.py` -- the N₂ application), `validation/n2/td_check.py` (the
harness's Group D wiring), `validation/n2/experiment.py` (Group D).
**Origin:** same local complex potential (LCP) model as
`docs/physics/n2-resonance.md`/`docs/physics/n2-cross-section.md`; this is a
second, independent (time-dependent) numerical route to the same physical
observable, `docs/superpowers/specs/2026-07-22-n2-td-cross-section-design.md`.
**Units:** atomic units throughout (energy in Hartree, length in Bohr, time in
`hbar`/Hartree, cross section in bohr²).

## Physical picture

`docs/physics/n2-cross-section.md` computes `sigma_{v_init->v'}(E)` by solving
the driven (resolvent) equation `(E_tot*I - H_res) @ xi = d_{v_init}` directly
for each collision energy `E` -- a **time-independent (TI)** approach. This
document is a **time-dependent (TD)** route to the identical observable:
prepare the doorway state as a wavepacket at `t=0`, propagate it forward under
the same non-Hermitian nuclear Hamiltonian `H_res`, record its correlation
with each doorway function `d_v'(R)` over time, and Fourier-transform that
correlation function into the energy domain. The two methods are related by

```
S_TD(E) = (1/i) * integral_0^inf exp(i*E_tot*t) * <d_v'|exp(-i*H_res*t)|d_v_init> dt
        = <d_v'|(E_tot*I - H_res)^-1|d_v_init>
        = S_TI(E)
```

in the long-propagation-time limit -- the Laplace transform of the propagator
is exactly the resolvent. TD and TI are not two different physical models;
they are two different numerical routes to the same `S`-matrix, and validating
one against the other is a genuine cross-check (see "Validation" below), not
a tautology, because the discretizations (finite `dt`, finite propagation
time `T`, vs. a single direct linear solve) are entirely different failure
modes.

## Method

1. **Crank-Nicolson propagator** (`qscat.evolution.make_cn_stepper(H, dt)`):
   the general, N₂-independent primitive. Advances `psi_{n+1}` from `psi_n`
   via the Cayley form `(I + i*H*dt/2) @ psi_{n+1} = (I - i*H*dt/2) @ psi_n`,
   exact to `O(dt^3)` per step, unconditionally stable, unitary for Hermitian
   `H`, and norm-decaying when `H` has a negative-imaginary-part
   (absorbing/optical) component. `H` may be a general complex,
   non-Hermitian matrix; the LHS matrix is LU-factored once per call, and
   each `stepper(psi)` call reuses that factorization
   (`scipy.linalg.lu_factor`/`lu_solve`).
2. **Initial wavepacket**: `psi(0) = d_{v_init}(R)`, the same doorway function
   `d_v(R) = sqrt(Gamma(R)/(2*pi)) * chi_v(R)` as the TI solver
   (`docs/physics/n2-cross-section.md`).
3. **Propagation**: `psi(t)` evolves under the *time-independent*, non-
   Hermitian `H_res = T_nuc(mu) + diag(V_d(R) - i*Gamma(R)/2)` -- identical
   to the TI solver's `H_res`, just propagated rather than resolvent-solved.
   The `-i*Gamma/2` term makes `||psi(t)||` decay: the resonance depletes
   into the dissociative-attachment/autodetachment channels as the ECS
   nuclear grid absorbs outgoing flux.
4. **Correlation function**: `c_{v'}(t_n) = <d_{v'}|psi(t_n)>` recorded at
   each time step -- the DVR **c-product** (a plain coefficient dot, no
   conjugate), matching the TI solver's S-matrix convention (the DVR basis
   is already `1/sqrt(weight)`-normalized, and `psi(t_n)` is a genuinely
   complex ECS-driven state, not a Hermitian-normalized eigenvector).
5. **Energy transform**:
   `S_{v'}(E) = (1/i) * sum_n w_n * exp(i*(E + eps[v_init])*t_n) * c_{v'}(t_n) * dt`,
   with composite Simpson weights `w_n` (trapezoidal fallback for an even
   sample count), then
   `sigma_{v_init->v'}(E) = 4*pi^3*|S|^2/(2*E)`, zero if `E<=0` or the final
   channel is energetically closed.

Because the correlation functions `c_v'(t_n)` are `E`-independent, ONE
propagation from `v_init=0` yields `sigma_TD` at every `(E, v')` pair of
interest -- not one propagation per anchor (`td_cross_section.td_ve_cross_section`
accepts `E` and `vprimes` as arrays/lists for exactly this reason).

## Propagation config: `dt`, `n_steps`, convergence

The resonance's own eigenmodes of `H_res` sit at `Re(E) ~ -0.7..-0.4` Ha
(`eps[0] ~ -0.745` Ha shifted by the ~2.3-2.4 eV Π_g resonance), and the
Crank-Nicolson Cayley-transform phase error per step grows as `~(E*dt)^3`,
accumulating over `n_steps` steps -- so accuracy requires `dt` small relative
to `1/|E|`, not just `n_steps` large.

- `T = n_steps*dt = 1500` a.u. is long enough to deplete the resonance
  (`Gamma(R0) ~ 0.017` Ha gives a decay time `~1/(Gamma/2) ~ 120` a.u., so
  `T=1500` is >10 decay times) -- confirmed empirically:
  `||psi(T)||/||psi(0)|| ~ 1e-2`, comfortably `< 0.1`.
- `dt = 0.025` a.u. (`n_steps = 60000`) keeps the per-step Cayley phase error
  for these modes under ~1e-4 rad/step.
- **Convergence check** (`test_td_cross_section.py`'s V2): halving `dt` (same
  total `T`, so `n_steps` doubles to 120000) changes `sigma_TD` at
  `(E=0.1, v'=1)` by `< 5%` -- the discretization is in the converged regime,
  not still drifting.

This same `(dt, n_steps) = (0.025, 60000)` config is reused, unchanged, for
Group D of the harness (`validation/n2/td_check.py`): one propagation costs
~9s (amortized across all 4 gated anchors), on top of the ~7s `vres_on_grid`
cost the harness already pays once for Group C5 (`_build_system` is
`functools.lru_cache`d and shared between `cross_section.py` and
`td_check.py`).

## Validation: TD vs. TI (an exact differential oracle, not a cross-model comparison)

Unlike the TI-vs-Houfek comparison in `docs/physics/n2-cross-section.md`
(genuinely different methods/dimensionality, so only loose, factor-of-3
agreement is expected), TD and TI solve for the *same* `S`-matrix by
construction (see "Physical picture" above) -- so `sigma_TD` is checked
against `sigma_TI` as an **exact differential oracle**, gated at `rtol <=
0.10` (the residual is entirely the finite-`dt`/finite-`T` discretization,
already shown converged by the V2 check above).

`projects/n2_td_cross_section/test_td_cross_section.py` (V1) at the tuned,
converged config:

| E (Ha) | v' | sigma_TD (bohr²) | sigma_TI (bohr²) | ratio |
|---|---|---|---|---|
| 0.1 | 1 | 6.223 | 6.182 | 1.007 |
| 0.2 | 2 | 9.203e-3 | 9.313e-3 | 0.988 |

`validation/n2/td_check.py` extends this to all 4 GATED C5 anchors (Group D
of the harness), also checking the same Houfek factor-3 bound C5 uses (since
TD shares the TI solver's `V_d(R)`/`Gamma(R)`/doorway machinery, it inherits
the same 2 documented LCP-model limitations for the elastic/near-threshold
anchors -- so only the 4 GATED anchors are re-checked under Group D, not the
2 DOCUMENTED-LIMITED ones):

| E (Ha) | v' | sigma_TD (bohr²) | sigma_TI (bohr²) | ratio TD/TI | Houfek (bohr²) | ratio TD/Houfek |
|---|---|---|---|---|---|---|
| 0.2 | 1 | 5.605e-2 | 5.593e-2 | 1.002 | 1.257e-1 | 0.446 |
| 0.2 | 2 | 9.203e-3 | 9.313e-3 | 0.988 | 1.203e-2 | 0.765 |
| 0.2 | 3 | 1.876e-3 | 1.812e-3 | 1.035 | 2.193e-3 | 0.855 |
| 0.1 | 1 | 6.223 | 6.182 | 1.007 | 6.121 | 1.017 |

All 4 within `rtol=0.10` of `sigma_TI` and within the factor-3 Houfek bound
-- **Group D is PASS**, harness now 0 PENDING.

## Promotion note

`make_cn_stepper` has no FEM-DVR-ECS or N₂-specific structure (it operates on
a generic complex matrix `H` and a scalar `dt`), so it was promoted directly
to `qscat.evolution.make_cn_stepper` rather than staying project-local;
`projects/n2_td_cross_section/propagator.py` is now a thin re-export
preserving the project's original import path.

## Validation summary

- `libs/qscat/tests/test_crank_nicolson.py`: exact-`expm` match for small
  `dt`, unitarity (norm conservation) for Hermitian `H`, norm decay for
  non-Hermitian (absorbing) `H` -- **PASS**.
- `projects/n2_td_cross_section/test_td_cross_section.py`: V1 (TD ~= TI at
  two anchors, `rtol<=0.10`, physical `sigma`), V2 (convergence in `dt`,
  `rtol<=0.05`; resonance depletion, `norm_ratio<0.1`) -- **PASS**.
- `validation/n2/experiment.py` Group D: all 4 GATED anchors, **PASS**;
  harness exit code `0`, 0 PENDING.
- `.superpowers/sdd/task-1-report.md`, `task-2-report.md`, `task-3-report.md`
  (this sub-project's own numbering): the full numeric record.
