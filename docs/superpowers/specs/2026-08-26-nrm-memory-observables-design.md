# TD-NRM memory observables — design

**Date:** 2026-08-26
**Status:** design, approved for planning
**Predecessors:** `2026-08-19-nrm-td-design.md` (the time-dependent nonlocal model,
merged). This is sub-project 3, the interpretive half the roadmap
(`2026-08-24-nrm-td-roadmap.md`, Stage 3) scoped and deferred.

## 1. What this is for

The nonlocal resonance model is an approximation between the LCP and the exact
2-D solver, and its strengths and weaknesses are currently known only as
cross-section ratios. A ratio says *how much* a model is wrong; it says nothing
about *what it is doing differently*.

The time-dependent route makes that answerable, because of an accident of how it
is implemented: **the memory is not a stored history, it is state.** The
resummation (`docs/physics/nrm-time-dependent.md` §2.2) turned Eq. (2.1)'s
convolution into auxiliary nuclear packets `φ_n`, one per projected electronic
state. The literal memory-integral form would give only the convolution; here the
memory is a set of wavepackets you can look at.

Three observables follow, all measured on the propagation we already run.

## 2. The observables

### 2.1 Partition — where the amplitude goes

`S_d(t) = ‖Ψ_d‖²` and `S_n(t) = ‖φ_n‖²`, over the real nuclear region.

The LCP's `−iΓ/2` is a sink with no destination: in the local model these
channel populations do not exist. In the nonlocal one they do, resolved by
channel and by `E_n(R)`.

Prototype (N₂, E = 0.10 Ha, 73 arms): the arms peak at **0.19 of `S_d(0)`** at
t ≈ 25 — a fifth of the norm in the electronic continuum at once.

### 2.2 Exchange rate — the observable that has no LCP counterpart

For `i ∂_t Ψ = H Ψ`, `d‖ψ‖²/dt = 2 Im⟨ψ|Hψ⟩`. Restricting to the coupling term:

```
exchange(t) = 2 Im ⟨Ψ_d | Σ_n V_dn φ_n⟩          nonlocal
exchange_markovian(t) = −⟨Ψ_d | Γ_loc | Ψ_d⟩      local limit, ≤ 0 by construction
```

The Markovian form **can only lose**. The nonlocal one can be positive, and in
the prototype it is: **positive at 85 of 4001 steps, first at t = 132, raw
maximum +8.776e-7** — which is `+2.420e-4` when divided by `S_d(0)`, the form
this document quoted unqualified through several revisions and the reason the
shipped `exchange` is unnormalized with all three forms in its docstring. That is
amplitude returning from the continuum into the discrete state — a process the
LCP cannot represent at all, rather than one it merely gets wrong.

This is the headline observable. It converts "the LCP is off by a factor" into a
rate, at a time, from a channel.

### 2.3 Decay law — how wrong a single rate constant is

`S(t)` against `exp(−Γ_eff t)` with the golden-rule `Γ_eff = ⟨Ψ_d(0)|Γ_loc|Ψ_d(0)⟩ / ‖Ψ_d(0)‖²`.

Prototype (N₂): faster than exponential at very short times, then **four orders
slower** by t = 517 (measured 1e-1 against `exp(−Γt)` = 2.2e-5). That reproduces
the behaviour Váňa 2017 Fig. 3.22 reports for N₂ — the packet decays only when it
returns to the autodetachment region — which is what makes it a validation of the
diagnostic as well as a result.

## 3. The Markovian reference comes from `F`, not from the LCP

`Γ_loc` is the local limit of `F(E)` itself — the **`√w`-weighted row sum**.
`nonlocal-resonance-model.md` §9 records it reproducing `Γ(ε_loc, R)` to median
0.977 on NO and 1.011 on F₂ — and **that comparison is against Eq. (68)**,
`2π|V_dk⁺(R; ε_loc)|²` via `coupling.v_dk_plus`, at the same LOCAL electron
energy `ε_loc(R) = E_tot − v₀(R)`:

```
Γ_loc(R_i) = −2 Im[ (F √w)_i / √w_i ]
```

(the `√w` factors are the coefficient-to-value conversion; `diag F` alone is
0.14× the local potential because the kernel spans ~10 nodes.)

**Not `local_complex_potential`'s `Γ`, and not as an oracle for `local_width`
either.** That `Γ` is the width at the RESONANCE energy `ε_res(R) = E_res(R) −
v₀(R)`, a different quantity with different R-dependence: measured on the N₂ gate
deck the ratio `Γ_loc/Γ_LCP` sweeps 0.12 → 8.9 monotonically, because `ε_loc`
crosses zero at R = 2.43 (the channel closes) while `ε_res` stays roughly flat.
No choice of `E` makes them the same comparison. It also carries a pole walk that
freezes silently and, on NO, does not converge in `r_max` at all (3.98e4 spread),
and it is a different construction — comparing against it would confound
nonlocality with the LCP's separate `V_d`/`Γ` build.

**A correction to this document's own first draft.** It quoted a freeze at
R = 1.5783 on the N₂ gate deck. That number came from a prototype call passing
the SAME electronic grid for both `elec_grid_a` and `elec_grid_b`;
`local_complex_potential` takes a pair at two different ECS angles and the
two-angle match is how the physical pole is selected, so a degenerate pair
discriminates nothing. The figure is struck rather than corrected — the route
adopted here needs no pole walk at all. Taking
both sides from the same `F` isolates the one thing under study.

The shipped Markovian propagation (`markovian=True`) keeps using the LCP's
`Vd`/`Gamma` — that is a different question ("how does this differ from the LCP
we ship") and is already answered.

## 4. The measurement that must come first

**The population interpretation is not free, and the prototype asserted it.**
The arms live on an ECS-rotated grid and `H_ext` is complex **symmetric**, not
Hermitian, so the conjugating norm is not conserved and `‖φ_n‖²` is not a
probability in the usual sense.

There is an exact identity to check instead. The amplitude the coupling removes
from `Ψ_d` must be the amplitude it adds to the arms. Writing the coupling block
as `C`, the two rates are `2 Im⟨Ψ_d|Cφ⟩` and `2 Im⟨φ|Cᵀ Ψ_d⟩`, and their sum is

```
4 Σ_n Re[ conj(Ψ_d) φ_n ] · Im[V_dn]
```

— **zero only where `V_dn` is real.** Under ECS it is not, so the imbalance is a
measurable quantity, and it is exactly the leakage of the electronic rotation
into nuclear-space bookkeeping.

**MEASURED 2026-08-27, and the population reading does not survive.** On the N₂
gate deck (73 arms, `n_states=None`, E = 0.10 Ha), the imbalance as a fraction of
the larger of the two one-sided rates is **median 0.822, max 1.057** — O(1), not
small. It is concentrated **entirely in the interaction region** (over 4000 steps: 57 % in R ∈ [2.0, 2.2], 32 % in [2.2, 2.5], 11 % in
[1.8, 2.0] and <0.5 % beyond 2.5, where `Γ_loc` peaks), and the ECS tail carries 2e-55 of it. That is the outcome this section
said would forbid the reading, not the one that confines it harmlessly.

The mechanism: `exchange_arm` is positive and 5–100× larger than `exchange_d`,
so the arms gain conjugating norm far faster than `Ψ_d` loses it. The two rates
are not one transfer — the electronic ECS rotation injects norm through the
coupling itself.

**Consequences, binding on everything downstream:**

- **§2.1 is a RELATIVE channel decomposition, not a population.** Every docstring,
  axis label and sentence must say so. `arm_norm` is where the amplitude goes
  *relative to itself over time and across channels*; it is not a probability and
  the two curves in a partition plot do not sum to anything conserved.
- **§2.2 is unaffected and remains the headline.** `exchange_d` was validated
  against a finite-difference `d‖Ψ_d‖²/dt` — **4.92e-5 (discrete block), 1.71e-5
  (arms)**, 3-point central difference at `h = 0.02` about `t = 20`, on
  **full-grid** block norms — so it is exactly the rate at which `Ψ_d` gains or
  loses through the coupling, well defined regardless of where the amplitude
  goes.

  The scheme and the *time* are both part of the number: the residual is `O(h²)`
  and strongly `t`-dependent, blowing up wherever a rate crosses zero (the arm
  residual runs 4.9e-5 at t = 5, 8.4e-6 at t = 20, 9.5e-4 at t = 50 near a
  crossing). Two independent measurements of "the same" quantity disagreed by an
  order of magnitude for exactly that reason, neither being wrong. Full-grid
  rather than real-region because only there is `d‖ψ‖²/dt = 2 Im⟨ψ|Hψ⟩` exact;
  at `t = 20` no flux has left, so the two agree to every digit. Its SIGN, which is the whole result, is
  untouched.
- **§2.3 is unaffected**, being a property of `S_d(t)` alone.
- **`imbalance` is a first-class output**, the same shape as `exchange`, reported
  per molecule rather than as a footnote.

**A mechanism this measurement handed over for free.** Checking the same identity
on the arms shows their own diagonal-block term is **−21×** their coupling term:
`E_n` is complex, so the arm blocks are strongly dissipative in their own right.
**That is where `H_ext`'s dissipation lives** — not in the discrete block, whose
own term is 4.7e-6 of its coupling. It explains mechanistically what
`nrm-time-dependent.md` §5 established only empirically: truncating the arm set
makes `H_ext` non-dissipative because the arms *are* the dissipation, and
removing channels removes it unevenly (hence non-monotonically).

## 5. API

Opt-in, on the propagation that already produces the state:

```python
propagate_nrm(h_ext, launch, nuclear_grid, *, dt, n_steps, order=3,
              memory: MemorySpec | None = None) -> TdNrmResult
```

`TdNrmResult` gains fields populated only when `memory` is given:

```python
arm_norm:            NDArray | None  # (n_steps+1, n_E) summed over channels
arm_norm_by_channel: NDArray | None  # (n_steps+1, n_arm|k, n_E) — see below
arm_peak:            NDArray | None  # (n_arm, n_E) running max per channel
exchange:            NDArray | None  # (n_steps+1, n_E) nonlocal
exchange_local:      NDArray | None  # (n_steps+1, n_E) Markovian limit of the same F
imbalance:           NDArray | None  # (n_steps+1, n_E) the §4 residual
```

`MemorySpec` carries `gamma_local` and `n_channels`. **`n_channels` defaults to `None`, which records
every channel** — correct wherever it fits, which is most of what this
sub-project runs (`n_arm × (n_steps+1) × n_E × 8` bytes: 2.3 MB for the N₂ gate
deck, ~113 MB for a production F₂ one). An integer `k` keeps the **first `k` in
block order**, i.e. `ing`'s adiabatic-tracking order — **not** the `k` largest by
norm, which is not knowable in one pass and must not be called that.

**`arm_peak` is what actually answers "which channels received the flux."** It is
the running maximum per channel: one pass, `O(n_arm × n_E)` memory, negligible.
With it a truncated run is honest rather than a guess — you can see whether the
first `k` were the right ones, and re-run pointed at the interesting ones if you
want their time series.

Deliberately NOT selecting by `|V_dn|` at t = 0: that sounds principled but
assumes coupling strength predicts which channel ends up carrying amplitude, and
that is a *finding* §6's campaign might test, not a rule to bake into the
recorder. Default off, so the existing hot loop is untouched by callers that do not ask —
every gate in `libs/qscat/tests` propagates.

**Measured cost: +0.33 %** — `record()` takes 0.097–0.100 ms against a 34.2 ms
order-3 Padé step, plus 0.017 ms to reconstruct the arm blocks. **Identical at
`n_channels=None` (73 series) and at `n_channels=4`**: the cost is the two sparse
mat-vecs and the real-region copy, not the per-channel write. **So there is no
cost argument for truncating** — `None` is free as well as correct.

Timed around the call rather than A/B'd across runs, deliberately. A paired
wall-clock A/B of two 4000-step propagations gave +2.8 %, and repeats of the same
configuration gave −0.6 %, +1.5 %, +0.1 %, with `memory=None` runs alone
spreading 1.8–3.0 % among themselves. The signal is about a tenth of that noise,
so no paired run of that size can resolve it, and the first A/B produced both a
wrong number and a plausible wrong explanation for it.

## 6. Scope

**Observables 1–3, on N₂, F₂ and NO.** The comparative result is the point: the
three molecules fail differently in the energy domain (N₂'s LCP is decent, F₂'s
sweeps through unity, NO's is undetermined), and whether that ordering survives
in the time domain is the question this sub-project exists to answer.

**Out of scope, named so it is not lost:** the projection of the exact 2-D
`Ψ(r,R,t)` onto `φ_d` — the time-resolved model error rather than the integrated
one. It needs both propagations on a shared deck and is a sadaharu workload; it
is the natural sequel and the roadmap keeps it.

## 7. What this will not establish

- **Not a validation of the model.** These are diagnostics of a model already
  gated against the exact solver. A striking return-flux signal does not make the
  NRM more correct; it explains what it is doing.
- **NO's `Γ_loc` inherits NO's ingredient problems.** Choice A's ingredient set is
  corrupted at `R = 2.2657` (`nonlocal-resonance-model.md` §11), so NO/A numbers
  stay unquotable here as everywhere else.
- **A frozen or non-converged `Γ` is not this sub-project's to fix.** §3 routes
  around it; NO's pole-walk non-convergence remains open and separately tracked.
