# NRM time-dependent — roadmap

Where the time-dependent nonlocal resonance model stands, and what comes next.
Supersedes the ordering (not the content) of the task list in
`2026-08-19-nrm-td.md`. Written 2026-08-24, after tasks 1–5.

## The premise, restated

The time-dependent route is **slower** than the resolvent — measured on F₂'s
production nuclear deck, 773 s of propagation against 246 s for the same
8-energy sweep time-independently. It is not justified on cost and no stage
below is motivated by speed.

It is justified because `Ψ_d(R,t)` carries `S(t)`, `⟨R⟩_t` and `⟨P⟩_t`, and
those turn a cross-section feature into a mechanism. Váňa & Houfek PRA 95,
022714 (2017) closes by proposing exactly this — time-dependent calculations
"within the LCP approximation or the nonlocal resonance model … and thus
interpret the results in the same way" (p. 022714-16) — as future work it does
not execute.

## Stage 0 — done, PR #31 (`nrm-td` → `main`)

Tasks 1–5: the extended-space resummation, the low-rank launch basis, the
propagation and half-Fourier transform, the vector gate, and the DA cross
section.

| Claim | Measured |
|---|---|
| Elimination identity ⇒ `nonlocal_operator`'s `F(E)` | 4.4e-14 relative |
| Half-Fourier transform vs closed form | 1.8e-13, `dt⁶` |
| `⟨P⟩_t` vs analytic Gaussian | 14 significant figures |
| `Ψ_d^TD` vs `Ψ_d^TI`, N₂, vector to vector | 1.73e-04 |
| `σ_TD/σ_TI`, F₂ DA, 0.010–0.060 Ha | 0.986 – 1.014 |

Six figures in `docs/physics/nrm-time-dependent.md`.

## Stage 1 — in flight, stacked PR (`nrm-td-ve-lcp` → `nrm-td`)

**Task 7 — the Markovian (LCP) limit.** Promoted from seventh to first.
PRA 47 Eq. (2.11) states the LCP *is* the Markovian limit, and the paper's
headline result is a nonlocal-versus-local **packet** comparison. This builds
the local half. It is also cheap: no arms means the propagated Hamiltonian is
`N_R` square (974) rather than `(1+n_states)·N_R` (53570), seconds instead of
~45 minutes. It must settle by measurement which `V_d` enters Eq. (2.15) —
`v_d_discrete` (PRA 77 Eq. 20) or `lcp`'s `E_res`, which differ by 0.0053 Ha on
F₂ and 0.0229 Ha on NO — against the shipped `lcp_da_cross_section`.

**Task 6 — vibrational excitation.** The channel where choice B is exact to
sub-0.7%, where PRA 77 publishes curves for all three molecules, and where the
packet decays **in place** by autodetachment rather than crossing the box —
cheaper and better conditioned than DA, as the N₂ gate already showed. `T^bg`
is energy-domain and static, so only `T^res` changes route.

Sequential, not parallel: both touch `td_cross_section.py`.

## Stage 2 — coverage and plumbing

- **Convergence study** (was task 8, now much smaller). `n_states` is gone as a
  knob — the complete arm set is a correctness constraint, not a tunable. What
  remains is `dt`, `t_max`, Padé order and `rank_tol`, per molecule and channel.
- **One production-electronic-deck run.** Every TD result so far uses a reduced
  electronic grid (55 points, `H_ext` = 53570). The production deck is 132
  points and `H_ext` = 128568, and **has never been propagated on any molecule**.
  One run tests whether the reduced-deck justification holds. sadaharu territory.
- **Three-molecule campaign** (task 9): N₂ VE, F₂ VE and DA, NO VE.
- **NO probe → VE, not DA** (task 10, redirected). NO's DA channel is a known
  open failure for the *time-independent* route, so a TD run there measures an
  already-broken oracle. PRA 77 publishes NO VE curves; that is the honest probe.
- **`qscat-run` wiring** (task 11) as `nrm-td`, and **docs consolidation** (12).

## Stage 2b — THE GATE: reproduce the published figures, TI and TD

**Nothing finalises until both routes reproduce the published figures.**

**The thesis figures are the primary target, and a better gate than PRA 77
alone.** Váňa 2017, **Figs. 3.14 (NO, p. 46), 3.20 (F₂, p. 52) and 3.23 (N₂,
p. 55)** plot, on our exact models, the time-independent 2-D reference together
with all three time-dependent S-matrix methods AND the LCP. They already ARE the
TI-vs-TD agreement gate, on a four-curve structure, and they include **DA panels
PRA 77 does not have**. Fig. 3.14 covers NO VE 0→0, 0→1, 0→8 over 0–0.10 Ha plus
DA over 0.170–0.200 Ha at `t_c = 35000`.

**PRA 77's Figs. 4/5/6/8 are the nonlocal overlay on top**, and its panel
inventory fixes exactly which curves:

| paper figure | discrete state | panels |
|---|---|---|
| Fig. 4 (N₂) | A, physical | VE 0→0, 0→1, 0→8 |
| Fig. 5 (NO) | A | VE 0→1, 0→2 |
| Fig. 6 (F₂) | A | VE 0→1, **DA** |
| Fig. 8 (all three) | B, asymptotic | N₂ 0→0, 0→8; NO 0→1, 0→2; F₂ 0→1, **DA** |

Four curves per panel: exact 2-D, LCP, nonlocal without background, nonlocal
with background. Fig. 4's N₂ panels are LINEAR-scale, and Fig. 8 omits N₂ 0→1
because "results of all calculations are practically the same in this particular
case" — a claim on that linear axis, which our numbers should reproduce rather
than contradict.

**Gaps against the target today.** NO VE has never been run by either route.
N₂ 0→8 has not. Choice A has not been swept for VE. The LCP's VE curve lives in
`projects/n2_ti_cross_section` and is not reachable from the unified surface, so
a four-curve panel cannot currently be assembled from one config.

**Fig. 3.14 already found a defect in our own exact 2-D NO DA** (2026-08-24).
Its bottom panel puts NO σ_DA on a `1e-9 a₀²` axis with the LCP scaled `×1e-5`
to fit. Our LCP matches the published LCP; our NRM choice B matches the published
EXACT to a factor ~1.5; **our exact 2-D is ~7 orders too large.** So
`nonlocal-resonance-model.md` §7.2's "unexplained NRM-B collapse" is not a
collapse — the nonlocal model was right and the oracle was wrong, which is why
auditing the NRM found nothing. Under repair, with the reproduction of Fig. 3.14
as its acceptance test.

**A methodological note worth keeping.** The campaign concluded "no published NO
DA reference exists" from PRA 77's panel inventory. PRA 77 publishes none — but
the thesis does, and it was in `reference/literature/` throughout. Its reference
note indexes the H₂⁺ chapter's figures and not Fig. 3.14, because the note was
scoped to what the repo already consumed. **A reference note that records only
what you have needed cannot tell you what you are missing.** When a result looks
anomalous, re-read the sources for what they contain, not for what the note
extracted.

## Stage 3 — the interpretive sub-project (separate spec + plan)

Its own brainstorm → spec → plan cycle, because it is a research question rather
than a capability. Scoped here only so the questions are not lost:

1. **Packet splitting.** PRA 47 explains the LCP's failure for H₂⁻ as a
   temporary splitting between ≈2 and ≈5 fs. Does the nonlocal packet split on
   N₂/NO/F₂ where the local one does not? PRA 47 had no exact oracle for its
   models; we do, on all three.
2. **Quasibound ↔ peak correspondence *within the nonlocal model*.** PRA 95
   established this for the LCP (p. 022714-15). Doing it nonlocally is new, and
   `lcp_resonance_levels` / `exact_resonance_states` already exist.
3. **Boomerang versus broadening.** PRA 47 (p. 1042) warns that oscillations
   from rapid packet broadening look like N₂-style boomerang oscillations and
   are not the same thing. The diagnostics can tell them apart.
4. **Formation times.** NO's first VE 0→1 peak forms at `t > 10 000` and its
   lowest quasibound state lives `> 30 000` a.u. (PRA 95, p. 022714-15).

## Carried risks and open questions

- **NO's long-lived states** mean very long propagations — a routine sadaharu
  workload, not a laptop one.
- **F₂ DA has a floor of ~2e-2** from the ≥24 near-real modes in its `V_d` well
  (5.08e-3 of the launch norm); `S(t)` plateaus and no absolute survival
  criterion is reachable there. Convergence is `σ` stationary in `T`.
- **The NO DA collapse is RESOLVED, in the opposite direction to the one
  assumed.** The nonlocal model reproduces the published exact result; our exact
  2-D NO DA is ~7 orders too large (Stage 2b). §7.2 and §9 of
  `nonlocal-resonance-model.md`, `diatomic-ve-cross-sections.md`, and the
  committed `no-2d-ti-da-cross-section.png` all state the inverted conclusion and
  must be corrected once the repair lands.
- **Memory, not time, is the binding constraint.** The order-3 Padé stepper
  holds three sparse LU factorizations of `H_ext` for the whole propagation.
  Prefer the MUMPS backend (~9× lower peak RSS) over SuperLU wherever the deck
  is large.
