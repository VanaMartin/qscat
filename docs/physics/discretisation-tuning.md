# FEM-DVR-ECS discretisation tuner: design, calibration, and gate

**Location:** `libs/qscat/qscat/tuning` (`analyze`, `mesh`, `ecs`, `probes`, `metrics`,
`propose`, `incident`), the `discretisation-tuner` skill (the supervised loop),
`validation/tuning/` (`calibrate.py`, `test_emoscat_decks.py` — this document's Task 8).
**Origin:** `docs/superpowers/specs/2026-07-28-discretisation-tuner-design.md`; realizes the
future-work note in `docs/physics/diatomic-ve-cross-sections.md`
("Future: automatic discretisation").
**Units:** atomic units throughout (energy in Hartree, length in Bohr).

## Key result

Calibrated once against F₂'s genuinely-open dissociative-attachment channel
($C = 0.10$), the tuner reproduces-and-beats the hand-tuned eMoScat F₂ deck
on the 1-D probes (37% fewer points, clean absolute convergence) and
correctly flags the coarse shared N₂-style grid's historical
under-resolution of the same K≈58–78 wave — the two things this
sub-project set out to prove — and the same clean result holds for H₂⁺'s
proxy deck. The resonance-aware `channel="dissociation"` nuclear path
converges F₂'s 2-D $\sigma_\mathrm{DA}$ on the FIRST a-priori pass
(1.6562 bohr², matching the eMoScat deck and finding #3's own refine²
value) at deck-parity size (1.027×), and sizes H₂⁺'s resonant grid ~4%
under its proxy deck. Two honest caveats carry the same weight as the
result: reaching convergence costs approximately deck-sized resolution —
the deliverable is convergence + automation at deck-competitive size, not
the "10–20% smaller" grid originally hoped for — and N₂/NO's proposed
nuclear grids cost more points than their decks (root-caused to a fixed
real-region extent default, not to $C$), a documented limitation for a
follow-on.

## Why

Discretisation errors have been the single most expensive class of bug in this repo: the
shared N₂-style nuclear grid silently under-resolved F₂'s K≈58 dissociative-attachment (DA)
wave (σ_DA off by ~36 orders of magnitude); the H₂⁺ Coulomb tail forced a 1300-bohr
electronic grid. Every grid before this sub-project was hand-tuned by a human's "good eye"
for element lengths. `qscat.tuning` replaces that with a computed, potential-adaptive
discretisation, and this document is where it is CALIBRATED and GATED against the
known-good eMoScat decks — the tool is only trustworthy once it reproduces the decks it is
meant to replace, and flags the specific failure that motivated it.

## The hybrid approach

A physics **prior** lays out a good grid a-priori (no eigensolve), then cheap **convergence
probes** validate it:

1. **Equidistribution mesh + h/p sweep** (`qscat.tuning.mesh`). Place real-region element
   boundaries so each element carries a constant de-Broglie phase
   $\int k(x)\,\mathrm{d}x \approx C\,(\mathrm{order} - 1)$ with $k(x) = \sqrt{2\,\mathrm{mass}\,\max(E_\mathrm{max} - V(x),\, 0)}$ — a GLL-DVR element
   of order `q` resolves ~`(q−1)/2` wavelengths, so this keeps every element "the same number
   of oscillations wide." In classically forbidden stretches, element length is instead
   capped by the local decay length `1/kappa`; elements adjacent to a turning point or
   singularity are halved. `optimal_real_mesh` sweeps `order ∈ {6, 8, 10, 14}` and keeps
   the `(mesh, order)` combination with the fewest total DVR points — the h/p optimum.
2. **ECS tail — a separate, exp-growth regime** (`qscat.tuning.ecs`). The rotated outgoing
   wave decays as $\exp[-K(x - R_0)\sin\theta]$, so the tail uses growing elements, not the
   oscillatory equidistribution mesh. `max_stable_angle` scans θ up to the ~35° double-ECS
   cap (the binding limit is the 2-D corner where both coordinates sit on their complex
   tails at once, not either potential alone); `tune_ecs_tail` sizes exponentially-growing
   elements to absorb the fastest channel wavenumber `K` down to `~1e-12`.
3. **`propose_grid`** (`qscat.tuning.propose`) wires a per-coordinate MODEL ADAPTER (which
   picks `V`/mass/extent/channel-`k` from a `ResonanceModel`) through analyze → mesh → ECS
   into one `FemDvrEcsGrid` — the one-shot a-priori half of the tuner.
4. **The decoupled 1-D probes** (`qscat.tuning.probes`) validate it empirically:
   - `probe_nuclear` — vibrational eigenvalues stable under one h-refinement.
   - `probe_electronic` — the anion bound electronic-state energy stable under refinement
     (a cheap proxy for the full two-angle resonance-pole match).
   - `probe_channel_representation` — THE cheapest and most diagnostic: no eigensolve at
     all. Compares the grid's own quadrature estimate of $\int |F|^2\,\mathrm{d}r$ (a channel function at
     wavenumber `k`, partial wave `l`) against a fine-uniform-grid reference over the same
     span. A grid whose elements are large compared to `1/k` aliases this badly — exactly
     the failure mode that cost F₂ its ~36 orders.
5. **`grid_cost`/`tensor_cost`** (`qscat.tuning.metrics`) give exact DVR point counts and
   rough, anchored 2-D sparse-LU cost estimates, for RELATIVE ranking of candidates.

The `discretisation-tuner` skill runs analyze → propose → probe at the energy extremes →
refine/coarsen → a final 2-D spot-check, as a supervised loop (see the skill for the worked
example); `qscat.tuning` itself is pure, deterministic primitives with no judgment baked in.

## Task 8: calibrating the phase constant `C`

The equidistribution mesh has exactly one free numerical knob, the phase-per-element
constant $C$ in $\mathrm{phase\ per\ element} = C\,(\mathrm{order} - 1)$ — smaller $C$ means a finer mesh.
`validation/tuning/calibrate.py` (`uv run python -m validation.tuning.calibrate`) measures
it by sweeping `C` and checking whether `propose_grid`'s nuclear grid for N₂/NO/F₂/H₂⁺
(proxy) reproduces-or-beats the corresponding committed eMoScat/proxy deck.

**F₂ is the decisive case.** F₂ has a genuinely OPEN dissociative-attachment channel within
its tested range (`(0.01, 0.05)` Ha): its anion bound electronic state
(`anion_electronic_states`, at the eMoScat deck's dissociation limit `R_inf = 10.7` bohr) sits
at $\varepsilon_e \approx -0.127$ Ha, an exothermic threshold, so
$E_\mathrm{DR} = E_\mathrm{max} - \varepsilon_e > E_\mathrm{max}$ and the outgoing
nuclear wavenumber $K = \sqrt{2\mu E_\mathrm{DR}} \approx 78$ at $E_\mathrm{max} = 0.05$ — the same wave
whose under-resolution on the coarse shared grid cost ~36 orders of magnitude (see
`docs/physics/diatomic-ve-cross-sections.md`). The calibrated `C` is the SMALLEST value at
which `propose_grid`'s F₂ nuclear mesh represents that wave to `rtol = 1e-3`
(`probe_channel_representation`) using FEWER points than the eMoScat F₂ DA deck.

A 40-candidate sweep, `C ∈ [0.05, 2.0]` step `0.05` (each candidate: `propose_grid` +
`probe_channel_representation` + `probe_nuclear`, all four molecules — a few minutes total),
found:

| Quantity | value |
|---|---|
| Calibrated `C` | **0.10** |
| F₂ eps_e (anion bound state) | −0.12694 Ha |
| F₂ K_DA at E_max=0.05 | 78.28 |

**Per-molecule result at $C = 0.10$:**

| Molecule | proposed n | deck n | ratio | channel rel_error | deck's own rel_error | vib converged |
|---|---|---|---|---|---|---|
| F₂ | 609 | 974 | **0.63×** | 2.89e-4 (**converged**) | 2.12e-4 | yes |
| N₂ | 614 | 428 | 1.44× | 1.66e-3 (not conv., but **beats deck 17×**) | 2.88e-2 | yes |
| NO | 604 | 597 | 1.01× | 5.32e-3 (not conv., but **beats deck 7×**) | 3.74e-2 | yes |
| H₂⁺ (proxy) | 584 | 510 | 1.15× | 7.55e-5 (**converged**) | 4.10e-4 | yes |

F₂ reproduces-and-beats the deck outright: fewer points (609 vs 974) AND absolute
convergence on the exact wave that used to fail catastrophically. H₂⁺'s proxy nuclear
deck is likewise a clean reproduce-and-beat (see "H₂⁺" below).

## Genuine finding: N₂/NO's floor-K channel check isn't a clean absolute bar

N₂ and NO do NOT have an open DA channel in their tested (VE-scale) energy ranges — N₂'s
is closed within the whole +0.5 Ha window; NO's opens at ~0.17 Ha, above the tested
`(0.004, 0.12)` range. Lacking an `eps_e` threshold for a closed channel, their
channel-representation check uses the conservative FLOOR $K = \sqrt{2\mu E_\mathrm{max}}$ (treating
the entire incident electron energy as if converted to nuclear translational energy — a
generous over-estimate, never actually reached since the true threshold is far more
negative). The sweep shows this floor is not met at `rtol = 1e-3` by ANY sane $C$ —
**not even by the eMoScat decks themselves** (N₂'s own deck: rel_error ≈ 0.029; NO's own
deck: rel_error ≈ 0.037, both `≫ rtol`). This is not a calibration failure: it means the
floor is a deliberately conservative bound these decks were never tuned to resolve, and an
absolute `rtol` gate there would fail on a bar the reference implementation itself never
cleared. `test_emoscat_decks.py` gates N₂/NO's channel-representation COMPARATIVELY instead
— rel_error no worse than the deck's own — exactly the design spec's stated criterion
("same-or-better probe precision... than the committed hand-tuned deck"), and both pass
comfortably (17× and 7× better than their decks respectively). Their REAL requirement, the
vibrational spectrum (`probe_nuclear`), converges cleanly at every `C` tried.

**H₂⁺ (proxy nuclear deck) is reported alongside N₂/NO, and is a clean case.** Like N₂/NO,
H₂⁺'s DR channel wavenumber uses the $\sqrt{2\mu E_\mathrm{max}}$ floor rather than a pinned $\varepsilon_e$
(its exit channel is a Rydberg SERIES, not one bound state, so pinning a single threshold
the way F₂'s DA `eps_e` is pinned is awkward). Unlike N₂/NO, though, H₂⁺ converges
absolutely: its much lighter reduced mass (918 vs 13000–17000 for N₂/NO/F₂) keeps even the
floor modest ($K \approx 9.6$ at $E_\mathrm{max} = 0.05$), and its proxy deck's real-region extent (14.0 bohr)
sits much closer to the fixed `_NUCLEAR_X_MAX_DEFAULT = 18.0` bohr than N₂'s (12.0) or NO's
(9.0) — so its point-count ratio (1.15×) fits the standard 1.3× margin without widening.
Both the proxy deck and the proposed grid represent the $K \approx 9.6$ wave to `rtol = 1e-3`.

**Second finding: N₂/NO's point counts exceed their decks' (1.0–1.5×).** Traced to
`qscat.tuning.propose`'s fixed `_NUCLEAR_X_MAX_DEFAULT = 18.0` bohr real-region default,
which is LARGER than N₂'s (12.0 bohr) and NO's (9.0 bohr) committed nuclear real regions
(though smaller than F₂'s 10.7-bohr real region plus its complex tail) — a
per-molecule-INDEPENDENT constant from Task 5's a-priori adapter, not derived from the
potential profile itself (e.g. "where has the interaction died"). More real-region span at
comparable density costs more points, independent of `C`. This is a Task-5 a-priori-adapter
limitation that Task 8's $C$-calibration cannot fix ($C$ controls density, not extent) — a
documented follow-on (deriving `x_max` from the potential rather than a fixed constant),
not addressed here.

## Genuine finding #3: RESOLVED (resonance-aware mesh)

The design spec's "final 2-D spot-check" (the ONE full observable solve confirming the
tensor-product grid delivers the claimed precision) is `test_f2_2d_da_cross_section_
spot_check` (`@pytest.mark.slow`, ~2.5 min on SuperLU, run on F₂ — the molecule that
reproduces-and-beats on the 1-D probes). It is NOT a rubber stamp: `propose_grid`'s F₂
nuclear grid (609 points — the SAME grid that passes both 1-D probes and "beats" the deck's
974 points) gives $\sigma_\mathrm{DA}(E=0.03) = 0.308$ bohr² — but ONE h-refinement of that same
nuclear grid (1189 points) gives `1.644` bohr², and a SECOND refinement (2369 points) gives
`1.658` bohr² (agreeing with the first refinement to 0.85%, and with the eMoScat deck's own
reference value, 1.66 bohr², to ~0.7%). Refining the ELECTRONIC grid instead (nuclear held
at the base 609) changes nothing (`0.30842 → 0.30842`) — isolating the gap squarely to the
nuclear grid's resolution.

**The 1-D probes (channel-representation + vibrational) are NECESSARY but NOT SUFFICIENT**
for this observable — that was the standing gap. The cause: eMoScat's own F₂ deck
hand-places extra-fine sub-0.1-bohr elements specifically around R≈2.5–2.7 bohr — a narrow
feature in the ELECTRON-NUCLEAR INTERACTION (`v_int`/`lambda(R)`, concretely the
adiabatic resonance curve `V_d(R)`), not in `v0` alone. The plain a-priori equidistribution
mesh is built purely from `v0`'s classical $k(x)$ profile (`_nuclear_adapter`/
`analyze_potential`), so it has no way to see a feature that lives only in the coupling
term.

### The fix: the resonance-aware nuclear mesh (`channel="dissociation"`)

`propose_grid(model, "nuclear", energy_range, channel="dissociation")` (nuclear-only;
`channel="ve"` remains the byte-identical default for every other path) makes the mesh
aware of exactly the structure the plain `v0`-only pass cannot see, via
`qscat.tuning.resonance.resonance_curve_arrays` (the two-angle ECS pole match already used by
`qscat.ecs.find_resonance_pole`, sampled densely inside `interaction_region` and once at
the asymptote):

1. **Exit-wave DVR order.** $K_\mathrm{exit} = \sqrt{2\mu\,\max(E_\mathrm{max} - \operatorname{Re} V_d^\mathrm{asym},\; E_\mathrm{max})}$ — the
   fast outgoing dissociation wavenumber the ECS tail must absorb — sizes the quadrature
   order directly via `order_for_wavenumber(K_exit, min_len, target_ppw=6)`, instead of
   letting the h/p sweep pick an order from `v0` alone.
2. **Crossing super-refinement.** `R* =` the outermost sign change of
   `Re(V_d(R)) − v0(R)` (`_outermost_crossing`) locates the resonance crossing (F₂:
   R*≈2.598, matching the eMoScat deck's own hand-placed fine region); a Γ-closing-width
   half-window around it (clamped to [0.15, 0.18] bohr) is then LOCALLY super-refined to
   ~0.03-bohr elements via `refine_elements_in_window` — overriding the global `min_len`
   only inside that window, which is the point (a floored global `min_len` is exactly why
   the earlier worst-case-`k`-merge design was inert).
3. **Trimmed real extent.** The resonant path uses a reduced real-region default (10.5
   bohr, vs the VE path's 18.0) — the ECS tail absorbs the outgoing wave, so the real
   region need only host the interaction region plus a few exit-wave wavelengths, close to
   the eMoScat F₂ deck's own ~10.7-bohr real extent.

**Verified numbers (controller-measured, 2026-07-28):**

| Molecule | Resonant grid | Deck | Ratio | 2-D observable |
|---|---|---|---|---|
| F₂ | 1000 pts, order 14 | 974 pts (eMoScat DA deck), order 14 | 1.027× (deck-parity) | σ_DA(E=0.03) = 1.6562 bohr², CONVERGED (deck 1.66, finding-#3 refine² 1.658) |
| H₂⁺ | 489 pts, order 8 | 510 pts (proxy deck), order 8 | 0.959× (~4% smaller) | not laptop-verifiable (full 2-D DR ~1.15M unknowns — Docker/MUMPS) |

(F₂'s old `channel="ve"` grid was smaller still — 609 points — but gave σ_DA≈0.31, the ~5×
gap this whole finding is about.) The gate is
`validation/tuning/test_resonance_aware.py`: `@pytest.mark.slow` SIZE tests assert
`F2_resonant.n <= 1.05·F2_deck.n` and `H2P_resonant.n <= H2P_proxy_deck.n` (both pass, plus
an order-floor sanity check), and an `@pytest.mark.slow` CONVERGENCE test reruns the 2-D DA
spot-check harness with `channel="dissociation"` and asserts
$|\sigma_\mathrm{base} - \sigma_\mathrm{refined}| / \sigma_\mathrm{refined} < 0.15$ with
$\sigma_\mathrm{base} > 1.0$ bohr².

**The honest finding (stated plainly, not spun): the "10–20% smaller than the hand deck"
expectation this sub-project set out with does NOT hold for F₂.** eMoScat's F₂ DA deck is
a near-optimal expert hand-tuning — it already hand-places exactly the fine crossing region
the resonance-aware mesh now finds automatically. The OLD `propose_grid` (`channel="ve"`)
was smaller than that deck (609 vs 974) ONLY because it was under-converged; that
under-convergence WAS finding #3. Reaching 2-D convergence costs approximately deck-sized
resolution, and the resonance-aware tuner reaches it AUTOMATICALLY — at deck-parity for F₂
(1.027×) and a few percent under for H₂⁺ (0.959×). The deliverable this closes is
**convergence + automation at deck-competitive size, not a point-count reduction.**

**`refine_to_2d_convergence` remains the general, model-agnostic fallback** (`qscat.tuning.
refine2d`, the skill's step 6) for any observable/coordinate combination the resonance-aware
adapter doesn't cover (a non-adiabatic channel, a structure in some other coupling term,
electronic-side sub-structure): it iteratively refines whichever of a caller-supplied
`(g_r, g_R)` pair moves the observable more, closing over ANY scalar cross-section, until
both relative moves fall under `rtol` or `max_iter` is hit. The resonance-aware nuclear path
above is the SPECIFIC, cheap, a-priori fix for the F₂/H₂⁺-style adiabatic-resonance case;
`refine_to_2d_convergence` is what to reach for when a new model's 2-D spot-check finds a
gap the a-priori adapters don't already know how to close.

**Practical consequence:** `propose_grid(..., channel="dissociation")` is now the
recommended nuclear grid for any resonant/dissociative-channel observable (DA, DR) — it
reaches 2-D convergence on the FIRST a-priori pass, without a probe/refine loop, at
deck-competitive size. The plain `channel="ve"` path is unchanged and remains correct for
non-resonant (VE-only) nuclear grids.

## The gate (`validation/tuning/test_emoscat_decks.py`, `test_resonance_aware.py`)

Three kinds of test:

1. **Reproduce-or-beat** (fast, probes only) — `propose_grid` at the calibrated `C` for
   N₂/NO/F₂/H₂⁺: F₂ and H₂⁺ gated strictly (absolute channel convergence + vibrational
   convergence + points ≤ 1.3× deck); N₂/NO gated on vibrational convergence (absolute) and
   channel-representation (comparative to their own deck, per finding #1 above), with N₂'s
   point-count margin widened to the measured ~1.44× ratio (documented, not silently
   loosened to 1.3×).
2. **Flag-the-failures** (fast) — the cheapest probe, `probe_channel_representation`, on the
   COARSE shared N₂-style grid (`qscat.core.grids.nuclear_grid()`), must report
   `converged=False` for F₂ DA's K≈58 wave (and across the whole K≈52–78 range the tested
   energies produce) — the regression guard for the exact bug that motivated this
   sub-project.
3. **The 2-D spot-check** (`@pytest.mark.slow`, F₂ only) — see finding #3 above: confirms
   the refined-grid FAMILY converges, and records the base-grid gap as an honest, actionable
   finding rather than asserting a false match.
4. **The resonance-aware re-tune** (`validation/tuning/test_resonance_aware.py`,
   `@pytest.mark.slow` throughout — each grid build pays a ~60–90s resonance scan): F₂'s and
   H₂⁺'s `channel="dissociation"` nuclear grids are SIZE-gated against their decks
   (deck-parity / no-larger-than-proxy) and, for F₂, CONVERGENCE-gated against a once-refined
   solve — see finding #3's resolution above.

An H₂⁺ Coulomb-incident coarse-grid FLAG check was explored but not included as a gate: at
H₂⁺ DR's low incident `k` (long de-Broglie wavelength), `probe_channel_representation`'s
failure mode on a truncated Coulomb electronic grid is dominated by real-region EXTENT, not
element density, and neither a 30-bohr nor the 60-bohr proxy grid gave a clean converged/
not-converged split — not a clean regression gate the way the F₂ DA case is.

See also: `docs/superpowers/specs/2026-07-28-discretisation-tuner-design.md` (the design),
the `discretisation-tuner` skill (the supervised loop + worked example),
`docs/physics/diatomic-ve-cross-sections.md` (the F₂ DA physics and the original K≈58 bug).
