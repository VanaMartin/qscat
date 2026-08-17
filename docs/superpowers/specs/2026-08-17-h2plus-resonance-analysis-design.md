# H₂⁺ dissociative recombination — exact resonances against the cross sections

Date: 2026-08-17

## Purpose

Produce, for the H₂⁺-like model, an **exact-resonance-annotated reading of the
DR cross sections**: `σ_DR0(E)` and `σ_DR1(E)` over three electron-energy
windows, with the positions of the *genuine* (non-Born-Oppenheimer) resonances
marked on them, measured against the Born-Oppenheimer quasi-bound levels, and
with the cross-section features that are **not** resonances explained through the
S-matrix.

Three questions, in order of what they establish:

1. **Do the exact resonance poles sit on the cross-section peaks?** The peaks are
   observable; the poles are computed independently. Their coincidence (or not)
   is the check that the resonance picture describes this system.
2. **How far are those poles shifted from the BO levels?** The BO levels are what
   the published treatment of this model marks on these curves; the shift is the
   non-adiabatic correction, and H₂⁺ is where it should be largest — a Rydberg
   electron is slow and its orbital extends far, so the separation of electronic
   and nuclear motion is at its weakest.
3. **What accounts for the features no resonance explains?** The published work
   identifies sharp local minima with no nearby quasi-bound level, and reads them
   as simple zeros of the S-matrix rather than resonance structure.

This is the H₂⁺ instance of the pattern the repository is built around: the exact
2-D solution is the oracle, the approximation (here BO/Rydberg quasi-bound
levels) is under test.

## Reference literature and notation

The primary reference is the thesis this model comes from:

- M. Váňa, *A model of resonant collisions of electrons with molecules and
  molecular ions*, doctoral thesis, Charles University, Prague, 2017.
  <https://dspace.cuni.cz/handle/20.500.11956/92902>

The figures this work is built to answer (printed page numbers as they appear in
the thesis):

| Source | Locator | What it gives |
|---|---|---|
| **Fig. 4.7** | p. 70 | The target figure: DR₀/DR₁ cross sections, TD and TI, over three energy windows, with the quasi-bound levels `ω_i^j` marked as dashed verticals and the ion vibrational thresholds shaded |
| **Fig. 4.3** | p. 64 | The Rydberg electron-energy curves `E_Ryn(R)` whose vibrational levels are the `ω_i^j` |
| **Fig. 4.10** | p. 73 | Zoomed `σ_DR1` minima paired with `Re S_DR1` / `Im S_DR1`, showing the minima are simple zero crossings of both parts |
| **Table 4.2** | p. 73 | the overlap `abs(<ω_i^j, Ψ⁺(R,r;E)>)` at two energies — what assigns a cross-section feature to a quasi-bound state |

**Notation.** `ω_i^j` denotes a vibrational level supported by a Rydberg
electronic curve; one index is the vibrational quantum number and the other the
Rydberg curve. The published panels group them so that consecutive panels are
offset by roughly one ion vibrational quantum (~0.0095 Ha), which fixes the
grouping empirically. **Stage 2 below must reproduce the published positions, and
that reproduction is what settles the index convention** — this spec deliberately
does not assert it in advance.

The three windows requested for this work are `[0, 0.008]`, `[0.010, 0.018]` and
`[0.020, 0.027]` Ha. They are the published panels trimmed to stop short of each
ion vibrational threshold, which is deliberate and worth keeping: the Rydberg
series accumulates at the threshold, and the accumulation region is where a
shift-invert pole search is least trustworthy.

## What already exists

- `qscat.core.dissociation.dr_cross_section(model, …, n_channels=3,
  return_wavefunction=…)` — the exact TI driven solve for a charged target:
  Coulomb incident channel, a loop over Rydberg exit channels, `σ` shaped
  `(len(E), n_channels)`. Returns `Ψ⁺` on request; **does not return the
  S-matrix**.
- `qscat.core.exact_resonance_states` — poles of the full 2-D S-matrix by
  two-angle ECS stability, validated on N₂ (docs/physics/exact-2d-resonances.md).
  Never run on H₂⁺.
- `qscat.core.resonance_levels` / `lcp_resonance_levels` — BO levels in the LCP
  complex curve. Built for the *anion resonance* curve; the H₂⁺ Rydberg curves
  are a different object (bound in the BO picture, no local width).
- `qscat.core.anion_electronic_states(grid, model, R, n_states)` — the neutral's
  bound electronic states at fixed `R`, i.e. the Rydberg series `E_Ryn(R)`.
- `qscat.model.H2P` with `charge = -1` and
  `max_nuclear_ecs_angle_deg = 22.5`.
- Execution: `apps/qscat-run/examples/h2p-dr-ti.yaml` + `docker/run.sh`, with
  `grid: {preset: emoscat}` the full ~1.15 M-unknown deck (electronic → 1300
  bohr, nuclear → 14 bohr) and `H2P:proxy` a laptop-sized reduction. MUMPS only.
- `qscat.core.plot_cross_sections`, `plot_resonance_levels`, `qscat.viz`.

## Design

### 1. The S-matrix element (library)

`dr_cross_section(..., return_smatrix=True)` returns the complex per-channel
amplitude alongside `σ`. The solver already forms this quantity and discards it;
exposing it is what Fig. 4.10's analogue needs, and the pole-fitting cross-check
below depends on it too.

The contract is pinned by a test asserting that `σ` recomputed from the returned
amplitude equals the `σ` the same call returns, so the relation between them
cannot drift. The docstring states the convention explicitly (which normalization
the returned amplitude carries, and how it relates to `σ_DR = C·|T|²/…`).

### 2. Born-Oppenheimer quasi-bound levels `ω_i^j`

For each Rydberg curve `n`, walk `anion_electronic_states` over the nuclear grid
to build `E_Ryn(R)`, then diagonalize the nuclear Hamiltonian `T(µ) + E_Ryn(R)`
in that curve. These levels are **real** in the BO picture — a Rydberg state is
bound in its own curve, and its width comes only from the coupling to the
dissociative channel that BO discards. That is precisely the approximation under
test.

Laptop-feasible (1-D nuclear diagonalizations on a 1-D electronic walk).

### Scope decision: a Rydberg-index cutoff

**States too close to a cation vibrational threshold are out of scope.** They are
omitted because they are not isolated resonances — the level spacing in a Rydberg
series falls as `n^-3`, so neighbouring peaks merge before their widths do, and a
"mark the resonance position on the peak" analysis has nothing to mark. This is a
physics scope decision, not a statement about what the method can reach.

The cut is at **`n_eff <= 12`**, equivalently: drop any level bound by less than
about **3.5 mHa** to the next threshold above it, where
`n_eff = 1/sqrt(2 * binding)`.

The measured cation vibrational thresholds
(`vibrational_states(full_grid().grids[1], H2P.mu, 6, H2P.v0)`) are

    -0.097604  -0.087802  -0.078519  -0.069754  -0.061507  -0.053780  Ha

spaced 9.80, 9.28, 8.76, 8.25, 7.73 mHa — anharmonically decreasing. Levels
located inside the published windows sit at `n_eff` ~ 9-11 (`<r>` ~ 133-176
bohr), so the cutoff keeps them all; the extreme admitted case is `<r>` ~ 216
bohr.

**Consequence for the electronic box:** a converged box wants ~3-4x the state's
extent, so ~650-900 bohr for the admitted levels. That is the range the
DR-window box probe must bracket, and it is what makes the pole campaign
tractable — without the cutoff, states bound by ~1 mHa reach `n_eff` ~ 22 and
`<r>` ~ 745 bohr, demanding a box several times larger for states that could not
be resolved as individual resonances anyway.

The published windows already stop short of each threshold, so this cutoff
formalises a restriction the chosen energy ranges half-imposed.

### 3. Exact resonances

`exact_resonance_states` seeded from the `ω_i^j`. The physically important
question is what electronic box the pole search needs. The cross-section deck
needs 1300 bohr because the **incident Coulomb wave** must be represented; a
Rydberg resonance state is a **closed-channel, localized** object whose extent
goes as `n²` (~120 bohr at `n = 11`). If the poles converge on a far smaller box,
the pole campaign is cheap and possibly laptop-feasible. Stage 3 measures this
rather than assuming it.

### 4. Independent pole extraction from `S(E)`

Fit the computed `S(E)` for poles (Lorentzian near isolated peaks; a Padé
continuation where they crowd). This is an independent determination from the
same physics, and it is what makes "the poles sit on the peaks" a check rather
than a tautology: the exact poles come from an eigenproblem that never sees the
cross section, the fitted poles come from the cross section and never see the
eigenproblem. Where the two disagree, that disagreement is a result.

### 5. Assignment

The Table 4.2 idea with exact states: `|⟨ψ̃_res | Ψ⁺(E)⟩|` at energies on and
off each feature, using `return_wavefunction=True` and the c-product. Large
overlap assigns a feature to a state. This also answers, in passing, the question
of how the resonance state relates to the scattering solution — near a pole,
`Ψ_sc ≈ ψ_res ⟨ψ̃_res|V|Ψ_i⟩/(E − E_res)`, so the overlap should peak sharply at
the pole and fall away from it.

### 6. Figures

- **A — the Fig. 4.7 analogue.** Three stacked panels, log `σ` against electron
  energy, `σ_DR0` and `σ_DR1`, with exact poles as solid verticals, BO `ω_i^j` as
  dashed verticals, thresholds shaded. Requires a `markers` facility on
  `plot_cross_sections` (generic: labelled verticals, no physics).
- **B — the Fig. 4.10 analogue.** For each identified non-resonant minimum: the
  zoomed `σ_DR1` above `Re S` / `Im S` on a linear scale, showing the double zero
  crossing.
- **C — the shift.** Exact pole vs nearest BO level per state, in meV, with the
  width alongside; the H₂⁺ counterpart of the N₂ table.

## Compute strategy

Everything full-deck runs on **sadaharu** (x86_64, Docker, system MUMPS) via
`docker/run.sh` with a qscat-run config derived from `h2p-dr-ti.yaml`. The image
is built natively there; emulated cross-builds from the ARM Mac are unreliable.

**Cost is the campaign's dominant unknown and is measured before anything is
committed to.** A single full-deck energy point at ~1.15 M unknowns has never
been timed in this repository. Stage 0 measures factor time, solve time and peak
RSS, and every downstream sampling decision follows from it. If a full-deck
energy point turns out to cost minutes rather than seconds, the response is to
narrow to the two windows with the richest structure rather than thin the
resolution everywhere — a coarse sweep through a Rydberg series is worse than no
sweep, because it aliases the peaks it is meant to locate.

## Validation

| Claim | How it is checked |
|---|---|
| `return_smatrix` is consistent | `σ` recomputed from the amplitude equals the returned `σ` (unit test) |
| BO levels are right | reproduce the published `ω_i^j` positions marked in Fig. 4.7 |
| Exact poles are converged | box-convergence in the electronic extent (Stage 3), plus the two-angle residuals the solver already reports |
| Poles are physical | coincidence with cross-section peaks, and agreement with the independent `S(E)` fits |
| Assignment is meaningful | overlap peaks at the pole and falls off it |
| Non-resonant minima | both `Re S` and `Im S` cross zero there, with no pole nearby |

## Milestones

0. Cost probe on sadaharu (full-deck energy point; one shift-invert solve).
1. `return_smatrix` + tests.
2. BO `ω_i^j` levels, checked against the published positions.
3. Pole box-convergence probe (electronic 150 / 300 / 600 bohr).
4. Pole campaign across the three windows.
5. `σ(E)` sweep on the full deck.
6. Independent pole extraction from `S(E)`; compare with 4.
7. Assignment overlaps.
8. Figures, physics note, reference-note rows for the ch. 4 sources.

## Risks and limits

- **Rydberg accumulation.** Poles crowd towards each threshold. Shift-invert
  returns a local window, so near-degenerate poles can merge or be missed. The
  trimmed windows help; residual risk is that the highest `n` in each window is
  unreliable, and that must be reported per level rather than hidden.
- **Widths below the noise floor.** These states are long-lived. The BO-levels
  work already records that widths below ~1e-6 Ha are noise; some exact widths
  will land there. Positions remain quotable when widths do not.
- **ECS angle ceiling.** Nuclear θ < 22.5° for this model (enforced), electronic
  < 45°: a pole is only exposed once `2θ` exceeds its argument, so sufficiently
  broad resonances are unreachable in principle. State it rather than discover it.
- **The BO comparison is not like-for-like in one respect**: BO levels here carry
  no width at all, so only positions can be compared for the shift. The exact
  widths have no BO counterpart to be measured against — they are new information,
  not a correction to something.
- **`exact_resonance_states` has never run on a Coulomb model.** Everything about
  Stage 3 and 4 is first-of-kind for this solver; a negative result there is a
  legitimate outcome and would redirect the work to the `S(E)`-fitting route
  alone.

## Out of scope

- The time-dependent DR route (thesis Fig. 4.7 also plots TD curves). This work
  is TI only; TD-DR is `td_da_cross_section`'s ionic generalization and a separate
  project.
- Any change to `dr_cross_section`'s physics — only the S-matrix return is added.
- Improving the H₂⁺ model, grids, or the `emoscat` deck.
- The `n > 3` Rydberg exit channels (the published work uses 3).
