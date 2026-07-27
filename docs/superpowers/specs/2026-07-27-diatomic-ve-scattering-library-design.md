# Diatomic VE-scattering library (`qscat.core` + `qscat.model`) + NO/F2 — Design Spec

**Date:** 2026-07-27
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — implementation pending
**Lifecycle:** `qm-method-lifecycle` **stage 5 (promote to `qscat`)** for the N2-validated
2-D scattering machinery, then stages 2–4 (toy → validate) for NO and F2 as thin
instantiations. Builds on sub-projects #6 (exact 2-D TI), #7 (2-D TD), #8 (MUMPS), #9
(sweep reuse), and the order-3 Padé + elastic free-reference fixes.

## Context

The exact 2-D electron–diatomic vibrational-excitation (VE) solver has been validated on N₂
to the point where the time-dependent (TD) route reproduces the time-independent (TI) oracle
— and Houfek's independent data — to ~1–2 % across the resonance. Surveying the eMoScat
decks (`reference/eMoScat/input/{N2,NO,F2}/2D_model.txt`) shows **N₂, NO and F₂ are the same
model and the same method — only parameters differ**:

| | μ (a.u.) | l | D₀ | α₀ | R₀ | α_c | (λ params: λ∞, λ₁, R_λ, λ_c, R_c) |
|---|---|---|---|---|---|---|---|
| N₂ | 12766.36 | 2 | 0.75102 | 1.15350 | 2.01943 | 0.40 | 6.21066, 1.05708, −27.9833, 5.38022, 2.40500 |
| NO | 13614.16 | 1 | 0.2363 | 1.5710 | 2.1570 | 1.00 | 6.3670, 5.0000, 2.0843, 6.0500, 2.2850 |
| F₂ | 17315.99 | 1 | 0.05980 | 1.51610 | 2.69060 | 3.00 | 18.8490, 3.21300, 1.8320, 18.1450, 2.5950 |

The potential *form* — Morse `v0(R)`, sigmoid `λ(R)`, Gaussian-in-r interaction
`V_int(r,R) = −λ(R) e^{−α_c r²}` — and the incident/outgoing wavepacket config are identical
across molecules. So the entire solver stack is molecule-agnostic; adding a molecule should
be **data + validation, not new solver code.** This spec promotes the validated machinery
into `qscat` and instantiates NO and F₂ on top of it. (The ion H₂⁺ is deferred, per Martin.)

## Guiding principle: `qscat.core` vs `qscat.model`

`qscat` already stands for *quantum scattering*, so the machinery does not live under a
redundant `qscat.scattering`. Instead the new work splits along the axis that matters:

- **`qscat.core` — the model-independent engine.** Everything that cannot be avoided and does
  not depend on which molecule is being solved: the TI driven solver, the TD Padé/Tannor-Weeks
  solver, channels, grids, vibrational states, wavepackets, correlation, plotting. It builds
  on the existing general primitives (`qscat.{units, linalg, dvr, ecs, special, evolution}`,
  which stay flat — they are already general and widely imported; renaming them into
  `qscat.core.*` would be pure churn).
- **`qscat.model` — everything tied to a specific model, gathered in one place.** The shared
  resonance-model *form*, the **model protocol** the `core` solvers consume, and the
  per-molecule **parameter data**. This is as far as generalization goes: a different
  functional form (or the future angular coupled-channel model) is a new `qscat.model` type
  implementing the same protocol.

**Hard boundary: `qscat.core` never imports `qscat.model`.** `core` depends only on the model
*protocol* (a structural type). Adding a molecule — or a wholly different model form — touches
`qscat.model` + validation only, never the solvers. This is the extensibility guarantee for
"a lot more later."

## Deliverable 1 — the generalized library

### `qscat.core` (new) — the model-independent VE-scattering engine

| module | promoted from | public API (molecule-agnostic) |
|---|---|---|
| `channels` | `n2_2d_cross_section/channels.py` | `channel_vector(tgrid, k, chi_v, l)`; open/closed-channel threshold helpers |
| `grids` | `electronic_grid.py`, `nuclear_grid.py` | `electronic_grid(r_max, order, n_complex, …)`, `nuclear_grid(r_max, quadrature, n_complex, …)` — the FEM-DVR-ECS element layout, parameterized (extents are config, not baked in) |
| `vibrational` | `n2_ti_cross_section/vibrational.py` | `vibrational_states(grid_R, mu, n_states, v0)` — neutral bound states on the nuclear grid |
| `driven` | `cross_section_2d.py` | `ve_cross_section(tgrid, model, eps, chi, v_init, vprimes, E, *, ordering, return_wavefunction)` — exact TI driven L-S; analyze-once/`SparseLU.refactor` sweep; σ = π\|S−δ\|²/2E |
| `wavepacket` | `td/wavepacket.py` | `gaussian_coeffs`, `initial_state`, `outgoing_channel` |
| `correlation` | `td/correlation.py` | `eta_incident`, `eta_outgoing` (Tannor-Weeks deconvolution factors) |
| `time_dependent` | `td/td_propagation.py` + `td_cross_section.py` | `td_ve_cross_section(tgrid, model, eps, chi, v_init, vprimes, E, *, dt, n_steps, wp_in, wp_out, order=3, subtract_free_reference=True)`, `sigma_from_correlations`, `propagate` |
| `plot` | `cross_section_plot.py` | `plot_cross_sections(E_grid, sigma, *, channels, reference, thresholds, …)` (already generic) |

The solvers take a **`model`** (a `qscat.model` protocol object), not a bare `H`: both the
Hamiltonian (`model.hamiltonian(tgrid)`) and the interaction diagonal
(`model.interaction_diag(tgrid)`) come from one object — the free-reference TD path needs
`H − diag(V_int)`, which the model supplies cleanly. `core` type-annotates against the
protocol, never the concrete `DiatomicResonanceModel`.

### `qscat.special` additions

The energy-normalized free radial functions (pure special functions, model-independent;
currently in `n2_2d_cross_section/channels.py`):
- `riccati_bessel_en(r, k, l)` = `sqrt(2k/π)·r·j_l(kr)` (regular, δ(E−E′)-normalized).
- `riccati_hankel_en(r, k, l)` = `sqrt(2k/π)·r·h_l^{(1)}(kr)` (outgoing; the TD `F_out` uses
  its half).

### `qscat.model` (new) — the model layer

- **`ResonanceModel` protocol** (structural type): `mu: float`, `ell: int`, `v0(R)`,
  `lam(R)`, `v_int(r, R)`, `surface(r, R)`, `hamiltonian(tgrid)`, `interaction_diag(tgrid)`.
  This is the entire contract `qscat.core` depends on.
- **`DiatomicResonanceModel`** (frozen dataclass) — the shared *form* implementing the
  protocol: Morse `v0`, sigmoid `λ`, Gaussian-in-r `V_int`; `hamiltonian`/`interaction_diag`
  are thin wrappers over `qscat.dvr.hamiltonian_nd` / `potential_nd`. Carries no grid state.
- **Parameter registry** (`qscat.model.library` or similar): the N₂/NO/F₂ `DiatomicResonanceModel`
  instances built from the table above, plus each molecule's default grid + wavepacket config
  (electronic/nuclear extents, dt/n_steps/order, wp_in/wp_out). One entry per molecule; the
  place a new molecule is added.

## Deliverable 2 — NO and F₂ (data + validation)

With `qscat.core` + `qscat.model` in place, each molecule is a **registry entry** (already in
`qscat.model`) plus a **validation module**:

```
validation/diatomic/
  <mol>/             per-molecule validation: exact-2D sigma(E) as the ORACLE; TD-vs-TI (and,
                     where the LCP is built, LCP-vs-TI) self-consistency at anchor energies;
                     the committed sigma(E) figure (elastic + first excitations)
```

No new solver code, no per-molecule `projects/` solver dirs. NO first (l=1, α_c=1.0 — nearest
N₂), then F₂ (l=1, α_c=3.0, weakly bound — the near-threshold stress case).

## Interface (key signatures)

```python
# qscat.model
class ResonanceModel(Protocol):
    mu: float
    ell: int
    def v0(self, R): ...
    def v_int(self, r, R): ...
    def surface(self, r, R): ...
    def hamiltonian(self, tgrid): ...        # H (sparse)
    def interaction_diag(self, tgrid): ...   # diag(V_int), for the TD free-reference

@dataclass(frozen=True)
class DiatomicResonanceModel:                # implements ResonanceModel
    mu: float; ell: int
    D0: float; alpha0: float; R0: float
    lambda_inf: float; lambda_1: float; R_lambda: float; lambda_c: float; R_c: float
    alpha_c: float
    # v0 / lam / v_int / surface / hamiltonian / interaction_diag ...

N2 = DiatomicResonanceModel(mu=12766.36, ell=2, D0=0.75102, ...)   # registry
NO = DiatomicResonanceModel(mu=13614.16, ell=1, ...)
F2 = DiatomicResonanceModel(mu=17315.99, ell=1, ...)

# qscat.core.driven
def ve_cross_section(tgrid, model, eps, chi, v_init, vprimes, E, *,
                     ordering="COLAMD", return_wavefunction=False): ...
# qscat.core.time_dependent
def td_ve_cross_section(tgrid, model, eps, chi, v_init, vprimes, E, *,
                        dt, n_steps, wp_in, wp_out, order=3,
                        subtract_free_reference=True): ...
```

## Validation

- **Regression net (Deliverable 1 is behavior-preserving):** the entire N₂ validation must
  stay bit-identical through the promotion — exact-2D vs Houfek (harness group E, `GATED_RTOL
  =1e-3`), TD-vs-TI (`test_td_cross_section.py` `@slow`, order-3 tolerances), the N₂ harness
  **23 PASS / 0 PENDING / 6 NOTE / 0 FAIL**. The N₂ projects are refactored to consume
  `qscat.core` + the `qscat.model` N₂ instance; no physics number changes.
- **NO/F₂ oracle:** no independent golden data ships for NO/F₂ (only N₂ has Houfek's
  `CSVE.V00.J00`). The **exact-2D TI solver is the oracle**; TD-vs-TI (target ~1–2 % as for
  N₂, order-3 Padé) and, where the LCP is built, LCP-vs-TI are the gates — the same
  self-consistency structure as N₂'s E1/F1 groups, minus the independent-data gate unique to
  N₂. An eMoScat cross-check for NO/F₂ is possible but deferred (eMoScat is read-only).
- Free-particle / first-Born / S-matrix reciprocity–unitarity limits (from #6) are
  molecule-agnostic and rerun per molecule as cheap sanity gates.

## Sub-project decomposition (execution order)

- **A — Promote & generalize (the crux).** Create `qscat.core` + `qscat.model` +
  `qscat.special` additions; move the machinery; refactor the N₂ projects to consume them;
  regression-gate on the full N₂ validation. Largest, highest-risk step — done first and
  behavior-preserving. **This spec's plan (Plan A) covers this sub-project.**
- **B — NO.** Registry entry + grid/wp config + validation; deliver the TD/TI σ(E) curves.
- **C — F₂.** Same; watch the near-threshold elastic limit (α_c=3.0, weakly bound).

B and C each get their own (short) spec → plan → execution → merge.

## Out of scope

- **H₂⁺ / ionic channels** (Coulomb tail; deferred per Martin).
- **Independent NO/F₂ golden data** (running eMoScat to generate it; a later cross-check).
- **The angular coupled-channel extension** (parked; `docs/physics/angular-coupled-channels.md`)
  — but the `ResonanceModel` protocol is designed so it can join as a new model type later.
- **Rust kernels / further optimization** (the sparse-LU hot path is #8's line).
- **Dissociative-attachment (DA) channels**, higher partial waves.

## Verification

- `uv run pytest -q -m "not slow"` → all pass, N₂ numbers unchanged (bit-identical where the
  promotion is a pure move).
- `uv run mypy libs/qscat` → 0; `uv run ruff check .` → clean; `qscat.core` importing
  `qscat.model` is forbidden (checked by inspection / an import guard).
- `uv run python -m validation.n2.experiment` → 23/0/6/0, unchanged.
- New NO/F₂ validation modules run and gate TD-vs-TI (and LCP-vs-TI where applicable) against
  the exact-2D oracle; committed σ(E) figures per molecule.
- `docs/physics/` notes for `qscat.core`/`qscat.model` and per-molecule results; `CLAUDE.md`
  updated (new `qscat.core`/`qscat.model` submodules, the `validation/diatomic/` layer).
