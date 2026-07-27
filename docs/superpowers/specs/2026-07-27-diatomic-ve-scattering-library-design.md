# Diatomic VE-scattering library (`qscat.scattering`) + NO/F2 — Design Spec

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

| | μ (a.u.) | l | D₀ | α₀ | R₀ | α_c | (λ params) |
|---|---|---|---|---|---|---|---|
| N₂ | 12766.36 | 2 | 0.75102 | 1.15350 | 2.01943 | 0.40 | 6.21066, 1.05708, −27.9833, 5.38022, 2.40500 |
| NO | 13614.16 | 1 | 0.2363 | 1.5710 | 2.1570 | 1.00 | 6.3670, 5.0000, 2.0843, 6.0500, 2.2850 |
| F₂ | 17315.99 | 1 | 0.05980 | 1.51610 | 2.69060 | 3.00 | 18.8490, 3.21300, 1.8320, 18.1450, 2.5950 |

The potential *form* — Morse `v0(R)`, sigmoid `λ(R)`, Gaussian-in-r interaction
`V_int(r,R) = −λ(R) e^{−α_c r²}` — and the incident/outgoing wavepacket config (r0=45,
σ=6, p0=−0.35; test fn r0=55, σ=4, p0=0.7) are identical across molecules. So the entire
solver stack is molecule-agnostic; adding a molecule should be **data + validation, not new
solver code.** This spec promotes the validated machinery into `qscat` and then instantiates
NO and F₂ on top of it. (The ion H₂⁺ is deferred, per Martin.)

## Guiding principle: methods → library, models → data

- **Methods (molecule-agnostic) graduate to `qscat`:** grids, Hamiltonian assembly, free
  radial functions, the TI driven Lippmann-Schwinger solver, the TD Padé/Tannor-Weeks
  solver, wavepackets, correlation factors, plotting.
- **Models (the given testbeds) stay as data** — never hardcoded in the library. The shared
  *form* (Morse + sigmoid + Gaussian) is a reusable `qscat.models` builder; each molecule is
  a parameter set. This keeps `qscat` about *how we compute*, and the research layer about
  *what we compute* (the potentials are testbeds for where approximations fail, per the
  research-program framing — not the library's concern except as a reusable form).

## Deliverable 1 — the generalized library

### `qscat.scattering` (new submodule) — electron–diatomic VE scattering

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

The solvers take a **`model`** (below) rather than a bare `H`, so both the Hamiltonian
(`model.hamiltonian(tgrid)`) and the interaction diagonal (`model.interaction_diag(tgrid)` /
`model.v_int`) come from one object — the free-reference TD path needs `H_2D − diag(V_int)`,
which the model supplies cleanly.

### `qscat.special` additions

The energy-normalized free radial functions (currently `n2_2d_cross_section/channels.py`),
pure special functions:
- `riccati_bessel_en(r, k, l)` = `sqrt(2k/π)·r·j_l(kr)` (regular, δ(E−E′)-normalized).
- `riccati_hankel_en(r, k, l)` = `sqrt(2k/π)·r·h_l^{(1)}(kr)` (outgoing; the TD `F_out` uses
  its half).

### `qscat.models` (new) — the shared model form

`DiatomicResonanceModel` (frozen dataclass): fields `mu`, `ell`, and the potential
parameters; methods `v0(R)`, `lam(R)`, `v_int(r, R)`, `surface(r, R)` (the full
`v0 + l(l+1)/2r² + V_int`), and grid-level `hamiltonian(tgrid)` / `interaction_diag(tgrid)`
(thin wrappers over `qscat.dvr.hamiltonian_nd` / `potential_nd`). Constructed from a
parameter set; carries no grid state. Per-molecule parameters live as **data** (a small
registry / JSON), not in the class.

## Deliverable 2 — NO and F₂ (data + validation)

A unified, extensible research layer (chosen for scale — "a lot more later"):

```
projects/diatomic_ve/
  models.py          registry: DiatomicResonanceModel param sets (N2, NO, F2, …) + per-
                     molecule grid/wavepacket config (from the eMoScat decks)
  curves.py          thin driver: TI + TD sigma(E) curves for a named molecule via
                     qscat.scattering (no solver code -- pure orchestration)
validation/diatomic/
  <mol>/             per-molecule validation: exact-2D sigma(E) as the ORACLE; TD-vs-TI and
                     (where built) LCP-vs-TI self-consistency at anchor energies; figure
```

Each new molecule = one registry entry (params + grid config) + one validation module.
NO first (l=1, α_c=1.0 — nearest N₂), then F₂ (l=1, α_c=3.0, weakly bound — the
near-threshold stress case).

## Interface (key signatures)

```python
# qscat.models
@dataclass(frozen=True)
class DiatomicResonanceModel:
    mu: float; ell: int
    D0: float; alpha0: float; R0: float
    lambda_inf: float; lambda_1: float; R_lambda: float; lambda_c: float; R_c: float
    alpha_c: float
    def v0(self, R): ...
    def lam(self, R): ...
    def v_int(self, r, R): ...            # -lam(R) exp(-alpha_c r^2)
    def surface(self, r, R): ...          # v0 + ell(ell+1)/2r^2 + v_int
    def hamiltonian(self, tgrid): ...     # H_2D (sparse)
    def interaction_diag(self, tgrid): ...# diag(V_int) for the TD free-reference

# qscat.scattering.driven
def ve_cross_section(tgrid, model, eps, chi, v_init, vprimes, E, *,
                     ordering="COLAMD", return_wavefunction=False): ...
# qscat.scattering.time_dependent
def td_ve_cross_section(tgrid, model, eps, chi, v_init, vprimes, E, *,
                        dt, n_steps, wp_in, wp_out, order=3,
                        subtract_free_reference=True): ...
```

## Validation

- **Regression net (Deliverable 1 is behavior-preserving):** the entire N₂ validation must
  stay bit-identical through the promotion — exact-2D vs Houfek (harness group E, `GATED_RTOL
  =1e-3`), TD-vs-TI (`test_td_cross_section.py` `@slow`, order-3 tolerances), the N₂ harness
  **23 PASS / 0 PENDING / 6 NOTE / 0 FAIL**. The N₂ project is refactored to consume
  `qscat.scattering` + a `DiatomicResonanceModel` N₂ instance; no physics number changes.
- **NO/F₂ oracle:** no independent golden data ships for NO/F₂ (only N₂ has Houfek's
  `CSVE.V00.J00`). The **exact-2D TI solver is the oracle**; TD-vs-TI (target ~1–2 % as for
  N₂, order-3 Padé) and, where the LCP is built, LCP-vs-TI are the gates — the same
  self-consistency structure as N₂'s E1/F1 groups, minus the independent-data gate unique to
  N₂. An eMoScat cross-check for NO/F₂ is possible but deferred (eMoScat is read-only, not
  built here).
- Free-particle / first-Born / S-matrix reciprocity–unitarity limits (from #6) are
  molecule-agnostic and rerun per molecule as cheap sanity gates.

## Sub-project decomposition (execution order)

- **A — Promote & generalize (the crux).** Move the machinery to `qscat.scattering` /
  `qscat.special` / `qscat.models`; refactor the N₂ projects to consume it; regression-gate
  on the full N₂ validation. Largest, highest-risk step — done first and behavior-preserving.
- **B — NO.** Registry entry + grid config + validation; deliver the TD/TI σ(E) curves.
- **C — F₂.** Same; watch the near-threshold elastic limit (α_c=3.0, weakly bound).

Each is its own spec → plan → subagent-driven execution → review → merge.

## Out of scope

- **H₂⁺ / ionic channels** (Coulomb tail — `sphHankel1En`'s Coulomb branch exists in eMoScat;
  deferred per Martin).
- **Independent NO/F₂ golden data** (running eMoScat to generate it; a later cross-check).
- **The angular coupled-channel extension** (parked; `docs/physics/angular-coupled-channels.md`).
- **Rust kernels / further optimization** (the sparse-LU hot path is #8's line; not this).
- **Dissociative-attachment (DA) channels**, higher partial waves.

## Verification

- `uv run pytest -q -m "not slow"` → all pass, N₂ numbers unchanged (bit-identical where the
  promotion is a pure move).
- `uv run mypy libs/qscat` → 0; `uv run ruff check .` → clean.
- `uv run python -m validation.n2.experiment` → 23/0/6/0, unchanged.
- New NO/F₂ validation modules run and gate TD-vs-TI (and LCP-vs-TI where applicable) against
  the exact-2D oracle; committed σ(E) figures per molecule.
- `docs/physics/` notes for `qscat.scattering` (the promoted method) and per-molecule results;
  `CLAUDE.md` updated (new `qscat.scattering`/`qscat.models` submodules, the `projects/
  diatomic_ve/` layer).
```
