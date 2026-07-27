# `qscat.core` + `qscat.model`: the electron–diatomic VE-scattering engine

**Location:** `qscat.core` (model-independent solvers) + `qscat.model` (the model form,
protocol, and per-molecule parameters) + energy-normalized radial functions in
`qscat.special`. **Origin:** promoted from the N₂ sub-projects (#6 exact 2-D TI, #7 2-D TD,
plus the order-3 Padé + elastic free-reference fixes) once validated to ~1–2 % TD-vs-TI on N₂
— see `docs/superpowers/specs/2026-07-27-diatomic-ve-scattering-library-design.md`. **Units:**
atomic units throughout.

## What this is

The exact 2-D electron–diatomic vibrational-excitation (VE) cross section, by two
interchangeable routes on the same FEM-DVR-ECS grid:

- **Time-independent (`qscat.core.driven.ve_cross_section`)** — the exact driven
  Lippmann-Schwinger solve `(E_tot·I − H)Ψ_sc = V_int·Ψ_i`, T-matrix projection, σ = π|S−δ|²/2E.
  One sparse factorization per collision energy, reused across channels; the energy sweep
  analyzes once and `SparseLU.refactor`s per energy (constant sparsity pattern).
- **Time-dependent (`qscat.core.time_dependent.td_ve_cross_section`)** — an incident Gaussian
  wavepacket propagated under `H` with the **order-3 diagonal-Padé** operator
  (`qscat.evolution.make_pade_stepper`, `O(dt⁷)`/step — order-1 Crank-Nicolson under-converges
  ~100 % over a multi-thousand-step run), then the Tannor-Weeks energy transform of the stored
  correlation `c_{v'}(t)`. The elastic (diagonal) channel subtracts the free-particle reference
  `S_free(E)` from a `V_int=0` propagation (not a literal 1 — the outgoing normalization makes
  `S_free ≈ 2π² ≠ 1`). Matches the TI oracle to ~1–2 % across the resonance for all channels
  (see `docs/physics/n2-2d-td-cross-section.md`).

Supporting modules: `qscat.core.channels` (`channel_vector` + threshold logic),
`qscat.core.grids` (`electronic_grid`/`nuclear_grid` — the parameterized FEM-DVR-ECS layout),
`qscat.core.vibrational` (`vibrational_states` — neutral bound states, `v0` passed in),
`qscat.core.wavepacket`, `qscat.core.correlation`, `qscat.core.plot` (`plot_cross_sections`).

## The `core` / `model` split

`qscat` already stands for *quantum scattering*, so the engine does not hide under a redundant
`qscat.scattering`. Instead:

- **`qscat.core`** — model-independent. Its solvers take a `model` object and know nothing
  about which molecule is being solved. They build on the general primitives
  (`qscat.{units, linalg, dvr, ecs, special, evolution}`).
- **`qscat.model`** — everything tied to a specific model, gathered in one place:
  - `ResonanceModel` — a `@runtime_checkable Protocol` (`mu`, `ell`, `v0`, `lam`, `v_int`,
    `surface`, `hamiltonian(tgrid)`, `interaction_diag(tgrid)`). **This is the entire contract
    `qscat.core` depends on.**
  - `DiatomicResonanceModel` — the shared Morse-`v0` + sigmoid-`λ(R)` + Gaussian-in-r-`V_int`
    form (N₂/NO/F₂ differ only in parameters).
  - `N2`, `NO`, `F2` — the per-molecule registry instances (from the eMoScat decks).

**Hard boundary — `qscat.core` never imports `qscat.model` (nor `projects.*`) at runtime.**
`core` type-annotates against the `ResonanceModel` protocol under `TYPE_CHECKING` only; it is
enforced by `libs/qscat/tests/test_core_no_model_import.py` (a fresh-interpreter subprocess
asserts neither `qscat.model` nor any `projects.*` lands in `sys.modules` after
`import qscat.core`). The payoff: adding a molecule — or a wholly different model form, e.g. the
parked angular coupled-channel model — is a `qscat.model` + validation change that never touches
the solvers.

## Usage

```python
from qscat.model import N2
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.core.vibrational import vibrational_states
from qscat.core.driven import ve_cross_section
from qscat.core.time_dependent import td_ve_cross_section
from qscat.dvr import TensorGrid

tgrid = TensorGrid([electronic_grid(r_max=16.0, order=8, n_complex=6),
                    nuclear_grid(r_max=22.0, quadrature=10, n_complex=5)])
eps, chi = vibrational_states(tgrid.grids[1], N2.mu, 4, N2.v0)
sigma_ti = ve_cross_section(tgrid, N2, eps, chi, 0, [0, 1, 2], E_grid)          # exact TI
sigma_td = td_ve_cross_section(tgrid, N2, eps, chi, 0, [0, 1, 2], E_grid,        # TD (order-3 Padé)
                               dt=1.0, n_steps=1500, wp_in=..., wp_out=...)
```

Swap `N2` for `NO`/`F2` (or any `ResonanceModel`) — nothing else changes. The N₂ research
projects (`projects/n2_*`) consume this engine (their old modules are thin shims binding `N2`);
`qscat.model.N2` is the single source of truth for the N₂ model.
