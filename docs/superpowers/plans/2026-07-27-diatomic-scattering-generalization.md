# Sub-project A: promote N2 VE-scattering machinery to `qscat.core` + `qscat.model` — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote the N2-validated exact-2-D electron–diatomic VE-scattering machinery into
the reusable library as `qscat.core` (the model-independent engine) + `qscat.model` (the
model form, protocol, and per-molecule parameters) + `qscat.special` radial functions, then
refactor the N2 projects/validation to consume it — **behavior-preserving, gated bit-identical
by the full N2 validation.**

**Architecture:** `qscat.core` solvers take a `model` object satisfying the `qscat.model.
ResonanceModel` protocol (`mu`, `ell`, `v0`, `v_int`, `surface`, `hamiltonian(tgrid)`,
`interaction_diag(tgrid)`); `core` NEVER imports `model` (depends only on the protocol). Each
promotion task moves validated code, changes the N2-specific signature to the model-taking
one, adapts the old N2 project module into a thin shim binding the N2 model (so the suite
stays green after every task), and gates the change differentially. The final task points the
validation at `qscat.core`/`qscat.model` directly and removes the shims.

**Tech Stack:** Python 3.12, scipy.sparse, `qscat.{dvr,linalg,ecs,special,evolution}`,
`typing.Protocol`, pytest, matplotlib.

**Design spec:** `docs/superpowers/specs/2026-07-27-diatomic-ve-scattering-library-design.md`

## Global Constraints

- **Behavior-preserving.** Every promoted piece must reproduce the current N2 result to
  round-off (differential test vs the pre-move code), and the N2 harness
  (`python -m validation.n2.experiment`) must stay **23 PASS / 0 PENDING / 6 NOTE / 0 FAIL**
  after every task. No physics number changes.
- **`qscat.core` must not import `qscat.model`** — solvers type-annotate against the
  `ResonanceModel` protocol only. Verified by inspection + a guard test.
- Existing general primitives (`qscat.{units, linalg, dvr, ecs, special, evolution}`) stay
  flat and are NOT renamed.
- Core `qscat` stays numpy/scipy-only (matplotlib is already a dev dep; `qscat.core.plot`
  imports it lazily/guarded, matching `observation.py`'s `matplotlib.use("Agg")` pattern).
- `uv run mypy libs/qscat` → **0**; `uv run ruff check .` → clean (line length 100).
- The TD physics fixes stay intact: order-3 Padé default (`make_pade_stepper`), the elastic
  free-reference subtraction, `SparseLU.refactor` energy-sweep reuse.
- After each task the suite stays green (old N2 project modules become thin shims until the
  final rewire); commit per task. Trailer:
  `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`

---

## File Structure

| File | Responsibility |
|---|---|
| `libs/qscat/qscat/special/…` | **Modify.** Add `riccati_bessel_en`, `riccati_hankel_en` |
| `libs/qscat/qscat/model/__init__.py`, `diatomic.py`, `library.py` | **Create.** `ResonanceModel` protocol; `DiatomicResonanceModel`; N2/NO/F2 registry |
| `libs/qscat/qscat/core/__init__.py` | **Create.** re-exports the engine API |
| `libs/qscat/qscat/core/{grids,vibrational,channels,driven,wavepacket,correlation,time_dependent,plot}.py` | **Create.** promoted, model-taking engine |
| `libs/qscat/tests/test_model.py`, `test_core_*.py` | **Create.** differential + protocol tests |
| `projects/n2_2d_cross_section/*`, `projects/n2_2d_td_cross_section/*`, `projects/n2_ti_cross_section/*` | **Modify.** thin shims → then consume `qscat.core`/`qscat.model` directly |
| `validation/n2/*` | **Modify (final task).** import the N2 model from `qscat.model`, call `qscat.core` |
| `docs/physics/qscat-core-scattering.md`, `CLAUDE.md` | **Create/Modify.** the promoted method + layout |

---

### Task 1: `qscat.special` radial functions

**Files:** Modify `libs/qscat/qscat/special/` (add module or extend); create
`libs/qscat/tests/test_special_radial.py`; modify `projects/n2_2d_cross_section/channels.py`
(re-export).

**Interfaces — Produces:** `riccati_bessel_en(r, k, l)` = `sqrt(2k/π)·r·j_l(kr)`;
`riccati_hankel_en(r, k, l)` = `sqrt(2k/π)·r·h_l^{(1)}(kr)` (`h_l^{(1)} = j_l + i y_l`). Both
accept real `r` arrays, `k>0`, integer `l`; return complex128.

- [ ] **Step 1: Failing test.** In `test_special_radial.py`: assert `riccati_bessel_en`
  matches the current `projects.n2_2d_cross_section.channels.riccati_bessel_en` elementwise to
  1e-14 on a sample `(r, k, l)`; assert `riccati_hankel_en` real part equals `riccati_bessel_en`
  (since `Re h^{(1)} = j_l`) and its imag part is `sqrt(2k/π)·r·y_l(kr)`; assert `k<=0` raises.
- [ ] **Step 2: Run → fail** (functions not in `qscat.special`).
- [ ] **Step 3: Implement.** Move `riccati_bessel_en`'s body into `qscat.special`; add
  `riccati_hankel_en` (the outgoing sibling, generalizing `_outgoing_coeffs`'s `riccati_h1`).
  Export from `qscat.special.__init__`.
- [ ] **Step 4: Rewire + gate.** `projects/n2_2d_cross_section/channels.py` re-exports
  `riccati_bessel_en` from `qscat.special` (keep `channel_vector` here for now); `td/correlation.py`'s
  `_outgoing_coeffs` uses `qscat.special.riccati_hankel_en`. Run `pytest libs/qscat projects/n2_2d_cross_section -q -m "not slow"` → pass; `mypy libs/qscat` → 0; ruff clean.
- [ ] **Step 5: Commit.**

---

### Task 2: `qscat.model` — protocol, form, registry

**Files:** Create `libs/qscat/qscat/model/{__init__.py, diatomic.py, library.py}`; create
`libs/qscat/tests/test_model.py`.

**Interfaces — Produces:** `ResonanceModel` (Protocol: `mu`, `ell`, `v0(R)`, `lam(R)`,
`v_int(r,R)`, `surface(r,R)`, `hamiltonian(tgrid)`, `interaction_diag(tgrid)`);
`DiatomicResonanceModel` (frozen dataclass implementing it — Morse/sigmoid/Gaussian, wrapping
`qscat.dvr.hamiltonian_nd`/`potential_nd`); `library.N2`/`NO`/`F2` instances (params from the
spec table).

**Background:** the form is `projects/n2_resonance/potential.py`'s `v0`/`lam`/`v_int` +
`hamiltonian2d.py`'s `build_h2d`/`interaction_diag`/`potential_2d`, parameterized. `mu`, `ell`
come from `hamiltonian2d.MU`/`ELL`. `hamiltonian(tgrid)` = `hamiltonian_nd(tgrid, [1.0, mu],
surface)`; `interaction_diag(tgrid)` = `potential_nd(tgrid, v_int)`.

- [ ] **Step 1: Failing test.** In `test_model.py`: build `library.N2`; assert `N2.v0(R)`,
  `N2.v_int(r,R)`, `N2.surface(r,R)` match the current `projects.n2_resonance.potential` /
  `hamiltonian2d` functions elementwise to 1e-14 on sample grids; assert `N2.mu == 12766.36`,
  `N2.ell == 2`; assert `N2.hamiltonian(small_tgrid)` equals `hamiltonian2d.build_h2d(small_tgrid)`
  to round-off (sparse `abs(A-B).max()`); assert `DiatomicResonanceModel` satisfies
  `isinstance(N2, ResonanceModel)` (runtime-checkable protocol).
- [ ] **Step 2: Run → fail.**
- [ ] **Step 3: Implement** `ResonanceModel` (runtime_checkable Protocol), `DiatomicResonanceModel`
  (fields per the spec; methods computing v0/lam/v_int/surface; `hamiltonian`/`interaction_diag`
  via `qscat.dvr`), and `library.py` with N2/NO/F2 instances (verbatim params from the spec table).
- [ ] **Step 4: Gate.** `pytest libs/qscat/tests/test_model.py -q` → pass; mypy 0; ruff clean.
- [ ] **Step 5: Commit.**

---

### Task 3: `qscat.core.grids` + `qscat.core.vibrational`

**Files:** Create `libs/qscat/qscat/core/{__init__.py, grids.py, vibrational.py}`; create
`libs/qscat/tests/test_core_grids.py`; modify `projects/n2_2d_cross_section/electronic_grid.py`,
`projects/n2_ti_cross_section/nuclear_grid.py`, `projects/n2_ti_cross_section/vibrational.py`
(→ thin shims).

**Interfaces — Produces:** `qscat.core.grids.electronic_grid(r_max, order, n_complex, …)`,
`nuclear_grid(r_max, quadrature, n_complex, …)` (the FEM-DVR-ECS element layout, extents
parameterized); `qscat.core.vibrational.vibrational_states(grid_R, mu, n_states, v0)`.

- [ ] **Step 1: Failing test.** Assert the promoted `electronic_grid`/`nuclear_grid` produce a
  `FemDvrEcsGrid` identical (points, weights, R0) to the current `n2_electronic_grid`/
  `n2_nuclear_grid` at the N2 defaults; assert `vibrational_states(grid_R, MU, 4, N2.v0)`
  reproduces the current `n2_ti_cross_section.vibrational.vibrational_states` `eps`/`chi` to
  round-off.
- [ ] **Step 2: Run → fail.**
- [ ] **Step 3: Implement** by moving the bodies (element layout, kinetic + diag(v0) eigensolve)
  into `qscat.core`; `vibrational_states` takes the neutral potential `v0` as a callable
  (was N2-hardcoded) — pass `N2.v0`.
- [ ] **Step 4: Rewire + gate.** Old `n2_electronic_grid`/`n2_nuclear_grid`/`vibrational_states`
  become thin wrappers calling `qscat.core` (nuclear/vibrational binding `N2.v0`). Run
  `pytest projects/n2_2d_cross_section projects/n2_ti_cross_section libs/qscat -q -m "not slow"`
  → pass; harness 23/0/6/0; mypy 0; ruff clean.
- [ ] **Step 5: Commit.**

---

### Task 4: `qscat.core.channels` + `qscat.core.driven` (TI solver, model-taking)

**Files:** Create `libs/qscat/qscat/core/{channels.py, driven.py}`; create
`libs/qscat/tests/test_core_driven.py`; modify `projects/n2_2d_cross_section/{channels.py,
cross_section_2d.py}` (→ shims binding `N2`).

**Interfaces — Produces:** `channels.channel_vector(tgrid, k, chi_v, l)`;
`driven.ve_cross_section(tgrid, model, eps, chi, v_init, vprimes, E, *, ordering="COLAMD",
return_wavefunction=False)` — the exact TI driven L-S, analyze-once/`SparseLU.refactor` sweep,
σ = π|S−δ|²/2E. Uses `model.hamiltonian(tgrid)` and `model.interaction_diag(tgrid)`.

**Background:** promote `cross_section_2d.ve_cross_section_2d`; replace its internal
`build_h2d(tgrid)` / `interaction_diag(tgrid)` / `ELL` with `model.hamiltonian(tgrid)` /
`model.interaction_diag(tgrid)` / `model.ell`. The `lam_scale` test lever becomes a `model`
built with a scaled `V_int` (or keep a `lam_scale` kwarg that scales `interaction_diag` — keep
whichever the free-particle/Born tests need; preserve their behavior).

- [ ] **Step 1: Failing test.** Assert `ve_cross_section(TG, N2, eps, chi, 0, [0,1,2], E_arr)`
  equals the current `cross_section_2d.ve_cross_section_2d(TG, eps, chi, 0, [0,1,2], E_arr)` to
  round-off on a small grid + the gated anchor energies; assert `channel_vector` matches the
  current one; assert the scalar/array/`return_wavefunction` contract is preserved.
- [ ] **Step 2: Run → fail.**
- [ ] **Step 3: Implement** the promotion (model-taking signature).
- [ ] **Step 4: Rewire + gate.** `cross_section_2d.ve_cross_section_2d` becomes a shim:
  `ve_cross_section(tgrid, N2, …)`. Run the exact-2D anchor tests + harness group E; harness
  23/0/6/0; the #6 free-particle/first-Born/reciprocity limits still pass; mypy 0; ruff clean.
- [ ] **Step 5: Commit.**

---

### Task 5: `qscat.core.wavepacket` + `correlation` + `time_dependent` + `plot` (TD solver, model-taking)

**Files:** Create `libs/qscat/qscat/core/{wavepacket.py, correlation.py, time_dependent.py,
plot.py}`; create `libs/qscat/tests/test_core_td.py`; modify `projects/n2_2d_td_cross_section/*`
+ `projects/n2_2d_cross_section/cross_section_plot.py` (→ shims binding `N2`).

**Interfaces — Produces:** `wavepacket.{gaussian_coeffs, initial_state, outgoing_channel}`;
`correlation.{eta_incident, eta_outgoing}`; `time_dependent.{propagate, sigma_from_correlations,
td_ve_cross_section}` (`td_ve_cross_section(tgrid, model, eps, chi, v_init, vprimes, E, *, dt,
n_steps, wp_in, wp_out, order=3, subtract_free_reference=True)`, using
`qscat.evolution.make_pade_stepper` and, for the free reference, `model.hamiltonian −
diag(model.interaction_diag)`); `plot.plot_cross_sections(…)`.

**Background:** promote `td_propagation.propagate` (replace `build_h2d` with
`model.hamiltonian`; the free path uses `model.interaction_diag`), `td_cross_section.*` (the
Tannor-Weeks transform, `_s_vector_one_energy`, `_sigma_one_energy` with the free-reference
elastic), `wavepacket.py`, `correlation.py`, `cross_section_plot.py` verbatim (these are
already model-agnostic apart from the Hamiltonian source).

- [ ] **Step 1: Failing test.** Assert `td_ve_cross_section(TG, N2, …, order=3)` reproduces the
  current `projects.n2_2d_td_cross_section.td_cross_section.td_ve_cross_section_2d` on a small
  grid / few-step propagation (public-API shape contract, cheap); assert the free-reference
  elastic path and `sigma_from_correlations` match; keep one `@slow` differential vs the exact
  TI at an anchor (v'=1, E=0.10, rel 0.06 — the same gate as `test_td_cross_section.py`).
- [ ] **Step 2: Run → fail.**
- [ ] **Step 3: Implement** the promotion (model-taking).
- [ ] **Step 4: Rewire + gate.** TD project modules become shims binding `N2`. Run the TD fast
  tests + the `@slow` v2a/v2c differential; harness 23/0/6/0; mypy 0; ruff clean.
- [ ] **Step 5: Commit.**

---

### Task 6: Point N2 validation at `qscat.core`/`qscat.model`; retire the shims

**Files:** Modify `validation/n2/*` (`exact2d.py`, `td_exact2d.py`, `cross_section.py`,
`model.py`, loaders) to import `N2` from `qscat.model` and call `qscat.core` directly; delete
the now-empty/shim N2 project modules whose sole content moved to `qscat` (keep any genuinely
N2-specific research code — e.g. resonance-pole studies — in `projects/`).

- [ ] **Step 1:** Update `validation/n2/exact2d.py`/`td_exact2d.py` to `from qscat.model import
  N2` + `from qscat.core.driven import ve_cross_section` / `qscat.core.time_dependent`, passing
  `N2`. Run the harness → **23/0/6/0**, byte-identical NOTE/PASS numbers.
- [ ] **Step 2:** Remove shim modules that are now pure re-exports (or leave a one-line
  re-export where an external test still imports the old path — prefer updating the importer).
  Grep for stale imports; fix.
- [ ] **Step 3: Gate.** `pytest -q -m "not slow"` (full repo) → pass; the `@slow` TD + exact2d
  differentials → pass; harness 23/0/6/0; mypy 0; ruff clean.
- [ ] **Step 4: Commit.**

---

### Task 7: Docs, import guard, final verification

**Files:** Create `docs/physics/qscat-core-scattering.md`; modify `CLAUDE.md`; create
`libs/qscat/tests/test_core_no_model_import.py`.

- [ ] **Step 1: Import guard.** `test_core_no_model_import.py`: import every `qscat.core.*`
  module and assert none of them (nor `qscat.core` transitively) pulls in `qscat.model` —
  inspect `sys.modules` after a fresh import, or parse `qscat.core`'s imports. Run → pass.
- [ ] **Step 2: Docs.** `docs/physics/qscat-core-scattering.md`: the promoted method (TI driven
  + TD Padé/Tannor-Weeks), the `core`/`model` split, the `ResonanceModel` protocol, the
  extensibility guarantee (`core` ↛ `model`), cross-referencing #6/#7 and the TD-fix notes.
- [ ] **Step 3: CLAUDE.md.** Add `qscat.core` + `qscat.model` to the repo map (submodule list);
  note the N2 projects now consume them; note `validation/diatomic/` is where NO/F2 land next.
- [ ] **Step 4: Full verification:**
```
uv run pytest -q -m "not slow"          # all pass, N2 numbers unchanged
uv run pytest -q -m slow                # TD/exact2d differentials pass
uv run mypy libs/qscat                  # 0
uv run ruff check .                     # clean
uv run python -m validation.n2.experiment   # 23/0/6/0, unchanged
```
- [ ] **Step 5: Commit.**

---

## Final verification

- [ ] `qscat.core` (engine) + `qscat.model` (protocol + `DiatomicResonanceModel` + N2/NO/F2
  registry) + `qscat.special` radial functions exist and are tested.
- [ ] Every N2 physics number is unchanged: harness 23/0/6/0; exact-2D-vs-Houfek and TD-vs-TI
  differentials pass at their existing tolerances.
- [ ] `qscat.core` does not import `qscat.model` (guard test passes).
- [ ] N2 projects/validation consume `qscat.core`/`qscat.model`; no dead shims left.
- [ ] mypy 0; ruff clean; docs + CLAUDE.md updated.
- [ ] NO/F2 registry entries exist (data only) — ready for sub-projects B and C, which add
  only `validation/diatomic/<mol>/` on top.
