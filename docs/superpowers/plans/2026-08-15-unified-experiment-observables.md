# Unified experiment observables — implementation plan

> **For agentic workers:** executed INLINE (subagent cap reached this session). Steps use
> checkbox (`- [ ]`) syntax for tracking. Design source: the harsh audit
> (`scratchpad/fable-audit.md`) the maintainer reviewed and approved ("run sequentially all").

**Goal:** make `apps/qscat-run` the single execution surface that turns one config into ALL common
observables (cross sections, energy levels + eigenstate wavefunctions, resonance states,
scattering/propagated wavefunctions) across all three methods (TI, TD, **LCP**) for any molecule or
ion, then retire the duplicative per-molecule `validation/*` curve drivers.

**Architecture:** grow `qscat-run` additively (R5, R2, R6, R3, R1 — each independently green), THEN
consolidate the `validation/*` drivers onto it (R4, R7, R8). `qscat.core`/`qscat.dvr` gain the small
shared primitives; `qscat_run` gains the observable dispatch + artifacts; `validation/*` keeps only
golden-data anchor gates.

**Tech stack:** Python 3.12, uv, pytest, numpy/scipy, `qscat` library, `qscat_run` app.

## Global constraints
- Atomic units throughout; no ad-hoc conversions.
- `qscat.core` never imports `qscat.model`/`projects` at runtime (protocol-only; enforced by
  `test_core_no_model_import.py`). New core code obeys this.
- Behavior-preserving refactors must keep every existing test green; capability adds ship with
  their own tests. `ruff check`, `ruff format`, `mypy libs/qscat` clean.
- Docker/MUMPS-only decks (h2plus full, ~1.15M unknowns) CANNOT be figure-verified on the laptop —
  migrations of those keep the driver until a container run confirms parity; do not delete blind.
- Commit per task with the standard trailers; open ONE PR at a coherent milestone (maintainer
  merges themselves).

---

## Task order & status

- [x] **R5** — hoist the real-region-DVR-index helper into `qscat.dvr` (kill 4 copies). *refactor*
- [x] **R2** — `eigenstates` observable: vibrational (eps+chi) + resonance poles in the result/artifacts. *capability*
- [x] **R6** — `return_wavefunction` flag on `da_cross_section`/`dr_cross_section` (parity with `ve`). *capability*
- [x] **R3** — full complex Ψ field in `WavefunctionSnapshot` + viz-ready npz artifact. *capability*
- [x] **R1** — `lcp` as a `qscat-run` method; LCP angles + nuclear deck fields on `MoleculePreset`. *capability*
- [ ] **R4** — migrate diatomic+h2plus curve drivers to committed `*.yaml` + generic figure step; delete superseded drivers/config schemas. *refactor/debt*
- [ ] **R7** — migrate `n2/ti_curve.py` curve-gen to qscat-run; keep anchor gates; drop `projects.*` reach-through. *refactor*
- [ ] **R8** — single `FIGURE_DIR`; wire `nuclear_density` as a real observable or delete. *refactor*

---

### R5 — shared real-region DVR index helper
**Files:** `libs/qscat/qscat/dvr/grid.py` (add method) or a small function in `qscat.dvr`; export in
`__init__`; test `libs/qscat/tests/test_dvr_index.py`. Callers: `apps/qscat-run/qscat_run/runner.py`
(`_index_near`), `validation/diatomic/td_da.py` (`_real_index_near`), `validation/h2plus/td_dr.py`
(`_nuclear_index_near`).
**Interface produced:** `FemDvrEcsGrid.real_index_near(r_value: float) -> int` — nearest unscaled
real-region DVR index (complex-tail points masked out). Byte-identical to the 4 existing copies.
**Acceptance:** new unit test passes; all callers delegate; full fast suite green.

### R2 — eigenstates observable
**Files:** `apps/qscat-run/qscat_run/config.py` (ArtifactSpec `eigenstates`), `runner.py` (attach
`EigenStates` to result), `artifacts.py` (write npz), tests under `apps/qscat-run/tests/`.
**Interface:** `EigenStates{kind: "vibrational"|"resonance", energies, states, axis}` on
`ExperimentResult`. Vibrational from the already-computed `eps, chi`; resonance via
`qscat.ecs.find_resonance_pole`.
**Acceptance:** a config requesting `artifacts: {eigenstates: {...}}` emits eps + chi columns; test
asserts shapes + that chi diagonalizes the nuclear Hamiltonian.

### R6 — da/dr return_wavefunction
**Files:** `libs/qscat/qscat/core/dissociation.py`, test in
`libs/qscat/tests/test_dissociation*.py`.
**Interface:** `da_cross_section(..., return_wavefunction=False)` /
`dr_cross_section(..., return_wavefunction=False)` mirroring `ve_cross_section`'s overloads.
**Acceptance:** flag returns Ψ⁺ (same object `ve` would); default path byte-identical.

### R3 — full-field wavefunction snapshot + viz npz
**Files:** `runner.py` (`WavefunctionSnapshot` optional complex `psi` on `tg.shape`), `artifacts.py`
(npz emit), maybe a `qscat.viz` convenience loader; tests.
**Acceptance:** opt-in full field is emitted and round-trips through
`qscat.viz.plot_wavefunction_2d`; marginals stay the default.

### R1 — lcp method
**Files:** `config.py` (`VALID_METHODS += lcp`), `presets.py` (`MoleculePreset` lcp fields +
`resolve_grid` lcp path), `runner.py` (`_run_lcp`), `artifacts.py`, examples, tests.
**Acceptance:** `methods: [ti, lcp]` on F2 overlays exact vs LCP σ_DA from one config, matching
`validation/diatomic/lcp_da_curves.py`'s numbers on a small grid.

### R4/R7/R8 — consolidation (BLOCKED on a shared-schema dependency — see finding)
Migrate each `validation/*` curve driver to a committed `*.yaml` + the generic figure step; delete
the driver + its rogue config schema once a run reproduces its committed figure (laptop-feasible
grids only; Docker-only decks stay until container-verified). Collapse the 3 `FIGURE_DIR` copies.

**FINDING (2026-08-15, discovered during R4 blast-radius mapping):** the two "rogue config
schemas" the audit flagged for deletion are NOT driver-internal — they are the eMoScat-deck source
of truth for a SEPARATE sub-project, the discretisation tuner:
- `validation/diatomic/config.py` (`MoleculeConfig`/`CONFIGS`) is imported by
  `validation/tuning/calibrate.py`, `test_resonance_aware.py`, `test_emoscat_decks.py` (the tuner's
  calibration + gate), plus `test_diatomic.py`/`test_da_grid.py`.
- `validation/h2plus/config.py` (the free functions) is imported by the same three tuning modules
  plus `test_config.py`.

So R4 is not "delete 7 drivers"; it is a cross-sub-project migration: the tuner's deck references
must first move onto `qscat_run.presets` (or the decks must be relocated to a shared, non-driver
home) BEFORE any config schema can be deleted. Combined with (a) figure/data-parity verification
needing slow full-preset runs, (b) the h2+ deck being Docker/MUMPS-only, and (c) the n2 curve's
bespoke Houfek-reference overlay (a golden-data gate the audit says KEEP), the deletion step is a
deliberate, reviewed refactor — NOT an autonomous sweep. **Recommended sequencing for R4:**
1. Relocate the eMoScat per-molecule deck definitions to a single shared module (e.g.
   `validation/decks.py` or fold into `qscat_run.presets`), repoint the tuner + drivers at it.
2. Verify qscat-run reproduces each committed figure's DATA on the preset deck (F2/NO on laptop;
   h2+ under Docker) — the acceptance gate.
3. Only then delete the superseded curve drivers, keeping the n2 Houfek gate + `validation/tuning`.
The capability half (R1/R2/R3/R5/R6) is independent and complete; R4/R7/R8 ship as a follow-on PR.
