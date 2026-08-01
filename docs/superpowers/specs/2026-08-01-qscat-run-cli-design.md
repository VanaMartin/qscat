# `qscat-run` — a config-driven CLI for 2-D model experiments — Design Spec

**Date:** 2026-08-01
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — spec for review
**Lifecycle:** a production/tooling layer on top of the validated `qscat` solvers (not a new
numerical method). It generalizes the ad-hoc per-molecule runners in `validation/*/` into ONE
config-driven, dockerized entry point.

## Goal

Give a newcomer an easy, reproducible, dockerized way to run the 2-D electron-diatomic model for a
chosen molecule in TI, TD, or both, and extract the relevant cross sections (and TD moment-resolved
cross sections, optional correlation functions, and wavefunction snapshots) — all declared in a
single YAML config. A `README.md` walks a newcomer from install to reproducing a known cross section.

## Non-goals (out of scope)

- Replacing the existing `validation/*` runners — they stay (they carry the golden gates); the CLI
  is a parallel production entry point that can later share the presets.
- The `auto`-tuned grid path (`qscat.tuning`) — explicit + preset grids only for now.
- New physics/solvers — the CLI orchestrates the existing `qscat.core` observables only.
- A GUI / web UI.

## Architecture

A new **uv workspace member `apps/qscat-run`** (so `click`/`pyyaml`/`matplotlib` stay out of the lean
`qscat` core lib), providing a `qscat-run` console script. Pipeline:

    YAML config -> parse+validate (schema) -> resolve (molecule->model, preset->grid per method,
    defaults) -> run (one driven solve / one propagation, feeding all requested observables) ->
    write artifacts + manifest -> output_dir

`apps/qscat-run` depends only on `qscat` (+ click/pyyaml/matplotlib/numpy) — **never on
`validation/`**. Modules:

- `config.py` — the dataclass schema for a parsed+validated config (`ExperimentConfig`) + the YAML
  loader + validation (raises actionable errors: unknown molecule, invalid (molecule, observable)
  combo, TD block missing when `td` requested, etc.).
- `presets.py` — the single source of the validated per-molecule PRODUCTION decks: per molecule, the
  TI (driven-solve) grid, the TD (launch-box) grid, default energy range, default incident/
  test-function params, and the valid observables. Consolidates what today lives in
  `validation/diatomic/config.py`, `validation/h2plus/config.py`, and the eMoScat decks.
- `runner.py` — resolves the config and drives the solve(s): for TI, one `driven` Ψ₊ solve projected
  onto every requested observable/channel (VE via `ve_cross_section`, DA via `da_cross_section`, DR
  via `dr_cross_section`); for TD, ONE `propagate` feeding the electronic-axis VE extractors AND the
  nuclear-axis DA/DR extractors together (propagate-once/record-all). Returns an in-memory result
  bundle.
- `artifacts.py` — writes the result bundle to `output_dir` (CSV/npz/PNG + manifest + resolved
  config).
- `cli.py` — the `click` command group.

## Config schema (YAML)

```yaml
molecule: F2                     # N2 | NO | F2 | H2P
methods: [ti, td]                # any subset of {ti, td}
observables:                     # a LIST -- VE and DA/DR can be observed in ONE experiment
  - {kind: ve, channels: [0, 1, 2]}   # v' final vibrational states (list) or a count (int)
  - {kind: da, channels: 1}           # dissociation channels (usually 1)
  # - {kind: dr, channels: 3}         # count of DR/Rydberg channels (H2P only)
energies: {min: 0.01, max: 0.05, step: 0.001}   # or {values: [0.03, 0.04]}
grid: {preset: emoscat}          # per-molecule; TI->driven deck, TD->launch box. or explicit (below)
v_init: 0                        # initial vibrational level
td:                              # required iff 'td' in methods
  dt: 1.0
  n_steps: 1800
  order: 3                       # diagonal-Pade order
  extractors: [flow, delta, tw]  # any subset (DA/DR default flow; VE default tw)
  incident: {r0: 45, p0: -0.35, sigma: 6}          # else preset default
  test_function: {r0_out: 8, p0_out: 72, sigma_out: 0.07}   # nuclear TW test packet (else preset)
artifacts:
  cross_section: true                                  # sigma(E) per observable/channel
  cross_section_vs_time: {moments: [500, 1000, 1500]}  # TD: sigma(E) re-extracted at truncation times
  correlations: false                                  # OPTIONAL: raw c(t)/flux(t) per observed channel (TD)
  wavefunction_snapshots:
    td_times: [200, 800, 1500]                         # TD: |psi(t)|^2 density (r- and R-projected)
    ti_energies: [0.03]                                # TI: Psi_+(E)
backend: auto                    # auto | mumps | scipy  (SparseLU default_backend)
output_dir: runs/f2-da
```

**Minimal config** (preset fills the rest): `molecule`, `methods`, `observables`, `output_dir`.

**Explicit grid** (override the preset):
```yaml
grid:
  electronic: {real: [[6, 1.5], [20, 90.0]], ecs: {angle: 40, elements: 16, quadrature: 10}}
  nuclear:    {real: [[40, 10.7]],            ecs: {angle: 35, elements: 10, quadrature: 14}}
```
(An explicit grid is used for BOTH TI and TD; the user then owns the launch-box requirement — the
README states the KEY LESSON. Presets handle the TI/TD split automatically.)

## Observables + validity

`observables` is a list of `{kind, channels}`. Channel counts are simple (no nesting): `ve` — a
count or explicit list of `v'`; `da` — a count (usually 1); `dr` — a count of Rydberg channels.
Validity is checked at resolve time against the preset's declared observables:

| molecule | ve | da | dr |
|---|---|---|---|
| N2 | ✅ | (closed in-range) | — |
| NO | ✅ | ✅ | — |
| F2 | ✅ | ✅ | — |
| H2P | — | — | ✅ |

An invalid combo (e.g. `dr` on `F2`, or `da` on `N2`) is a validation error with an actionable
message. TD-DR/DA channel counts are bounded by what the electronic grid resolves (the Rydberg-count
= electronic-extent constraint from the H2P work) — a too-large count surfaces the existing
`anion_electronic_states` ValueError, wrapped with a hint.

## The shared-work runner

- **TI**: one driven Ψ₊(E) sweep (`SparseLU.refactor` per energy) — projected onto every requested
  observable's channels. VE, DA, DR each reuse the same Ψ₊ (DA/DR already do internally). One solve
  serves all TI observables.
- **TD**: ONE `propagate(psi0, H, extractors=[...])` — the extractor list assembled from ALL
  requested observables (electronic-axis TW/Dirac/Flux for VE + nuclear-axis for DA/DR), plus any
  wavefunction-snapshot recorders. Each extractor then yields its σ(E). The `cross_section_vs_time`
  moments re-run each extractor's transform on the correlation truncated at each `t_i`.

## Artifacts (in `output_dir`)

- `cross_section.csv` + `.npz` + `.png` — σ_{v→v'}(E) per observable/channel; TI and TD (all
  requested extractors) overlaid on the plot, tabulated in the CSV.
- `cross_section_vs_time.npz` + `.png` (TD) — σ(E) at the truncation moments (the converging view).
- `correlations.npz` (TD, opt-in) — raw c(t)/flux(t) per observed channel per extractor.
- `wavefunction/psi_t{t}.npz`/`.png` (TD density) and `psi_E{e}.npz`/`.png` (TI Ψ₊).
- `manifest.json` — resolved config, `qscat` version, git SHA, timestamp, backend, per-stage
  timings, host info — the reproducibility record.
- `config.resolved.yaml` — the fully-expanded config (presets inlined) for exact re-runs.

## CLI commands (`click` group `qscat-run`)

- `qscat-run run CONFIG [--output DIR] [--backend auto|mumps|scipy] [--dry-run]` — parse → resolve →
  run → write. `--dry-run` resolves + prints the plan (grids, solve sizes, estimated cost) without
  solving.
- `qscat-run validate CONFIG` — schema + preset resolution + validity, no solve (fast; for CI/edit
  loops).
- `qscat-run init MOLECULE [--observables ve,da --methods ti,td] -o config.yaml` — scaffold a
  fully-commented starter config seeded from the molecule's preset.
- `qscat-run list` — molecules, presets, and the valid (molecule, observable) matrix.

## Docker + README

- `docker/run.sh CONFIG [OUTPUT_DIR]` — builds the `test` image (has MUMPS) via `docker/build.sh
  test`, then `docker run` mounting `CONFIG` in (read-only) and `OUTPUT_DIR` out, invoking
  `qscat-run run /config.yaml --output /out`. MUMPS is auto-selected for large decks (H2P). Replaces
  the molecule-specific `docker/run-n2.sh` with a general one (keep `run-n2.sh` or fold it into a
  sample config).
- Top-level **`README.md`** — the newcomer walkthrough: (1) what qModeling/QSCAT is, one paragraph;
  (2) install (`uv sync --all-packages`); (3) `qscat-run init F2 -o f2.yaml`, edit; (4) run locally
  (`qscat-run run f2.yaml`) on a small/proxy grid; (5) run in Docker (`docker/run.sh f2.yaml
  runs/f2`) for the full decks; (6) the output artifacts explained one by one; (7) reproduce a
  validated cross section (F2 DA, N2 VE) and compare to the committed figures; (8) the KEY LESSON
  (TD needs a launch-box grid — presets handle it) and where the physics docs live.

## Validation / tests

- `apps/qscat-run/tests/`: schema validation (good/bad configs), preset resolution (per molecule/
  method), the validity matrix, `init`/`list`/`validate` command behavior, and a FAST end-to-end
  smoke — `qscat-run run` on a tiny-grid config (small molecule, few energies / few TD steps)
  producing all artifact types (cross_section, vs_time, snapshots, manifest) — asserting files exist
  + are finite/well-formed. Heavier full-deck runs are `@slow`/Docker, not the fast gate.
- The physics is already validated by the underlying solvers; the CLI tests cover orchestration +
  I/O + config correctness.

## Deliverables

- **D1** `apps/qscat-run` workspace member: `config.py` (schema+loader+validation), `presets.py`
  (per-molecule decks), `cli.py` (the 4 commands) — with `validate`/`init`/`list` fully working.
- **D2** `runner.py` (TI shared-solve + TD shared-propagation, all observables) + `artifacts.py`
  (all artifact writers + manifest + resolved config).
- **D3** `docker/run.sh` + the top-level `README.md` newcomer walkthrough.
- **D4** tests (schema/preset/validity/commands + the fast end-to-end smoke) + a couple of committed
  sample configs (`examples/*.yaml`: F2 DA, N2 VE, H2P DR).

## Verification

- `uv run pytest -q -m "not slow"` incl. the CLI tests pass; `uv run mypy` clean over the new
  package; `uv run ruff check .` clean.
- `qscat-run run examples/f2-da.yaml` on a proxy grid produces every declared artifact + a manifest;
  `qscat-run validate` rejects the documented bad configs.
- `docker/run.sh examples/f2-da.yaml runs/f2` runs end-to-end in the container (documented; the full
  H2P deck is the MUMPS path).
- The README reproduces a committed validated cross section.
