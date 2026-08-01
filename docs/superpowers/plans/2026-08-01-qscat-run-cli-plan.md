# `qscat-run` config-driven experiment CLI — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A YAML-config-driven, dockerized CLI (`qscat-run`) on top of the `qscat` solvers: pick a molecule, run TI/TD/both, extract a LIST of observables (VE + DA/DR in one experiment), and write cross sections + TD moment-resolved σ(E) + optional per-channel correlations + wavefunction snapshots, with a reproducibility manifest.

**Architecture:** A new uv workspace member `apps/qscat-run` (keeps `click`/`pyyaml`/`matplotlib` out of core `qscat`). Pipeline: YAML → parse+validate → resolve (molecule→model, preset→grid per method, defaults) → run (ONE TI driven-Ψ₊ sweep / ONE TD propagation feeding all observables) → write artifacts + manifest. Depends only on `qscat`, never on `validation/`.

**Tech Stack:** Python ≥3.12, `click`, `pyyaml`, `matplotlib`, numpy; `qscat.model` (N2/NO/F2/H2P), `qscat.core` (`ve/da/dr_cross_section`, `td_ve/td_da_cross_sections_all`, the `Extractor` classes, `propagate`, `initial_state`, `vibrational_states`), `qscat.core.grids` (`segmented_grid`, `fem_grid_exp_tail`, `electronic_grid`, `nuclear_grid`), `qscat.dvr.TensorGrid`, `qscat.linalg.SparseLU` (`set_default_backend`). pytest; mypy; ruff.

## Global Constraints

- **Atomic units** throughout (the CLI passes floats to `qscat`; no unit conversion).
- **`apps/qscat-run` depends ONLY on `qscat`** (+ click/pyyaml/matplotlib/numpy) — NEVER `import validation.*` or `projects.*`. Per-molecule decks are re-declared in `presets.py` (consolidated from `validation/diatomic/config.py` + `validation/h2plus/config.py` values, copied as data, not imported).
- **The physics is unchanged** — the CLI only orchestrates existing `qscat.core` observables; it adds no numerics.
- **Config format is YAML.** DR channel counts are FLAT integers (no nesting). `observables` is a LIST.
- **Grid presets are TI/TD-aware**: a preset resolves to the driven-solve deck for TI and the launch-box deck for TD (the KEY LESSON: an off-box TD incident diverges).
- mypy clean over `apps/qscat-run/qscat_run`; ruff clean; new deps added to `apps/qscat-run/pyproject.toml` only.
- The root `pyproject.toml` `[tool.uv.workspace].members` gains `"apps/*"`; setup stays `uv sync --all-packages`.

## File Structure

- `apps/qscat-run/pyproject.toml` — the member (name `qscat-run`, `[project.scripts] qscat-run = "qscat_run.cli:main"`, deps click/pyyaml/matplotlib/qscat, hatchling).
- `apps/qscat-run/qscat_run/__init__.py`
- `apps/qscat-run/qscat_run/config.py` — `ExperimentConfig` (+ nested dataclasses) + `load_config(path)` + `validate_config` (actionable errors).
- `apps/qscat-run/qscat_run/presets.py` — `MODELS` (name→model), `PRESETS` (per-molecule decks: TI grid, TD launch grid, defaults, valid observables), `resolve_grid(cfg, method)`, `resolve_defaults(cfg)`.
- `apps/qscat-run/qscat_run/runner.py` — `run_experiment(cfg) -> ExperimentResult` (TI + TD).
- `apps/qscat-run/qscat_run/artifacts.py` — `write_artifacts(result, cfg, output_dir)`.
- `apps/qscat-run/qscat_run/cli.py` — the `click` group (`run`/`validate`/`init`/`list`) + `main`.
- `apps/qscat-run/tests/` — schema/preset/validity/command tests + the fast end-to-end smoke.
- `apps/qscat-run/examples/{f2-da,n2-ve,h2p-dr}.yaml` — sample configs.
- `docker/run.sh` — general dockerized runner. `README.md` (top-level) — newcomer walkthrough.

---

### Task 1: Package scaffold + config schema + presets + no-solve commands (`validate`/`list`/`init`)

**Files:** create `apps/qscat-run/pyproject.toml`, `qscat_run/{__init__,config,presets,cli}.py`; modify root `pyproject.toml` (workspace members); test `apps/qscat-run/tests/test_config.py`, `test_cli_nosolve.py`.

**Interfaces (Produces):**
- `config.ExperimentConfig` — frozen dataclass: `molecule: str`, `methods: tuple[str,...]` (subset of `{"ti","td"}`), `observables: tuple[Observable,...]` (`Observable(kind: str, channels: int | tuple[int,...])`), `energies: EnergySpec` (`min/max/step` or `values`), `grid: GridSpec` (`preset: str | None`, `electronic/nuclear: SegmentSpec | None`), `v_init: int`, `td: TdSpec | None` (`dt/n_steps/order/extractors/incident/test_function`), `artifacts: ArtifactSpec`, `backend: str`, `output_dir: str`.
- `config.load_config(path: str | Path) -> ExperimentConfig` — YAML → dataclass; `config.validate_config(cfg) -> None` — raises `ConfigError` (a `click.ClickException` subclass) with actionable messages: unknown molecule, invalid `(molecule, observable.kind)` per the validity matrix, `td` block required when `"td" in methods`, `dr` requires `H2P`, explicit grid missing electronic/nuclear, unknown extractor, etc.
- `presets.MODELS: dict[str, ResonanceModel]` (`{"N2": N2, ...}`); `presets.PRESETS: dict[str, MoleculePreset]` where `MoleculePreset` carries `ti_grid()`, `td_grid()` (TensorGrid builders), `default_energies`, `default_incident`, `default_test_function`, `valid_observables: frozenset[str]`, `n_vib`. Values copied from `validation/diatomic/config.py` (F2/NO/N2) + `validation/h2plus/config.py` (H2P full/proxy) — as literals, no import.
- `presets.resolve_grid(cfg, method) -> TensorGrid`; `presets.resolve_defaults(cfg) -> ExperimentConfig` (fills omitted energies/incident/test_function/channels from the preset).

**Design notes:** the validity matrix — `{"N2": {"ve"}, "NO": {"ve","da"}, "F2": {"ve","da"}, "H2P": {"dr"}}` (N2 `da` is closed in-range → allow with a warning, per the spec; treat `da` as valid for N2 but `validate` emits a warning). `MoleculePreset.td_grid()` uses the launch-box electronic grid (large `r_max`) + the fine nuclear deck; `ti_grid()` uses the driven-solve deck (`da_grid`-style). For H2P, `ti_grid`/`td_grid` = `full_grid` (Docker) with a `proxy` size variant (`preset: proxy` → the reduced grids). Preset names: `emoscat` (full/validated), `proxy` (laptop). Keep `presets.py` pure data + builder functions.

- [ ] **Step 1: Write failing tests** (`test_config.py`, `test_cli_nosolve.py`):
```python
# test_config.py
from qscat_run.config import load_config, validate_config, ConfigError
import pytest, textwrap

def _write(tmp_path, text):
    p = tmp_path / "c.yaml"; p.write_text(textwrap.dedent(text)); return p

def test_minimal_config_loads_and_resolves(tmp_path):
    cfg = load_config(_write(tmp_path, """
        molecule: F2
        methods: [ti]
        observables: [{kind: ve, channels: 2}, {kind: da, channels: 1}]
        output_dir: out
    """))
    validate_config(cfg)                       # no raise
    assert cfg.molecule == "F2"
    assert [o.kind for o in cfg.observables] == ["ve", "da"]

def test_dr_on_neutral_rejected(tmp_path):
    cfg = load_config(_write(tmp_path, """
        molecule: F2
        methods: [ti]
        observables: [{kind: dr, channels: 3}]
        output_dir: out
    """))
    with pytest.raises(ConfigError, match="dr.*H2P|not valid for F2"):
        validate_config(cfg)

def test_td_method_requires_td_block(tmp_path):
    cfg = load_config(_write(tmp_path, """
        molecule: F2
        methods: [td]
        observables: [{kind: da, channels: 1}]
        output_dir: out
    """))
    with pytest.raises(ConfigError, match="td"):
        validate_config(cfg)
```
```python
# test_cli_nosolve.py -- click CliRunner over validate/list/init (no solve)
from click.testing import CliRunner
from qscat_run.cli import main

def test_list_shows_molecules():
    r = CliRunner().invoke(main, ["list"])
    assert r.exit_code == 0 and "F2" in r.output and "H2P" in r.output

def test_init_scaffolds_valid_config(tmp_path):
    out = tmp_path / "f2.yaml"
    r = CliRunner().invoke(main, ["init", "F2", "--observables", "ve,da", "-o", str(out)])
    assert r.exit_code == 0 and out.exists()
    r2 = CliRunner().invoke(main, ["validate", str(out)])
    assert r2.exit_code == 0

def test_validate_rejects_bad(tmp_path):
    bad = tmp_path / "bad.yaml"; bad.write_text("molecule: XX\nmethods: [ti]\nobservables: []\noutput_dir: o\n")
    assert CliRunner().invoke(main, ["validate", str(bad)]).exit_code != 0
```
- [ ] **Steps 2–5:** run→fail; create the package (pyproject + `apps/*` in workspace + `uv sync --all-packages`), implement `config`/`presets`/`cli` (validate/init/list; `run` is a stub raising "implemented in Task 2/3" for now); run→pass; `uv run mypy apps/qscat-run/qscat_run` + `uv run ruff check apps/qscat-run`; commit `feat(cli): qscat-run scaffold + config schema + presets + validate/init/list`.

---

### Task 2: TI runner + cross-section artifacts + `run` (method=ti)

**Files:** create `qscat_run/runner.py`, `qscat_run/artifacts.py`; modify `qscat_run/cli.py`; test `tests/test_runner_ti.py`.

**Interfaces:**
- `runner.run_experiment(cfg) -> ExperimentResult` (TI path here; TD in Task 3). `ExperimentResult` holds: `cross_sections: dict[str, CrossSectionSeries]` keyed by `"{method}:{observable}:{extractor?}"`, the energy grid, `wavefunctions: list[WavefunctionSnapshot]`, `timings: dict`, resolved grids. TI path: set `SparseLU` backend from `cfg.backend`; build `eps, chi = vibrational_states(nuc, model.mu, n_vib, model.v0)`; for each observable dispatch `ve_cross_section`/`da_cross_section`/`dr_cross_section` on the SAME `resolve_grid(cfg,"ti")` (they internally reuse the driven Ψ₊). TI wavefunction snapshots: `ve_cross_section(..., return_wavefunction=True)` at `artifacts.wavefunction_snapshots.ti_energies`.
- `artifacts.write_artifacts(result, cfg, output_dir)` — writes (this task: TI): `cross_section.csv` (columns: energy, then one per `method:observable:channel`), `cross_section.npz`, `cross_section.png` (matplotlib, log-y σ vs E, one line per series), `wavefunction/psi_E{e}.npz`+`.png` (TI Ψ₊ density projected onto r and R), `config.resolved.yaml` (dump the resolved cfg), `manifest.json` (`qscat` version via `importlib.metadata`, git SHA via `subprocess` best-effort, timestamp passed IN — do not call `datetime.now()` at import; the CLI stamps it, timings, backend, host `platform.platform()`).
- `cli.run`: parse→validate→resolve→`run_experiment`→`write_artifacts`; `--output` overrides `output_dir`; `--backend`; `--dry-run` prints the resolved plan (grids, sizes, observables) and exits.

**Design notes:** `manifest.json` timestamp — the CLI captures `time.time()`/an ISO string once and passes it into `write_artifacts` (keeps `artifacts.py` import-time-pure). Plot helper is shared with Task 3. Keep each series' provenance (method/observable/channel) explicit in the CSV header + npz keys.

- [ ] **Step 1: Write the failing test** (`tests/test_runner_ti.py`) — a TINY grid TI config (small F2, 2 energies, ve+da), run end-to-end, assert artifacts exist + finite:
```python
def test_ti_run_writes_cross_section_and_manifest(tmp_path):
    # explicit tiny grid (fast); ve+da; 2 energies
    cfg = load_config(_write(tmp_path, tiny_f2_ti_yaml(output=str(tmp_path/"out"))))
    validate_config(cfg)
    result = run_experiment(cfg)
    write_artifacts(result, cfg, tmp_path/"out", timestamp="2026-01-01T00:00:00")
    assert (tmp_path/"out"/"cross_section.csv").exists()
    assert (tmp_path/"out"/"manifest.json").exists()
    assert (tmp_path/"out"/"config.resolved.yaml").exists()
    arr = np.load(tmp_path/"out"/"cross_section.npz")
    assert np.all(np.isfinite(arr["energy"]))
```
Use an EXPLICIT tiny grid in the yaml (electronic r_max~8, nuclear r_max~10, low order) so it runs in seconds.
- [ ] **Steps 2–5:** run→fail; implement `runner` (TI) + `artifacts` + wire `cli.run`; run→pass; also `CliRunner().invoke(main, ["run", cfg, "--dry-run"])` prints a plan; mypy+ruff; commit `feat(cli): TI runner + cross-section/manifest artifacts + run command`.

---

### Task 3: TD runner + TD artifacts (moments, correlations, snapshots) + `run` (td/both)

**Files:** modify `qscat_run/runner.py`, `qscat_run/artifacts.py`, `qscat_run/cli.py`; test `tests/test_runner_td.py`.

**Interfaces:**
- `runner.run_experiment` TD path: `resolve_grid(cfg,"td")`; `psi0 = initial_state(tg, chi[v_init], **incident)`; build the extractor set for ALL requested observables — electronic-axis `TannorWeeks/Dirac/Flux` (VE, with `vprimes`/`wp_out`) and nuclear-axis ones (DA/DR, with `surface`/`position`/`wp_out`/`n_channels`, from the test_function + nuclear-surface-near-R) — for each `extractor` in `cfg.td.extractors`; ONE `propagate(tg, psi0, [], dt, n_steps, hamiltonian=model.hamiltonian(tg), extractors=[...])`; then `ext.sigma(E)` per extractor → cross sections keyed `td:{observable}:{extractor}`. (Reuse `td_ve_cross_sections_all`/`td_da_cross_sections_all` where the observable set is exactly VE-only or DA-only; use the raw extractors when mixing VE+DA in one propagation.)
- `cross_section_vs_time`: propagate once to `max(moments)`; for each `t_i` in `artifacts.cross_section_vs_time.moments`, re-extract σ(E) from each extractor's correlation TRUNCATED at `n_i = round(t_i/dt)` samples. **The one non-trivial piece:** the extractors accumulate a full series but expose no "σ up to step n". Add a minimal truncation path — either (a) a `truncate(n)` on the extractor returning a copy transformed on `series[:n]`, or (b) a runner-level helper that slices the extractor's recorded rows and calls the same transform. Prefer (b) if it avoids touching `qscat.core`; if the extractor internals aren't sliceable from outside, add a small `sigma(E, *, n_steps=None)` kwarg to the extractors (a `qscat.core` change — then VE byte-identical must hold: `n_steps=None` == today).
- `correlations` (opt-in): dump each extractor's recorded series (`ext` exposes its `t`/`c` via `.result` or `._arrays`) to `correlations.npz` keyed by observable/channel/extractor.
- TD wavefunction snapshots: a lightweight recorder passed to `propagate` (or reuse the `snapshots`/`sample_period` path) capturing `|psi(t)|^2` projected onto r and R at `artifacts.wavefunction_snapshots.td_times`; written to `wavefunction/psi_t{t}.npz`+`.png`.
- `artifacts`: add `cross_section_vs_time.npz`+`.png` (σ(E) curves per moment), `correlations.npz`, the TD `psi_t{t}` writers. `methods: [ti, td]` overlays both on `cross_section.png`.

**Design notes:** flag the moment-truncation as the piece needing care; keep any `qscat.core` change (the `n_steps=` truncation kwarg) minimal + guarded by the existing VE golden/differential tests (`n_steps=None` == today, byte-identical). The nuclear surface/position DVR index comes from a `_nuclear_index_near(tg, R)` helper (as in `validation/h2plus/td_dr.py`).

- [ ] **Step 1: Write the failing test** (`tests/test_runner_td.py`) — a TINY-grid TD config (small molecule, few steps, ve+da, one moment, one snapshot), run end-to-end:
```python
def test_td_run_writes_all_td_artifacts(tmp_path):
    cfg = load_config(_write(tmp_path, tiny_td_yaml(output=str(tmp_path/"out"))))  # dt, n_steps~10, moments:[..], td_times:[..]
    result = run_experiment(cfg)
    write_artifacts(result, cfg, tmp_path/"out", timestamp="t")
    assert (tmp_path/"out"/"cross_section_vs_time.npz").exists()
    assert list((tmp_path/"out"/"wavefunction").glob("psi_t*.npz"))
    # correlations opt-in:
    # (a second config with correlations: true asserts correlations.npz exists)
```
- [ ] **Steps 2–5:** run→fail; implement the TD runner path + moment truncation + correlations + snapshots + the TD artifact writers + wire `run` for `td`/`both`; if a `qscat.core` extractor `n_steps=` kwarg was needed, run the VE golden/differential tests to confirm byte-identical; run→pass; mypy+ruff; commit `feat(cli): TD runner (moments, correlations, snapshots) + full run command`.

---

### Task 4: Docker + README + examples + end-to-end

**Files:** create `docker/run.sh`, top-level `README.md` (or update if present), `apps/qscat-run/examples/{f2-da,n2-ve,h2p-dr}.yaml`; modify `docker/README.md`; test `apps/qscat-run/tests/test_examples.py`.

**Interfaces / content:**
- `docker/run.sh CONFIG [OUTPUT_DIR]`: `docker/build.sh test` then `docker run --rm -v $(realpath CONFIG):/config.yaml:ro -v $(realpath OUTPUT_DIR):/out qmodeling:test uv run --no-sync qscat-run run /config.yaml --output /out`. MUMPS auto-selected in the `test` image (H2P). `set -euo pipefail`; default `OUTPUT_DIR=runs/<config-stem>`.
- Top-level `README.md`: the newcomer walkthrough (spec §Docker+README): what QSCAT is (1 para); install (`uv sync --all-packages`); `qscat-run init F2 -o f2.yaml` → edit; run locally on a proxy grid; `docker/run.sh f2.yaml runs/f2` for full decks; each artifact explained; reproduce a validated cross section (F2 DA / N2 VE) vs the committed figure; the KEY LESSON (presets handle the TD launch box) + links to `docs/physics/`.
- `examples/*.yaml`: F2 DA (ti+td), N2 VE (ti), H2P DR (td, proxy + a note the full deck is Docker/MUMPS) — each `qscat-run validate`-clean and using a proxy/small grid so the examples double as runnable smokes.
- Test `test_examples.py`: every `examples/*.yaml` passes `validate_config`; ONE fast end-to-end (`n2-ve` on a proxy/tiny grid) runs + produces artifacts. Mark heavier ones `@slow`.

- [ ] **Steps:** write `test_examples.py` (validate all examples + one fast e2e) → fail; create examples + `docker/run.sh` + README + docker/README update; run→pass (fast e2e); shellcheck `docker/run.sh` mentally + `bash -n`; mypy+ruff; commit `feat(cli): docker/run.sh + newcomer README + example configs + e2e smoke`.

---

## Verification (whole sub-project)

- `uv run pytest -q -m "not slow"` (incl. the CLI tests) passes; the VE golden/differential tests still pass if a `qscat.core` truncation kwarg was added (byte-identical). `uv run mypy apps/qscat-run/qscat_run` + `uv run mypy libs/qscat/qscat` clean; `uv run ruff check .` clean.
- `qscat-run validate examples/*.yaml` all clean; `qscat-run run examples/n2-ve.yaml` (proxy) writes cross_section + manifest + resolved config; `--dry-run` prints a plan; `list`/`init` work.
- `apps/qscat-run` imports no `validation`/`projects` (a test asserting this, mirroring `test_core_no_model_import.py`).
- `bash -n docker/run.sh` clean; the README reproduces a committed validated cross section (documented; full decks are the Docker/MUMPS path).

## Out of scope (this plan)

- Replacing the `validation/*` runners (they stay; may later import the CLI presets).
- The `auto`-tuned grid path (`qscat.tuning`).
- Parallel/distributed runs, checkpoint/resume of a propagation, a web UI.
