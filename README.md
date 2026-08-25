# qModeling

[![CI](https://github.com/VanaMartin/qscat/actions/workflows/ci.yml/badge.svg)](https://github.com/VanaMartin/qscat/actions/workflows/ci.yml)
[![Docs](https://github.com/VanaMartin/qscat/actions/workflows/docs.yml/badge.svg)](https://github.com/VanaMartin/qscat/actions/workflows/docs.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD--3--Clause-blue.svg)](LICENSE)
[![Python 3.12+](https://img.shields.io/badge/python-3.12%2B-blue.svg)](https://www.python.org/)

qModeling is a Python-first quantum-mechanics research monorepo, home to **QSCAT**
(`libs/qscat`) — a CPU-first electron-diatomic-molecule scattering toolkit built
around FEM-DVR-ECS (finite-element discrete-variable-representation, exterior-
complex-scaled) grids. It computes electron-scattering cross sections
(vibrational excitation, dissociative attachment/recombination) for a model
diatomic molecule, both time-independently (a driven Lippmann-Schwinger solve)
and time-dependently (wavepacket propagation), and treats the exact numerical
solution as an oracle against which cheaper approximations (e.g. the
local-complex-potential method) are validated. `qscat-run` (`apps/qscat-run`) is
the config-driven CLI on top of it: point it at a YAML config and a molecule
preset, and it runs the experiment and writes cross sections + plots + a
reproducibility manifest. See `CLAUDE.md` for the full operating manual and repo
map.

## Install

```bash
uv sync --all-packages   # installs qscat + qscat-run + builds the Rust qscat_kernels
uv run pytest -q -m "not slow"
```

`uv sync --all-packages` is required (not plain `uv sync`, which prunes the
workspace members) — see `CLAUDE.md`'s tech-decisions section.

## Quickstart: `qscat-run`

Scaffold a starter config for a molecule, seeded from its validated preset deck:

```bash
uv run qscat-run init F2 --observables ve,da --methods ti,td -o f2.yaml
```

This writes a fully-commented YAML file (energies, grid preset, TD incident/test-
function defaults all filled in) — open it and edit to taste. `qscat-run list`
shows every molecule, its preset variants, and which observables are valid for it:

```bash
$ uv run qscat-run list
Molecules (preset variants -- valid observables):
  F2   presets=['emoscat'] observables=['da', 'resonance_levels', 've']
  H2P  presets=['emoscat', 'proxy'] observables=['dr']
  N2   presets=['emoscat'] observables=['da', 've']  (closed-in-range: ['da'])
  NO   presets=['emoscat'] observables=['da', 'resonance_levels', 've']
```

Check a config is well-formed (schema + molecule/observable validity + grid
resolution — fast, no solve) before running it:

```bash
uv run qscat-run validate f2.yaml
```

### Run locally (small/fast grid)

`apps/qscat-run/examples/` has seventeen committed, `validate`-clean example
configs (twelve at the top level, five under `examples/figures/` that drive the
dense committed curves). `n2-ve.yaml` uses a tiny hand-written grid and runs in a
fraction of a second — good for a first local run:

```bash
uv run qscat-run run apps/qscat-run/examples/n2-ve.yaml --output runs/n2-ve
```

`--dry-run` resolves the config and prints the plan (grid sizes, energy count,
output location) without solving anything — useful for sanity-checking a bigger
config before committing to a real solve:

```bash
uv run qscat-run run f2.yaml --dry-run
```

### Run in Docker (full production decks)

The molecule presets' `emoscat` variant is the real, validated deck (up to
~1.15M unknowns for H2P) — solving it, especially with the MUMPS sparse-LU
backend, wants a Docker image with MUMPS built in (the `runtime` image
deliberately has none, to keep it lean). `docker/run.sh` builds one and runs a
config inside it in one step:

```bash
docker/run.sh apps/qscat-run/examples/h2p-dr.yaml runs/h2p-dr
```

This builds `qmodeling-base` then the `test-deps` app image — deliberately not
`test`, which would run both test tiers (the slow one is measured in minutes)
before every compute invocation. It then runs `qscat-run run` inside the
container with your config mounted read-only and the output directory mounted
out. Default output directory is `runs/<config-stem>` if you omit it. Use
`docker/build.sh test` when you actually want the suite; see `docker/README.md`
for the image layering.

### Reproducing a validated cross section

`apps/qscat-run/examples/f2-da.yaml` and `n2-ve.yaml` ship with a tiny explicit
grid for fast local iteration (their own YAML comments explain why); swap in
`grid: {preset: emoscat}` (as `qscat-run init` would generate) to run the real,
convergence-tested deck and reproduce the committed reference curves:

- **F2 DA/VE**: `docs/physics/figures/f2-2d-ti-da-cross-section.png` /
  `f2-2d-ti-cross-section.png`, computed from the dense configs in
  `apps/qscat-run/examples/figures/` (`f2-da-dense.yaml` / `f2-ve-dense.yaml`);
  the per-molecule curve drivers they replaced were retired into this CLI — see
  `docs/physics/diatomic-ve-cross-sections.md`.
- **N2 VE**: `docs/physics/figures/n2-2d-ti-cross-section.png`, validated against
  Karel Houfek's independent `CSVE.V00.J00` data — see
  `docs/physics/n2-2d-cross-section.md` / `docs/physics/ti-energy-sweep-reuse.md`.

`qscat-run`'s preset decks are literal, by-value transcriptions of the same
numbers those validated `validation/*` runners use (see
`apps/qscat-run/qscat_run/presets.py`'s module docstring), so a preset-grid CLI
run and the committed figure are the same physics, different entry point.

## Output artifacts

Every `qscat-run run` writes into `output_dir` (or `--output DIR`):

| file | when | contents |
|---|---|---|
| `cross_section.csv` / `.npz` / `.png` | `artifacts.cross_section` (default on) | σ(E) per requested observable/channel; TI and TD series overlay on one plot (disjoint `ti:`/`td:` key prefixes) |
| `cross_section_vs_time.npz` / `.png` | TD + `artifacts.cross_section_vs_time.moments` | σ(E) re-extracted at each requested truncation time — the converging view of a TD extraction |
| `correlations.npz` | TD + `artifacts.correlations: true` (opt-in) | the raw per-step correlation/flux series behind each extractor's transform |
| `wavefunction/psi_*.npz` / `.png` | `artifacts.wavefunction_snapshots` | electronic- and nuclear-projected `\|Ψ\|²` density: `psi_E{e}.*` (TI, at requested energies) or `psi_t{t}.*` (TD, at requested times) |
| `manifest.json` | always | `qscat` version, git SHA, timestamp, backend, per-stage timings, host platform — the reproducibility record |
| `config.resolved.yaml` | always | the fully-expanded config (presets inlined: energies, grid, TD incident/test-function) — enough to re-run exactly |

## Config schema (summary)

```yaml
molecule: F2                     # N2 | NO | F2 | H2P
methods: [ti, td]                # any subset of {ti, td, lcp, nrm}
observables:                     # a LIST -- VE and DA/DR can be requested together
  - {kind: ve, channels: [0, 1, 2]}   # final vibrational states (list) or a count (int)
  - {kind: da, channels: 1}           # dissociation channels (usually 1; dr for H2P)
energies: {min: 0.01, max: 0.05, step: 0.001}   # or {values: [0.03, 0.04]}
grid: {preset: emoscat}          # per-molecule preset, or an explicit {electronic, nuclear} grid
v_init: 0                        # initial vibrational level
td:                               # required iff 'td' in methods
  dt: 1.0
  n_steps: 1800
  order: 3                       # diagonal-Pade order (order-1 == Crank-Nicolson)
  extractors: [flow, delta, tw]
  incident: {r0: 45, p0: -0.35, sigma: 6}          # else the preset default
  test_function: {r0_out: 8, p0_out: 72, sigma_out: 0.07}   # or a per-kind {ve: {...}, da: {...}}
artifacts: {cross_section: true}
backend: auto                    # auto | mumps | scipy
output_dir: runs/f2-da
```

Only `molecule`, `methods`, `observables`, `output_dir` are required — everything
else is filled from the molecule's preset. `apps/qscat-run/qscat_run/config.py`
is the source of truth for the full schema — the frozen dataclasses and
`validate_config`'s checks — and `apps/qscat-run/README.md` walks through the
methods and the observable-to-artifact mapping in prose. `qscat-run list` gives
the per-molecule observable validity matrix at the command line.

## The key lesson: TD needs a launch-box grid

A time-dependent run's incident wavepacket needs somewhere to *launch from* —
the electronic real region has to be large enough to hold the incident packet,
its interaction with the target, and (for DA/DR) the outgoing test packet,
entirely inside the real (unscaled) grid. An incident placed in or past the ECS
tail diverges instead of propagating. This is why every molecule preset defines
**two** grids, not one: `ti_grid()` (a compact driven-solve deck) and `td_grid()`
(a wider "launch box", same nuclear grid, bigger electronic real region) — see
`apps/qscat-run/qscat_run/presets.py`. Using `grid: {preset: ...}` handles this
automatically; if you supply an **explicit** grid instead, it is used for both TI
and TD, so *you* own making sure it's TD-launch-box-sized if you request `td`.
See `docs/physics/n2-2d-td-cross-section.md` for the underlying physics (order-3
Padé propagation, the elastic free-reference subtraction).

## Further reading

- `CLAUDE.md` — the full repo map, lifecycle, and tech decisions.
- `docs/physics/` — the physics notes behind every validated method (FEM-DVR-ECS,
  MUMPS backend, N2/NO/F2/H2+ cross sections, the discretisation tuner, ...).
- `apps/qscat-run/README.md` — the `qscat-run` methods, the mapping from each
  observable to its config knob and artifact, and the architecture. Its dated
  design record is
  `docs/superpowers/specs/2026-08-01-qscat-run-cli-design.md`, which captures the
  original intent rather than the current schema.
- `docker/README.md` — the Docker image layering
  (base/build/test-deps/test/runtime).
