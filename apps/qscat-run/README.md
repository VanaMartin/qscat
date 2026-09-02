# qscat-run

**The single execution surface for qModeling.** One YAML config turns any
molecule or ion into all the common observables across all three methods —
no per-experiment Python. Adding a molecule is a `MoleculePreset` entry
(`qscat_run/presets.py`), never solver code.

qscat-run is **repo-only**: it is not published to PyPI (and will not be
until the qscat citation article is out). Install it from a clone with
`uv sync --all-packages`.

```bash
qscat-run list                                                      # molecules, presets, valid observables
qscat-run init F2 --observables ve,da --methods ti,lcp -o f2.yaml   # scaffold
qscat-run validate f2.yaml                                          # actionable schema check
qscat-run run f2.yaml --output runs/f2                              # solve + write artifacts
qscat-run fetch validation/factory/results/o2-ve                    # download published results
```

The production preset decks are `O(10^4–10^6)` unknowns — run those under
Docker (`docker/run.sh <config> <out>`), which provides MUMPS. The committed
`examples/*.yaml` use tiny explicit grids for fast local iteration.

A run directory records the commit it came from in `manifest.json`'s
`git_sha`. It is a record of what produced the numbers, not an address for
them — published artifacts are addressed by content digest — so a run that
cannot determine it warns and carries on rather than failing. Inside a
container, where the build context excludes `.git`, the host passes it in
(`docker/run.sh` and `docker/build.sh` do this via `--build-arg GIT_SHA`);
`QSCAT_ALLOW_UNKNOWN_SHA=1` silences the warning where there is genuinely no
provenance to report. Publishing is stricter: the publisher refuses a
manifest that cannot name its commit.

## Results that are not in the repository

Converged sweeps are a few hundred kilobytes of CSV costing minutes to hours
of MUMPS solves, so they are published to object storage rather than
committed, and the run directory carries a small `artifacts.json` naming
them. `qscat-run fetch DIR` downloads what it names and verifies every byte
against the digest recorded at publication; a directory with no
`artifacts.json` keeps its results in git and needs no fetching.

Only outputs are published. The run's `config.resolved.yaml` and
`manifest.json` stay committed beside the pointer, so a clone with no network
still has the input and the provenance for every published number.

Reads are anonymous HTTPS — no account, no credentials. Publishing is
maintainer-only and lives outside this repository. What stays in git is
everything a test or a note depends on; see
`docs/adr/0008-computed-artifacts-live-in-public-object-storage.md` for where
the line is drawn, and `docs/artifacts.md` for the reader's side.

## Registry

| molecule | charge | valid observables | preset variants |
|---|---|---|---|
| `N2` | 0 | `ve` (`da` closed-in-range) | `emoscat` |
| `NO` | 0 | `ve`, `da`, `resonance_levels` | `emoscat` |
| `F2` | 0 | `ve`, `da`, `resonance_levels` | `emoscat` |
| `H2P` | −1 | `dr` | `emoscat`, `proxy` |
| `O2` | 0 | `ve` (DA closed below 3.7 eV) | `tuner` — the fitted model (`qscat.model.O2`, potential factory) on the discretisation tuner's deck; no eMoScat deck exists |
| `O2_SO12`, `O2_SO32` | 0 | `ve` | `tuner` — O₂'s two spin–orbit components (²Π_{1/2}, ²Π_{3/2}; statistical factor ⅓ each, summed by `validation/factory/o2_ve_figure.py`) on the same deck |

## Methods (`methods: [...]`)

- **`ti`** — exact time-independent driven Lippmann–Schwinger solve.
- **`td`** — time-dependent wavepacket propagation (order-3 Padé), with the
  Tannor–Weeks / Dirac / Flux energy extractors (`td.extractors: [tw,delta,flow]`).
- **`lcp`** — the local-complex-potential *approximation* of DA **and VE**
  (N2/F2/NO). `methods: [ti, lcp]` overlays the exact oracle and the
  approximation on one `cross_section.png` (keys `ti:da:ch0` vs `lcp:da:ch0`,
  or `ti:ve:v0->1` vs `lcp:ve:v0->1`) — the "where does the approximation
  fail?" comparison, from one config.
- **`nrm`** — the nonlocal resonance model (N2/NO/F2): the rung above `lcp`,
  keeping the energy dependence and the nonlocality the LCP discards. It serves
  **both** `ve` (PRA 77 Eq. 28/31/37) and `da` (Eq. 52–54). Takes an optional
  `nrm:` block — `choices: [a, b]` (PRA 77's two discrete-state choices,
  default `[b]`), `n_states` (the Eq. 60 state-sum truncation, default 100, a
  measured value), and `include_background` (default `true`, the paper's
  "nonlocal + background" curve; `false` is its bare "nonlocal" one — a `ve`
  knob, ignored by `da`). One series per choice, keyed `nrm-a:da:ch0` /
  `nrm-b:ve:v0->1`, so `methods: [ti, lcp, nrm]` puts the exact oracle and all
  three approximations on one axis — see
  `examples/f2-da-nrm-vs-lcp-vs-exact.yaml`,
  `examples/n2-ve-nrm-vs-exact.yaml` and
  `docs/physics/nonlocal-resonance-model.md`.

`lcp` and `nrm` both need the preset's grids (neither has an explicit-grid
form). `lcp` needs a `ve`, `da` or `resonance_levels` observable; `nrm` needs a
`ve` or `da` one. `nrm` runs its electronic Hamiltonian on `ti_grid()`'s own
factors, so a `methods: [ti, nrm]` ratio measures the model reduction rather
than two discretisations.

## Observables → config knob → artifact

| observable | how to request | artifact |
|---|---|---|
| **cross sections** | `observables: [{kind: ve\|da\|dr, channels: N}]` | `cross_section.{csv,npz,png}` (keys `{method}:{kind}:{chan}` — e.g. `ti:ve:v0->1`, `nrm-b:da:ch0`. TD adds the extractor as a fourth part, `td:{kind}:{extractor}:{chan}`, so `td:ve:tw:v0->1`, since one run can carry all three extractors) |
| **σ vs time** (TD) | `artifacts.cross_section_vs_time.moments: [t1,...]` | `cross_section_vs_time.{npz,png}` |
| **correlations** (TD) | `artifacts.correlations: true` | `correlations.npz` (raw per-step series) |
| **wavefunctions** (Ψ⁺/Ψ(t)) | `artifacts.wavefunction_snapshots: {ti_energies:[...], td_times:[...]}` | `wavefunction/psi_*.{npz,png}` (per-axis density) |
| **full complex Ψ field** | add `full_field: true` to `wavefunction_snapshots` | `psi=` in the npz + a domain-coloured `*_field.png` (feeds `qscat.viz`) |
| **energy levels + eigenstate wavefunctions** | `artifacts.eigenstates: true` | `eigenstates/eigenstates_{method}_vibrational.{npz,png}` |
| **resonance state** (LCP) | `artifacts.eigenstates: true` (LCP-capable molecule) | `resonance/resonance_lcp_resonance.{npz,png}` — complex pole `E_r−iΓ/2` + electronic eigenfunction at the width peak |
| **LCP scattering states** `ψ_sc(R)` | `wavefunction_snapshots: {ti_energies:[...], full_field: true}` on an `lcp` run | `eigenstates/eigenstates_lcp_scattering.{npz,png}` |
| **quasi-bound resonance levels** (BO/LCP) | `observables: [{kind: resonance_levels, channels: N}]` on `methods: [lcp]` — needs no `energies:` block — or `artifacts.resonance_levels: true` on an existing LCP run | `resonance_levels_{label}.{csv,npz,png}` — complex levels `E_v−iΓ_v/2`, DVR eigenvectors, the `V_d(R)`/`Γ(R)` curve, and a golden-rule comparator column |
| **published reference overlay** | `reference: [{path: ..., format: houfek, channels: [...]}]` | `reference.{csv,npz}` + dashed overlay on `cross_section.png` (keys `ref:...`) |

A `reference:` entry keeps its own energy axis: it is never interpolated onto
the run's energies, so it is written to `reference.{csv,npz}` rather than as
extra columns in `cross_section.csv` (whose rows are the run's own energies).
`path` resolves relative to the config file's own directory, not the current
working directory. `label`, when given, is what the PNG legend shows instead
of the raw `ref:...` series key; with multiple `channels` from one reference,
each gets its own `(chN)` suffix so the legend stays unambiguous. See
`examples/n2-ve-vs-houfek.yaml` for the flagship case — N2 vibrational
excitation overlaid on Karel Houfek's independent published data.

Every run also writes `config.resolved.yaml` (the fully default-filled config)
and `manifest.json` (qscat version, git SHA, timestamp, backend, timings) for
reproducibility. `methods: [ti, td, lcp, nrm]` merges everything into one result
under disjoint `ti:`/`td:`/`lcp:`/`nrm-a:`/`nrm-b:` key prefixes.

## Architecture

- `config.py` — YAML → frozen dataclasses + `validate_config` (actionable errors).
- `presets.py` — the `MoleculePreset` registry: per-molecule grid builders,
  default energies/incident/test-functions, LCP grids. **The one place a new
  molecule is added.**
- `runner.py` — `run_experiment(cfg)` → `ExperimentResult`; `_run_ti`/`_run_td`/
  `_run_lcp`/`_run_nrm` dispatch observables against `qscat.core` on their own grid.
- `artifacts.py` — `write_artifacts(result, cfg, out_dir, timestamp=...)`.

`qscat_run` depends on `qscat` (the library) and **must not import
`validation`** (enforced by `tests/test_no_validation_import.py`). The
per-molecule eMoScat deck definitions are duplicated in `validation/{diatomic,
h2plus}/config.py` (the discretisation tuner's reference) and locked
byte-identical to the presets by `test_diatomic_decks_match_presets` /
`test_h2p_decks_match_presets` — a guarded duplication, deliberate because the
layering forbids a shared import.

## History

The per-molecule *curve/figure drivers* that used to live in
`validation/{diatomic,h2plus}/` were retired into this tool (see
`docs/superpowers/plans/2026-08-15-unified-experiment-observables.md` and the
audit under `docs/superpowers/reviews/`). The plain σ(E) curves genuinely are
configs now — see `examples/figures/`. Four drivers are deliberately kept, and
*not* because qscat-run lacks a reference overlay: it has one (`reference:`, see
the table above, and `examples/n2-ve-vs-houfek.yaml` overlays the very Houfek
data in question). They are kept because:

- `validation/n2/ti_curve.py` is the driver behind the Houfek anchor **gate**
  (`validation/n2/test_ti_curve.py`), and it drives the `projects/`-side exact-2D
  solver — sharing `_build_system` verbatim with `validation/n2/exact2d.py`, so
  the same grid and vibrational states are exercised. qscat-run runs the
  `qscat.core` engine instead, a different code path, so it cannot stand in for
  that gate.
- `validation/diatomic/ve_nrm_figure.py`'s figure needs both `include_background`
  settings (one flag per `nrm` block, so two runs) — not reachable from a
  single config.
- `validation/diatomic/da_figure.py`'s LCP/exact ratio panel and its
  published-reference overlay on a `da` observable are not reachable from
  `qscat_run.artifacts`.
- `validation/diatomic/td_nrm_figures.py`'s seven panels need packet diagnostics
  no config exposes.

`examples/n2-ve-nrm-vs-exact.yaml` is the config form of the rest of it.
