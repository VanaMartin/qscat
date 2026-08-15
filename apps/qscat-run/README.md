# qscat-run

**The single execution surface for qModeling.** One YAML config turns any
molecule or ion into all the common observables across all three methods —
no per-experiment Python. Adding a molecule is a `MoleculePreset` entry
(`qscat_run/presets.py`), never solver code.

```bash
qscat-run init F2 --observables ve,da --methods ti,lcp -o f2.yaml   # scaffold
qscat-run validate f2.yaml                                          # actionable schema check
qscat-run run f2.yaml --output runs/f2                              # solve + write artifacts
```

The production preset decks are `O(10^4–10^6)` unknowns — run those under
Docker (`docker/run.sh <config> <out>`), which provides MUMPS. The committed
`examples/*.yaml` use tiny explicit grids for fast local iteration.

## Registry

| molecule | charge | valid observables | preset variants |
|---|---|---|---|
| `N2` | 0 | `ve` (`da` closed-in-range) | `emoscat` |
| `NO` | 0 | `ve`, `da` | `emoscat` |
| `F2` | 0 | `ve`, `da` | `emoscat` |
| `H2P` | −1 | `dr` | `emoscat`, `proxy` |

## Methods (`methods: [...]`)

- **`ti`** — exact time-independent driven Lippmann–Schwinger solve.
- **`td`** — time-dependent wavepacket propagation (order-3 Padé), with the
  Tannor–Weeks / Dirac / Flux energy extractors (`td.extractors: [tw,delta,flow]`).
- **`lcp`** — the local-complex-potential *approximation* of DA (F2/NO only).
  `methods: [ti, lcp]` overlays the exact oracle and the approximation on one
  `cross_section.png` (keys `ti:da:ch0` vs `lcp:da:ch0`) — the "where does the
  approximation fail?" comparison, from one config.

## Observables → config knob → artifact

| observable | how to request | artifact |
|---|---|---|
| **cross sections** | `observables: [{kind: ve\|da\|dr, channels: N}]` | `cross_section.{csv,npz,png}` (keys `{method}:{kind}:{chan}`) |
| **σ vs time** (TD) | `artifacts.cross_section_vs_time.moments: [t1,...]` | `cross_section_vs_time.{npz,png}` |
| **correlations** (TD) | `artifacts.correlations: true` | `correlations.npz` (raw per-step series) |
| **wavefunctions** (Ψ⁺/Ψ(t)) | `artifacts.wavefunction_snapshots: {ti_energies:[...], td_times:[...]}` | `wavefunction/psi_*.{npz,png}` (per-axis density) |
| **full complex Ψ field** | add `full_field: true` to `wavefunction_snapshots` | `psi=` in the npz + a domain-coloured `*_field.png` (feeds `qscat.viz`) |
| **energy levels + eigenstate wavefunctions** | `artifacts.eigenstates: true` | `eigenstates/eigenstates_{method}_vibrational.{npz,png}` |
| **resonance state** (LCP) | `artifacts.eigenstates: true` (LCP-capable molecule) | `resonance/resonance_lcp_resonance.{npz,png}` — complex pole `E_r−iΓ/2` + electronic eigenfunction at the width peak |
| **LCP scattering states** `ψ_sc(R)` | `wavefunction_snapshots: {ti_energies:[...], full_field: true}` on an `lcp` run | `eigenstates/eigenstates_lcp_scattering.{npz,png}` |

Every run also writes `config.resolved.yaml` (the fully default-filled config)
and `manifest.json` (qscat version, git SHA, timestamp, backend, timings) for
reproducibility. `methods: [ti, td, lcp]` merges everything into one result
under disjoint `ti:`/`td:`/`lcp:` key prefixes.

## Architecture

- `config.py` — YAML → frozen dataclasses + `validate_config` (actionable errors).
- `presets.py` — the `MoleculePreset` registry: per-molecule grid builders,
  default energies/incident/test-functions, LCP grids. **The one place a new
  molecule is added.**
- `runner.py` — `run_experiment(cfg)` → `ExperimentResult`; `_run_ti`/`_run_td`/
  `_run_lcp` dispatch observables against `qscat.core` on their own grid.
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
audit under `docs/superpowers/reviews/`). `validation/n2/ti_curve.py` is kept
— it is the Houfek golden-data gate qscat-run has no reference overlay for.
