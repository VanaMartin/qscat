# Execution-surface audit (2026-08-15)

> The code-quality audit that drove the qscat-run consolidation (PRs #8/#9/#10).
> Findings R1-R8 are tracked + resolved in
> `docs/superpowers/plans/2026-08-15-unified-experiment-observables.md`. Preserved
> here as the reasoning record. (Requested as a Fable-model review; the subagent
> cap was hit, so it was produced inline against the same rubric.)

# Harsh audit: per-molecule/ion execution sprawl & the path to "one config → all observables"

Scope: `validation/diatomic/`, `validation/h2plus/`, `validation/n2/`, cross-checked against
`apps/qscat-run/` and the `qscat.core` return surface. Code-quality only (compactness, readability,
generalization) — physics correctness is out of scope.

## TL;DR (the one finding that reframes everything)

**The generalization you're asking for already exists — it's `apps/qscat-run`.** It is a
config-driven runner with a molecule/ion registry that already unifies neutrals *and* the ion
(`presets.MODELS = {"N2": N2, "NO": NO, "F2": F2, "H2P": H2P}`, presets.py:67), one shared
TI+TD solve/propagation (`runner.run_experiment`), an `ExperimentResult` carrying
cross sections + Ψ⁺ density snapshots + correlations + timings + grids, and per-molecule code
reduced to a `MoleculePreset` data entry (presets.py:259-330). Adding a molecule there is
genuinely a registry entry.

So the real problem is not "we lack the abstraction." It is: **there are two parallel execution
surfaces, and the mature one (`qscat-run`) is starved while the hand-rolled one (`validation/*`)
is fed.** Every per-molecule driver under `validation/` re-implements — badly and inconsistently —
a subset of what `run_experiment` already does. The cleanup is *consolidation into qscat-run*, plus
three genuine capability gaps in qscat-run (LCP, eigenstates, full-field wavefunctions).

Grade for the recent additions (`dr_curves.py`, `td_dr.py`, the diatomic cluster): **works,
readable in isolation, architecturally redundant.** They should mostly not exist.

---

## 1. Duplication & sprawl map

Every one of these drivers is the same skeleton: *build grid → `vibrational_states(grid.grids[1],
model.mu, n_vib, model.v0)` → call one `qscat.core` fn → `FIGURE_DIR.mkdir` + `np.savez` +
`plot_cross_sections`*. The only genuinely molecule-specific bytes are the grid extents and the
energy window — both of which already live in `presets.MoleculePreset`.

| Driver | Lines | The `vibrational_states` boilerplate | The savez+plot boilerplate | Genuinely unique |
|---|---|---|---|---|
| `diatomic/curves.py` (TI VE) | 82 | :56 | :64-78 | grid extents, VPRIMES |
| `diatomic/da_curves.py` (TI DA) | 53 | :27 | :34-49 | E windows, stems |
| `diatomic/lcp_da_curves.py` (LCP) | 62 | :32 | :41-58 | LCP grids + overlay |
| `diatomic/td_da.py` (TD DA) | 138 | :95,:128 | main() prints only | launch grid, WP packets |
| `h2plus/dr.py` (TI DR) | 64 | :43 | :57 | full_grid, N_CHANNELS |
| `h2plus/dr_curves.py` (TI DR fig) | 147 | (via compute_dr) | :69-85,:104-133 | log-log fig |
| `h2plus/td_dr.py` (TD DR) | 153 | :101 | :140-146 | incident r0=800 |
| `n2/ti_curve.py` | 128 | :79 | :112-122 | Houfek overlay |
| `n2/exact2d.py`,`td_exact2d.py`,`experiment.py` | 173+152+208 | each | each | anchor tables |

The `vibrational_states(grid.grids[1], model.mu, n_vib, model.v0)` line is copy-pasted **≥9 times**
verbatim. The `FIGURE_DIR = Path(__file__).resolve().parents[2] / "docs"/"physics"/"figures"`
constant is redefined in `diatomic/curves.py:32`, `h2plus/dr_curves.py:43`, and
`n2/ti_curve.py:50` (three copies of the same path). `_nuclear_index_near` / `_real_index_near` /
`_index_near` — the "nearest real-region DVR index, complex tail masked to inf" helper — exists in
**four** places: `td_da.py:65`, `h2plus/td_dr.py:74`, and twice inside `qscat_run/runner.py`
(:287 as the shared body, :299/:307 as axis wrappers). runner.py:293 even *cites* `td_dr.py` as
the source it was copied from — the duplication is documented, not fixed.

`compute_dr` (h2plus/dr.py:25) and the h2plus `config.py` **free functions** (`full_grid`,
`proxy_grid`, `energy_grid`, module-level `N_CHANNELS`/`E_LO`/`E_HI`/`E_STEP`) are a *different
config shape entirely* from the diatomic `MoleculeConfig` frozen dataclass (config.py:29-102) —
which is itself a different shape from `qscat_run.MoleculePreset`. **Three incompatible
per-molecule config schemas for the same concept.**

---

## 2. Compactness / readability findings (ranked)

**F1 — SEVERE — the whole `validation/{diatomic,h2plus}` driver layer duplicates `qscat-run`.**
`da_curves.compute_da_curve`, `curves.compute_ti_curve`, `dr.compute_dr`, `td_da.compute_td_da_three_way`,
`td_dr.compute_td_dr` are each a 1:1 re-derivation of a `run_experiment` observable dispatch
(runner.py:225-248 TI, :445-513 TD) with a bespoke grid getter. They exist only because the
`validation/` figures predate qscat-run and were never migrated. This is the sprawl. Everything
below is a symptom of it.

**F2 — SEVERE — three config schemas for one concept.** diatomic `MoleculeConfig` dataclass
(config.py:29) vs h2plus free functions (config.py:32-79) vs `qscat_run.MoleculePreset`
(presets.py:85-115). The ion diverging into free functions is not justified by anything physical —
`H2P` is already a first-class `MoleculePreset` (presets.py:306). The diatomic dataclass's
`da_grid`/`lcp_elec_grids`/`lcp_nuclear_grid` grid-builder methods (config.py:61-102) are grid
construction masquerading as config; they belong next to `presets.resolve_grid`.

**F3 — HIGH — the observables you asked for cannot be extracted at all from these drivers.**
Grep is unambiguous: the `validation/*` drivers emit **cross sections only**. Not one of them can
dump an eigenstate wavefunction, a resonance state, or a scattering field:
  - Energy levels: `eps` is computed everywhere and *saved once* (curves.py:68 stuffs it in the
    npz); the eigenfunctions `chi` are computed and **discarded** everywhere.
  - Scattering wavefunction Ψ⁺: `qscat.core.ve_cross_section` grew a `return_wavefunction=True`
    path (driven.py:159,225-226) — but **only `qscat_run/runner.py:255` ever calls it.** No
    `validation/` driver does. `da_cross_section`/`dr_cross_section` have no such flag at all.
  - Resonance states: `qscat.ecs.find_resonance_pole` exists but no driver surfaces the pole
    eigenpair as an artifact.
  - Nuclear density: a real observable exists (`projects/n2_2d_cross_section/nuclear_density.py:69`)
    and is wired into **nothing**.

**F4 — HIGH — `qscat-run`'s own wavefunction artifact is density-marginals only, so it can't feed
`qscat.viz`.** `WavefunctionSnapshot` (runner.py:110-125) stores `rho_r`/`rho_R` (|Ψ|² summed onto
each axis, `_project_density` runner.py:189-198) — 1-D marginals. The `qscat.viz` domain-coloring /
animation stack you just built needs the **full complex 2-D field** on the tensor grid. So today
"produce the wavefunction" and "animate the wavefunction" are two disconnected code paths. The
snapshot throws away exactly the phase information viz needs.

**F5 — MEDIUM — figure/IO boilerplate is copy-pasted, not a helper.** The
`mkdir → np.savez → plot_cross_sections(..., title=f"{name} ...", path=stem.png)` block recurs in
curves.py:64-78, da_curves.py:34-49, lcp_da_curves.py:41-58, dr_curves.py:69-85. `qscat-run`
already centralizes this in `artifacts.write_artifacts` — the validation copies are the redundant
ones. Three redefinitions of `FIGURE_DIR` (see §1).

**F6 — MEDIUM — `dr_curves.py` (just added) is a fourth copy of the pattern for a third config
schema.** It's clean *locally* (reuses `compute_dr`, closes over one `_sweep`), but it exists only
because h2plus has no artifact pipeline. Under qscat-run it is a `h2p-dr-ti.yaml` example plus a
generic figure step — ~140 lines → ~20 lines of YAML.

**F7 — LOW — `n2/` is the oldest and most bespoke, and reaches into `projects/`.** `ti_curve.py`
imports four `projects.n2_2d_cross_section.*` internals + a `projects.n2_ti_cross_section`
vibrational solver (ti_curve.py:38-42) rather than the promoted `qscat.core`. `exact2d.py`,
`td_exact2d.py`, `experiment.py` are anchor-table harnesses — legitimately test infrastructure, but
their *curve/figure* halves overlap qscat-run. Keep the anchor gates; migrate the curve generation.

---

## 3. The missing generalization — concrete target

You do **not** need a new abstraction. You need to (a) fill three gaps in `qscat-run` and (b) delete
the `validation/*` curve drivers in favor of committed example configs + a generic figure step.

### 3a. Add LCP as a method (the one real architectural gap)
`VALID_METHODS = {"ti","td"}` (config.py:68) — add `"lcp"`. A `_run_lcp` in runner.py mirroring
`_run_ti`, calling `qscat.core.lcp.local_complex_potential` + `lcp_da_cross_section`
(the exact bodies now in `validation/diatomic/lcp_da_curves.py:30-34`). LCP needs two extra
preset fields — the two ECS pole-matching angles + the fine nuclear deck — which already exist as
`MoleculeConfig.lcp_angle_a/b` + `lcp_nuclear_grid` (config.py:55-56,95). Fold them into
`MoleculePreset`. Then `methods: [ti, lcp]` overlays exact vs approximation from one config, which
is precisely what `lcp_da_curves.main` hand-wires today (lcp_da_curves.py:46-58).

### 3b. Promote "energy levels + eigenstates" and "resonance states" to observables
These are cheap and already computed. Proposed additions to `ExperimentResult`:
```python
@dataclass(frozen=True)
class EigenStates:                 # vibrational OR resonance
    kind: str                      # "vibrational" | "resonance"
    energies: NDArray[float64]     # eps (real) or complex poles
    states: NDArray                # chi columns on grids[1], or Ψ_res on the tensor grid
    grid_axis: NDArray[float64]    # real_points of the relevant axis
```
Driver work is nil — `_run_ti` already holds `eps, chi` (runner.py:221); just attach them when
`artifacts.eigenstates` is requested. Resonance states come from the existing
`qscat.ecs.find_resonance_pole` two-spectrum matcher.

### 3c. Store the full-field wavefunction, not just marginals (unblocks viz)
Widen `WavefunctionSnapshot` to optionally carry the complex `psi` reshaped to `tg.shape` (TI Ψ⁺ is
already in hand at runner.py:265; TD `result.snapshots` likely carries the field pre-projection).
Then `artifacts` can emit an `.npz` that `qscat.viz.animate_wavefunction` /
`plot_wavefunction_2d` consumes directly — closing F4. Keep the cheap marginals as the default;
full field is opt-in (memory).

### 3d. The unified surface, end state
`run_experiment(load_config(yaml))` → `ExperimentResult{cross_sections, eigenstates, wavefunctions
(full field), cross_section_vs_time, correlations}` for any `methods ⊆ {ti, td, lcp}` and any
`molecule ∈ MODELS` (neutral or ion). Per-molecule dedicated code = one `MoleculePreset` +
committed `*.yaml`. `validation/*` keeps ONLY the golden-data gates (n2 anchors vs Houfek,
`test_*`), which assert numbers and legitimately need bespoke reference wiring.

**Collapses into qscat-run (delete after migration):** `diatomic/curves.py`, `diatomic/da_curves.py`,
`diatomic/lcp_da_curves.py`, `diatomic/td_da.py`, `h2plus/dr.py`, `h2plus/dr_curves.py`,
`h2plus/td_dr.py`, `n2/ti_curve.py`, and the diatomic `MoleculeConfig` + h2plus config free
functions. **Stays bespoke:** `n2/exact2d.py`/`td_exact2d.py`/`experiment.py`/`test_anchors.py`
(golden-data anchor gates), `validation/tuning/*` (calibration), all `test_*`.

---

## 4. Prioritized recommendations

| # | Item | Sev | Effort | Files | Type |
|---|---|---|---|---|---|
| R1 | Add `lcp` as a `qscat-run` method + fold LCP angles/nuclear-deck into `MoleculePreset` | HIGH | M | `qscat_run/{config,presets,runner,artifacts}.py` | capability |
| R2 | Add `eigenstates` observable (vibrational eps+chi, resonance poles) to `ExperimentResult`/artifacts | HIGH | S–M | `qscat_run/{runner,artifacts,config}.py`, maybe `qscat.core` resonance helper | capability |
| R3 | Store full complex Ψ field in `WavefunctionSnapshot`; emit viz-ready npz; wire to `qscat.viz` | HIGH | M | `qscat_run/{runner,artifacts}.py`, `qscat.viz` | capability |
| R4 | Migrate diatomic + h2plus curve drivers to committed `*.yaml` + a generic figure step; **delete** the 7 drivers + 2 rogue config schemas | SEVERE (debt) | M–L | `validation/diatomic/*`, `validation/h2plus/*`, `examples/*.yaml` | refactor |
| R5 | Hoist the `_index_near` real-region-DVR helper into `qscat.dvr` (kill 4 copies) | MED | S | `qscat.dvr`, remove from `runner.py`/`td_da.py`/`td_dr.py` | refactor |
| R6 | Add `da_cross_section`/`dr_cross_section` a `return_wavefunction` flag to match `ve_cross_section` (enables R3 for DA/DR) | MED | S | `qscat.core.dissociation` | capability |
| R7 | Migrate `n2/ti_curve.py` curve-gen to qscat-run; keep anchor gates; drop the `projects.*` reach-through | LOW | M | `validation/n2/*` | refactor |
| R8 | Single `FIGURE_DIR` constant; wire `projects/n2_2d/nuclear_density.py` as a real observable or delete it | LOW | S | validation, `qscat_run/artifacts.py` | refactor |

**Sequencing:** R5→R2→R6→R3→R1 are additive and independently shippable (each a small PR that grows
qscat-run without touching validation). R4/R7 are the payoff — do them *after* R1/R3 land, so the
deleted drivers have a full-featured replacement. Do not delete a driver until its config +
figure step reproduces its committed figure byte-for-byte (the h2plus/n2 committed PNGs are the
acceptance test).

**The harsh version:** you already built the right thing once (`qscat-run`) and then kept writing the
wrong thing (`validation/*` scripts) beside it — most recently `dr_curves.py`. Stop feeding the
hand-rolled layer. Every future "run molecule X, get observable Y" should be a YAML, and the only
per-molecule code anyone writes should be a `MoleculePreset` entry.
