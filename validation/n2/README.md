# N₂ LCP benchmark

Validation harness for the electron–N₂ local complex potential (LCP) model: elastic and
vibrationally-inelastic scattering off the N₂ ²Π_g shape resonance (the well-known ~2.3 eV
"boomerang" resonance seen in electron–N₂ collisions).

## Physics

`model.py` provides the closed-form pieces of the model, extracted from
`reference/eMoScat/input/experimental/N2-model.json` (see `config.json`):

- `v0(R)` — neutral N₂ ground-state potential (Morse form), minimum `-D_0` at `R_0`.
- `lam(R)` — R-dependent interaction strength λ(R), with `λ(R_c) == λ_c`.
- `v_int(r, R)` — electron–molecule interaction potential, a Gaussian well in the electron
  coordinate `r` scaled by `λ(R)`.
- `v_eff_el(r, R)` — fixed-R electronic effective potential, `v_int` plus the
  `l(l+1)/2r²` centrifugal term (`l` = 2, ²Π_g symmetry).

The resonance position `E_res(R)` and width `Γ(R)` are **not** closed-form — they are poles
of the fixed-nuclei electronic scattering problem, found via an exterior-complex-scaling
(ECS) eigensolver (`resonance.py`, built on `qscat.dvr`/`qscat.ecs`; see
`docs/physics/n2-resonance.md` for the method and the `projects/n2_resonance/` toy model it
was validated against). Likewise the vibrational-excitation cross sections `σ(v=0→v')`
require a time-independent (TI) nuclear scattering solve driven by `E_res(R)`/`Γ(R)`, and
time-dependent (TD) cross sections require propagating the nuclear wavepacket through the
same complex potential surface — those two are still pending. All atomic units throughout
(Hartree, bohr, bohr²).

## Check groups

`experiment.run_checks()` returns a flat list of `(group, name, status, detail)` tuples,
`status ∈ {PASS, PENDING, FAIL}`:

- **Group A — model** (`model.model_checks()`): closed-form identities on `v0`/`lam`/`v_int`
  (minimum location/depth, asymptote, λ at `R_c`, well sign). **Green now.**
- **Group C1–C4 — data integrity** (`loader.integrity_checks()`): shape, monotonic energy
  grid, non-negativity, ordered channel thresholds of the golden cross-section data.
  **Green now.**
- **Group B — resonance position** (`resonance.e_res_at_R0()`): computes `E_res(R_0)` via
  two-angle ECS pole matching (`qscat.ecs.find_resonance_pole`) and checks it against the
  literature window `reference.LITERATURE["E_res_eV"]` (2.3–2.5 eV, Schulz; Berman/Domcke).
  **Green now** (`E_res(R_0) ≈ 2.445 eV`).
- **Group C5 — cross-section anchors** (`cross_section.compute_anchor_results()`): computes
  `σ(E, v=0→v')` at the 6 `reference.ANCHOR_COORDS` via the validated time-independent (TI)
  solver (`projects/n2_ti_cross_section`) and compares against the Houfek golden data. Each
  anchor is classified GENERALLY, not by hardcoded coordinate: **GATED** (a VE channel,
  `channel>=1`, clear of its own threshold by `reference.ANCHOR_MARGIN_HA`) gets a real
  **PASS**/FAIL at `1/reference.ANCHOR_FACTOR <= ratio <= reference.ANCHOR_FACTOR`; elastic
  (`channel==0`) or near-threshold anchors are **DOCUMENTED-LIMITED** — reported as `NOTE`
  with the LCP mechanism named, never failing the harness (see `cross_section.py`'s module
  docstring and `docs/physics/n2-cross-section.md`). **4/4 GATED anchors PASS, 2 anchors
  NOTE** with the current model.
- **Group D — time-dependent**: cross sections from a time-dependent wavepacket
  propagation, cross-checked against the TI solver. **PENDING** — needs the TD model.

`experiment.main()` prints the table and returns exit code `0` unless any check is `FAIL`
(`PENDING`/`NOTE` checks never fail the run).

## Data provenance

`data/CSVE.V00.J00` — 400×32 table of electron–N₂ vibrational-excitation cross sections
(Hartree energy grid; columns = elastic + v=0→1..30, atomic units bohr²), computed by Karel
Houfek via a time-independent method, external to this repo. See `data/MANIFEST.md` for the
full column layout and units. Anchor reference values in `reference.anchors()` are always
resolved by looking up the golden data file — never hardcoded — so they track the data file
if it is ever regenerated.

## Running

```bash
# Tests (loader + model + harness)
uv run pytest validation/n2

# Print the benchmark table locally
uv run python -m validation.n2.experiment

# Same, inside the CPU runtime Docker image
docker/run-n2.sh
```
