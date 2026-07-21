# N₂ electron-impact ²Π_g shape resonance: fixed-R pole search

**Location:** `projects/n2_resonance/` (toy model: `potential.py`, `grid_n2.py`,
`pole.py`), `qscat.ecs.find_resonance_pole` (the promoted general matcher),
`validation/n2/` (`model.py`, `resonance.py`, `experiment.py` — the benchmark
harness, Group B).
**Origin:** LCP (local complex potential) model extracted from
`reference/eMoScat/input/experimental/N2-model.json`; see
`.superpowers/sdd/n2-lcp-model-extraction.md`. Builds on the FEM-DVR-ECS grid
(`docs/physics/femdvr-ecs.md`).
**Units:** atomic units throughout (energy in Hartree, length in Bohr).

## Physical picture

Electron scattering off N₂ near ~2 eV is dominated by the ²Π_g shape
resonance — a metastable N₂⁻ state (the well-known "boomerang" resonance)
that decays back to N₂ + e⁻. In the local complex potential model, the
electron sees a fixed-nuclei ("clamped-R") effective potential
`V_eff_el(r, R) = v_int(r, R) + l(l+1)/(2r²)` (`l = 2`), and the resonance
shows up as a complex pole `E(R) = E_res(R) - i*Γ(R)/2` of the associated
scattering problem: `E_res(R)` is the resonance position (vertical
attachment energy) and `Γ(R)` the autodetachment width, both functions of
the N₂ bond length `R`.

`v0(R)`, `lam(R)`, `v_int(r, R)`, `v_eff_el(r, R)` are closed-form (Morse
neutral curve, Gaussian well in `r` scaled by an R-dependent strength — see
`validation/n2/model.py`, the validated "single source" these formulas are
kept in lockstep with). `E_res(R)`/`Γ(R)` are *not* closed-form: they require
diagonalizing the fixed-R electronic Hamiltonian
`H_el(R) = T + diag(V_eff_el(r, R))` on an ECS-rotated radial grid and
locating the pole in its spectrum.

## Method: two-angle ECS pole matching

Exterior complex scaling (`docs/physics/femdvr-ecs.md`) turns the divergent
continuum wavefunctions of the unscaled problem into decaying ones by
rotating the radial coordinate beyond a pivot `R0` into the complex plane at
angle `theta`. A resonance pole becomes (nearly) `theta`-independent once
`theta` is large enough to fully expose it, while the *discretized* continuum
eigenvalues keep rotating with `theta`.

This gives a direct way to find the pole without a `theta`-stabilization
scan: diagonalize the *same* `H_el(R)` on two grids built at two different
(but both "large enough") ECS angles, restrict both spectra to a search
window, and find the pair of eigenvalues — one from each spectrum — that
agree most closely. That pair is the angle-stable pole; nearby discretized
continuum eigenvalues, having rotated by different amounts on the two grids,
do not match this closely. Formalized in
`qscat.ecs.find_resonance_pole(eigs_a, eigs_b, window)`:

```
E_pole  = 0.5 * (ea + eb)      # ea, eb: the closest-matching pair
residual = |ea - eb|            # small (<< Gamma) => genuine angle-stable pole
E_res   = Re(E_pole)
Gamma   = max(0, -2 * Im(E_pole))
```

`projects/n2_resonance/pole.find_pole(R, grid_a, grid_b, window)` assembles
`H_el(R)` on each grid (`electronic_hamiltonian`, `qscat.dvr.hamiltonian` +
`eigen`) and calls `find_resonance_pole` on the resulting spectra; it is the
thin, physics-specific caller of the promoted general matcher.
`validation/n2/resonance.py` mirrors the same grid construction and calls
`find_pole` equivalent logic against `model.v_eff_el` independently, so the
harness (Group B) does not import from `projects/`.

## R0 result (validated)

At the N₂ equilibrium bond length `R0 = 2.01943` Bohr, with grids built at
`theta = 35°` and `44°` (`grid_n2.n2_electronic_grid`, `r_pivot=10`,
`n_real=n_complex=8`, `quadrature=8`) and search window
`Re in [0.04, 0.16]` Ha, `Im in [-0.05, 0]` Ha:

```
E_res(R0) = 2.445 eV   Gamma(R0) = 0.455 eV   residual ~ 3e-6 Ha
```

Both lie inside the literature plausibility bands
(`E_res ∈ [2.3, 2.5]` eV, `Gamma ∈ [0.35, 0.55]` eV — Schulz; Berman/Domcke;
see `validation/n2/reference.py`), and the match residual (~1e-6 to 1e-5 Ha)
is far smaller than `Gamma` (~0.017 Ha), the signature of a genuine
angle-stable pole rather than an accidental near-miss between rotating
continuum states. The pole is also resolution-stable: a coarser grid
(`n_real=n_complex=6`) reproduces `E_pole` to within a few percent
(`projects/n2_resonance/test_pole.py`, V1/V2).

## R-scan: `E_res(R)`, `Gamma(R)`, `V_d(R)`

`projects/n2_resonance/pole.resonance_curve(R_grid, grid_a, grid_b)` traces
the pole across a range of bond lengths (e.g. `R ∈ [1.6, 3.0]` Bohr,
`np.linspace(1.6, 3.0, 15)`), returning `(E_res, Gamma, V_d)` with
`V_d(R) = v0(R) + E_res(R)` — the anion/resonance curve measured from the
neutral Morse potential.

A single fixed search window cannot span this whole range: at short R the
pole is a genuine complex shape resonance (as at `R0`), but as `R` stretches
past roughly 2.3–2.4 Bohr, `Gamma(R)` narrows to zero and the state crosses
into a real, angle-independent **bound** state — the standard
dissociative-electron-attachment picture, where the anion curve dips below
the neutral one at large R. `resonance_curve` handles this with a
**continuation walk**: starting from the `R_grid` index nearest the R0
region (where the fixed seed window `Re ∈ [0.04, 0.16]` Ha is valid), it
walks outward in both directions, each step recentering the search window on
a `±0.05` Ha box around the *previous* step's matched pole. This tracks the
pole smoothly through the resonance-to-bound crossing without mode-hopping
onto a neighboring (continuum or other-pole) eigenvalue.

Over `R ∈ [1.6, 3.0]` Bohr (validated in `projects/n2_resonance/test_curve.py`,
V3):

- `E_res(R)` decreases smoothly and monotonically from `~5.7 eV` (R=1.6) to
  `~-2.9 eV` (R=3.0) — a well-behaved, mode-hop-free curve (curvature
  `|ΔΔE_res|` stays a small multiple of the local step throughout).
- `Gamma(R)` decreases monotonically from `~3.2 eV` (R=1.6) to exactly 0 for
  `R ≳ 2.4` Bohr, matching the expected resonance-to-bound-state closure.
- `Gamma(R) >= 0` everywhere by construction (`Gamma = max(0, -2*Im(E_pole))`).

## Model caveat: neutral N₂ Morse is a model potential, not a spectroscopic fit

`v0(R)` (`validation/n2/model.py`, `projects/n2_resonance/potential.py`) is
eMoScat's neutral N₂ Morse curve, extracted as-is from
`reference/eMoScat/input/experimental/N2-model.json`. Its dissociation
energy `D_0 = 0.75102` Ha (=~20.4 eV) is =~2x real N₂'s actual dissociation
energy (=~9.8 eV). Consequently the model's neutral vibrational spacing
(`omega_e =~ 0.0125` Ha analytic; FEM-DVR `eps1-eps0 =~ 0.0124` Ha, see
`projects/n2_ti_cross_section/test_vibrational.py`) is =~16% larger than
real N₂'s (0.01074 Ha / 2358 cm⁻¹).

This is a deliberate, accepted property of the model (maintainer decision),
not a bug: the FEM-DVR-ECS vibrational solver is verified correct against
the *analytic* Morse spectrum of this same potential (residuals ~1e-14 Ha),
so the gap lives entirely in the model's `D_0`/`alpha_0` choice, not in the
numerics. Meanwhile the model's resonance parameters `E_res(R0)`/`Γ(R0)`
(=~2.44 eV / 0.46 eV, see the "R0 result" section above) DO match real N₂
electron-scattering data — only the *neutral* vibrational ladder departs
from real N₂ spectroscopy. This model-vs-reality gap is folded into the
LCP-vs-Houfek-2D differences seen in the time-independent (TI) cross-section
benchmark (`projects/n2_ti_cross_section/`; method, resonance-peak agreement,
and the two structural LCP limitations are documented in
`docs/physics/n2-cross-section.md`).

## Simplifications / out of scope

- **Diagonal-potential DVR approximation** (inherited from `qscat.dvr`): `V`
  is represented by its values on the grid points, not exact quadrature of
  `⟨i|V|j⟩`. Adequate here because `v_eff_el` is smooth away from `r=0`
  (the centrifugal term diverges there, but the grid does not place a point
  at the origin).
- **Fixed, single ECS angle per grid, hand-picked window/angle pair**: the
  two-angle match at `(35°, 44°)` with the `Re ∈ [0.04, 0.16]`,
  `Im ∈ [-0.05, 0]` Ha window is validated at `R0` and used as the
  continuation seed for the R-scan; it is not an automated
  `theta`-stabilization scan (which would locate the pole without any prior
  window guess). Automating that is future work, as noted in
  `docs/physics/femdvr-ecs.md`'s limitations.
- **Fixed reduced-mass / no nuclear dynamics here**: this is the *electronic*
  (fixed-R) problem only. `E_res(R)`/`Γ(R)`/`V_d(R)` are inputs to a
  nuclear-motion (vibrational excitation, dissociative attachment) solve; the
  time-independent version of that solve is now implemented and validated
  (`docs/physics/n2-cross-section.md`, `validation/n2/README.md` Group C5 —
  4/4 gated anchors PASS); the time-dependent (wavepacket propagation)
  version remains out of scope / PENDING (Group D).
- **Continuation window widths are hand-tuned** (`re_half_width =
  im_half_width = 0.05` Ha in `resonance_curve`), chosen empirically to
  track the pole across the observed drift rate without catching a
  neighboring branch (branches sit ~0.02–0.05 Ha apart on this grid); a
  different `R_grid` spacing or a much wider R range could require
  re-tuning.

## Validation

- `projects/n2_resonance/test_pole.py`:
  - **V1** — `E_res(R0) ∈ [2.3, 2.5]` eV, `Gamma(R0) ∈ [0.35, 0.55]` eV.
  - **V2** — match residual `< 1e-3` Ha (angle-stability); resolution
    stability (coarser grid agrees to `<5%`).
- `projects/n2_resonance/test_curve.py`:
  - **V3** — `Gamma(R) >= 0` everywhere; `E_res(R)` monotonically
    non-increasing (no mode-hop) across the scan; curvature of `E_res(R)`
    bounded by a small multiple of the local step; `Gamma(R)` monotonically
    non-increasing to 0.
- `libs/qscat/tests/test_find_resonance_pole.py`: the promoted matcher on
  synthetic spectra (a shared "pole" eigenvalue plus differing "continuum"
  eigenvalues; empty-window `ValueError`).
- `validation/n2/experiment.py` Group B (`B1`): `E_res(R0)` computed
  independently via `validation/n2/resonance.py` against the literature
  window — **PASS**.
- The nuclear-motion time-independent cross-section solve built on top of
  `E_res(R)`/`Γ(R)` (Group C5) is documented separately in
  `docs/physics/n2-cross-section.md`.
