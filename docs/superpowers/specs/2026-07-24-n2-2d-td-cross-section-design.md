# N₂ Exact 2-D Time-Dependent VE Cross-Section — Design Spec (sub-project #7)

**Date:** 2026-07-24
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — implementation pending
**Lifecycle:** `qm-method-lifecycle` stages 1–3. Builds on #5 (`qscat.dvr`/`qscat.linalg`,
the N-dimensional sparse Hamiltonian), #6 (the exact 2-D `H_2D`, grids, channel functions,
and the exact σ_TI that is the differential oracle here), and #4 (the 1-D TD route this is
the 2-D twin of).

## Context

Sub-project #6 solved the exact 2-D electron–N₂ problem in the **energy domain** (driven
Lippmann–Schwinger, one sparse solve per energy) and certified it against Houfek's data to
5–6 figures — it is now the oracle. Sub-project #7 solves the **same** exact problem in the
**time domain**: launch an electron wavepacket, propagate it under `H_2D`, and read the
cross section off the correlation functions.

The physics parallel to #3→#4 (1-D TI→TD) is exact, and the differential oracle is the same
identity: by the Laplace/Fourier relation between the resolvent and the propagator,

```
S_TD(E) = (1/i) ∫₀^∞ e^{i E_tot t} ⟨Φ_{v'}| e^{-i H_res t} |Ψ_i⟩ dt  =  ⟨Φ_{v'}|(E_tot − H)⁻¹|Ψ_i⟩  =  S_TI(E)
```

so **TD must converge to #6's exact TI** — an *exact* differential oracle, not a loose
cross-model bound. (The `(1/i)` single-state form above states the *principle* — resolvent =
time integral of the propagator. In practice the incident state is a wavepacket, a
superposition over energies, so the single-energy `S(E)` is recovered by the Tannor–Weeks
deconvolution below, which divides out the wavepacket's own spectral content; the two are the
same object, one written for a fixed energy and one for a wavepacket.)

**Martin's stated goal:** *observe the wavefunction and cross-section formation from the
correlation functions.* This makes #7 a **numeric-output-first** project: the deliverable is
the sampled time-series data of the dynamics, from which the cross section forms; figures are
the visual layer on top, and animation is an explicit later extension.

**What TD buys over TI:** one propagation yields the **entire σ(E) curve** (the energy
transform is cheap once `c(t)` is stored), so the boomerang oscillation structure across the
resonance becomes directly observable.

## Method (atomic units, bohr²)

Reuses #6 verbatim: the tensor FEM-DVR-ECS grid, `H_2D = build_h2d(tgrid)`, the neutral
vibrational states `(ε_v, χ_v)`, and the energy-normalized radial channel function
`F_{E,l}` (`riccati_bessel_en`). Then:

- **Initial state** (`eMoScat TimeDependentModel2d.cpp:338-384`):
  `Ψ(r,R,0) = g(r) ⊗ χ_{v_init}(R)`, with `g` an inward-moving L²-normalized Gaussian
  wavepacket `g(r) = (2π σ²)^{-1/4} e^{-(r-r₀)²/4σ²} e^{i p₀ r}` (reference: `r₀≈45`,
  `p₀≈−0.35`, `σ≈6`). **No driving term** — the wavepacket carries the incident flux, and
  the ECS tail absorbs everything that leaves. `g` is masked to the unscaled region and
  converted to DVR coefficients (`·sqrt(w_r)`, same convention as #6).

- **Crank-Nicolson propagation** under the time-independent, complex-symmetric non-Hermitian
  `H_2D`: `(I + i H_2D dt/2) Ψ_{n+1} = (I − i H_2D dt/2) Ψ_n`. `H_2D` is constant ⇒ **factor
  `A = I + i H_2D dt/2` once** (sparse LU), reuse each step (one back-substitution per step).
  The `−2θ`-rotated ECS continuum makes `‖Ψ‖` decay as outgoing flux is absorbed; propagate
  until the norm is negligible.

- **Correlation functions** (fine cadence — the core numeric time series):
  `c_{v'}(t_n) = c_product(Φ_{v'}, Ψ(t_n))`, where `Φ_{v'} = g_out(r) ⊗ χ_{v'}(R)`
  (`g_out` an outgoing Gaussian; reference `r₀≈75`, `p₀≈+0.5`, `σ≈4`), masked, c-product
  (no conjugate — ECS). Recorded every step (or a small fixed period) for all tracked `v'`.

- **Tannor–Weeks energy transform → σ(E)** (`eMoScat TestFunction2d.cpp:298-307`):
  `S_{v→v'}(E) = [2π · η*_{v'}(E) · η_v^{in}(E)]^{-1} ∫₀^∞ e^{i E_tot t} c_{v'}(t) dt`,
  `E_tot = E + ε_{v_init}`, with the deconvolution factors
  `η_v^{in}(E) = ⟨g | F_{E,l}⟩` (incident Gaussian on the energy-normalized regular free
  function) and `η_{v'}(E) = ⟨g_out | outgoing radial fn at k'⟩`; a composite quadrature in
  time. Then `σ_{v→v'}(E) = π|S − δ|²/(2E)` (elastic subtracts 1). This is the **S-matrix**
  form of #6's **T-matrix** formula: with `S = 1 − 2πiT`, `|S − δ|² = 4π²|T|²`, so
  `π|S − δ|²/(2E) = 4π³|T|²/(2E)` — identical to #6's `σ = 4π³|T|²/(2E)`, inheriting the same
  reference-fixed partial-wave convention noted in `docs/physics/n2-2d-cross-section.md`.
  Closed channels (`E_tot − ε_{v'} ≤ 0`) give 0.

## Numeric outputs (the primary deliverable), on two cadences

- **Fine:** `c_{v'}(t)` — complex array `[n_t × n_channels]`. The raw material of the
  transform and the literal "formation from the correlation functions."
- **Coarse (static points or a period):** at a configurable schedule of times, the nuclear
  density `ρ(R,t_k) = Σ_r |Ψ(r,R,t_k)|²` (unscaled region), the electronic density
  `ρ(r,t_k)`, and the norm `‖Ψ(t_k)‖`. The full `Ψ(r,R,t_k)` retained only at a few key
  times (it is N≈143k complex).
- **Derived:** `σ(E)` on a dense energy grid; the norm-decay curve `‖Ψ(t)‖` vs `t` (the
  resonance depleting).

All are inspectable/exportable numeric arrays returned from functions (and/or saved as
`.npz`), independent of any plotting.

## Interface

```
projects/n2_2d_td_cross_section/
├── propagator.py          thin re-export of qscat.evolution.make_sparse_cn_stepper
├── wavepacket.py          gaussian_wavepacket(grid, r0, p0, sigma) -> masked DVR coeffs
├── td_propagation.py      propagate H_2D from Psi(0); sample c(t) (fine) + densities (coarse)
├── correlation.py         correlation functions c_{v'}(t) and the Tannor-Weeks eta factors
├── td_cross_section.py    energy transform c(t) -> S(E) -> sigma(E)  (the crux)
├── observation.py         numeric-output layer + PNG snapshots (rho(R,t), rho(r,t), norm, c, sigma)
└── convergence.py         box/dt/T sizing; the TD working grid
```

The generic **sparse Crank-Nicolson propagator** promotes to
`qscat.evolution.make_sparse_cn_stepper(H_sparse, dt)` — the sparse sibling of the existing
dense `make_cn_stepper`, factoring `A = I + iH dt/2` once with `qscat.linalg.SparseLU` and
returning a `stepper(psi)` that does one back-substitution. Reuses
`projects.n2_2d_cross_section.{hamiltonian2d, channels, convergence}`,
`projects.n2_ti_cross_section.vibrational`, and #6's `validation.n2.exact2d` (the σ_TI oracle).

## Validation

**V1 — sparse CN vs dense.** `make_sparse_cn_stepper` matches the dense `make_cn_stepper` (and
`scipy.expm` for small Hermitian H) to round-off; norm-conserving for Hermitian H,
norm-decaying for an absorbing (ECS) H; O(dt³) local error.

**V2 (crux) — TD converges to the exact TI oracle.** At the gated anchors and a set of extra
(E, v'), σ_TD from the energy transform agrees with #6's exact σ_TI within a tight tolerance
(limited only by dt / finite-T / η-deconvolution discretization). This is the exact
differential oracle.

**V3 — the full σ(E) curve.** From one propagation, σ_TD(E) on a dense grid overlaid on the
exact σ_TI(E) across the **usable** window; report where the two agree and where the
Tannor–Weeks η-deconvolution noise (small `η_in(E)` near the window edges) makes the TD
result unreliable — plot the honest usable window, not the noise.

**V4 — propagation convergence.** Halving dt and/or lengthening T changes σ negligibly;
`‖Ψ(T)‖` has decayed to a small fraction of `‖Ψ(0)‖` (resonance depleted, wavepacket gone).

**V5 — the numeric outputs are sound.** `ρ(R,t)`, `ρ(r,t)` real and ≥0; `‖Ψ(t)‖` monotone
non-increasing under the absorbing contour; `c_{v'}(t)` finite; the energy transform of the
stored `c(t)` reproduces σ(E) (i.e. the saved numeric time series is self-consistent with the
reported curve).

**Harness group F.** Add a group **F** to `validation/n2/experiment.py`, guarded like C5/D/E:
σ_TD at the gated anchors vs #6's σ_TI (tight rtol) and vs Houfek (factor-3), degrading to a
labeled FAIL on error. Because TD needs the large box, group F likely runs a **reduced**
configuration (shorter box/T, documented looser tolerance) to stay CPU/Docker-affordable,
while the converged result lives in the test suite / a one-off study — decided by measurement,
and stated in the docs.

## Out of scope

- The DA channel; higher partial waves; rotation.
- **Optimizing the sparse LU / back-substitution.** This is explicitly the eventual next
  step, NOT part of #7. eMoScat's Intel **MKL PARDISO** did all of N2/NO/F2 in **under an
  hour**; our `scipy` SuperLU stack is the validated toy, and the optimize-in-Rust stage (a
  PyO3 kernel and/or a PARDISO/MKL binding, plus the deferred MKL x86 base image) is where
  that bar gets met. #7 validates the method; performance is a separate lifecycle stage.
- Chebyshev / higher-order Padé propagation (order-1 CN is the validated choice, as in #4).
- The animation (a later extension; #7 stores the frames so it is a render pass over existing
  data).
- Promoting anything beyond the sparse CN stepper into `qscat`.

## Verification

- `uv run pytest projects/n2_2d_td_cross_section libs/qscat validation/n2 -q` → all pass.
- `uv run mypy libs/qscat` → 0 errors. `uv run ruff check .` → clean.
- TD σ matches the exact TI σ at the gated anchors within the stated tolerance; the full
  σ(E) curve overlays σ_TI across the usable window; propagation converged; norm decayed.
- `uv run python -m validation.n2.experiment` → group F PASS/FAIL/NOTE rows; exit 0; no
  regression of the existing 23 PASS / 0 PENDING / 4 NOTE / 0 FAIL.
- `docker/build.sh test` → passes (group F at its reduced config if needed).
- `qscat.evolution.make_sparse_cn_stepper` promoted, mypy-clean; the numeric outputs,
  convergence table, the σ(E)-vs-σ_TI overlay figure, and the snapshots recorded in
  `docs/physics/n2-2d-td-cross-section.md`; `CLAUDE.md` updated.
