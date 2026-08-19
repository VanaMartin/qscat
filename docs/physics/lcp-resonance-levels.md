# BO/LCP resonance levels — quasi-bound vibrational states of the anion

**Location:** `qscat.core.lcp.ResonanceLevels`/`lcp_resonance_levels`/`resonance_levels`
(`libs/qscat/qscat/core/lcp.py`), `qscat.ecs.match_angle_stable`
(`libs/qscat/qscat/ecs/pole.py`), `IonicResonanceModel.max_nuclear_ecs_angle_deg`
(`libs/qscat/qscat/model/ionic.py`), `apps/qscat-run` (`resonance_levels` observable
kind + `artifacts.resonance_levels` flag, `runner.ResonanceLevelsRun`,
`artifacts._write_resonance_levels`), `libs/qscat/tests/test_lcp_resonance_levels.py`,
`apps/qscat-run/tests/test_runner_resonance_levels.py`,
`apps/qscat-run/tests/test_artifacts_resonance_levels.py`,
`apps/qscat-run/examples/f2-resonance-levels.yaml`.
**Origin:** step 2 of the Born-Oppenheimer/LCP approximation to the 2-D model's
resonance energies. Step 1 (the fixed-`R` electronic pole, giving the complex curve
`V_d(R) - i*Gamma(R)/2`) already existed as `local_complex_potential` and is unchanged
by this work; this note covers step 2, the nuclear eigenvalue problem in that curve.
**Units:** atomic units throughout (energy in Hartree, length in bohr).

## 1. What this computes

The thesis writes the LCP curve and the resulting nuclear Hamiltonian as

```text
V_res(R) = E_res(R) - (i/2) Gamma(R)          (Vana & Houfek 2017, PRA 95, Eq. 41)
H_LCP    = -(1/2 mu) d^2/dR^2 + E_res(R) - (i/2) Gamma(R)      (ibid., Sec. IV)
```

qscat's `Vd` **is** the thesis's `E_res` (both names appear in docstrings so the code
reads against the thesis). `lcp_resonance_levels` builds exactly this operator on a
nuclear FEM-DVR-ECS grid,

```text
H_N = T(mu) + diag(W),        W(R) = Vd(R) - i*Gamma(R)/2
```

using the same `qscat.dvr.kinetic(grid, mu)` kinetic-energy assembly as every other
nuclear problem in the library, and diagonalizes it (complex-symmetric, not Hermitian).
The eigenvalues are the quasi-bound levels

```text
E_v - i*Gamma_v/2,     Gamma_v = max(0, -2 Im E_v)
```

`resonance_levels(model, nuclear_grid_a, nuclear_grid_b, elec_grid_a, elec_grid_b)`
is the model-facing convenience: it runs `resonance_pole_walk` once to get `Vd`/`Gamma`,
lays that curve onto two nuclear grids, and calls `lcp_resonance_levels`.
`lcp_resonance_levels(nuclear_grid_a, nuclear_grid_b, mu, Vd_a, Vd_b, Gamma, ...)` is
the array-taking numeric core (no `ResonanceModel` dependency — `qscat.core` never
imports `qscat.model`).

In the ECS tail, `local_complex_potential` already sets `Vd = model.v0(z) + s_asym`
(the analytic continuation of the anion curve) and `Gamma = 0`. Consequently levels
above the anion dissociation limit pick up a genuine **nuclear (dissociative)** width
from the ECS rotation, while levels below it carry only the **electronic
autodetachment** width baked into `Gamma(R)`. Both come out of one diagonalization —
the reason for doing this complex rather than perturbatively (see Section 7 for when
the perturbative version fails).

## 2. What eMoScat and the thesis actually did

There is no step-2 analog anywhere in the reference implementation. At
`reference/eMoScat/source/Model2d/TimeDependentModel2d.cpp:58-79` (and identically
`source/module_LCP.cpp:295-308`), the code executes `vRes[i] <- Re(vRes[i])` —
**it discards the imaginary part of the curve** — diagonalizes `T(mu) + Re(V_res)`,
and keeps the 15 lowest-`Re E` states blindly, with no angle-stability selection. Those
real-only states are the thesis's `omega_j`: they are used as vertical markers against
cross-section peaks (thesis Fig. 3.26, for NO) and as a projection basis for the
populations `|<omega_i|Psi_d(t)>|^2` (thesis Fig. 3.27). The thesis reports their
*widths* only qualitatively — "the lower the energy of the state, the smaller its
width, the larger its lifetime and the narrower the corresponding structure" — with
lifetimes **estimated from time-dependent peak-formation times**, never computed
directly.

`ResonanceLevels.golden_rule` reproduces exactly this: it is the `Gamma = 0`
diagonalization (`Re(Vd)` implicitly, since setting `Gamma` to the literal array `0`
rather than taking `Vd.real` preserves the ECS-tail analytic continuation) plus the
first-order perturbative width `<chi_v|Gamma|chi_v>_c`. Where `golden_rule` agrees
with the full complex `energies`, this run reproduces what eMoScat/the thesis actually
computed; the complex `energies`/`widths` columns are the genuine extension — the
widths and lifetimes the thesis could only estimate from propagation.

## 3. Terminology: these are NOT Siegert pseudostates

Hvizdoš et al. (Phys. Rev. A **97**, 022704 (2018), App. A) build a *Siegert
pseudostate* basis for the closely related H₂⁺ problem, and the two constructions must
not be conflated:

| | Siegert pseudostates (Hvizdoš et al., App. A) | This work |
|---|---|---|
| Operator | `H_N = -1/(2M) d^2/dR^2 + V0(R)` — the **ion/neutral** curve, real | `H_N = T(mu) + V_res(R)` — the **anion** curve, complex |
| Boundary | outgoing-wave at a **finite** radius `a`: `(d/dR - i K_j) phi_j\|_a = 0` (their Eq. A4) | ECS rotation of the tail; Dirichlet at both ends |
| Orthogonality | bilinear **plus a surface term**, `int_0^a phi_j phi_j' dR + i phi_j(a) phi_j'(a)/(K_j + K_j') = delta` (their Eq. A5) | plain bilinear c-product over the whole rotated grid |
| Purpose | a complete basis for an MQDT frame transformation | the physical quasi-bound levels themselves |

Docstrings and this note call these **complex-scaled (ECS) resonance eigenstates** or,
following the thesis, **quasi-bound vibrational states `omega_j`** — "Siegert
pseudostate" is reserved for the Hvizdoš et al. construction.

## 4. Normalization

`ResonanceLevels.states` are c-product-normalized: `sum_i c_i^2 = 1` under
`qscat.linalg.c_product` (the bilinear, non-conjugated inner product), computed over
the **whole rotated grid**, never `np.vdot`/`.conj()`/`np.linalg.norm`.

This is the direct consequence of Section 3's comparison. The Siegert pseudostates'
surface term (their Eq. A5) exists *precisely because* that construction truncates the
domain at a finite radius `a` with an outgoing-wave boundary condition there — the
surface term is what restores completeness/orthogonality across that truncation. ECS
does not truncate; it rotates the tail into the complex plane and lets the wavefunction
decay there. So the plain bilinear c-product over the full (real + rotated-tail) grid
is already the correct and complete inner product, with no surface correction needed.

Note this differs deliberately from the existing `resonance_eigenstate`
(the *electronic* resonance eigenstate at fixed `R`), which normalizes over the real
region only — that function reports a real-region-restricted quantity by design; this
one reports the genuine Siegert-type nuclear state, and the two docstrings say so.

## 5. Selection of the physical levels

`Vd(R)` at real `R` does not depend on the nuclear ECS tail angle, so the expensive
electronic pole walk (`resonance_pole_walk`) runs **once**. A second nuclear grid,
sharing every real node and quadrature with the first but rotated at a different tail
angle, costs one extra nuclear diagonalization and buys the two-angle stability test —
the same criterion eMoScat used electronically
(`reference/eMoScat/source/DiscreteStates.cpp:204-309`), generalized from a single pole
(`qscat.ecs.find_resonance_pole`) to many simultaneous states
(`qscat.ecs.match_angle_stable`).

`match_angle_stable(eigs_a, eigs_b, window, rel_tol=1e-4, atol=1e-8)` keeps every
eigenvalue of grid `a` inside `window` whose nearest grid-`b` partner satisfies
`|E_a - E_b| < max(rel_tol*|E_a|, atol)`, returning the midpoint `(E_a+E_b)/2` and the
residual `|E_a - E_b|`. The rotated discretized continuum fails this test (its
eigenvalues genuinely move with the rotation angle) and drops out; the physical,
angle-invariant levels survive. Sharing the real nodes between the two grids is what
makes the comparison meaningful — a mismatch there raises `ValueError` rather than
silently comparing two different discretizations.

Because the nuclear ECS tail's `Vd` continuation is angle-independent at real `R`
(`local_complex_potential` assembles it once from the pole walk's output), the second
grid needs no second electronic solve — only a second nuclear diagonalization of an
already-cheap 1-D problem.

**One walk per run, and the curve you plot is the curve the levels came from.** The
electronic pole walk dominates the cost (~25 s of a ~27 s F2 run), so a caller that
wants both the levels *and* the LCP curve — to solve `lcp_da_cross_section` in it, or to
draw it under the level bars — must take `resonance_levels(..., return_curve=True)`,
which returns `(levels, Vd_a, Gamma)` from the single walk it already ran. Calling
`local_complex_potential` separately would double the run time *and* — if any setting
differed, as it did before this was fixed — report a `V_d`/`Gamma` in the artifacts that
is **not** the curve the levels were diagonalized in. `apps/qscat-run`'s LCP path takes
the `return_curve` route and uses one half-width setting (`runner._LCP_WALK_HALF_WIDTH`)
for every LCP observable it produces.

**The default window's `Im` band is sized for autodetachment only.** With no explicit
`window`, `lcp_resonance_levels` spans `Re` over the anion curve at the real nodes and
`Im` down to `-max_R Gamma(R)`. That floor is correct for a level *below* the anion
dissociation limit (`Gamma_v = <chi_v|Gamma|chi_v> <= max Gamma`), but it is **not** a
bound on the nuclear (dissociative) width of a level *above* it: that width is generated
by the ECS rotation of the tail and bears no relation to `Gamma(R)` — for a barrierless
curve it is `~1e-3` Ha, orders of magnitude below the default floor. Such levels fall
outside the default window and are simply absent from the result. Pass an explicit
`window` with a low enough `im_lo` to look for them. The degenerate case — `Gamma ~ 0`
over the whole grid, so the band collapses to `+-atol` and *no* dissociative level
whatsoever can be represented — now emits a `UserWarning` rather than silently returning
a bound-states-only spectrum.

**The real/ECS-tail junction in `Gamma`.** `_assemble_lcp` force-zeroes `Gamma` in the
rotated tail. If the walk's `Gamma` is not already ~0 at the outermost *real* node, the
local complex potential `W = V_d - i Gamma/2` **steps** at the junction and the tail
reflects the outgoing dissociative wave instead of absorbing it. That condition now
warns (`_JUNCTION_GAMMA_TOL = 1e-8`); the fix is to extend the real region outward until
the autodetachment width has died off. It does not fire on the F2 preset deck.

## 6. The ECS angle bound

Hvizdoš et al. (§II) show that for the H₂⁺ model potential the ECS rotation angles
must satisfy `theta_nuclear < pi/8` (22.5 deg) and `theta_electronic < pi/4` (45 deg);
otherwise `V(R,r)` **diverges** at large `R`/`r` under the rotation. The `a3*R^4` term
in their potential (their Eq. 8) needs `4*theta < pi/2` to stay bounded; the electronic
`exp(-r^2/3)` factor needs `2*theta < pi/2`. `IonicResonanceModel.max_nuclear_ecs_angle_deg
= 22.5` records the nuclear bound (eMoScat's own H₂⁺ nuclear deck uses 22.0 degrees —
just inside it); `_check_angle_bound` rejects any nuclear grid whose worst element
angle **reaches or exceeds** that bound (strict rejection at the boundary itself, since
`4*theta < pi/2` is a strict inequality and `theta == bound` is already the marginal,
non-decaying case). Neutral diatomics (`DiatomicResonanceModel`) carry no such bound —
`getattr(model, "max_nuclear_ecs_angle_deg", None)` is `None` and the check is a no-op.

## 7. Validation

Two exact analytic oracles, gated in `libs/qscat/tests/test_lcp_resonance_levels.py`:

- **Zero-width limit** (`test_zero_width_reproduces_the_analytic_morse_spectrum`):
  `Gamma = 0`, `Vd` a bare Morse curve. The bound-state eigenvalues reproduce the
  closed-form Morse spectrum to `rtol=1e-5` (measured 1e-12 to 1e-7, growing mildly
  with `n` — textbook FEM-DVR spectral truncation error, nothing anomalous), `Im(E)`
  and `widths` at round-off, `real_weight == 1.0`.
- **Constant-width rigid shift** (`test_constant_width_shifts_the_spectrum_rigidly`):
  `H(Gamma_0) = H(0) - i*(Gamma_0/2)*I` is an algebraic identity; the complex machinery
  reproduces it to round-off (`shifted.energies == zero.energies - 0.5j*g0`,
  `shifted.widths == g0`).

Plus a convergence study, a golden-rule consistency test, and a genuine two-angle
selection test:

- **h-refinement** (`test_levels_converge_under_h_refinement`): `energies` at a
  coarser and a finer real-region discretization agree to `rtol=1e-6` (fixture chosen
  as `n_real=30` vs `n_real=60`, both cross-checked against the analytic Morse
  spectrum: errors 7.03e-8 and 8.29e-12 respectively — the discretization is genuinely
  converging, not coincidentally agreeing).
- **Angle independence** (`test_levels_are_independent_of_the_ecs_angle_pair`): the
  selected levels don't depend on which pair of (safe) angles is used — gated at
  `rtol=1e-6, atol=1e-9` on the complex `energies`. Note the `atol`: this test does
  **not** constrain the imaginary part below `1e-9`, which is exactly the regime
  Section 7.2 shows to be noise.
- **Golden-rule agreement in the perturbative regime**
  (`test_golden_rule_matches_the_complex_result_for_a_weak_constant_width`) and
  **golden-rule divergence in the non-perturbative regime**
  (`test_energies_diverge_from_golden_rule_for_a_broad_r_dependent_level`, an `R`-
  localized Gaussian `Gamma(R)` bump straddling the outer well): the ground level's
  `|energies - golden_rule|` stays near-perturbative (`< 5e-5`) while the most extended
  level's divergence exceeds its own comparator-predicted width — a qualitative, not
  just quantitative, breakdown, confirming `golden_rule`'s divergence from `energies`
  is real physics (the wavefunction rearranging away from a strong, localized loss
  region — a second-order-and-higher effect the linear comparator cannot see), not a
  numerical artifact (`residuals < 1e-10` throughout, still genuinely angle-stable).
- **The two-angle selection doing real work, not passing vacuously**
  (`test_two_angle_selection_isolates_the_resonance_from_the_rotated_continuum`): at a
  near-threshold window with real tail amplitude, the *raw* windowed spectrum on grid
  `a` alone changes count with the rotation angle (6 vs 7, across three angle pairs —
  the discretized continuum literally rotates with `theta`), while the *selected*
  count is angle-invariant (always 1) and the surviving level's `real_weight < 0.96`
  (genuine, non-trivial tail amplitude, not a compact bound state sitting entirely in
  the real region).

**Honest limitation.** No prior computed width `Gamma_v` exists anywhere — in eMoScat,
the thesis, or Hvizdoš et al. — to check the imaginary parts against; the thesis's own
widths are qualitative estimates from propagation timing, not numbers. Until the
milestone-3 comparisons (thesis `omega_j` reproduction, elastic-peak correspondence,
NO lifetime bounds — see Section 9) land, the imaginary parts are gated only by
**internal consistency**: the two analytic oracles above (both of which are `Gamma=0`
or rigid-shift limits, so neither exercises a genuinely R-dependent resonance width),
the h-refinement/angle-independence checks, and the golden-rule agreement/divergence
pair. This is real, useful evidence that the machinery is self-consistent — it is not
yet an external check on whether `Gamma_v` is *correct* for a real molecule.

Two consequences of that gap are load-bearing enough to have their own sections below,
and **both must be read before quoting any number out of this capability**: what
angle-stability does and does not certify (7.1), and the width scale below which the
reported `Gamma_v` is noise (7.2). The single check that would close most of the gap —
`lcp_da_cross_section` in the *same* curve, whose `sigma(E)` peaks must sit at `Re E_v`
with FWHM `Gamma_v` — is not built yet and is the top item in Section 9.

### 7.1 Angle-stability is NOT sufficient to trust a level

**`residuals` measures ECS-tail angle stability only, not real-region convergence.**
`nuclear_grid_a` and `nuclear_grid_b` are required to share every real node and
quadrature by construction (`lcp_resonance_levels` raises `ValueError` otherwise), so
real-region discretization error is common to both spectra and cancels out of their
difference — `residuals` stays pinned near machine precision (`~1e-14` to `~1e-15`)
regardless of how coarse or fine the shared real grid is. It is a genuine, useful
diagnostic for whether a level is contaminated by the rotated continuum, but it says
nothing about whether the shared real-region discretization itself is converged.

**This is not a theoretical caveat — a pure grid artifact passes the test.** A state
localized in the real region is *bit-identical* on both grids (they share every real
node), so its two eigenvalues agree to round-off and it is accepted at residual `~1e-16`
whether or not it corresponds to anything physical. Demonstrated directly on the shipped
Morse fixture: an accepted "level" at

```text
Re E = +2.615e-2 Ha,  Gamma = 0,  residual = 3.8e-16,  real_weight = 1.000
```

sits **above the dissociation limit of a barrierless Morse curve**, where no bound or
quasi-bound state can exist. It is a discretization artifact, and it moves with the
discretization exactly as an artifact does and a physical level does not:

| discretization | accepted level above the limit |
|---|---|
| `n_real = 30`, `r_max = 6` (the fixture) | `2.615e-2` Ha |
| `r_max = 8` | `3.457e-2` Ha |
| `n_real = 45` | `1.134e-2` Ha |
| `n_real = 60` | **none — it is gone** |

Note in particular that a `real_weight` floor would **not** have rejected it —
`real_weight` is `1.000`, the most "trustworthy-looking" value the diagnostic can take.

**Trusting a reported level therefore requires three checks, and this API only performs
the first:**

1. **Angle-pair invariance** — provided (`residuals`, `match_angle_stable`). Rejects
   contamination by the rotated continuum. Nothing else.
2. **Invariance of `energies` under real-node h-refinement** — *the caller must do
   this.* Re-run with the shared real region refined (more points per element, or more
   elements) and confirm the level does not move and does not disappear. This is what
   catches the artifact above.
3. **Invariance under `r_max` / tail extent** — *the caller must do this.* Re-run with
   the real region extended outward and the ECS tail moved with it. A physical level is
   insensitive to where the box ends; an artifact is not.

Checks 2 and 3 are deliberately **not** automated here (see Section 9); they are
required user practice, and the h-refinement test in Section 7 is the worked example of
what check 2 looks like.

### 7.2 Widths below ~1e-6 Ha are inside the method's own noise floor

The acceptance tolerance is `atol = 1e-8` — *the same scale as the widths being
reported* when those widths are small. Below roughly `1e-6` Ha a reported `Gamma_v` is
not a measurement, it is round-off that survived selection.

Measured on the `test_lcp_resonance_levels.py` fixture (`grid_pair()`,
`window = (-0.05, -1e-4, ±1e-6)`), at the near-threshold level
`Re E = -2.22597e-2` Ha, whose true width is exactly **zero**:

| angle pair (deg) | `Im E_v` (Ha) | reported `Gamma_v` (Ha) | residual |
|---|---|---|---|
| 35 / 25 | −6.69e-10 | 1.338e-9 | 1.31e-9 |
| 45 / 20 | **+**6.34e-10 | 0 (clamped) | 4.98e-9 |
| 40 / 30 | −4.50e-10 | 9.00e-10 | 2.21e-9 |

Every row passes the acceptance test. Yet the same level is reported with three
different widths spanning an order of magnitude — and on one angle pair `Im E` comes out
with the **wrong sign** entirely, so `Gamma_v = max(0, -2 Im E)` clamps it to zero. At
this scale even the sign of the imaginary part is noise.

**Do not quote widths below ~1e-6 Ha**, in this note or anywhere downstream — including
the `~1e-15` values in the F2 table below, which are shown only to demonstrate they are
at the floor, not as measurements.

**`golden_rule` can legitimately be `nan`.** All-`nan` when the comparator is
unavailable (the `Gamma=0` problem has no angle-stable state in the tight
`[-atol, atol]` diagnostic window — this can happen even though the *primary* solve
correctly keeps the corresponding complex level as physical) or disabled
(`golden_rule=False`). Per-level `nan` when a distance guard rejects the nearest
candidate pairing: either the Gamma-induced real-part shift exceeds half the local
level spacing in the primary spectrum (a strongly non-perturbative shift — `nan` is
the honest answer, not a wrong number), or output levels are near-degenerate. Both are
tested directly (`test_golden_rule_returns_nan_when_comparator_window_is_empty`,
`test_golden_rule_nans_a_level_with_no_plausible_comparator`).

## 8. First results: F2

Config: `apps/qscat-run/examples/f2-resonance-levels.yaml` — `molecule: F2`,
`methods: [lcp]`, `observables: [{kind: resonance_levels, channels: 6}]`,
`grid: {preset: emoscat}` (the production F2 deck; nuclear angle pair 35/25 degrees,
via `presets.resolve_lcp_grids`).

```bash
uv run qscat-run run apps/qscat-run/examples/f2-resonance-levels.yaml --output runs/f2-levels
```

Real run output (`resonance_levels_lcp_resonance_levels.csv`), the first numbers this
capability has produced:

| v | Re E (Ha) | Gamma_v (Ha) | residual | real_weight | Re E₀ (golden rule) | Gamma_v⁽¹⁾ |
|---|---|---|---|---|---|---|
| 0 | −0.148167188 | 2.94e-15 | 4.60e-16 | 1.0000000000 | −0.148167188 | 3.70e-15 |
| 1 | −0.145851066 | 4.37e-15 | 2.68e-15 | 1.0000000000 | −0.145851066 | 4.00e-15 |
| 2 | −0.143668367 | 5.68e-15 | 8.96e-16 | 1.0000000000 | −0.143668367 | 5.26e-15 |
| 3 | −0.141619191 | 5.31e-15 | 1.26e-15 | 0.9999999999 | −0.141619191 | 5.44e-15 |
| 4 | −0.139703632 | 5.67e-15 | 2.83e-15 | 0.9999999998 | −0.139703632 | 6.14e-15 |
| 5 | −0.137921779 | 4.16e-15 | 1.11e-15 | 1.0000000000 | −0.137921779 | 4.51e-15 |

All six levels are effectively real-valued (`Gamma_v ~ 1e-15`, the solver's own
round-off floor — well inside the `~1e-6` Ha noise floor of Section 7.2, so those
figures say "zero", not "1e-15") — deeply bound anion vibrational states, `energies` and
`golden_rule` agreeing to `~1e-9` (fully perturbative), `real_weight ~ 1.0`
(negligible ECS-tail amplitude), `residuals` at the machine-precision floor
(angle-stable). Density-weighted mean position `<R>` ranges 3.38 to 3.89 bohr across
v=0..5 (F2's neutral Morse minimum is at `R0 = 2.6906` bohr).

**Why these six are all zero-width, and what that does and doesn't demonstrate.**
This is not a null result about the imaginary machinery — it is a statement about
which part of the F2 anion spectrum the shipped example shows. F2's anion is bound
at large `R`, so its lowest vibrational levels sit entirely in the region where
`Gamma(R) ~ 0` (autodetachment closed, `Vd < v0`) and are genuinely bound against
autodetachment, not merely numerically narrow. The complex machinery itself is
validated by the two analytic oracles in Section 7 (in particular the constant-`Gamma`
rigid shift, which exercises a *nonzero*, exactly-known width), not by this example.

A scratch run (`channels=40` on the same config, not shipped) shows where non-zero
widths appear: `v=0..17` stay bound (`Gamma_v ~ 1e-15`, `real_weight > 0.999`),
then at `v=18` the width jumps abruptly to `~5e-8` and climbs smoothly through
`~2.8e-5` by `v=39`, while `real_weight` simultaneously collapses from `0.9997` to
`0.001-0.05`.

**These `v>=18` states are numerical tail states and must not be quoted as physics.**
The `real_weight` collapse is not a side effect of the onset — it is what rules the
onset out, and the argument does not depend on the R≈2.597 pole-walk artifact at all:

- **As autodetaching resonances they are impossible.** A Siegert state's outgoing tail
  goes as `exp(i k R)` with `k = k_r + i k_i` and `k_i/k_r ~ Gamma_v/4E`. At
  `Gamma_v ~ 5e-8` and `E ~ 0.13 Ha` that ratio is `~1e-7`, so under the tail's 25°
  rotation the state's amplitude is killed essentially immediately past `R0`: a genuine
  resonance this narrow **must** have `real_weight ~ 1`. A measured `real_weight` of
  `1e-3` — three orders of magnitude of the norm living in the rotated tail — is
  flatly incompatible with a width of `5e-8`.
- **As dissociative (nuclear-width) levels they are equally impossible.** The benign
  alternative would be a level above the anion dissociation limit picking up a nuclear
  width from the ECS rotation. But F2's anion curve is barrierless there, and a
  barrierless dissociative width is `~1e-3` Ha — **five orders of magnitude larger**
  than the `5e-8` reported. A level cannot be delocalized into the rotated tail and
  simultaneously that narrow.

What is left is the discretized rotated continuum: states that live mostly in the tail,
whose apparent "width" is a discretization scale rather than a decay rate. They pass the
two-angle test only in the sense Section 7.1 describes — nothing has been checked against
h-refinement or `r_max`. Independently of this argument, their width scale (`1e-8` to
`~3e-5`) also overlaps the known R≈2.597 pole-walk artifact discussed below and their
`Re E` (`-0.1269...`) sits at the crossing region, so even taken at face value they would
be inseparable from that defect. For both reasons the shipped example stays at
`channels: 6` — an unambiguous slice of the spectrum — and no `v>=18` number from this
scratch run should be reported, plotted, or cited. Whether genuinely autodetaching F2
anion levels exist above the crossing is open, and the way to settle it is the
`sigma(E)` cross-check in Section 9, not a wider `channels`.

**A known rough edge, stated plainly rather than smoothed over.** On this runner's
real F2 preset grid, the electronic pole walk leaves a spurious `Gamma` ~2.3e-5 at
R≈2.597, right where the anion curve crosses the neutral (`Vd = v0`) and `Gamma` must
be exactly zero — the imaginary part is nonzero only where `v0(R) < E_res(R)`
(Vana & Houfek 2017, PRA 95, Sec. IV). This persists at **both** the default 0.05 and
the tightened 0.01 electronic-walk half-widths — it is grid-dependent, traced to the
outer real-segment count, not simply a half-width setting (the library's own toy test
grid *does* clear it at 0.01: 2.17e-5 -> 3.0e-14 there). Consequently **every F2
`resonance_levels` run through this runner emits a `UserWarning`**, including the
config above. The reported levels are unaffected: for the eight lowest F2 anion
levels, the fraction of `|c|^2` falling inside the crossing's own R∈[2.5976, 2.608]
window ranges from 5.4e-18 (v=0) to 4.8e-13 (v=7) of the total density — these
low-lying states are radially localized 0.7-1.2 bohr outward of the crossing, in the
steeply rising part of `Vd(R)`, so their overlap with the artifact-affected region is
many orders of magnitude below anything that could matter. Both halves of this are
real: the warning fires and is honest signal about a genuine curve defect at R≈2.6;
the six/eight reported levels are demonstrably fine. This is recorded as a known rough
edge for whoever next tunes F2's discretisation or the pole walk's step size on this
crossing (the same R≈2.6 feature the discretisation tuner's own findings already flag
as under-resolved by the a-priori mesh) — not fixed here, since it is
`resonance_pole_walk`/`local_complex_potential` library work, out of this plan's
scope.

## 9. Limits and what is next

- **The highest-value gap: no oracle with a NONZERO width, and the repo already owns
  one.** Every `Gamma` assertion in Section 7 is exactly zero (`test_zero_width_...`),
  exactly constant (`test_constant_width_...`, a rigid algebraic shift), or
  self-consistency (angle/h-refinement/golden-rule). **Nothing checks a genuinely
  `R`-dependent width against a known answer.** The right differential oracle is already
  in this library and is not being used: `lcp_da_cross_section` solves the resolvent in
  the **same** `(V_d, Gamma)` curve these levels are eigenvalues of, so
  `sigma_DA(E)` **must** show peaks at `Re E_v` with FWHM `Gamma_v` — an independent,
  fully internal check of both parts of every complex eigenvalue, on a real molecule,
  with no external data needed. It would also have settled Section 8's `v>=18` question
  immediately: those states, if real, would appear as structure in `sigma(E)`, and they
  do not. Building this is milestone 3 and is **the top follow-up**, ahead of everything
  else in this list.
- **Refinement-invariance (Section 7.1, checks 2 and 3) is user practice, not API.**
  Nothing in `ResonanceLevels` surfaces whether a level survives h-refinement or a
  larger `r_max`; a caller who skips those checks can and will be handed an artifact at
  `residual ~ 1e-16` and `real_weight = 1.000`. Automating them (re-solving on a refined
  grid and reporting per-level movement) is a deliberate non-goal here, not an oversight
  — but it is the obvious next API-level improvement.
- **No comparison to a real target yet.** These six F2 levels are internally
  consistent (Section 7) but have not been checked against anything external. That is
  milestone 3, deferred because it needs data the thesis author holds (the actual
  `omega_j` values, Figs. 3.26/3.27) and is not a prerequisite for anything in this
  plan: reproducing the thesis's `omega_j` (via `golden_rule`, which is built to match
  exactly what eMoScat computed), the elastic-peak correspondence (which cross-section
  feature belongs to which level), NO lifetime bounds from the reported qualitative
  widths, and running both the eMoScat-deck and thesis-grid preset variants as a
  convergence cross-check.
- **`Psi_d` projection / level populations are deferred.** The thesis's second use of
  `omega_j` — projecting a propagated discrete-state amplitude `Psi_d(t)` onto these
  levels to get populations `|<omega_i|Psi_d(t)>|^2` (its Fig. 3.27) — needs the
  discrete-state formalism this plan does not build; `ResonanceLevels.states` already
  carries what a future projection would need (c-product-normalized DVR coefficients).
- **A full 2-D complex-scaled diagonalization** (the direct Siegert states of the 2-D
  model, without the BO reduction) is the natural next comparison and is not attempted
  here.
- **`H2P.mu` already corrected.** The design spec for this plan flagged a discrepancy —
  eMoScat's JSON deck (and `qscat.model.library`, at the time) carried `H2P.mu =
  918.25`, while thesis Table 1.2 and Hvizdoš et al. §II A both give `M = 918.076`
  (`m_p/2` for the modern proton mass) — and deferred fixing it. That fix landed ahead
  of this plan's own commits (`fix(model): H2+ reduced mass 918.25 -> 918.076`); see
  the Changelog. It is noted here only because no H₂⁺ `resonance_levels` run has yet
  been validated against a published number to confirm the corrected value resolves
  the ~1e-4 relative spacing error the spec predicted.
