# BO/LCP resonance levels — quasi-bound states of the 2-D model

Date: 2026-08-15

## Purpose

Compute the **Born-Oppenheimer approximation to the resonance energies** of the
2-D electron-diatomic model: solve the electronic problem at fixed `R` to get the
complex potential curve `V_d(R) - i Gamma(R)/2` (the LCP curve), then diagonalize
the *nuclear* problem in that curve to get complex vibrational eigenvalues
`E_v - i Gamma_v/2` — the quasi-bound (Siegert) states of the anion.

Three uses, staged:

1. an accuracy statement — BO/LCP levels vs the structure of the exact 2-D
   solver's `sigma(E)`;
2. an interpretive tool — which feature of `sigma_DA`/`sigma_VE` belongs to which
   quasi-bound level;
3. a cheap resonance-energy predictor — a per-molecule table computable without
   any energy sweep.

This is the research program's standard shape: the exact 2-D solver
(`qscat.core.driven`/`dissociation`) is the oracle, and the BO/LCP reduction is
the approximation under test.

## Reference literature and notation

The primary reference is the thesis these models were studied in:

- M. Váňa, *A model of resonant collisions of electrons with molecules and
  molecular ions*, doctoral thesis, Charles University, Prague, 2017.
  <https://dspace.cuni.cz/handle/20.500.11956/92902>

and the paper that treats the closely-related H2+ problem with an explicitly
Siegert basis:

- D. Hvizdoš, M. Váňa, K. Houfek, C. H. Greene, T. N. Rescigno, C. W. McCurdy,
  R. Čurík, *Dissociative recombination by frame transformation to Siegert
  pseudostates: a comparison with a numerically solvable model*, Phys. Rev. A
  **97**, 022704 (2018); arXiv:1710.10333.

Supporting citations, all as given in that bibliography:

- K. Houfek, T. N. Rescigno, C. W. McCurdy, *Numerically solvable model for
  resonant collisions of electrons with diatomic molecules*, Phys. Rev. A **73**,
  032721 (2006) — the 2-D model itself.
- K. Houfek, T. N. Rescigno, C. W. McCurdy, *Probing the nonlocal approximation
  to resonant collisions of electrons with diatomic molecules*, Phys. Rev. A
  **77**, 012710 (2008) — the discrete-state choices.
- W. Domcke, *Theory of resonance and threshold effects in electron-molecule
  collisions: the projection-operator approach*, Phys. Rep. **208**, 97 (1991) —
  the LCP approximation.
- K. Houfek, M. Čížek, J. Horáček, *On irregular oscillatory structures in
  resonant vibrational excitation cross-sections in diatomic molecules*, Chem.
  Phys. **347**, 250 (2008).
- D. Birtwistle, A. Herzenberg, J. Phys. B **4**, 53 (1971); L. Dubé,
  A. Herzenberg, Phys. Rev. A **20**, 194 (1979) — the boomerang model.

**Notation.** The thesis writes the LCP curve as

```
V_res(R) = E_res(R) - (i/2) Gamma(R)                         (thesis Eq. 1.63)
H_LCP    = -(1/2 mu) d^2/dR^2 + E_res(R) - (i/2) Gamma(R)    (thesis Eq. 1.65)
```

and calls the quasi-bound vibrational levels of the molecular anion `omega_j`.
qscat's `Vd` **is** the thesis's `E_res`; `H_LCP` is exactly the `H_N` below.
Docstrings and the physics note should give both names so the code is readable
against the thesis.

### Terminology: these are NOT Siegert pseudostates

Hvizdoš et al. use a *Siegert pseudostate* basis, and the two constructions must
not be conflated:

| | Siegert pseudostates (Hvizdoš et al., App. A) | This work |
|---|---|---|
| Operator | `H_N = -1/(2M) d^2/dR^2 + V0(R)` — the **ion/neutral** curve, real | `H_N = T(mu) + V_res(R)` — the **anion** curve, complex |
| Boundary | outgoing-wave at a **finite** `a`: `(d/dR - i K_j) phi_j|_a = 0` (Eq. A4) | ECS rotation of the tail; Dirichlet at both ends |
| Orthogonality | bilinear **plus a surface term**, `int_0^a phi_j phi_j' dR + i phi_j(a) phi_j'(a)/(K_j + K_j') = delta` (Eq. A5) | plain bilinear c-product over the rotated grid |
| Purpose | a complete basis for an MQDT frame transformation | the physical quasi-bound levels themselves |

**This settles the normalization question.** The surface term in Eq. (A5) exists
precisely *because* the pseudostate construction truncates at finite `a` with an
outgoing boundary condition. ECS rotates the tail instead of truncating it, so the
plain bilinear c-product over the whole grid is the correct and complete inner
product — no surface correction. The spec's choice stands, and now has a citation
for why it differs from the pseudostate literature.

The docs should therefore call these **complex-scaled (ECS) resonance
eigenstates** or, following the thesis, **quasi-bound vibrational states
`omega_j`** — reserving "Siegert pseudostate" for the Hvizdoš et al. construction.

### ECS angle constraints (hard, model-dependent)

Hvizdoš et al. (§II) state that for the H2+ model the bending angle must satisfy
`theta_nuclear < pi/8` (22.5 degrees) and `theta_electronic < pi/4` (45 degrees),
otherwise `V(R,r)` **diverges** at large `R`/`r` — the quartic `alpha_3 R^4` term
in Eq. (8) needs `4*theta < pi/2`, and the `exp(-r^2/3)` factor needs
`2*theta < pi/2`.

This binds the second-angle grid. The default `angle_b = angle_a - 10` degrees
moves *downward* and is therefore always safe, but the constraint must be asserted
rather than assumed: `resonance_levels` validates both nuclear grids' tail angles
against a model-supplied bound, and the ionic model reports 22.5 degrees.
(eMoScat's H2+ nuclear deck uses 22.0 degrees — just inside it.)

### Consistency checks already performed

- Every N2/NO/F2 constant in `qscat.model.library` (`mu`, `l`, `D0`, `alpha0`,
  `R0`, `lambda_inf`, `lambda_1`, `R_lambda`, `lambda_c`, `R_c`, `alpha_c`)
  matches thesis Table 1.1 exactly.
- H2+: `V0 = 0.1027`, `R0 = 2.0`, `alpha = 0.69`, `a1 = 1.6435`, `a2 = 6.2`,
  `a3 = 0.0125`, `a4 = 1.15`, `l = 1` all match thesis Table 1.2 and Hvizdoš
  et al. Eqs. (7)-(8).
- **Discrepancy found: `H2P.mu`.** `qscat.model.library` has `mu = 918.25`, but
  thesis Table 1.2 and Hvizdoš et al. §II A both give `M = 918.076`, which is
  `m_p/2` for the modern proton mass (1836.15267/2). The 918.25 value came in
  from eMoScat's JSON deck. The error is 0.019%, so it shifts H2+ vibrational
  spacings by ~1e-4 relative — harmless qualitatively, but wrong for reproducing
  published numbers. **Fix separately, before any H2+ level comparison is
  quoted**; it is not part of this spec's scope.

## What already exists, and what eMoScat actually did

**Step 1 is already ported and is not changed by this work.**
`qscat.core.lcp.local_complex_potential` reproduces eMoScat's
`ModelLCP::make_vres` (`reference/eMoScat/source/ModelLCP/ModelLCP.cpp:152-286`):
the two-electronic-ECS-angle pole of `H_el(R) = -1/2 d^2/dr^2 + model.surface(r,R)`,
seeded from the asymptotic anion bound state, walked inward over real `R`,
electronic shift frozen on breakdown, `Gamma = max(0, -2 Im E)`, `Gamma = 0` in the
nuclear ECS tail. The qModeling version re-centres the search window on every
accepted pole, which the reference did not (`ModelLCP.cpp` updates its bound only
on failure) — keep that improvement.

**Step 2 does not exist in eMoScat.** At
`reference/eMoScat/source/Model2d/TimeDependentModel2d.cpp:58-79` (and identically
`source/module_LCP.cpp:295-308`) the code executes `vRes[i] <- Re(vRes[i])` — it
**discards the imaginary part** — diagonalizes `T(mu) + Re V_res`, and keeps the 15
lowest-`Re E` states blindly with no stability selection. Those states are used
only as a projection basis for `|<v_k|psi(t)>|^2` populations. There are no complex
anion vibrational eigenvalues anywhere in the reference, no 2-D complex-scaled
diagonalization, and no committed numeric outputs to gate against. The nearest
genuine BO-nuclear-eigenvalue code is `source/time_independent_model.cpp:344-415`
(H2+ `compute_vibration_thresholds`), which diagonalizes `T(mu) + U_k(R)` on
**real** adiabatic curves.

**These real-part-only levels are exactly the thesis's `omega_j`.** The thesis
describes them as "eigenstates in the potential energy `V_res(R)`" and uses them
in two places: as vertical markers against the cross-section peaks (Fig. 3.26 for
NO), and as a projection basis for the populations `|<omega_i|Psi_d(t)>|^2`
(Fig. 3.27). It reports their *widths* only qualitatively — "the lower the energy
of the state, the smaller its width, the larger its lifetime and the narrower the
corresponding structure" — with lifetimes **estimated from time-dependent peak
formation times**, never computed.

That reconciliation sets the terms of this work precisely:

- the **golden-rule comparator** column `E_v^(0)` reproduces the thesis's
  `omega_j` and is therefore directly checkable against data already in hand;
- the **full complex** `E_v - i Gamma_v/2` is the genuine extension — it computes
  the widths and lifetimes the thesis could only estimate.

So: step 1 is a completed port, step 2's real part is a reproduction target, and
step 2's imaginary part is new physics output.

## Design

### The nuclear Siegert problem

```
H_N = T(mu) + diag(W),      W(R) = V_d(R) - i Gamma(R)/2
```

complex-symmetric FEM-DVR-ECS on the nuclear grid, Dirichlet at both ends —
`qscat.dvr.kinetic(grid, mu)` builds exactly this operator. Eigenvalues
`E_v - i Gamma_v/2` are the BO/LCP approximation to the resonance energies;
`Gamma_v = max(0, -2 Im E_v)`.

In the ECS tail `local_complex_potential` already sets
`V_d = model.v0(z) + s_asym` (the analytic continuation of the anion curve, with
`s_asym` the shift at the outermost real `R`) and `Gamma = 0`. Consequently levels
above the anion dissociation limit acquire a genuine **nuclear (dissociative)**
width, while levels below it are nuclear-bound and carry only the **electronic
autodetachment** width. Both emerge from one diagonalization — the reason for doing
it complex rather than perturbatively.

Two corrections relative to the reference, both required for ECS correctness:

- eigenvectors are normalized with the **bilinear c-product** (`c . c = 1`,
  `qscat.linalg.c_product`), not LAPACK's `||c||_2` (`EigenSystem.hpp:68-75`);
- no conjugated inner product (`zdotc`, `Vector.hpp:299-306`) anywhere in the
  rotated region.

Note the existing `resonance_eigenstate` normalizes its *electronic* eigenvector
over the real region only. For a Siegert state the full-grid bilinear norm is the
correct choice and is what this function uses; the real-region weight is reported
separately as a diagnostic rather than used as a normalization. The two functions
therefore differ deliberately, and the docstrings say so.

### Selection of the physical levels

`V_d(R)` at real `R` is angle-independent, so the expensive electronic pole walk
runs **once**. A second nuclear grid with identical `real_segments` and
`quadrature` but a different tail `angle_deg` shares every real node and differs
only in the tail; one extra nuclear diagonalization then buys the two-angle
stability test — the same criterion eMoScat used electronically
(`DiscreteStates.cpp:204-309`) and the same one `qscat.ecs.find_resonance_pole`
already implements for a single pole.

New vectorized primitive in `qscat/ecs/pole.py`:

```python
def match_angle_stable(
    eigs_a, eigs_b, window, *, rel_tol=1e-4, atol=1e-8
) -> tuple[NDArray[complex128], NDArray[float64], NDArray[intp]]:
    """Every eigenvalue of `eigs_a` inside `window` whose nearest `eigs_b`
    partner satisfies |E_a - E_b| < max(rel_tol*|E_a|, atol).
    Returns (energies, residuals, indices_into_a)."""
```

`energies` are the midpoints `(E_a + E_b)/2`, matching `find_resonance_pole`'s
convention; `indices_into_a` lets the caller recover the eigenvector on grid *a*.

`find_resonance_pole` and `match_angle_stable` share one private pairing helper
(window filter plus nearest-partner distances), but keep their distinct accept
rules: `find_resonance_pole` returns the globally smallest-residual pair
**whatever its residual** and lets the caller judge, whereas `match_angle_stable`
returns every pair under the tolerance. `find_resonance_pole`'s observable
behaviour is unchanged — its existing tests must pass untouched.

Rejected alternative: a single grid with a hand-set `|Im E|` threshold. Cheaper,
but the accept/reject boundary would be a tuned constant rather than a physical
invariance — unusable exactly where the interesting broad levels live.

### The golden-rule comparator

The same routine also diagonalizes `T(mu) + diag(Re V_d)` (with the same two-angle
selection) and forms `Gamma_v^(1) = <chi_v|Gamma(R)|chi_v>_c` using the c-product.
This is eMoScat's actual step 2, plus the first-order width it never computed. It
costs one extra pair of diagonalizations on an already-cheap 1-D problem.

Where `E_v^(0) - i Gamma_v^(1)/2` agrees with the full complex `E_v`, the level is
perturbative; where they diverge, the level is genuinely broad and the
non-perturbative treatment is load-bearing. That divergence is a reported result,
not an error.

### Library API — `qscat.core.lcp`

Naming follows the module's existing split: array-taking functions are named for
the observable (`lcp_da_cross_section`), model-taking functions for the quantity
(`local_complex_potential`).

```python
@dataclass(frozen=True)
class ResonanceLevels:
    energies:    NDArray[complex128]   # E_v - i Gamma_v/2, ascending Re
    widths:      NDArray[float64]      # Gamma_v = max(0, -2 Im E_v)
    states:      NDArray[complex128]   # (n_levels, grid.n) DVR coefficients,
                                       # c-normalized: sum(c_i^2) = 1
    residuals:   NDArray[float64]      # two-angle stability residual per level
    real_weight: NDArray[float64]      # fraction of |c|^2 inside the real region
    golden_rule: NDArray[complex128]   # E_v^(0) - i Gamma_v^(1)/2 comparator


def lcp_resonance_levels(
    nuclear_grid_a, nuclear_grid_b, mu, Vd_a, Vd_b, Gamma, *,
    window=None, n_levels=None, rel_tol=1e-4, atol=1e-8, golden_rule=True,
) -> ResonanceLevels: ...


def resonance_levels(
    model, nuclear_grid_a, nuclear_grid_b, elec_grid_a, elec_grid_b, *,
    re_half_width=0.05, im_half_width=0.05, resid_tol=1e-3, **kwargs,
) -> ResonanceLevels: ...
```

`resonance_levels` calls the already-public `resonance_pole_walk` **once** on the
shared real points and assembles `V_d` onto both nuclear grids. The ~10-line
assembly currently at the tail of `local_complex_potential` is extracted into a
private helper used by both; `local_complex_potential`'s behaviour is unchanged.

Validation and errors:

- `mu > 0`; `Vd`/`Gamma` shapes match their grids; `nuclear_grid_a` and
  `nuclear_grid_b` must agree on the real nodes (asserted) or the eigenvalue
  comparison is meaningless.
- `window` defaults to `Re in [min Re V_d, max Re V_d]` over the **real** nodes
  (the full span of the anion curve, so neither the well bottom nor levels lying
  above the neutral dissociation limit `v0(inf) = 0` are cut) and
  `Im in [-max(Gamma), atol]`. Levels with `Im E > atol` are unphysical and are
  dropped with a warning.
- An empty window raises `ValueError` with the same message shape as
  `find_resonance_pole`.
- `n_levels` truncates to the lowest-`Re E` levels after selection; `None` returns
  all selected.

### `qscat-run` surface

Both entry points call `resonance_levels`.

**(a) Standalone observable kind.** `observables: [{kind: resonance_levels,
channels: N}]` on `methods: [lcp]`, `channels` being the number of levels to
report (omit → all selected). This run needs **no `energies` block and no incident
wavepacket**. Two targeted changes in `config.py`:

- `energies` becomes optional when every requested observable is
  `resonance_levels` (it is currently required);
- the `(molecule, observable.kind)` validity gate accepts `resonance_levels` for
  the LCP-capable molecules, using the same capability check the existing `lcp`
  method path applies.

**(b) Opt-in artifact.** `artifacts.resonance_levels: true`, alongside the existing
`eigenstates` flag. On an LCP `da` run this reuses the `V_d`/`Gamma` already
computed and adds only the nuclear diagonalizations.

**Second-angle grid.** The tail angle of grid *b* defaults to `angle_a - 10` degrees
(eMoScat's decks pair 44/35 and 40/30), overridable as `grid.nuclear_angle_b`. Real
segments and quadrature are copied verbatim from grid *a* — the shared real region
is what makes the comparison valid, so it is derived, not configurable.

**Artifacts.** A new `ResonanceLevelsRun` result dataclass in `runner.py` —
`(label, levels: ResonanceLevels, R_axis, Vd, Gamma)`, wrapping the library
dataclass rather than restating its fields, and *not* an overload of
`EigenStates`, whose `energies` field is `float64` and stays so:

- `resonance_levels.csv` — `v, Re_E, Gamma_v, residual, real_weight, Re_E0,
  Gamma_v_1`;
- `resonance_levels.npz` — complex energies, nuclear eigenvectors, the real `R`
  axis, and the `V_d`/`Gamma` curves;
- `resonance_levels.png` — the levels as horizontal bars on the `V_d(R)` curve,
  bar thickness proportional to `Gamma_v`.

## Validation

Two exact analytic oracles, one convergence study, one consistency check, and one
cross-model comparison.

1. **Zero-width limit.** `V_d` = a bare Morse curve, `Gamma = 0`. The bound levels
   must reproduce the analytic Morse spectrum
   `E_n = -D (1 - alpha (n + 1/2) / sqrt(2 mu D))^2` with `Im E_v ~ 0`. Gates the
   kinetic term, the mass, the boundary conditions and the level ordering at once.
2. **Constant-width exactness.** With `Gamma = Gamma_0` constant, a constant
   imaginary term commutes with everything, so `E_v = E_v^real - i Gamma_0/2` must
   hold to round-off. Gates the complex machinery independently of grid
   convergence.
3. **Angle and grid independence.** Energies stable under the theta pair and under
   h/p refinement. The `residuals` field is the per-level metric; a level whose
   residual does not fall under refinement is reported as unconverged rather than
   silently returned.
4. **Golden-rule consistency.** Full complex vs perturbative comparator must agree
   in the narrow-level limit; the divergence for broad levels is recorded, not
   asserted away.
5. **`Gamma(R)` support condition.** The thesis states (§1.5) that the imaginary
   part is nonzero only where `V0(R) < E_res(R)` — i.e. `Gamma = 0` outside the
   autodetachment region, to the right of the crossing `R_c`. This is a property
   of the *existing* `local_complex_potential` output; assert it as a precondition
   of the level solver, since a spurious `Gamma` tail beyond `R_c` would
   contaminate every width.
6. **Reproduction of the thesis `omega_j`, milestone 3.** The golden-rule
   comparator's `E_v^(0)` are the thesis's `omega_j` (see above), so they are
   checkable directly against the values behind Figs. 3.23/3.26 — data available
   from the author. This is the tightest available gate on the real parts.
7. **Cross-section peak correspondence, milestone 3.** Per the thesis §3.4.2, in
   the **elastic** cross section the boomerang maxima sit at approximately the
   `omega_j` energies, whereas in **VE 0→1** they are displaced, which is what
   makes those structures asymmetric. So the comparison is: `Re E_v` against the
   elastic-channel peak positions of the exact-2D `sigma(E)` already gated to
   Houfek — and, as a documented *non*-match, against the displaced VE 0→1 peaks.
   A cross-model comparison with a stated tolerance, not a tight gate.
8. **Lifetimes, milestone 3.** `tau_v = 1/Gamma_v` against the thesis's
   TD-estimated formation times for NO: the lowest state has a lifetime above
   30 000 a.u. (so `Gamma_0 < 3.3e-5` Ha), and the level behind the first VE 0→1
   peak — identified there as the *second* vibrational state — forms fully at
   `t > 10 000` (so `Gamma ~ 1e-4` Ha). Order-of-magnitude bounds, but genuine
   numbers from the reference and the only external check that touches the
   imaginary parts at all.

**Stated limitation:** eMoScat itself commits no numeric output — no `.dat`/`.asc`
files — and never computed a complex nuclear eigenvalue. Items 6-8 lean on the
thesis and on data held by the author, and item 7 compares against qModeling's own
gated exact-2D solver. No prior computed width exists anywhere to check
`Gamma_v` against directly; items 8's bounds are inferred from reported
time-dependent behaviour. The docs must say this plainly.

## Milestones

- **M1** — `match_angle_stable`, `find_resonance_pole` rewritten over it,
  `lcp_resonance_levels`/`resonance_levels`, validations 1-4. Deliverable: complex
  levels and nuclear wavefunctions per molecule.
- **M2** — both `qscat-run` entry points, the three artifacts, an example YAML, and
  `docs/physics/lcp-resonance-levels.md`.
- **M3** — the thesis grids as `*:thesis` preset variants; validations 6-8
  (reproduce the thesis `omega_j`, the elastic-peak correspondence against the
  exact 2-D, the NO lifetime bounds), each run on both decks per above.

Each milestone gets a `physics-reviewer` pass before merge, per the repo's
promotion rule.

## The two grid parametrizations — run both

Thesis Tables 2.1/2.2 disagree with the eMoScat JSON decks that
`qscat_run.presets` and `validation/diatomic/test_da_grid.py` lock. F2 nuclear:
thesis `nq = 20`, `theta = 35`, real elements from 1.8 — deck: order 12,
`theta = 25`, from 2.0. N2 nuclear: thesis `nq = 12`, real region to 5.0 — deck:
order 14, to 12.0. Both are plausible discretizations of the same model.

For H2+ there are **three**: eMoScat's JSON deck, thesis Table 2.2 (electronic
`nq = 8`, `theta = 5`; nuclear `nq = 8`, `theta = 22`), and Hvizdoš et al. Table I
(`nq = 6`, `theta = 20` on both coordinates).

M3 therefore computes the levels on **both** decks (all three for H2+) and reports
the difference.
The `omega_j` are grid-convergent quantities, so agreement turns the discrepancy
into positive evidence that the reported levels are converged; disagreement is
itself a finding and identifies which deck under-resolves. This needs the thesis
grids added as preset variants (`*:thesis`) — additive, leaving the existing
locked presets and their byte-identical test untouched.

Not a blocker for M1/M2: the library function takes whatever grids it is given.

## Out of scope

- **The diabatic-state projection and level populations.** The thesis's other use
  of these levels — `Psi_d(R,t) = int dr Psi(R,r,t) phi_d*(R,r)` (§1.6, Eq. 1.70)
  and the populations `|<omega_i|Psi_d(t)>|^2` (Fig. 3.27) — needs the discrete
  state `phi_d` with its suppression function `f(r)` (Eq. 1.69) and eMoScat's
  phase-locking across `R` (`ModelLCP.cpp:361-395`). That is the natural follow-on
  once the levels exist, and is deliberately a separate spec.
- Any change to `local_complex_potential`'s physics or to the exact 2-D solvers.
- A full 2-D complex-scaled diagonalization (the direct Siegert states of the 2-D
  model without the BO reduction). That is the natural next comparison and is
  deliberately deferred; nothing here forecloses it.
- The nonlocal (Feshbach `QHQ`/`PHP`) route, `reference/eMoScat/source/module_NRM.cpp`,
  which carries its own unresolved `FIXME` about a failed Fortran comparison.
- GPU/CUDA and Rust optimization. The nuclear problem is 1-D and dense-diagonalizes
  in about a second; it is not a hot path.
