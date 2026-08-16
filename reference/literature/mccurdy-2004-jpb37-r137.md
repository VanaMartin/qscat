# C. W. McCurdy, M. Baertschy, T. N. Rescigno, J. Phys. B 37, R137 (2004) — Solving the three-body Coulomb breakup problem using exterior complex scaling

**Source:** `reference/literature/mccurdy-2004-jpb37-r137.pdf` (gitignored) · DOI 10.1088/0953-4075/37/17/R01
**Pagination:** per-article, R137–R187; extractor page N = printed page `R(136+N)` (offset +136,
checked at extractor page 1 → R137 and extractor page 5 → R141).

## Why this repository cares

This is the standard topical-review reference for exterior complex scaling (ECS) — broader in
scope than `reference/literature/rescigno-2000-pra62-032706.md` (already noted), which this
review itself cites as ref. [38] and explicitly defers to for the FEM-DVR kinetic-energy matrix
elements ("the kinetic energy ... has a blocked structure that can be deduced from the
definitions of the FEM-DVR basis functions [38]", p. R152). What this review adds *beyond* that
note is the surrounding formal theory: the general complex-contour formalism `q(r)` that the
`qscat.ecs.ecs_map` sharp-ECS map is a special case of, the explicit condition under which a
resonance wavefunction decays under rotation, the driven-equation derivation with a
finite-range-truncated potential and the "order of limits" argument that justifies solving on a
*finite* `R0` at all, and a Wronskian amplitude-extraction formula that is the conceptual
ancestor of `qscat.core.td_extractors.Flux`. It does not add anything to the FEM-DVR mechanics
itself — for that, see the Rescigno 2000 note.

## What this repository uses

| Fact | Locator | Used by |
|---|---|---|
| General complex-contour formalism `R(r) = integral_0^r q(r') dr'`, with `q(r) -> 1` as `r->0` and `q(r) -> e^{i eta}` as `r->infinity` (covers both smooth and sharp ECS) | p. R144, Eqs. (17)-(18) | the general form `qscat.ecs.ecs_map` specializes to sharp ECS (`q(r)=1` for `r<=R0`, `q(r)=e^{i theta}` for `r>R0`) |
| Resonance-wave decay condition: `exp(i k r)` on the rotated contour decays iff the scaling angle `eta` exceeds `alpha = arg(k_res)` | p. R142, prose above Eq. (10) / Fig. 1(c) discussion | the formal reason ECS renders a resonance wavefunction square-integrable; underlies `docs/physics/femdvr-ecs.md`'s "exposes resonances as isolated, theta-stationary complex poles" claim |
| Sharp ECS contour `r -> r` for `r<=R0`, `R0+(r-R0)e^{i eta}` for `r>R0` | p. R143, Eq. (10) | identical to `qscat.ecs.ecs_map` (same functional form as Rescigno-McCurdy 2000 Eq. 21) |
| Driven scattered-wave equation solved with the potential truncated at `R0`, `V_{R0}(r) = V(r)` for `r<=R0`, `0` for `r>R0` | p. R147, Eqs. (30)-(31) | the formal justification for solving the driven equation on a *finite* ECS-truncated grid — `qscat.core.driven.ve_cross_section`, `qscat.core.dissociation.da_cross_section` |
| "Order of limits" theorem: `r->infinity` first (imposing the ECS outgoing condition), then `R0->infinity`, gives the correct scattered wave; larger `R0` needed nearer threshold | p. R148, prose above Eq. (34) | the formal correctness argument behind convergence-in-`R0` grid practice (per-molecule grids; `discretisation-tuner`'s tail sizing) |
| Wronskian amplitude-extraction formula, `A = -(1/k) W(j_l(kr), psi_sc(r))|_{r=R0-eps}` | p. R148, Eq. (34c) | conceptual ancestor of the fixed-surface Wronskian-flux extractor `qscat.core.td_extractors.Flux` (ported from eMoScat's `FluxTestFunction2d`, not derived directly from this paper, but the same formal object) |
| Sharp ECS (grid / compact-support basis) handles the derivative discontinuity at `R0` exactly; smooth ECS complicates matrix elements and is numerically worse for analytic basis sets | p. R145-R146 Sec. 2.3, R148-R149 Sec. 3 | restates, more generally, the same conclusion as Rescigno-McCurdy 2000 that motivates `qscat.dvr`'s sharp, node-aligned ECS design |
| Two-potential formalism note: a Coulomb potential can serve as the reference potential in electron-ion scattering | p. R148, prose after Eq. (34c) | background relevant to a possible future Coulomb-reference treatment of `IonicResonanceModel` (H2+); **not implemented** — `qscat.special.coulomb` instead evaluates Coulomb wavefunctions directly (mpmath) |

## Equations

```
R(r) = integral_0^r q(r') dr'                                        p. R144, Eq. (17)
q(r) -> 1 as r -> 0,  q(r) -> e^{i eta} as r -> infinity              p. R144, Eq. (18)   [general ECS contour condition]
r -> r,                       r <= R0                                p. R143, Eq. (10)   [sharp ECS map]
  -> R0 + (r - R0) e^{i eta},  r > R0
psi_+ = psi_sc + psi_0                                                p. R147, Eq. (26)
(E - H) psi_sc = (H - E) psi_0                                        p. R147, Eq. (27a)  [driven equation]
V_{R0}(r) = V(r) for r<=R0, 0 for r>R0                                p. R147, Eq. (30)   [finite-range truncation]
(E - H) psi_sc_{R0}(r) = V_{R0}(r) psi_0(r)                           p. R147, Eq. (31)
psi_sc_{R0}(R(r)) -> 0 as r -> infinity                               p. R147, Eq. (32)   [ECS outgoing condition]
psi_sc_{R0}(R(r)) -> psi_sc(R(r)) as R0 -> infinity                   p. R147, Eq. (33)   [order-of-limits convergence]
A = -2*mu/k <j_l(kr)|V(r)|psi_+(r)>                                   p. R148, Eq. (34a)
  = (2*mu/k) <j_l(kr)|T-E|psi_sc(r)>_{R0}                             p. R148, Eq. (34b)
  = -(1/k) W(j_l(kr), psi_sc(r))|_{r=R0-eps}                          p. R148, Eq. (34c)  [Wronskian amplitude extraction]
```

Atomic units throughout (implicit in the mass `mu`/momentum `k` notation, p. R146 Eq. 14) — the
same convention `libs/qscat/qscat/units.py` uses.

## Parameters and numeric values

None — this is a methods review, not a molecular-model paper. No constants checked against the
repo.

## Findings and limits

- The paper's central methodological claim (restated from Rescigno-McCurdy 2000, generalized to
  arbitrary contours): sharp ECS on a grid or compact-support basis has no numerical
  difficulties, while smooth ECS complicates the matrix elements and introduces numerical error
  that has been "difficult to eliminate" for analytic basis sets — p. R146, end of Sec. 2.3;
  p. R148-R149, Sec. 3. This is the same conclusion `docs/physics/femdvr-ecs.md`'s "Method"
  section relies on.
- The "order of limits" argument (p. R148) is the formal correctness proof for the practice of
  solving the driven equation on a finite `R0` and increasing it until convergence — it explains
  *why* a finite-grid ECS calculation converges to the exact scattering solution, rather than
  merely asserting it.
- This review explicitly treats the FEM-DVR construction (Sec. 3.2, Eqs. 41-49) as a summary of
  Rescigno-McCurdy 2000 (its ref. [38]) and defers the kinetic-energy matrix-element derivation
  to that paper — confirming the existing `rescigno-2000-pra62-032706.md` note is the primary
  FEM-DVR source, not this one.
- The review's own worked problem — three-body Coulomb breakup (electron-impact ionization of
  hydrogen, double photoionization of helium) — is a genuinely different physical problem from
  qscat's 2-electron/nuclear-coordinate diatomic model: it involves a singular, non-diagonal
  `1/|r1-r2|` electron-electron potential and multi-electron amplitude-extraction machinery
  (surface-integral methods, the volume-dependent overall phase) that qscat does not use. See
  "Not used here."

## Terminology map

| Paper symbol | qModeling name |
|---|---|
| `R(r)` / `q(r)` (general complex contour / its derivative) | `qscat.ecs.ecs_map`'s sharp special case (`z(x)`) |
| `eta` (ECS rotation angle) | `theta_deg` in `qscat.ecs.ecs_map` |
| `R0` (ECS pivot radius) | `GridSpec.R0` / `FemDvrEcsGrid.R0` |
| `psi_0`, `psi_sc`, `psi_+` (incident/scattered/total wave) | same driven-equation split used by `qscat.core.driven.ve_cross_section` |
| `j_l(kr)`, `h_l^+(kr)` (Riccati-Bessel / outgoing Riccati-Hankel) | `qscat.special.radial` |
| `W(a,b)` (Wronskian amplitude extraction) | conceptual ancestor of `qscat.core.td_extractors.Flux`'s fixed-surface flux transform |

## Not used here

- Sec. 3.1 (finite differences), 3.3 (B-splines), 3.4 (direct numerical integration) — alternative
  ECS numerical implementations; `qscat.dvr` uses FEM-DVR exclusively (Sec. 3.2, already covered
  by the Rescigno 2000 note).
- Sec. 4 (amplitude and cross-section extraction for multi-electron breakup: surface-integral
  methods, the volume-dependent overall phase in the formal theory of ionization, Peterkop's
  asymptotic form) — specific to the three-body Coulomb breakup problem, not the 2-electron/
  nuclear-coordinate diatomic resonance model this repo implements.
- Sec. 5-7 (numerical implementation details for e-H and e-He2 double photoionization, comparison
  with experiment and other theories) — application results, not method or formula sources for
  this repo.
- The multi-electron singular potential `1/|r1-r2|` treatment (single-centre expansion, Sec. 3.2
  continued) — qscat's 2-D model has no analogous singular interaction term.
