# C. W. McCurdy, C. K. Stroud, Comput. Phys. Commun. 63, 323 (1991) — Eliminating wavepacket reflection from grid boundaries using complex coordinate contours

**Source:** `reference/literature/mccurdy-1991-cpc63-323.pdf` (gitignored) · DOI 10.1016/0010-4655(91)90259-N
**Pagination:** printed 323–330; extractor page N = printed page 322+N (offset zero, checked at page 1 and page 8).

## Why this repository cares

This is the paper that first demonstrated **exterior complex scaling (ECS) eliminates
wavepacket reflection from grid boundaries in time-dependent propagation** — analytically,
not by tuning an absorbing potential. It is a precursor, not a literal dependency: the
repo's actual FEM-DVR-ECS grid construction and coordinate map come from Rescigno &
McCurdy (2000) (`reference/literature/rescigno-2000-pra62-032706.md`), and no code here
implements this paper's specific tent-function/Peaceman–Rachford finite-element scheme.
But the *architectural claim* this paper establishes — propagate a wavepacket under a
Hamiltonian continued onto a complex exterior-scaling contour, and outgoing flux vanishes
before it can reflect off the outer boundary, without truncating the physical region — is
exactly the mechanism `qscat.evolution`'s time-dependent route relies on: a Crank-
Nicolson/Padé stepper (`make_cn_stepper`/`make_pade_stepper`) propagating `H` built on
`qscat.dvr`'s FEM-DVR-**ECS** grid, where (per `docs/physics/n2-2d-td-cross-section.md`)
"`H_2D`'s absorbing ECS tail makes `||psi(t)||` genuinely decay" as the wavepacket leaves
the interaction region — the same effect this paper analytically demonstrates 15 years
earlier, in a simpler one-dimensional finite-element setting. Without this source, the
repo's ECS-grid-plus-CN/Padé-propagator TD architecture would be an empirically-observed
convenience rather than a method with a proven foundation.

## What this repository uses

| Fact | Locator | Used by |
|---|---|---|
| Complex-scaling a coordinate near a grid boundary analytically prevents wavepacket reflection from that boundary | p. 323, Abstract; p. 329, Sec. 5 | conceptual foundation for `qscat.dvr`/`qscat.ecs`'s complex tail combined with `qscat.evolution`'s TD steppers — see `docs/physics/n2-2d-td-cross-section.md` |
| Exterior scaling contour: identity for `x < x_M`, `x_M + (x-x_M)e^{iθ}` for `x >= x_M` — a sharp (piecewise-linear) bend, not a smooth turn | p. 326, Eq. (13), Fig. 1 | the same functional form as `qscat.ecs.ecs_map`'s `R0 + (r-R0)e^{iθ}` tail (formalized later, with the ECS radius pinned to a FEM node, by Rescigno & McCurdy 2000) |
| Mechanism: the outgoing-packet phase factor decays like `exp(-p0 sin(θ) x)` beyond `x_M` on the complex contour, vanishing before reaching the outer grid edge if the contour is long enough | p. 329, prose around Eq. (17) | the qualitative reason `qscat.dvr`'s complex tail (element angle + length) must be long/steep enough to absorb the outgoing flux — the repo's `qscat.tuning.ecs` module (double-ECS-capped angle + exponential tail growth) automates this sizing, not derived from this paper's formulas |
| Complex-symmetric (non-Hermitian-conjugate) propagator for a complex-scaled `H`: `exp(-iHΔt) ≈ (S + iHΔt/2)^{-1}(S - iHΔt/2)`, with the complex conjugate replaced by transpose-without-conjugation | p. 327, Eq. (25) | the same Cayley/Crank-Nicolson form `qscat.evolution.make_cn_stepper`/`make_sparse_cn_stepper` implement for a general non-Hermitian complex `H` (order-1 Padé) |
| Numerical demonstration: 700 complex finite elements on a 30-bohr ECS-scaled grid (scaling starting at 17.5 bohr, θ=15°) reproduce a 2400-element real-coordinate calculation on a 100-bohr grid to 4 significant figures, even though most of the wavepacket has crossed the scaling boundary | p. 328–329, Figs. 2–3 | independent, much-earlier confirmation of the size/accuracy tradeoff the repo's own ECS+CN/Padé TD convergence studies (`n2_2d_td_cross_section`) reproduce on the diatomic model; not itself reused as a repo benchmark |

## Equations

```
H Psi(x,t) = i hbar dPsi(x,t)/dt                                      p. 324, Eq. (1)
lim_{eta->0} <f|(E-H+i eta)^-1|g> = sum_n <f|n><n|g>/(E-E_n)          p. 324, Eq. (4)   [discretized Green's function, does not converge as N grows]
exp(-iHt) ~ U_N(t) = sum_n |n> exp(-iE_n t) <n|                        p. 324, Eq. (5)
C(x) = x,                    0 <= x < x_M                              p. 326, Eq. (13)  [exterior scaling contour]
     = x_M + (x-x_M) e^{i theta},   x_M <= x < infinity
Psi(C(x),t) = exp{ i a_t [C(x)-x_t]^2 + i p_t [C(x)-x_t] + i gamma_t }  p. 326, Eq. (14)  [free Gaussian on the contour]
exp(i p_0 [C(x)-x_M]) -> exp(-p_0 sin(theta) x)   for x > x_M           p. 329, Eq. (17)  [decay mechanism beyond the scaling point]
a_n(z) : tent-function finite-element basis, tridiagonal H and S        p. 327, Eq. (18)
z_n = x_n,                      x_n < x_M                              p. 327, Eq. (22)  [nodes moved onto the contour]
    = x_M + (x_n - x_M) e^{i theta},   x_n >= x_M
<n|H(z)|m> = (1/2) integral_C a_n(z) H(z) a_m(z) dz                     p. 327, Eq. (23)  [complex contour integral for matrix elements]
exp(-iHdt) ~ (S + (i/2)H dt)^{-1} (S - (i/2)H dt)                       p. 327, Eq. (20)  [Peaceman-Rachford / Crank-Nicolson]
exp(-iHdt) ~ (S + (i/2)H dt)^{-1} (S - (i/2)H dt), complex-symmetric    p. 327, Eq. (25)  [modified for complex-scaled H, conjugate -> transpose]
```

Atomic units throughout (`ħ = m_e = e = 1`, p. 324) — matches `libs/qscat/qscat/units.py`;
no conversion needed.

## Parameters and numeric values

No model parameters (masses, potential constants) are published here — this is a methods
paper with a single illustrative 1-D test problem (`V(x) = -1/x`, `0 < x < ∞`; incident
packet `a=10, k=-2, d=1`, p. 328, Eq. (28)), not a molecular model the repo's constants
derive from. Not applicable to compare against repo code.

## Findings and limits

- Central claim (p. 329, Sec. 5): exterior complex scaling formally requires the potential
  be zero beyond the scaling radius `x_M` for arbitrary accuracy, but in practice it is
  sufficient that the outgoing wavepacket "not be appreciably scattered back towards
  smaller x by the action of the potential beyond `x_M`" — i.e. the scaling radius need
  only be placed past the range of the interaction, not past the whole potential's support.
- The scaling angle `θ` trades off against the range of outgoing momenta the calculation
  must resolve: smaller momenta need larger `θ` (p. 329) — the same qualitative tradeoff
  `qscat.tuning.ecs`'s automated angle/tail sizing addresses for the repo's grids, though
  via its own calibrated procedure, not this paper's formulas.
- Stated scope (p. 330, Sec. 5 "Discussion"): the paper works a 1-D, single-particle
  example for illustration only; the authors state multidimensional applications "have
  proven to be equally successful and will be reported in a forthcoming publication" — the
  seed of what becomes the 2-D FEM-DVR-ECS diatomic-scattering model this repository
  implements (via Rescigno & McCurdy 2000 and Houfek, Rescigno & McCurdy 2006).
- The paper explicitly frames itself as connecting an established *time-independent*
  technique (complex coordinates / complex basis functions for Green's-function matrix
  elements, p. 325) to the *time-dependent* propagation problem — the repo's TD route
  (`qscat.core.time_dependent`) is exactly this connection exploited on the diatomic model,
  cross-validated against the TI route (`qscat.core.driven`) as an independent oracle.

## Terminology map

| Paper symbol | qModeling name |
|---|---|
| `C(x)` (exterior scaling contour) | `qscat.ecs.ecs_map` (same piecewise form, later pinned to a FEM node by Rescigno & McCurdy 2000) |
| `x_M` (scaling radius) | `R0` in `qscat.ecs.ecs_map` / `GridSpec.R0` |
| `θ` (scaling angle) | `theta_deg` in `qscat.ecs.ecs_map`; `ElementSpec.angle_deg` |
| `a_n(x)` (real tent-function finite element) | superseded in the repo by the Lobatto/GLL shape functions of `qscat.dvr` (Rescigno & McCurdy 2000), not the simple tent basis used here |
| Eq. (25)'s complex-symmetric Crank-Nicolson propagator | `qscat.evolution.make_cn_stepper` (order 1) / `make_pade_stepper` (higher order, via van Dijk & Toyama 2007 and Formánek, Váňa & Houfek 2010 — see `reference/literature/vandijk-2007-pre75-036707.md` and `reference/literature/formanek-2010-aip1281-667.md`) |

## Not used here

- Sec. 2's detailed derivation connecting discretized-spectrum Green's-function matrix
  elements (Eqs. 4–12) to the complex-coordinates literature — read for the physical
  motivation, not itself a formula the repo implements.
- The specific tent-function finite-element scheme (Eqs. 18–24) and its tridiagonal
  kinetic/overlap matrix elements — the repo's `qscat.dvr` uses the Gauss-Lobatto
  DVR/Lobatto-shape-function construction of Rescigno & McCurdy (2000) instead, not this
  paper's simpler piecewise-linear elements.
- The Bottcher electron-impact-ionization application context (refs [22]–[23], p. 328) and
  the paper's own numerical example (`V(x) = -1/x`, Figs. 2–3) — read for qualitative
  confirmation only; not reproduced or reused as a repo benchmark.
- The references list's earlier complex-coordinates/complex-basis-function literature
  (refs [4]–[20]) — background lineage, not sources this note or the repo separately track.
