# D. Hvizdoš, master's thesis, Charles University, Prague (2016) — Two-dimensional model of dissociative recombination

**Source:** `reference/literature/hvizdos-2016-thesis.pdf` (gitignored) ·
<https://dspace.cuni.cz/handle/20.500.11956/96080>
**Pagination:** the PDF carries 7 unnumbered front-matter pages plus a blank
page before the body begins; extractor page **E** = printed page **P + 8**
throughout (checked at the table of contents, extractor p. 9, against every
section/page number it lists — Introduction p.3/extractor 11 through
Bibliography p.31/extractor 39 — and again at the last body page,
Conclusion p.29/extractor 37: no drift). All locators below are **printed**
page numbers.

## Why this repository cares

This is the **first time-independent solution of the e⁻+H₂⁺ model** —
the source `qscat.core.dissociation.dr_cross_section` descends from (see
`docs/physics/h2plus-dr.md`). It supplies three things this repo leans on
directly: the model parameters for H₂⁺ (§1.2.1, Table 1.1, the third of
three independent published sources agreeing on the reduced mass μ=918.076
against eMoScat's contradicted 918.25); the statement that the **driven
Schrödinger equation** is solved in place of Lippmann-Schwinger on the
ECS grid (§1.3.3, Eq. 1.62) — exactly what `qscat.core.driven` implements;
and the Coulomb-function machinery (§1.1) that `qscat.special.coulomb`
provides for the ionic channel. Without this source, `qscat.model.H2P`'s
parameters and the "why a driven equation, not Lippmann-Schwinger" design
choice in `qscat.core` would be unfounded.

## What this repository uses

| Fact | Locator | Used by |
|---|---|---|
| Energy-normalized spherical Coulomb functions `F_l`, `G_l`, combined `u±_l = G_l ± iF_l` | p. 5, Eq. (1.1)-(1.5) | `qscat.special.coulomb.coulomb_f_en`/`g_en`/`h1_en` |
| Model Hamiltonian `H = H0 + V_int(R,r) = H_ion_0 + H_el_0 + V_int(R,r)` | p. 5-6, Eq. (1.6) | `qscat.model.IonicResonanceModel.hamiltonian` (via `qscat.core`) |
| Ion vibrational Hamiltonian `H_ion_0 = -1/(2μ) d²/dR² + V0(R)` | p. 6, Eq. (1.7) | `IonicResonanceModel.v0`, nuclear kinetic term |
| Electronic Hamiltonian with Coulomb + centrifugal terms `H_el_0 = -1/2 d²/dr² - 1/r + l(l+1)/2r²` | p. 6, Eq. (1.8) | `IonicResonanceModel.surface`'s `charge/r` + centrifugal terms |
| Interaction potential `V_int(R,r) = -λ1(R) e^{-λ2(R) r²}/r` | p. 6, Eq. (1.9) | `IonicResonanceModel.v_int` |
| Lippmann-Schwinger form of `ψ+_E` | p. 7, Eq. (1.13)-(1.15) | conceptual basis for `qscat.core.driven` |
| Incident state `ψ_in = χ_vi(R) φ_ki,l(r)`, energy-normalized Coulomb `φ_ki,l` | p. 7, Eq. (1.16)-(1.18) | `qscat.core.channels.channel_vector(charge=-1)` |
| VE/DR outgoing boundary conditions | p. 7-8, Eq. (1.21)-(1.22) | `qscat.core.dissociation.dr_cross_section`'s exit-channel construction |
| Rydberg bound-state equation `[H_el_0 + V_int^∞(r)] ρ_n(r) = E_n ρ_n(r)` | p. 8, Eq. (1.23) | `qscat.core.dissociation.anion_electronic_states` |
| T-matrix definitions `T_VE`, `T_DR`, DR channel potential `V_DR = V_int - V_int^∞ + V0` | p. 8, Eq. (1.25)-(1.27) | `qscat.core.dr_cross_section`'s `v_dr_diag` |
| Cross-section normalization `σ_VE = σ_DR = (4π³/k_i²) \|T\|²` | p. 8, Eq. (1.30)-(1.31) | `qscat.core.dissociation.da_cross_section`/`dr_cross_section`'s `4π³\|T\|²/2E` |
| e⁻+H₂⁺ model form: Morse `V0(R)`, σ-capture `λ1(R)`, `λ2(R)=1/3`, `V_int^∞(r)` limit | p. 9, Eq. (1.32)-(1.35) | `qscat.model.IonicResonanceModel.v0`/`.v_int` |
| Table 1.1 — e⁻+H₂⁺ model parameters | p. 9, Table 1.1 | `qscat.model.library.H2P` |
| Driven Schrödinger equation solved instead of Lippmann-Schwinger: `(E-H)\|ψ_sc⟩ = V_int\|ψ_in⟩` | p. 15, Eq. (1.62) | `qscat.core.driven.ve_cross_section`/`dissociation.dr_cross_section` (the driven-equation solve) |
| ECS coordinate map | p. 15, Eq. (1.60) | `qscat.ecs.ecs_map` |
| ECS bending-angle bound `θ < π/8` from the `V_int` `tanh` argument | p. 14 (prose) | the nuclear-ECS-angle bound in `qscat.model.ionic.IonicResonanceModel` (`max_nuclear_ecs_angle_deg`), cross-checked against the tighter π/8 bound from Hvizdoš et al. 2018 |
| Final FEM-DVR-ECS grid (Table 2.1) | p. 17, Table 2.1 | not used directly — a third, independent H₂⁺ grid parametrization; see Not used here |
| σ_DR/σ_VE curves over [0, 0.05] Ha, three open DR + five open VE channels | p. 17-18, Fig. 2.1-2.2 | qualitative cross-check only; no golden numeric data extracted (see Findings) |
| Convergence tests (ECS bending point, electronic/nuclear grid quadrature order) | p. 19-21, §2.2, Fig. 2.3-2.5 | motivates but is not itself gated against — this repo's own convergence gate is `validation/h2plus/` |
| Rydberg potential curves `V_n(R) = V0(R) + E_el_n(R)` and structure interpretation | p. 9, 22-27, Eq. (1.36), §2.3, Fig. 1.3, 2.6-2.12 | narrative basis for `docs/physics/h2plus-dr.md`'s Rydberg-series discussion |

## Equations

```
d^2 w_l(eta,rho)/d rho^2 + (1 - 2 eta/rho - l(l+1)/rho^2) w_l = 0          p. 5, Eq. (1.2)
F_l(eta,rho) --rho->inf--> sin[rho - eta log(2rho) - l pi/2 + arg Gamma(l+1+i eta)]   p. 5, Eq. (1.3)
G_l(eta,rho) --rho->inf--> cos[rho - eta log(2rho) - l pi/2 + arg Gamma(l+1+i eta)]   p. 5, Eq. (1.4)
u_l^+/- = G_l +/- i F_l                                                    p. 5, Eq. (1.5)
H = H0 + V_int(R,r) = H_ion_0 + H_el_0 + V_int(R,r)                       p. 5-6, Eq. (1.6)
H_ion_0 = -1/(2 mu) d^2/dR^2 + V0(R)                                       p. 6, Eq. (1.7)
H_el_0 = -1/2 d^2/dr^2 - 1/r + l(l+1)/(2 r^2)                              p. 6, Eq. (1.8)
V_int(R,r) = -lambda1(R) e^{-lambda2(R) r^2} / r                          p. 6, Eq. (1.9)
V(R,r) = V0(R) - 1/r + l(l+1)/(2 r^2) + V_int(R,r)                        p. 6, Eq. (1.10)
V_int^inf(r) = lim_{R->inf} V_int(R,r)                                     p. 6, Eq. (1.11)
|psi+_E> = |psi_in> + G+(E) V_int |psi_in> = |psi_in> + |psi_sc>          p. 7, Eq. (1.13)
psi_in(R,r) = chi_vi(R) phi_{ki,l}(r)                                      p. 7, Eq. (1.16)
phi_{ki,l}(r) = sqrt(2/(pi ki)) F_l(-1/ki, ki r)                           p. 7, Eq. (1.18)
E = E_vi + ki^2/2                                                          p. 7, Eq. (1.20)
psi_sc --r->inf--> sqrt(2/(pi ki)) sum_vf f_VE_{vi->vf} chi_vf(R) u_l^+(kef r)   p. 7, Eq. (1.21)
psi_sc --R->inf--> sqrt(2/(pi ki)) sum_n f_DR_{vi->n} rho_n(r) Kn R h_l^+(Kn R)  p. 7, Eq. (1.22)
[H_el_0 + V_int^inf(r)] rho_n(r) = E_n rho_n(r)                            p. 8, Eq. (1.23)
E = E_n + Kn^2/(2 mu)                                                      p. 8, Eq. (1.24)
T_VE_{vi->vf}(E) = <chi_vf phi_{kf,l}|V_int|psi+_E> = sqrt(ki/kef) f_VE/pi p. 8, Eq. (1.25)
T_DR_{vi->n}(E) = <psi_DR_n|V_DR|psi+_E> = sqrt(ki/(mu Kn)) f_DR/pi        p. 8, Eq. (1.26)
V_DR(R,r) = V_int(R,r) - V_int^inf(r) + V0(R)                              p. 8, Eq. (1.27)
psi_DR_n(R,r) = sqrt(2 mu/(pi Kn)) sin(Kn R) rho_n(r)                      p. 8, Eq. (1.28)
sigma_VE_{vi->vf}(E) = (4 pi^3 / ki^2) |T_VE_{vi->vf}(E)|^2                p. 8, Eq. (1.30)
sigma_DR_{vi->vf}(E) = (4 pi^3 / ki^2) |T_DR_{vi->vf}(E)|^2                p. 8, Eq. (1.31)
V0(R) = beta1 (e^{-2 beta2(R-R0)} - 2 e^{-beta2(R-R0)})                    p. 9, Eq. (1.32)
lambda1(R) = alpha1 (1 - tanh[(alpha2-R-alpha3 R^4)/7]) (tanh[R/alpha4])^4 p. 9, Eq. (1.33)
lambda2(R) = 1/3                                                           p. 9, Eq. (1.34)
V_int^inf(r) = -2 alpha1 e^{-r^2/3} / r                                    p. 9, Eq. (1.35)
[H_el_0 + V_int(R,r)] rho_n(R,r) = E_n(R) rho_n(R,r)   [fixed-R eigenproblem]   p. 9, Eq. (1.36)
r'(r) = r for r<R0; R0 + (r-R0) e^{i theta} for r>=R0   [ECS map]          p. 13, Eq. (1.60)
(E-H)|psi_sc> = V_int|psi_in>              [driven Schrodinger equation]  p. 15, Eq. (1.62)
```

Atomic units throughout (energies in hartree, lengths in bohr), matching
`libs/qscat/qscat/units.py` — no conversion needed.

## Parameters and numeric values

Table 1.1, p. 9 (atomic units), the e⁻+H₂⁺ model:

| Parameter | Value | Parameter | Value |
|---|---|---|---|
| μ | 918.076 | α1 | 1.6435 |
| l | 1 | α2 | 6.2 |
| β1 | 0.1027 | α3 | 0.0125 |
| β2 | 0.69 | α4 | 1.15 |
| R0 | 2.0 | | |

Checked against the repo:

```
grep -n "918.076\|0.1027\|1.6435" libs/qscat/qscat/model/library.py
```

`qscat.model.library.H2P` (`mu=918.076, ell=1, charge=-1, V0=0.1027, R0=2.0,
alpha=0.69, a1=1.6435, a2=6.2, a3=0.0125, a4=1.15`) matches Table 1.1
**exactly, field by field** (verified 2026-08-16). `IonicResonanceModel.v0`
and `.v_int` (`libs/qscat/qscat/model/ionic.py:54-83`) implement Eq. (1.32)
and Eq. (1.33) literally, including the fixed `λ2(R)=1/3` (Eq. 1.34) baked in
as the `exp(-r^2/3)` factor rather than a named field. This is the third of
three independent published sources agreeing on μ=918.076 (with Váňa 2017
Table 1.2 and Hvizdoš et al. 2018 §II A) against eMoScat's JSON deck value of
918.25, which all three contradict — the correction landed in
`qscat.model.library.H2P` on 2026-08-15 (see `docs/physics/h2plus-dr.md`).

Table 2.1, p. 17 — the thesis's own final FEM-DVR-ECS grid (`nq=6`, θ=20° on
both coordinates; electronic real region to 1300 bohr, nuclear real region to
12 bohr) is a **third, independent** H₂⁺ grid parametrization, distinct from
both the eMoScat deck (`validation/h2plus/config.py`, electronic real region
also to 1300 bohr but θ=5°/exp-growth tail) and the coarser grid in Hvizdoš
et al. 2018 Table I (`nq=6`, θ=20° — matching this thesis's angle exactly).
Not checked against the repo's grids field-by-field; noted as a convergence
cross-check, not a locked reference (see Not used here).

## Findings and limits

- The DR cross sections of the three open channels differ by several orders
  of magnitude; DR0 is negligible at all tested energies, DR1 dominates the
  sum almost everywhere — p. 18, Fig. 2.1.
- A finite grid can only resolve finitely many of the infinitely many Rydberg
  states accumulating below each vibrational-excitation threshold, so a band
  of energies just below each threshold is genuinely unconverged regardless
  of grid refinement (it shrinks but never vanishes) — p. 18-19.
- Convergence with respect to the electronic ECS bending point `R0` is good
  away from these threshold bands; the DR1 channel is well converged in
  electronic-grid quadrature order except in the same bands, though a
  lower-density grid starts diverging at higher energies (so a wider energy
  range would need a denser grid); the DR0 channel is the most sensitive to
  nuclear grid quadrature order because its cross section is smallest — p.
  19-21, §2.2, Fig. 2.3-2.5.
- Rydberg-state vibrational-level energies (from the potential curves
  `V_n(R)`) coincide closely with most cross-section peaks; agreement
  degrades for higher vibrational levels of the *lowest* Rydberg states,
  where the electron-nuclear coupling is strongest — p. 27, §2.3. Some
  features (e.g. a peak present only in DR1, absent from DR0 and elastic VE0,
  at 0.00034 Ha) are not explained and are left as future work (possibly a
  Ramsauer-Townsend effect) — p. 27.
- The thesis states its own scope limit explicitly: this is a first
  time-independent solution, with comparison against approximative methods
  (local complex potential, frame-transformation/quantum-defect theory) and
  extension to a two-channel (direct+indirect) DR model both left as future
  work — p. 29, Conclusion. `qscat.core.lcp` and the frame-transformation
  comparison in Hvizdoš et al. 2018 are the eventual continuations of that
  stated plan.
- No cross-section values are tabulated in the PDF (Fig. 2.1-2.12 are plots
  only); this repo's `dr_cross_section` oracle is validated standalone (see
  `docs/physics/h2plus-dr.md`'s Validation section) rather than gated against
  digitized figures from this thesis.

## Terminology map

| Thesis symbol | qModeling name |
|---|---|
| `mu` | `qscat.model.library.H2P.mu` |
| `l` | `.ell` |
| `beta1` (V0(R) depth) | `.V0` |
| `beta2` (V0(R) rate) | `.alpha` |
| `R0` | `.R0` |
| `alpha1` (lambda1 depth) | `.a1` |
| `alpha2` | `.a2` |
| `alpha3` | `.a3` |
| `alpha4` | `.a4` |
| `lambda2(R) = 1/3` | not a named field — baked into `.v_int`'s `exp(-r^2/3)` |
| `Z=1` (Coulomb charge) | `.charge = -1` (sign convention: `charge/r` term is `-Z/r`) |
| `V0(R)` | `.v0(R)` |
| `V_int(R,r)` | `.v_int(r, R)` |
| `V(R,r)` | `.surface(r, R)` |
| `H_el_0`, `H_ion_0` | assembled via `qscat.dvr`/`qscat.core` kinetic + potential terms, no single-named equivalent |
| `psi_sc`, driven equation `(E-H)\|psi_sc> = V_int\|psi_in>` | `qscat.core.driven.ve_cross_section`'s scattered-wave solve |
| `F_l`, `G_l`, `u_l^+` | `qscat.special.coulomb.coulomb_f_en`, `.g_en`, `.h1_en` |
| `eta = -Z/k` | Sommerfeld parameter, same role as `qscat.special.coulomb`'s `eta` |
| `rho_n(r)`, `E_n` (Rydberg bound states) | `qscat.core.dissociation.anion_electronic_states` |
| `V_DR(R,r)` | `v_dr_diag` in `qscat.core.dissociation`/`.dr_cross_section` |
| `T_VE`, `T_DR` | T-matrix elements in `qscat.core.driven`/`dissociation.dr_cross_section` |
| `V_n(R) = V0(R) + E_el_n(R)` (Rydberg potential curves) | discussed in `docs/physics/h2plus-dr.md`, no single named function |

## Not used here

- Table 2.1's FEM-DVR-ECS grid parameters are read only as an independent
  convergence cross-check; `qscat`'s H₂⁺ grids (the eMoScat-derived deck in
  `validation/h2plus/config.py`, and `qscat.tuning`'s a-priori proposal) are
  built and validated on their own terms, not ported from this table.
- Fig. 2.1-2.12 (cross-section curves, potential-curve/vibrational-level
  overlays) are read for their qualitative conclusions only; no pixel/curve
  data was digitized or extracted from them.
- The two-channel (direct+indirect DR mechanism) extension and the local
  complex potential / frame-transformation comparisons the thesis proposes
  as future work (p. 29) are out of scope for this note — they are covered,
  where implemented, by `qscat.core.lcp` and by the separate
  Hvizdoš et al. 2018 reference note, not by this thesis itself.
- Front matter (declaration, abstract, acknowledgments), List of Figures/
  Tables/Abbreviations, and the Bibliography carry no technical content
  transcribed here beyond the citations already reflected in
  `reference/literature/README.md`.
