# M. Váňa, doctoral thesis, Charles University, Prague (2017) — A model of resonant collisions of electrons with molecules and molecular ions

**Source:** `reference/literature/vana-2017-thesis.pdf` (gitignored) ·
<https://dspace.cuni.cz/handle/20.500.11956/92902>
**Pagination:** Charles University thesis, unnumbered front matter. Printed page =
extractor page − 8; verified against the printed page numbers embedded in the
extracted text at extractor pages 12–18 (→ printed 4–10, matching the table of
contents' own page numbers) and spot-checked again at extractor pages 90 and 98
(→ printed 82, 90 — the final page). The offset is constant across the whole
document; all locators below use the printed number.

## Why this repository cares

This is the repo owner's own thesis (supervisor K. Houfek) and the single
broadest source behind qModeling: it is the only place all four models in
`qscat.model` (N₂-like, NO-like, F₂-like, H₂⁺-like) are parametrized together
(Tables 1.1–1.2), the only source for the FEM-DVR-ECS grids actually used in the
author's own calculations (Tables 2.1–2.2), and the source of the discrete-state
projection (§1.6) that is a named future direction for `qscat.core`. Chapter 3's
electron-molecule content was subsequently published as Váňa & Houfek, Phys. Rev.
A 95, 022714 (2017) (`vana-houfek-2017-pra95-022714.md`) — cited there in
preference to this thesis for anything that paper also covers; this note points
to that overlap rather than duplicating it. The H₂⁺ dissociative-recombination
model (ch. 4) has no separate journal-article note in this collection; Hvizdoš
(2016) and Hvizdoš et al. (2018) (see `README.md`) cover its time-independent
side.

## What this repository uses

| Fact | Locator | Used by |
|---|---|---|
| 2-D Hamiltonian, potential form, Morse `V0(R)`, interaction `Vint(R,r)` | p. 8, Eq. (1.11)–(1.15) | `qscat.model.DiatomicResonanceModel` (same equations as Houfek et al. 2006, renumbered) |
| `α(R) = αc`, sigmoid `λ(R)` | p. 9, Eq. (1.16)–(1.18) | `DiatomicResonanceModel.alpha_c`, `.lam` |
| **Table 1.1** — N₂-like/NO-like/F₂-like constants | p. 9 | `qscat.model.library.N2`/`NO`/`F2` |
| **Table 1.2** — H₂⁺-like model parameters (incl. µ = 918.076) | p. 18 | `qscat.model.library.H2P` |
| H₂⁺ Hamiltonian split, Coulomb `Hel0(r)`, interaction `Vint(R,r)` | p. 18, Eq. (1.72)–(1.75) | `qscat.model.IonicResonanceModel` |
| Cross-section normalization `σ = 4π³/k² |T|²` | p. 15, Eq. (1.60) | `qscat.core.driven.ve_cross_section`, `qscat.core.dissociation.da_cross_section`/`dr_cross_section` |
| T-matrix from S-matrix, `T = (S−δ)/2πi` | p. 15, Eq. (1.61) | `qscat.core` T-matrix bookkeeping generally |
| Fixed-R electronic Hamiltonian, LCP curve `Vres(R) = Eres(R) − iΓ(R)/2` | p. 16, Eq. (1.62)–(1.63) | `qscat.ecs.find_resonance_pole`, `qscat.core.lcp.local_complex_potential` |
| LCP initial wave packet `√(Γ/2π) χvi(R)`, `HLCP(R)` | p. 16, Eq. (1.64)–(1.65) | conceptual basis for `qscat.core.lcp`'s TI resolvent route (a different but equivalent formulation — see Not used here) |
| Discrete-state suppression function `f(r)` | p. 17, Eq. (1.69) | named future direction — not yet implemented (see Not used here) |
| Projection onto the discrete state, `Ψd(R,t) = ∫dr Ψ(R,r,t) φd*(R,r)` | p. 17, Eq. (1.70) | named future direction — not yet implemented |
| **Table 2.1** — N₂-like/NO-like FEM-DVR-ECS grids | p. 28 | reference grid, differs from the eMoScat deck (see Findings) |
| **Table 2.2** — F₂-like/H₂⁺-like FEM-DVR-ECS grids | p. 29 | reference grid, differs from the eMoScat deck (see Findings) |
| §3.4 boomerang / quasi-bound-state interpretation | p. 52–56 | `docs/physics/n2-resonance.md`'s narrative; conceptual basis for `qscat.ecs`'s resonance-eigenstate work |

## Equations

```
H = -1/(2mu) d^2/dR^2 - 1/2 d^2/dr^2 + V(R,r)                    p. 8, Eq. (1.11)
V(R,r) = V0(R) + l(l+1)/(2r^2) + Vint(R,r)                       p. 8, Eq. (1.12)
V0(R) = lim_{r->inf} V(R,r)                                      p. 8, Eq. (1.13)
Vint(R,r) = lambda(R) exp(-alpha(R) r^2)                         p. 8, Eq. (1.14)
V0(R) = D0 ( exp(-2 alpha0(R-R0)) - 2 exp(-alpha0(R-R0)) )       p. 8, Eq. (1.15)
alpha(R) = alpha_c                                               p. 9, Eq. (1.16)
lambda(R) = lambda_inf + lambda0 / (1 + exp(lambda1(R-Rlambda))) p. 9, Eq. (1.17)
lambda0 = (lambda_c - lambda_inf)(1 + exp(lambda1(Rc-Rlambda)))  p. 9, Eq. (1.18)

sigma_{vi,beta}(E) = (4 pi^3 / k_i^2) |T_{vi->beta}(E)|^2        p. 15, Eq. (1.60)
T_{vi->beta}(E) = (S_{vi->beta}(E) - delta_{vi,beta)) / (2 pi i) p. 15, Eq. (1.61)

Hel(r;R) = -1/2 d^2/dr^2 + V(R,r)                                p. 16, Eq. (1.62)
Vres(R) = Eres(R) - (i/2) Gamma(R)                               p. 16, Eq. (1.63)
Psi_LCP(R, t=0) = sqrt(Gamma(R)/2pi) chi_vi(R)                   p. 16, Eq. (1.64)
H_LCP(R) = -1/(2mu) d^2/dR^2 + Eres(R) - (i/2) Gamma(R)          p. 16, Eq. (1.65)

f(r) = 1 - 1/(1 - exp(-(r-rd)))            [discrete-state suppression]  p. 17, Eq. (1.69)
Psi_d(R,t) = Int_0^inf dr Psi(R,r,t) phi_d*(R,r)                 p. 17, Eq. (1.70)

H(R,r) = Hion0(R) + Hel0(r) + Vint(R,r)                          p. 17, Eq. (1.71)
Hion0(R) = -1/(2mu) d^2/dR^2 + V0(R)      [V0 = Morse, Eq. 1.15] p. 18, Eq. (1.72)
Hel0(r) = -1/2 d^2/dr^2 + l(l+1)/(2r^2) - 1/r  [Coulomb]         p. 18, Eq. (1.73)
Vint(R,r) = -a1 (1 - tanh(a2 - R - a3 R^4)/7) (tanh(R/a4))^4 exp(-r^2/3)/r   p. 18, Eq. (1.74)
Vinf_int(r) = lim_{R->inf} Vint(R,r) = -2 a1 exp(-r^2/3)/r       p. 18, Eq. (1.75)
```

Atomic units throughout — matches `libs/qscat/qscat/units.py`, no conversion
needed. Eq. (1.11)–(1.18) are the same equations as Houfek et al. (2006) Eq.
(3)–(6), (45)–(48) with the thesis's own numbering (see
`houfek-2006-pra73-032721.md` for that source, transcribed independently there).

## Parameters and numeric values

**Table 1.1**, p. 9 (atomic units):

| Parameter | N₂-like | NO-like | F₂-like |
|---|---|---|---|
| µ | 12 766.36 | 13 614.16 | 17 315.99 |
| l | 2 (d-wave) | 1 (p-wave) | 1 (p-wave) |
| D0 | 0.75102 | 0.2363 | 0.0598 |
| α0 | 1.15350 | 1.5710 | 1.5161 |
| R0 | 2.01943 | 2.1570 | 2.6906 |
| λ∞ | 6.21066 | 6.3670 | 18.8490 |
| λ1 | 1.05708 | 5.0000 | 3.2130 |
| Rλ | −27.9833 | 2.0843 | 1.8320 |
| λc | 5.38022 | 6.0500 | 18.1450 |
| Rc | 2.40500 | 2.2850 | 2.5950 |
| αc | 0.40000 | 1.0000 | 3.0000 |

Checked against the repo:

```
grep -n "mu=\|ell=\|D0=\|alpha0=\|R0=\|lambda_inf=\|lambda_1=\|R_lambda=\|lambda_c=\|R_c=\|alpha_c=" \
  libs/qscat/qscat/model/library.py
```

`qscat.model.library.N2` and `.NO` match Table 1.1 exactly (verified
2026-08-16 — same values already checked against Houfek et al. (2006)'s
identical Table I). `qscat.model.library.F2` (`mu=17315.99, ell=1,
D0=0.05980, alpha0=1.51610, R0=2.69060, lambda_inf=18.8490,
lambda_1=3.21300, R_lambda=1.8320, lambda_c=18.1450, R_c=2.5950,
alpha_c=3.0`) also matches Table 1.1 exactly (verified 2026-08-16). F₂ is
**not** in Houfek et al. (2006), but it *is* in Houfek, Rescigno & McCurdy,
Phys. Rev. A **77**, 012710 (2008), Table I, p. 012710-6 — nine years before
this thesis, with identical values. **That 2008 paper, not this thesis, is
F₂'s primary published source**; see
`reference/literature/houfek-2008-pra77-012710.md`.

**Table 1.2**, p. 18 (atomic units):

| Parameter | Value |
|---|---|
| µ | 918.076 |
| l | 1 (p-wave) |
| a0 | 0.1027 |
| α | 0.69 |
| R0 | 2.0 |
| a1 | 1.6435 |
| a2 | 6.2 |
| a3 | 0.0125 |
| a4 | 1.15 |

Note the thesis's own table header uses `a0`/`α` for the Morse well-depth/width
(elsewhere `D0`/`α0` in Eq. 1.15) — a naming collision with the unrelated
per-electron-molecule-model `α(R) = αc` of Eq. (1.16); not a repo issue, just a
thesis notation quirk worth flagging so a reader isn't confused chasing `α`
across sections.

Checked against the repo:

```
grep -n "mu=918.076\|R0=2.0\|alpha=0.69\|a1=1.6435\|a2=6.2\|a3=0.0125\|a4=1.15" \
  libs/qscat/qscat/model/library.py
```

`qscat.model.library.H2P` (`mu=918.076, ell=1, charge=-1, V0=0.1027, R0=2.0,
alpha=0.69, a1=1.6435, a2=6.2, a3=0.0125, a4=1.15`) matches Table 1.2 exactly,
field for field, including `µ = 918.076` (verified 2026-08-16). This is a
third independent published source (alongside Hvizdoš (2016) Table 1.1 and
Hvizdoš et al. (2018) §II A) agreeing on `µ = 918.076`, against eMoScat's own
JSON deck value of `918.25` — see `libs/qscat/qscat/model/library.py:63-66`'s
own comment on this discrepancy, already recorded in the Hvizdoš et al.
(2018) note.

## Findings and limits

- **Tables 2.1/2.2 (FEM-DVR-ECS grids, p. 28–29) differ from the eMoScat
  decks** locked by `validation/diatomic/test_da_grid.py`. Example: F₂-like
  nuclear coordinate here uses `nq = 20`, θ = 35°, `nc = 15` (p. 29); the
  eMoScat deck (`reference/eMoScat/input/F2/`) uses a coarser quadrature order
  and θ = 25°. Both are valid discretizations of the same model and
  quantities computed on them converge to the same answer (per
  `qscat.tuning`'s own reproduce-and-beat checks) — the two are a
  convergence cross-check on each other, not one superseding the other. Not
  independently re-verified numerically as part of writing this note beyond
  the locators above; see the `discretisation-tuner` skill's own gates for
  live convergence evidence.
- **§3.4** (p. 52–56): the boomerang mechanism — interference between the
  prompt autodetachment reflection and the time-delayed reflection after a
  full vibrational period — explains the *regular* oscillatory structure in
  the VE cross sections, but not the *asymmetric* peaks; those require at
  least one further vibrational-period contribution, and in general
  longer-lived quasi-bound states in the well outside the autodetachment
  region (p. 52–53). For the N₂-like model the nuclear-motion period is
  `T ≈ 655 a.u. ≈ 16 fs`, read off the spacing between successive maxima of
  the mean internuclear distance `⟨R⟩(t)` (Fig. 3.22, p. 54). The cross
  section decomposes cleanly into a first-reflection contribution (2-D model
  at t=650, LCP at t=320) and a second, boomerang contribution (2-D at
  t=1270, LCP at t=940) that modulates it with a series of symmetric peaks
  (Fig. 3.23, p. 55). The propagation is cut off once the wave-packet
  normalization has dropped ~8 orders of magnitude (`tc = 4000` for N₂-like,
  p. 56), beyond which further contributions are negligible.
- The peer-reviewed Váňa & Houfek (2017) §VIII gives the same
  quasibound-state interpretation in terms of eigenstates of `Vres(R)`
  directly (elastic boomerang maxima ≈ those eigenvalues, VE 0→1 maxima
  displaced) and adds NO-like lifetime estimates not repeated in this
  section of the thesis in the same form — see
  `vana-houfek-2017-pra95-022714.md` rather than duplicating here.
- The thesis states its own scope limit for §1.6/Eq. (1.69)-(1.70): the
  discrete-state choice `φd` used is "the physical choice... described in
  Houfek et al. [2008a]" (p. 17) — one of several possible choices discussed
  in that 2008 paper, not derived independently here.

## Terminology map

| Thesis symbol | qModeling name |
|---|---|
| `mu`, `l`, `D0`, `alpha0`, `R0`, `lambda_inf`, `lambda1`, `Rlambda`, `lambda_c`, `Rc`, `alpha_c` | same names in `qscat.model.library.N2`/`NO`/`F2` (Table 1.1 fields) |
| H₂⁺ table's `a0` (well depth) | `qscat.model.library.H2P`'s `V0` field |
| H₂⁺ table's `α` (Morse width) | `qscat.model.library.H2P`'s `alpha` field |
| `Eres(R)` (pole of `Hel(r;R)`, already includes `V0(R)`) | `V_d(R)` in `qscat.core.lcp` (also already includes `v0(R)` — see that module's docstring) |
| `Gamma(R)` | `Gamma(R)` (same name), `qscat.ecs.find_resonance_pole` |
| `Vres(R)`, `H_LCP(R)` | `qscat.core.lcp.local_complex_potential`'s complex `E_pole(R) = V_d - i*Gamma/2` (equivalent quantity, different code route — see Not used here) |
| `phi_d(R,r)`, `Psi_d(R,t)` (discrete-state projection) | not implemented; see Not used here |
| `sigma_{vi,beta} = 4pi^3/k^2 |T|^2` | `qscat.core.driven`/`dissociation`'s cross-section normalization |

## Not used here

- §1.6's discrete-state projection `Ψd(R,t)` (Eq. 1.69–1.70) is **not**
  implemented in `qscat.core` — it remains a named future direction (see the
  `td-alternative-extractors` memory's SP2/SP3 follow-ons), distinct from the
  three extractors (`tw`/`delta`/`flow`) that are implemented.
- The literal time-dependent LCP propagation of Eq. (1.64)–(1.65)
  (`Ψ_LCP(R,0) = √(Γ/2π) χvi`, evolved under `H_LCP`) is not the route
  `qscat.core.lcp` takes: `lcp_da_cross_section` solves the equivalent
  time-*independent* resolvent problem on `Vres(R)` directly (boundary
  wavefunction-value flux), not a propagated wave packet — a different but
  physically equivalent formulation of the same LCP approximation.
- Chapters 1–2's general multichannel-scattering formalism (§1.1) and the
  FEM-DVR/Padé numerical-method exposition (§1.3, §2.1–2.4) are read for
  context only; the corresponding qscat implementations (`qscat.dvr`,
  `qscat.evolution.make_pade_stepper`) are validated and cited against their
  own dedicated sources (see `docs/physics/femdvr-ecs.md`,
  `docs/physics/n2-2d-td-cross-section.md`), not against this thesis.
- Chapter 3's N₂-like/NO-like/F₂-like cross-section results and figures
  (§3.1–3.3, most of §3.4's NO-like/F₂-like subsections) substantially
  overlap Váňa & Houfek (2017)'s published content — see that paper's note
  rather than this one for anything beyond the N₂-like boomerang timing
  extracted above.
- Chapter 4's dissociative-recombination results and Attachments (animation
  descriptions, p. 89–90) are read for context only; no numeric values are
  taken from them.
