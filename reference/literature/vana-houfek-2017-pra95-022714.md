# M. Váňa, K. Houfek, Phys. Rev. A 95, 022714 (2017) — Time-dependent formulation of the two-dimensional model of resonant electron collisions with diatomic molecules and interpretation of the vibrational excitation cross sections

**Source:** `reference/literature/vana-houfek-2017-pra95-022714.pdf` (gitignored) · DOI 10.1103/PhysRevA.95.022714
**Pagination:** per-article, 022714-1 … 022714-17 (extractor page N = printed page 022714-N;
offset zero, checked at page 1 and page 17).

## Why this repository cares

This is **the peer-reviewed publication of this repository's own antecedent work** —
the author of qModeling is a co-author (M. Váňa) of this paper, and it is the published
counterpart of `reference/literature/vana-2017-thesis.pdf`. It takes the exact 2-D
(one nuclear + one electronic coordinate) model introduced by
`houfek-2006-pra73-032721.md` and gives it a **time-dependent** formulation: prepare a
wavepacket, propagate it under the FEM-DVR-ECS Hamiltonian with a generalized
Crank-Nicolson (diagonal Padé) evolution operator, and extract the S-matrix from the
propagated wave function by one of three methods. This is the direct published
ancestor of `qscat.evolution` (the Padé/CN steppers), `qscat.core.time_dependent` /
`qscat.core.td_extractors` (the Tannor-Weeks/delta/flux extractors), and
`qscat.core.lcp` (the LCP nuclear eigenvalue problem and its quasibound-state
interpretation `docs/physics/lcp-resonance-levels.md` still owes a milestone-3
external comparison to). Without this source, the order-3 Padé stepper's `N=3, dt=1.0`
default, the three TD extractors' shared propagate-once design, and the
quasibound/boomerang narrative used to sanity-check `resonance_levels` output would be
unfounded choices rather than a reproduction of a published method.

## What this repository uses

| Fact | Locator | Used by |
|---|---|---|
| 2-D Hamiltonian `H`, potential `V(R,r)`, Morse/sigmoid parametrization | p. 022714-2, Eq. (3)-(10) | same as `houfek-2006-pra73-032721.md`; here confirms the identical model is reused unchanged |
| Table I now also tabulates **F₂-like** parameters (absent from the 2006 paper) | p. 022714-2, Table I | `qscat.model.library.F2` (see Parameters below — new parity check) |
| Channel Hamiltonians `H_0^VE`, `H_0^DA` and asymptotic eigenstates | p. 022714-2 – 022714-3, Eq. (11)-(19) | `qscat.core.channels`, `qscat.core.driven` (conceptual — same asymptotics as the TI formulation) |
| Initial Gaussian wavepacket `Ψ_in_vi(R,r)` | p. 022714-3, Eq. (20) | `projects/n2_2d_td_cross_section/wavepacket.py`, `qscat.core.time_dependent`'s incident-wavepacket construction |
| Unitary evolution operator `ψ(t+t0) = e^{-iHt}ψ(t0)` | p. 022714-3, Eq. (23) | the operator `qscat.evolution` approximates |
| **Method 1 — Tannor-Weeks correlation function + S-matrix** | p. 022714-3, Eq. (24), (26)-(29) | `qscat.core.td_extractors.TannorWeeks`, `time_dependent.sigma_from_correlations` |
| **Method 2 — δ-function (line-projection) variant** | p. 022714-4, Eq. (32)-(35) | `qscat.core.td_extractors.Dirac` (eMoScat's `DiracTestFunction2d`) |
| **Method 3 — probability-flux (Wronskian) variant** | p. 022714-4, Eq. (36)-(37) | `qscat.core.td_extractors.Flux` (eMoScat's `FluxTestFunction2d`) |
| Cross section `σ = 4π³/k_i² |T|²`, `T = (S-δ)/2πi` | p. 022714-4, Eq. (38)-(39) | **the literal source** of `sigma_{v->v'}(E) = π|S-S_ref|^2/(2E)` used throughout `qscat.core.time_dependent`/`td_extractors` (unlike PRA73's Eq. 25-26, which the companion note flags as a *different*-looking normalization — see Equations) |
| LCP: `H_el(r;R)`, `V_res(R) = E_res(R) - (i/2)Γ(R)`, and the statement that `Im V_res` is nonzero only where `V0(R) < E_res(R)` | p. 022714-4, Eq. (40)-(41) + following sentence | `qscat.core.lcp.local_complex_potential`, `lcp_resonance_levels`; `docs/physics/lcp-resonance-levels.md`'s Γ-support condition and `test_lcp_resonance_levels.py`'s docstring citation |
| LCP nuclear Hamiltonian `H_LCP(R)` and initial `Ψ_LCP(R,0) = sqrt(Γ(R)/2π) χ_vi(R)` | p. 022714-4 – 022714-5, Eq. (42)-(43) | `qscat.core.lcp.lcp_resonance_levels` (`H_N = T(mu) + diag(Vd - i*Gamma/2)`) |
| ECS coordinate maps for `r` and `R` | p. 022714-5, Eq. (49)-(50) | `qscat.ecs.ecs_map` (same map as PRA73 Eq. 27-28, renumbered here) |
| **Generalized Crank-Nicolson (diagonal Padé) evolution operator** | p. 022714-5, Eq. (51) | `qscat.evolution.make_pade_stepper` (`make_cn_stepper`/`make_sparse_cn_stepper` are its order-1 special case) |
| Padé order `N=3`, time step `dt=1.0` used for **all three** models | p. 022714-6, p. 022714-16 (Appendix) | `docs/physics/n2-2d-td-cross-section.md`'s "eMoScat's setting" claim for `make_pade_stepper(H, dt=1.0, order=3)` |
| Quasibound-state interpretation of cross-section peaks as `V_res(R)` eigenstates | p. 022714-15 | `docs/physics/lcp-resonance-levels.md`'s milestone-3 comparison target |
| NO-like lifetime estimates: first VE 0→1 peak forms at `t > 10 000`; lowest (elastic) quasibound state has lifetime `> 30 000` a.u. | p. 022714-15 | `docs/physics/lcp-resonance-levels.md` milestone-3 "NO lifetime bounds from the reported qualitative widths" |
| NO-like `V_res(R)` has its minimum *behind* the crossing point (an outer well outside the autodetachment region), explaining the long-lived states | p. 022714-8 | narrative context for the NO lifetime numbers above; not itself a numeric input to any module |

## Equations

```
H = -(1/2mu) d^2/dR^2 - (1/2) d^2/dr^2 + V(R,r)                  p. 022714-2, Eq. (3)
V(R,r) = V0(R) + l(l+1)/(2r^2) + Vint(R,r)                       p. 022714-2, Eq. (4)
V0(R) = lim_{r->inf} V(R,r)                                      p. 022714-2, Eq. (5)
Vint(R,r) = lambda(R) exp(-alpha(R) r^2)   [see note below]      p. 022714-2, Eq. (6)
H_0^VE = -(1/2) d^2/dr^2 + l(l+1)/(2r^2) - (1/2mu) d^2/dR^2 + V0(R)   p. 022714-2, Eq. (11)
H_0^DA = -(1/2) d^2/dr^2 - (1/2mu) d^2/dR^2 + Vb(r)              p. 022714-2, Eq. (12)
Vb(r) = lim_{R->inf} V(R,r)                                      p. 022714-2, Eq. (13)
[-(1/2mu) d^2/dR^2 + V0(R)] chi_vi(R) = E_vi chi_vi(R)           p. 022714-3, Eq. (17)
Psi_in_vi(R,r) = (pi sigma^2)^{-1/4} chi_vi(R) exp(-(r-r0)^2/2sigma^2 - i p0 r)  p. 022714-3, Eq. (20)
psi(t+t0) = U(t) psi(t0) = exp(-iHt) psi(t0)                     p. 022714-3, Eq. (23)
Phi_out_vf(R,r) = (pi sigma^2)^{-1/4} chi_vf(R) exp(-(r-r0)^2/2sigma^2 + i q0 r)  p. 022714-3, Eq. (24)
Phi_out_DA(R,r) = (pi sigma^2)^{-1/4} phi_b(r) exp(-(R-R0)^2/2sigma^2 + i Q0 R)   p. 022714-3, Eq. (25)
C_beta(t) = Int Int [Phi_out_beta(R,r)]^* psi(R,r,t) dr dR      p. 022714-3, Eq. (26)
S^{T&W}_{vi->beta}(E) = (2pi)^{-1} [eta_out_beta(E)]^* eta_in_vi(E) Int C_beta(t) e^{iEt} dt   p. 022714-3, Eq. (27)
E = k_i^2/2 + E_vi = k_f^2/2 + E_vf = K^2/2mu + E_b              p. 022714-3, Eq. (28)
eta_out_beta(E) = Int Int [phi_tilde_out_beta(E,R,r)]^* Phi_out_beta(R,r) dr dR   p. 022714-3, Eq. (29)
phi_tilde_out_vf(E,R,r) = chi_vf(R) sqrt(kf/2pi) r h_l^+(kf r)  p. 022714-4, Eq. (30)
phi_tilde_out_DA(E,R,r) = sqrt(mu/2piK) e^{iKR} phi_b(r)         p. 022714-4, Eq. (31)
S^delta_{vi->vf}(E) = (1/2pi)[eta'_out_vf(E)]^* eta_in_vi(E) Int Int e^{iEt} chi_vf(R)^* psi(R,r0,t) dt dR   p. 022714-4, Eq. (34)
S^delta_{vi->DA}(E)  = (1/2pi)[eta'_out_DA(E)]^* eta_in_vi(E) Int Int e^{iEt} phi_b(r)^* psi(R0,r,t) dt dr   p. 022714-4, Eq. (35)
S^F_{vi->vf}(E) = [2 eta_in_vi(E)]^{-1} (1/2i) Int Int e^{iEt} [(phi_tilde_out_vf)^* d(psi)/dr - psi^* d(phi_tilde_out_vf)/dr]_{r=r0} dt dR   p. 022714-4, Eq. (36)
S^F_{vi->DA}(E) = [2 eta_in_vi(E)]^{-1} (1/2i mu) Int Int e^{iEt} [(phi_tilde_out_DA)^* d(psi)/dR - psi^* d(phi_tilde_out_DA)/dR]_{R=R0} dt dr   p. 022714-4, Eq. (37)
sigma_{vi->beta}(E) = (4 pi^3/k_i^2) |T_{vi->beta}(E)|^2         p. 022714-4, Eq. (38)
T_{vi->beta}(E) = [S_{vi->beta}(E) - delta_{vi,beta}] / (2 pi i)  p. 022714-4, Eq. (39)
H_el(r;R) = -(1/2) d^2/dr^2 + V(R,r)                             p. 022714-4, Eq. (40)
V_res(R) = E_res(R) - (i/2) Gamma(R)                             p. 022714-4, Eq. (41)
   [Im part nonzero only where V0(R) < E_res(R) -- stated in prose immediately after Eq. (41)]
Psi_LCP(R,t=0) = sqrt(Gamma(R)/2pi) chi_vi(R)                    p. 022714-4, Eq. (42)
H_LCP(R) = -(1/2mu) d^2/dR^2 + E_res(R) - (i/2) Gamma(R)         p. 022714-4, Eq. (43)
sigma^LCP_{vi->beta}(E) = (4 pi^3/k_i^2) |T^LCP_{vi->beta}(E)|^2  p. 022714-5, Eq. (46)
r'(r) = r for r<r0; r0+(r-r0)e^{i theta_r} for r>=r0             p. 022714-5, Eq. (49)
R'(R) = R for R<R0; R0+(R-R0)e^{i theta_R} for R>=R0             p. 022714-5, Eq. (50)
psi(t+dt) = e^{-iH dt} psi(t) ~= prod_{j=1}^{N} (1+c_j H dt)/(1-c_j^* H dt) psi(t)   p. 022714-5, Eq. (51)
```

Atomic units throughout (energies in hartree, lengths in bohr, `hbar = m_e = 1`;
1 a.u. of time `= 2.418884e-17` s) — the same convention `libs/qscat/qscat/units.py`
uses; no conversion needed (unit definitions stated p. 022714-2).

**Eq. (6) sign note.** As printed (and cross-checked against the extractor's raw text
twice), this paper's Eq. (6) reads `Vint(R,r) = λ(R)e^{-α(R)r²}` — **no leading minus
sign**. `vana-2017-thesis.pdf`'s equivalent equation (its Eq. 1.12/1.13-adjacent form,
line-checked in the extracted thesis text) agrees: also no leading minus. By contrast,
`houfek-2006-pra73-032721.pdf`'s own Eq. (6), re-extracted and checked directly against
its PDF page image for this note, unambiguously carries the minus:
`Vint(R,r) = −λ(R)e^{−α(R)r²}` (matching `houfek-2006-pra73-032721.md`'s existing
transcription). `qscat.model.diatomic.DiatomicResonanceModel.v_int` implements the
**minus-sign** form (`libs/qscat/qscat/model/diatomic.py:82`,
`V_int(r, R) = -lambda(R) exp(-alpha_c r^2)`), i.e. it follows PRA73/the physical
attractive-interaction convention, not the literal text of this 2017 paper (or its
companion thesis). This reads as a dropped minus sign carried from the thesis into this
paper's typesetting, not a physical sign change — the models, Table I values, and all
downstream Vres/LCP results are otherwise identical to PRA73 — but it is recorded here
as a checked, real textual discrepancy rather than silently normalized away.

The physics settles it independently of the typesetting: `λ(R)` is **positive**
everywhere in these models (`λ∞ = 6.21066` for N₂, `18.8490` for F₂ — Table I,
p. 022714-2), so `Vint = +λ e^{−αr²}` would be a purely **repulsive** interaction,
which cannot bind the electron into the resonance whose cross sections this paper
computes. The minus is required, and the repo's implementation is correct.

## Parameters and numeric values

Table I, p. 022714-2 (atomic units) — identical N₂/NO columns to
`houfek-2006-pra73-032721.md`'s Table I (already verified against
`qscat.model.library.N2`/`NO`), **plus a new F₂ column absent from the 2006 paper**:

| Parameter | F₂ |
|---|---|
| μ | 17 315.99 |
| l | 1 (p wave) |
| D0 | 0.0598 |
| α0 | 1.5161 |
| R0 | 2.6906 |
| λ∞ | 18.8490 |
| λ1 | 3.2130 |
| R_λ | 1.8320 |
| λc | 18.1450 |
| Rc | 2.5950 |
| αc | 3.0000 |

Checked against the repo:

```
grep -n "17315.99\|0.05980\|1.51610\|2.69060\|18.8490\|3.21300\|1.8320\|18.1450\|2.5950" \
  libs/qscat/qscat/model/library.py
```

`qscat.model.library.F2` (`mu=17315.99, ell=1, D0=0.05980, alpha0=1.51610,
R0=2.69060, lambda_inf=18.8490, lambda_1=3.21300, R_lambda=1.8320,
lambda_c=18.1450, R_c=2.5950, alpha_c=3.0`) matches this table exactly
(verified 2026-08-16). This is a stronger provenance for F₂'s electronic-surface
constants than previously recorded: `qscat-core-model-architecture`/CLAUDE.md
attribute `F2`'s parameters to the eMoScat deck (`reference/eMoScat/input/F2/`)
because `houfek-2006-pra73-032721.pdf` does not tabulate F₂ — but this 2017 paper
does, independently confirming the same numbers. (F₂'s nuclear *discretisation*
deck, `MoleculeConfig.da_grid()`, is still sourced from eMoScat/the thesis, not from
this paper's Table II grid parametrization — see Not used here.)

The Appendix (p. 022714-16) states the diagonal Padé evolution operator, Eq. (51), was
used at **order `N=3`, time step `dt = 1.0` a.u. ≈ 0.024 fs, for all three models**
(N₂-like, NO-like, F₂-like) — matches CLAUDE.md's `qscat.evolution` description of
`make_pade_stepper(H, dt=1.0, order=3)` as "eMoScat's setting" (verified: p. 022714-16
text, not independently re-derived).

## Findings and limits

- **Three S-matrix extraction methods agree at the converged cutoff time but differ
  markedly at early times**: the original Tannor-Weeks method (Eq. 27) gives the
  smoothest early-time energy transform; the δ-function (Eq. 34-35) and flux (Eq. 36-37)
  variants oscillate at early times before converging — p. 022714-13. This motivates why
  `qscat.core.td_extractors` implements all three as a shared propagate-once family
  rather than picking one.
- **Boomerang oscillations explain only regular structure, not asymmetric peaks.**
  Asymmetric cross-section peaks require considering interference of *more than two*
  time-separated contributions (initial reflection + repeated vibrational periods, in
  some cases plus a long-lived quasibound state) — the paper's own conclusion states
  "boomerang motion"/"boomerang oscillations" are not quite accurate terms for what it
  observes — p. 022714-16 (Conclusion).
- **NO-like model: long-lived quasibound states arise because `V_res(R)`'s minimum sits
  behind (outward of) the crossing point**, forming an outer potential well outside the
  autodetachment region — p. 022714-8. This directly produces the long formation times
  (first VE 0→1 peak at `t > 10 000`; lowest elastic quasibound state, lifetime
  `> 30 000` a.u.) reported at p. 022714-15.
- **N₂-like model**: nuclear motion period `T ≈ 655` a.u. `≈ 16` fs (from the spacing
  between mean-internuclear-distance maxima); normalization cutoff `t_c = 4000` a.u.
  (8-orders-of-magnitude decay) used to truncate the propagation for the reported cross
  sections — p. 022714-13.
- **Quasibound-state correspondence with peaks**: each narrow cross-section peak
  corresponds to an eigenstate of `V_res(R)` (states of the complex nuclear Hamiltonian,
  Eq. 43); elastic-channel boomerang maxima line up roughly with these quasibound-state
  energies, while VE 0→1 maxima are systematically displaced from them, producing highly
  asymmetric structure — p. 022714-15.
- **All three S-matrix methods (Tannor-Weeks, δ, flux) reproduce the time-independent
  cross sections "in perfect agreement" for the energies of interest** at the converged
  cutoff time, across all three models (N₂/NO/F₂-like) — p. 022714-16 (Conclusion). No
  numeric agreement percentage is stated in the text; this is a qualitative claim, unlike
  this repo's own quantitative TD-vs-TI gates (~1-2% median, see
  `docs/physics/n2-2d-td-cross-section.md`).
- The paper's own scope statement: for real systems, similar time-dependent
  calculations could be performed within the LCP approximation or the nonlocal
  resonance model "and thus interpret the results in the same way" — offered as a
  future direction, not executed here — p. 022714-16.

## Terminology map

| Paper symbol | qModeling name |
|---|---|
| `psi(t)` (propagated 2-D wave function) | `Psi(t)` in `qscat.core.time_dependent`/`td_extractors` |
| `C_beta(t)` (Eq. 26) | `TannorWeeks`'s recorded `c_{v'}(t_n)` series |
| `S^{T&W}` (Eq. 27) | `sigma_from_correlations`'s Tannor-Weeks transform |
| `S^delta` (Eq. 34-35) | `qscat.core.td_extractors.Dirac` |
| `S^F` (Eq. 36-37) | `qscat.core.td_extractors.Flux` |
| `eta_in`, `eta_out` | `sigma_from_correlations`'s deconvolution factors (same names) |
| diagonal Padé approximant, Eq. (51) | `qscat.evolution.make_pade_stepper` (`pade_roots` for `c_j`) |
| `V_res(R)`, Eq. (41) | `qscat.core.lcp`'s `Vd(R) - i*Gamma(R)/2`; also `docs/physics/n2-resonance.md`'s `V_d(R)` |
| `Psi_LCP(R,t)` (Eq. 42-43) | `qscat.core.lcp.lcp_resonance_levels`'s complex-symmetric nuclear eigenproblem (eigenstates, not a propagated wavepacket, in the current implementation — see Not used here) |
| `phi_d(R,r)` discrete state (Sec. V) | not named in this repo (see Not used here) |
| ECS `r0`, `R0`, `theta_r`, `theta_R` | `qscat.ecs.ecs_map`'s `R0`, `theta` |

## Not used here

- **Section V, the discrete-state projection `Psi_d(R,t)` (Eq. 47-48) and the
  nonlocal-theory framing it is drawn from.** This repo's `qscat.core.lcp` solves the
  LCP nuclear eigenvalue problem directly (Eq. 41/43) rather than propagating a 2-D wave
  packet and projecting it onto a fixed-`R` diabatic electronic state; no
  `phi_d(R,r)`-style projection or suppressing function (Eq. 47) is implemented.
- **Table II's FEM-DVR-ECS grid parametrization** (quadrature order `nq`, ECS angles
  `theta`, per-element counts for all three models, p. 022714-16) — this repo's N₂/NO/F₂
  discretisations are sourced from the eMoScat decks and `qscat.tuning`'s own
  discretisation tuner, not transcribed from this table.
- **Figs. 1-25** (potential surfaces, wave-packet/normalization evolution plots,
  cross-section-vs-cutoff-time comparisons) are read for their qualitative conclusions
  (recorded above) only; no pixel/curve data was extracted from them.
- The Supplemental Material videos ([18]) are not available in this snapshot and are
  not referenced by anything in this repo.
- Acknowledgments and the reference list ([1]-[25], p. 022714-17) carry no technical
  content for this repo beyond the two citations already covered by their own notes
  ([1] = `houfek-2006-pra73-032721.md`, [2] = `houfek-2008-pra77-012710.pdf`).
