# W. Domcke, Phys. Rep. 208, 97 (1991) — Theory of resonance and threshold effects in electron-molecule collisions: the projection-operator approach

**Source:** `reference/literature/domcke-1991-physrep208-97.pdf` (gitignored) ·
Elsevier PII `0370-1573(91)90125-6` (read from the PDF's own document title;
no DOI is printed anywhere on the offprint — the Elsevier DOI
<https://doi.org/10.1016/0370-1573(91)90125-6> is the mechanical PII mapping,
not something visible in the source)
**Pagination:** per-journal-volume, printed pages 97–188. Extractor page N =
printed page `96 + N` for all body text. Extractor page 1 is the unnumbered
offprint cover; extractor page 2 carries the article opening (printed 97–98:
running head `PHYSICS REPORTS … 208, No. 2 (1991) 97—188`, contents, abstract);
extractor page 3 is printed 99. Offset checked at extractor 3 (printed 99),
38 (134), 41 (137), 69 (165), 90 (186) and 92 (188) — stable throughout.
**OCR quality:** this PDF is a scanned Elsevier offprint with a TextBridge OCR
layer. The `.txt` extraction is usable for *locating* equations by number but
badly garbles their content (`Vdk.(R)xV.(R)`, `~1’d,E(R)`, …). **Every equation
transcribed below was read off a rendered page image** (`pdftoppm -r 160 -png`),
not from the `.txt`. Do not trust the `.txt` for symbols, superscripts or
complex conjugates.

## Why this repository cares

This is the canonical, book-length statement of the **nonlocal resonance model**
— the Feshbach projection-operator (PO) reduction of resonant electron-molecule
scattering to a one-dimensional nuclear equation with a complex, energy-
dependent, *nonlocal* effective potential. `qscat.core.nrm` implements that
model. It does **not** implement it in Domcke's form: it follows Houfek,
Rescigno & McCurdy, Phys. Rev. A **77**, 012710 (2008)
(`houfek-2008-pra77-012710.md`, "PRA 77"), which differs from this review on two
points that matter to code — the coupling that appears in the resonant T matrix
(the disagreement recorded below), and how the nonlocal potential's singular
energy integral is evaluated. This note exists so both differences have a
tracked home with page locators, instead of a code comment asserting one side.

The review is also where the **local-complex-potential (LCP) limit** of the
nonlocal model is derived, step by step, from the nonlocal kernel. That
derivation (Eq. 4.30c → 4.31 → 4.33 → 4.36) is the published justification for
the bridge test between `qscat.core.nrm` and `qscat.core.lcp`: LCP is what the
nonlocal model *becomes* when the energy dependence of the discrete-state–
continuum coupling is frozen at the fixed-nuclei resonance energy.

## What this repository uses

| Fact | Locator | Used by |
|---|---|---|
| Feshbach projectors `Q = \|φ_d⟩⟨φ_d\|`, `P = 1 − Q`, with `P` spanned by energy-normalized continuum states orthogonal to `φ_d` | p. 105, Eq. (2.8)-(2.10) | the `Q`/`P` split `qscat.core.nrm` builds; same construction as PRA 77 Eq. (15) |
| Discrete-state energy `ε_d`, discrete–continuum coupling `V_dk`, complex level shift `F(k)` | p. 106, Eq. (2.21)-(2.23) | the definitions behind `nrm`'s `V_d(R)`, `V_dn(R)` ingredients |
| Fixed-nuclei width `Γ(E) = 2π ∫dΩ_k \|V_dk\|²` and level shift `Δ(E) = Re F(k)` | p. 107, Eq. (2.26)-(2.29) | the `Γ = 2π\|V_dk\|²` relation `nrm`'s coupling gate rests on; same as PRA 77 Eq. (68) |
| **The nonlocal nuclear wave equation** `[E − T_N − V_d(R)]Ψ_{d,E}(R) − ∫dR' F(R,R';E)Ψ_{d,E}(R') = V_{dk_i}(R)χ_{v_i}(R)` | p. 133, Eq. (4.10) | the equation `qscat.core.nrm` solves; identical to PRA 77 Eq. (52) |
| The nonlocal kernel `F(R,R';E) = ∫k dk dΩ_k V_dk(R) G̃₀⁺(R,R';E − k²/2) V*_dk(R')` | p. 133, Eq. (4.11) | PRA 77 Eq. (53) is the same object; `nrm.nonlocal_potential` builds it |
| **The VE T-matrix element, Eq. (4.14)** — the one PRA 77 disputes | p. 134, Eq. (4.14) | recorded verbatim below; qscat implements PRA 77's replacement |
| Alternative resolvent form of the same T matrix, `T = ⟨v_f\|V*_{dk_f} G(E) V_{dk_i}\|v_i⟩` | p. 134, Eq. (4.16)-(4.17) | shows the disputed conjugate propagates into the resolvent form too |
| DA amplitude from the asymptotics, `lim_{R→∞} Ψ_{d,E}⁺(R) = T_DA e^{iKR}` | p. 134, Eq. (4.18) | the asymptotic read-out `nrm`'s σ_DA uses |
| **VE and DA integral cross sections** `σ = (4π³/k_i²) ν \|…\|²`, and the equivalent `σ_DA = (2π²/k_i²) ν (K/μ) \|lim Ψ_d\|²` | p. 136, Eq. (4.25)-(4.28) | the prefactor parity check below |
| Nonlocal kernel split into level shift and width, `F = Δ − iΓ/2`, both as vibrational-state sums | p. 136, Eq. (4.30a)-(4.30c) | the "state sum" `nrm`'s truncation-convergence work has to converge |
| **The LCP limit**: freeze `E − ε_v` at `E_res(R)` in Eq. (4.30c) and use vibrational completeness → local `Γ_L(R) = 2π\|V_{d,E_res(R)}(R)\|²`, `Δ_L(R)`, and `V_eff^(L)(R) = V_d(R) + Δ_L(R) − iΓ_L(R)/2` | p. 137, Eq. (4.31)-(4.33); p. 138, Eq. (4.34)-(4.36) | the published nonlocal→LCP bridge; the basis for comparing `qscat.core.nrm` against `qscat.core.lcp` |
| `Re V_eff^(L)(R)` **is** the fixed-nuclei resonance potential-energy curve | p. 138 | why qscat's pole-derived `V_d(R) = Re E_pole(R)` is the same object as `V_d + Δ_L` — see Parity below |
| LCP's threshold failure: it ignores the threshold condition `V_{d,E−ε_v} = 0` for `ε_v > E`, so multi-channel S-matrix unitarity is violated near threshold | p. 138 | published explanation for qscat's own documented near-threshold LCP departure (`docs/physics/diatomic-ve-cross-sections.md`) |
| LCP's elastic failure in e+N₂: agreement is excellent for inelastic channels but the *elastic* channel deviates significantly (Dubé & Herzenberg; Hazi et al.), largely repaired by the semilocal approximation | p. 165 | independent confirmation of qscat's "LCP misses the non-resonant elastic background" finding, for the same molecule |

## Equations

Transcribed from rendered page images (see the OCR caveat above). Atomic units
throughout (`ħ = 1`), the same convention as `libs/qscat/qscat/units.py`.

Projection-operator setup:

```
Q = |phi_d><phi_d| ,  P = 1 - Q                                  p. 105, Eq. (2.8)
P = Int k dk dOmega_k |k^(±)><k^(±)|                             p. 105, Eq. (2.9)
eps_d = <phi_d|H_el|phi_d>                                       p. 106, Eq. (2.21)
V_dk  = <phi_d|H_el|phi_k^(+)>                                   p. 106, Eq. (2.22)
F(k)  = <phi_d|H_el G_bg^(+)(k) H_el|phi_d>                      p. 106, Eq. (2.23)
Im F(k) = -Gamma(E)/2 ,  Delta(E) = Re F(k)                      p. 107, Eq. (2.26)-(2.27)
Gamma(E) = 2 pi Int dOmega_k |V_dk|^2                            p. 107, Eq. (2.28)
Delta(E) = (1/2pi) P Int dE' Gamma(E')/(E - E')                  p. 107, Eq. (2.29b)
```

The nonlocal nuclear problem:

```
[E - T_N - V_d(R)] Psi_dE(R) - Int dR' F(R,R';E) Psi_dE(R')
    = V_dk_i(R) chi_v_i(R)                                       p. 133, Eq. (4.10)
F(R,R';E) = Int k dk dOmega_k V_dk(R) G0~^(+)(R,R';E - k^2/2) V_dk*(R')
                                                                 p. 133, Eq. (4.11)
G0~^(+)(R,R';E) = <R|(E - T_N - V_0 + i eta)^-1|R'>              p. 133, Eq. (4.12)
V_eff(R,R';E) = V_d(R) delta(R-R') + F(R,R';E)                   p. 134, Eq. (4.13)
G(E) = [E - T_N - V_d(R) - F]^-1                                 p. 134, Eq. (4.16)
lim_{R->inf} Psi_dE^(+)(R) = T_DA(k_i,v_i) e^{iKR}               p. 134, Eq. (4.18)
```

**Eq. (4.14) — the disputed T-matrix element**, transcribed exactly as printed
on p. 134:

```
T(k_f, v_f; k_i, v_i) = <v_f|<phi_k_f^(-)|V|Psi_E^(+)>>
                      = Int dR chi_v_f*(R) V_dk_f*(R) Psi_dE^(+)(R)
                                                                 p. 134, Eq. (4.14)
```

and its resolvent restatement, carrying the same conjugated coupling:

```
T(k_f, v_f; k_i, v_i) = <v_f| V_dk_f* G(E) V_dk_i |v_i>          p. 134, Eq. (4.17)
```

Cross sections (`ν` counts the spatial degeneracy of the electronic resonance
state — `ν = 1` for Σ, `ν = 2` for Π):

```
sigma_vf<-vi(E) = (4 pi^3 / k_i^2) nu |<v_f|V_dE_f* G(E) V_dE_i|v_i>|^2
                                                                 p. 136, Eq. (4.25)
sigma_DA(E)     = (4 pi^3 / k_i^2) nu |<K^(-)|V_dE_i|v_i>|^2     p. 136, Eq. (4.26)
sigma_DA        = (2 pi^2 / k_i^2) nu (K/mu) |lim_{R->inf} Psi_dE(R)|^2
                                                                 p. 136, Eq. (4.28)
```

The nonlocal kernel as a vibrational-state sum (the sum runs over bound target
levels *and* the target dissociation continuum):

```
F(R,R';E)     = Delta(R,R';E) - i Gamma(R,R';E)/2                p. 136, Eq. (4.30a)
Delta(R,R';E) = Sum_v Int dE' V_dE'(R) chi_v(R) [E - E' - eps_v]^-1
                       chi_v*(R') V_dE'*(R')                     p. 136, Eq. (4.30b)
Gamma(R,R';E) = 2 pi Sum_v V_{d,E-eps_v}(R) chi_v(R) chi_v*(R')
                       V_{d,E-eps_v}*(R')                        p. 136, Eq. (4.30c)
```

**The LCP limit.** Approximate `E − ε_v` by the fixed-nuclei resonance energy
`E_res(R)` in Eq. (4.30c) and use vibrational completeness
`Sum_v chi_v(R) chi_v*(R') = delta(R−R')` (Eq. 4.32):

```
Gamma(R,R';E) ~= 2 pi |V_{d,E_res(R)}(R)|^2 delta(R - R')        p. 137, Eq. (4.31)
Gamma_L(R)     = 2 pi |V_{d,E_res(R)}(R)|^2                      p. 137, Eq. (4.33)
Delta(R,R';E) ~= Delta_L(R) delta(R - R')                        p. 138, Eq. (4.34)
Delta_L(R)     = P Int dE' |V_dE'(R)|^2 / [E_res(R) - E']        p. 138, Eq. (4.35)
V_eff^(L)(R)   = V_d(R) + Delta_L(R) - i Gamma_L(R)/2            p. 138, Eq. (4.36)
```

`E_res(R)` here is defined either as the pole of the K matrix (Eq. 3.4) or as
the real part of the pole of the S matrix (Eq. 3.1) — p. 137.

Standard LCP additionally replaces the energy-dependent entry/exit amplitudes by
`R`-dependent ones:

```
V_dE(R) -> V_{d,E_res(R)}(R) E_res(R)^{-1/4}                     p. 138, Eq. (4.37)
```

Keeping the energy dependence of the entry/exit amplitudes while applying the
local approximation only inside `G` is called the **semilocal** approximation
(p. 138).

## Parameters and numeric values

The review is a formalism review; it contains no parameter table for the
N₂/NO/F₂ two-dimensional models this repository runs. Those constants come from
`houfek-2006-pra73-032721.md` and `houfek-2008-pra77-012710.md`. Nothing in
`qscat.model.library` is sourced here, so no constant-by-constant parity check
applies.

Two *conventions* were checked against the code:

- **DA cross-section prefactor.** Eq. (4.26) gives
  `σ_DA = (4π³/k_i²) ν |T_DA|²`; with `ν = 1` and `k_i² = 2E` this is
  `4π³|T|²/(2E)`. Verified against the repo:

  ```
  grep -rn "np.pi\*\*3" libs/qscat/qscat/core/
  ```

  `libs/qscat/qscat/core/dissociation.py:230` and `:352`,
  `libs/qscat/qscat/core/lcp.py:491` and `libs/qscat/qscat/core/driven.py:105`
  all compute `4.0 * np.pi**3 * abs(t)**2 / (2.0 * float(e))` — **matches
  Eq. (4.26) exactly (verified 2026-08-17)**. The same prefactor appears in
  Gertitschke & Domcke 1993 Eq. (2.10) (`gertitschke-1993-pra47-1031.md`) and in
  PRA 77 Eq. (13)-(14), so all three published statements and the code agree.
  Note this is the *time-independent* convention; `qscat.core.td_extractors`
  deliberately uses `C_DA = π` for its propagated route
  (`libs/qscat/qscat/core/td_extractors.py:957`), reconciled there via
  `S = 1 − 2πiT` — a different convention, not a disagreement with this paper.

- **Where qscat's LCP curve sits in Eq. (4.36).** `qscat.core.lcp` builds
  `V_d(R) = Re(E_pole(R))`, `Γ(R) = max(0, −2 Im(E_pole(R)))` from the fixed-`R`
  electronic ECS pole (`libs/qscat/qscat/core/lcp.py:14`), and propagates on
  `H_N = T(μ) + V_d(R) − iΓ(R)/2` (`lcp.py:676`). Checked against p. 138:
  the review states there that `Re V_eff^(L)(R)` *is* the fixed-nuclei resonance
  potential-energy curve, so qscat's single pole-derived real part corresponds to
  Domcke's `V_d(R) + Δ_L(R)` **combined** — qscat never forms `Δ_L(R)` (Eq. 4.35)
  separately. **Not a numerical parity check** (the two are different routes to
  the same curve, not the same expression), but it fixes which of Domcke's
  symbols qscat's `V_d` maps onto; see the Terminology map.

## Findings and limits

- **The disagreement with PRA 77 — the load-bearing content of this note.**
  Domcke's Eq. (4.14) (p. 134, transcribed above) writes the VE T-matrix element
  with `V*_{dk_f}`, i.e. the complex conjugate of the coupling `V_dk` defined in
  Eq. (2.22) with the **outgoing** (`+`) background continuum state. PRA 77
  instead derives (its Eq. 31, p. 012710-4)

  ```
  T_res^{v_i -> v_f} = <chi_v_f| V_dk_f^{-*} |Phi_d^+>            PRA 77, p. 012710-4, Eq. (31)
  ```

  with the **incoming** (`−`) continuum state, and states immediately after it,
  in its own words: the matrix element in "Domcke [Ref. [6], Eq. (4.14)] where
  the matrix `V_dk` without a superscript, which corresponds to the matrix
  element `V_dk^+` defined by Eq. (21), was, **in our opinion, used
  incorrectly**" (PRA 77, p. 012710-4, quoted). PRA 77's reasoning, on the same
  page: `⟨φ⁻_{k_f}|φ⁺_k⟩ ≠ δ(k_f²/2 − k²/2)` (its Eq. 32), so
  `⟨φ⁻_{k_f}|PH_el Q|φ_d⟩ ≠ ⟨φ⁺_{k_f}|H_el|φ_d⟩` (its Eq. 33). Only in the
  special case of a **real** discrete state, where `φ⁻_k(r) = (φ⁺_k(r))*` in the
  radial case (its Eq. 34), do the two collapse — and then to `V⁺_{d k_f}`
  **without complex conjugation** (its Eq. 35). In three dimensions Eq. (34)
  becomes `φ⁻_k⃗ = (φ⁺_{−k⃗})*` (its Eq. 36), "and thus `V*_{d k⃗}` in Eq. (4.14)
  of [6] should be replaced by `V_{d,−k⃗}` under the assumption that `φ_d` is
  real, otherwise `V^{−*}_{d k⃗}` must be used" (PRA 77, p. 012710-4, quoted).
  PRA 77 records that the difference is small on its own but "becomes important
  when the background terms defined below are added to the resonant T matrix,
  since the coupling matrix elements `V^±_dk` are in general complex even when
  the discrete state is real" (p. 012710-4).

  **`qscat.core.nrm` implements PRA 77's form**, not Eq. (4.14)'s. This note
  records both sides; it does not adjudicate the physics beyond reporting what
  each source states.

- **The LCP is a limit of this model, not an independent theory.** The review
  derives it in one direction only — freeze the coupling's energy dependence at
  `E_res(R)`, invoke vibrational completeness, and the nonlocal kernel collapses
  to a local complex potential (Eq. 4.31-4.36, pp. 137-138). Domcke states the
  price explicitly: the nonlocality of `Γ` and `Δ` "would disappear if `V_dE`
  were strictly independent of energy" (p. 137), so everything the LCP loses is
  the coupling's energy dependence.

- **Stated limits of the LCP, from the review itself:**
  - The definition is not unique for broad resonances, because the resonance
    energy itself is ambiguous there — "There can thus be no unique local
    effective potential for very short-lived collision complexes" (p. 138).
  - `V_eff^(L)(R)` can be *multiple-valued* over a range of `R`, in which case
    the LCP does not exist as a potential-energy function at all (p. 138).
  - The LCP ignores the threshold condition `V_{d,E−ε_v} = 0` for `ε_v > E`,
    so multi-channel S-matrix unitarity is violated for resonances near
    threshold (p. 138).
  - The additional replacement of entry/exit amplitudes (Eq. 4.37) "generally
    destroys" the proper threshold behavior of inelastic cross sections
    (p. 138). The semilocal approximation exists precisely to avoid that.
- **Where the LCP nevertheless works, in the review's own applications.** For
  the ²Π_g shape resonance in e+N₂ the LCP describes the dynamics "perfectly"
  for inelastic channels; the *elastic* channel deviates significantly, a
  failure "largely eliminated when the semilocal approximation is employed"
  (p. 165). The review attributes the N₂ agreement to non-Born-Oppenheimer
  effects being minor in that particular system, and warns this is "no longer
  the case" for broader resonances, resonances closer to threshold, virtual
  states, and polar molecules (p. 165).

- **Not a settled question in the literature, and the repo should not present it
  as one.** PRA 77 (p. 012710-5) records an explicit disagreement with this
  review over whether the nonlocal theory goes *beyond* the Born-Oppenheimer
  approximation: PRA 77 considers its own smoothness conditions on `φ_d(r;R)`
  and `φ⁺_k(r;R)` to be equivalent to invoking Born-Oppenheimer, and contrasts
  that "with the view expressed by Domcke et al. [6,18], who claim that the
  nonlocal resonance theory goes beyond the Born-Oppenheimer approximation".
  Recorded here for completeness; nothing in qscat depends on the outcome.

## Terminology map

| Paper symbol | qModeling name | Note |
|---|---|---|
| `φ_d(r;R)` | `qscat.core.nrm`'s discrete state | same object; PRA 77's `φ_d` too |
| `V_d(R)` (p. 134, Eq. (4.13); p. 140, Eq. (4.43) — the discrete-state potential) | **not** `qscat.core.lcp`'s `V_d` | a genuine collision — `lcp.py:41` already flags that its `V_d` is not Houfek's; the same warning applies here. Domcke's `V_d` is the *unshifted* discrete-state curve; qscat's is `Re E_pole(R)`, which corresponds to Domcke's `V_d + Δ_L` |
| `Δ_L(R)` (Eq. 4.35) | never formed explicitly in qscat | absorbed into `Re E_pole(R)` |
| `Γ_L(R)` (Eq. 4.33) | `Gamma(R)` in `qscat.core.lcp` | same quantity, different route: Domcke's is `2π\|V_{d,E_res}\|²`, qscat's is `−2 Im E_pole(R)` |
| `V_eff^(L)(R)` (Eq. 4.36) | `V_d(R) − iΓ(R)/2` in `qscat.core.lcp` | same object |
| `F(R,R';E)` (Eq. 4.11) | `qscat.core.nrm`'s nonlocal potential | PRA 77 Eq. (53); qscat builds it via PRA 77's Eq. (60) |
| `V_dk`, `V*_dk` (Eq. 2.22, Eq. 4.14) | `V_dk⁺` / `V_dk^{−*}` following PRA 77 | **the disagreement** — see Findings |
| `Ψ_{d,E}⁺(R)` (Eq. 4.9) | the resonant nuclear wavefunction `nrm` solves for | PRA 77 writes `Φ_d⁺(R)` |
| `ν` (spatial degeneracy, Eq. 4.25-4.26) | not carried as a variable in qscat | qscat's models are single-channel with `ν = 1` implied by the `4π³/2E` prefactor |
| `T_N`, `V_0(R)`, `χ_v(R)` | `qscat`'s nuclear kinetic operator, `v0`, vibrational levels | same |

## Not used here

- **The analytic route to the singular energy integral.** Expanding the nuclear
  Green's function (Eq. 4.12) in target vibrational eigenstates (Eq. 4.29) turns
  the nonlocal kernel into the principal-value energy integral of Eq. (4.30b)
  and, for the dissociative case, the coupled integral equations of Eq. (4.57)
  (p. 142). The review makes that tractable by a **separable ansatz** for the
  coupling's energy dependence,
  `V_dE(R) = Σ_{i=1..N} f_i(E) g_i(R)` (p. 143, Eq. 4.60-4.61), with the `f_i(E)`
  "chosen such that threshold laws are fulfilled and that the Hilbert transform
  of `|f_i(E)|²` can be calculated analytically" (p. 143). **qscat does not take
  this route.** PRA 77 says why, at p. 012710-5: such an expansion "leads to a
  singular integral over electron energies, which is difficult to treat unless
  one assumes a particular energy dependence of the discrete-state–continuum
  coupling `V_dk⁺(R)` to be able to evaluate this singular integral analytically
  (see, e.g., Ref. [6])" — Ref. [6] being this review. PRA 77's ECS + DVR
  discretization of the `P` space (its Eq. 55-61, pp. 012710-5–6) evaluates
  `F(E,R,R')` without any such assumption, because the discretized electronic
  energies `E_n(R)` are complex while the total energy `E` is real, so
  `M(n) = E − T_R − V_0 − E_n` is never singular (PRA 77, p. 012710-6). That is
  what `qscat.core.nrm` implements.
- The whole of section 3 (pp. 112-131) — analytic continuation of fixed-nuclei
  potential-energy curves, threshold behavior of S-matrix singularities, virtual
  states, polar molecules, electron-ion resonances. Read only for the definition
  of `E_res(R)` used at p. 137; none of the threshold-law analysis is ported.
- Section 4.3's computational machinery beyond the point above: the Lanczos
  basis for the target dissociation continuum (p. 143), the harmonic-oscillator
  and shifted-oscillator model reductions (p. 140, Eq. 4.45-4.47), and the
  semiclassical Green's-function shortcuts (p. 143).
- **Section 5, the time-dependent formulation** (pp. 147-161) — equation of
  motion, time-dependent LCP and beyond, friction and memory effects. Not used
  by this repository's time-independent NRM work. The repository's source for
  the time-dependent nonlocal treatment of DA is
  `gertitschke-1993-pra47-1031.md`, which builds directly on this section.
- Section 6's applications (pp. 162-179) other than the two conclusions recorded
  above (N₂ elastic-vs-inelastic LCP behavior, p. 165). No curve or figure data
  was extracted; §6.2's F₂ discussion (pp. 165-167) is a survey of other
  authors' calculations, not a parameter source for `qscat.model.library.F2`.
- The review's bibliography (pp. 181-188) is not mined; where a claim here is
  attributed to a third author (Bardsley, Dubé & Herzenberg, Hazi et al.), the
  attribution is the review's own and was not chased to the primary source.
