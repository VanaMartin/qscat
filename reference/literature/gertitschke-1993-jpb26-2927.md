# P. L. Gertitschke, W. Domcke, J. Phys. B 26, 2927 (1993) — Systematically improved local complex potential approximation for the dynamics of electron-molecule collision complexes

**Source:** `reference/literature/gertitschke-1993-jpb26-2927.pdf` (gitignored) ·
[DOI 10.1088/0953-4075/26/17/024](https://doi.org/10.1088/0953-4075/26/17/024)
(confirmed via the IOPscience volume-26/issue-17 table of contents; the article
page itself is behind bot-protection and could not be fetched directly to
cross-check title/authors a second time)
**Pagination:** per-journal-volume, 2927–2942. Extractor page N = printed page
`2926 + N` (offset checked at both ends: extractor page 1 shows the printed
footer "2927", extractor page 16 shows the printed header "2942").
**OCR quality:** this PDF is a scanned capture (`Acrobat Capture 3.0`, no text
layer beyond OCR); many equations extract with substituted/garbled characters
(e.g. `AL`/`rL` for `Δ_L`/`Γ_L`, stray digits, dropped operators). Equations
below are transcribed only where the OCR is unambiguous; where it is not, the
content is described in prose instead of presented as a literal transcription.

## Scope caveat — this is NOT the paper qModeling's docs cite as "Domcke 1991"

`docs/superpowers/specs/2026-08-15-bo-lcp-resonance-levels-design.md` (and, by
inheritance, other qModeling docs that cite "Domcke 1991" for the LCP
formalism) refer to **W. Domcke, *Theory of resonance and threshold effects in
electron-molecule collisions: the projection-operator approach*, Phys. Rep.
**208**, 97 (1991)** — a review article. That review is a different paper, and
it **is** now held in this collection: see
[`domcke-1991-physrep208-97.md`](domcke-1991-physrep208-97.md). Cite that note
for the formalism. This note is for the *different*, later Gertitschke & Domcke
*J. Phys. B* paper, which is a numerical companion: it repeatedly cites "Domcke
1991" itself (e.g. p. 2928, p. 2930, p. 2932) as the source of the
time-dependent projection-operator formalism and the LCP/FOCLCP derivation,
stating that formalism "has been outlined in detail previously (Domcke 1991)"
(p. 2930) rather than re-deriving it here. Do not use this note as a substitute
citation for Phys. Rep. 208, 97 (1991); the two are related but distinct
sources.

Note also that a **third** related paper is held: Gertitschke & Domcke,
Phys. Rev. A **47**, 1031 (1993) — the *time-dependent* wave-packet treatment of
DA ([`gertitschke-1993-pra47-1031.md`](gertitschke-1993-pra47-1031.md)). Same
authors, same year, different journal and different content.

## Why this repository cares

qModeling's `qscat.core.lcp` implements the local-complex-potential (LCP)
reduction of dissociative attachment / vibrational excitation as the
"approximation under test" against the exact 2-D solver (see
`docs/physics/diatomic-ve-cross-sections.md`,
`docs/physics/lcp-resonance-levels.md`). This paper is the clearest published
account of exactly what the "standard" LCP approximation *is* as a
time-dependent wavepacket equation (§2, Eq. 1–13) and of the physical
mechanism the repo's own documented LCP departures echo: LCP fails when a
resonance is broad and/or close to threshold, because it assumes
exponential (memoryless) decay from t=0, which breaks down at short times.
The repo's own finding — LCP under-predicts F₂'s near-threshold DA spike
(`docs/physics/diatomic-ve-cross-sections.md`, "LCP under-predicts the
near-threshold spike") and misses the elastic non-resonant background — is
the same class of failure this paper documents quantitatively for e⁻+H₂
(§4): the standard LCP DA peak is too large by a factor of 14, and the
elastic cross section has the wrong shape, both traced to the wrong (purely
exponential) short-time behaviour of the survival probability `S(t)`. This
paper does not itself supply any parameter, formula or numeric result
presently implemented in `qscat`; its value here is corroborating,
citable evidence for the LCP validity limits the repo already documents from
its own oracle comparison.

## What this repository uses

| Fact | Locator | Used by |
|---|---|---|
| WP equation of motion under the projection-operator (PO) theory, `i d/dt Psi_d(R,t) = [T_N + V_d(R)]Psi_d(R,t) + (memory-kernel term)` | p. 2928, Eq. (1) | conceptual background for `qscat.core.lcp`'s "reduces the full non-local problem to a local potential" framing |
| The LCP approximation is the t→short-memory-time (delta-function memory kernel) limit of the exact non-local PO equation | p. 2930, Eq. (12)–(13) | corroborates why `local_complex_potential`/`lcp_da_cross_section` is documented as an *approximation under test*, not the exact theory |
| Separable discrete-continuum coupling `V_dE(R) = f(E) g(R)` and the width relation `Gamma(E) = 2*pi*\|f(E)\|^2` | p. 2930 Eq. (15); p. 2932 Eq. (27a) | see "Entry/exit amplitude" below — the `sqrt(Gamma/2pi)` structure `qscat.core.lcp.lcp_da_cross_section`'s doorway uses, in a related but not identical parametrization |
| Standard LCP fails for broad, near-threshold resonances (elastic cross-section shape, DA magnitude) | p. 2941, §4 discussion and Fig. 7 caption | corroborates (does not itself establish) the repo's own documented near-threshold LCP departure for F₂/NO DA (`docs/physics/diatomic-ve-cross-sections.md`) |
| Standard LCP performs comparably well for a narrow, not-too-near-threshold resonance (the N₂ ²Π_g model) | p. 2935–2936, §3 | corroborates why LCP is a *good* approximation for N₂-like cases in the repo's own comparison |

## Equations

Legibly transcribed (atomic units, `hbar = 1`, stated explicitly p. 2928):

```
i d/dt Psi_d(R,t) = [T_N + V_d(R)] Psi_d(R,t)
                     + Int dt' Int dR' F(R,R',t-t') Psi_d(R',t')      p. 2928, Eq. (1)
T_N = -(2 mu)^-1 d^2/dR^2                                             p. 2928, Eq. (2)
Psi_d(R, t=0) = V_{d,k_i}(R) chi_vi(R)                                p. 2929, Eq. (5)
[E - T_N - V_d(R)] Psi_{d,E}(R)
    - Int dR' F(R,R',E) Psi_{d,E}(R') = V_{d,E_i}(R) chi_vi(R)        p. 2929-2930, Eq. (10)
V_dE(R) = f(E) g(R)                                                   p. 2930, Eq. (15)
Gamma(E) = 2*pi*|f(E)|^2                                              p. 2932, Eq. (27a)
S(t) = <Psi_d(t)|Psi_d(t)>                                            p. 2934, Eq. (34)
g(R) = [1 + C(R - R1)^2]^-1     (H2 model discrete-continuum coupling) p. 2937, Eq. (37)
```

Not reliably transcribable from the OCR (garbled beyond safe reconstruction),
described in prose instead:

- **Eq. (3)–(4), p. 2928**: the memory kernel `F(R,R',t)` as a continuum sum
  over electron momentum `k` of `V_dk(R) exp(-i k^2 t/2) K(R,R',t) V_kd(R')`,
  with `K` the Feynman propagator of the target-molecule nuclear motion under
  `T_N + V0(R)`.
- **Eq. (6)–(9), p. 2929**: the VE and DA T-matrix elements as Laplace
  transforms in time of `Psi_d(R,t)` projected onto final vibrational/exit-
  channel states, and the resulting integral cross sections. The cross-section
  formulas (Eq. 8–9) did not survive OCR extraction at all (blank in the
  `.txt`) — **no cross-section formula from this paper is transcribed here**;
  only their existence and role (T-matrix → cross section via a Laplace
  transform of the propagated wavepacket) is used.
- **Eq. (12)–(13), p. 2930**: the standard LCP approximation, obtained by
  replacing the memory kernel with an instantaneous (delta-function-in-time,
  delta-function-in-`R`) term carrying a real level-shift function `Delta_L(R)`
  and a width function `Gamma_L(R)`, reducing Eq. (1) to a local, energy-
  independent complex potential `V_d(R) + Delta_L(R) - (i/2)*Gamma_L(R)`.
- **Eq. (28), p. 2932**: the implicit equation fixing the fixed-nuclei
  resonance energy `E_res(R)` self-consistently from `V_d(R)` and the
  (energy-dependent) level-shift function — OCR too corrupted to state the
  precise algebraic form; not needed by the repo, which gets `E_res(R)`
  directly from `qscat.ecs.find_resonance_pole` rather than this implicit
  equation.
- **Eq. (29)–(32), p. 2932–2933**: the paper's own contribution, the
  "first-order corrected LCP" (FOCLCP) — a systematic first-order-in-memory-
  time correction to Eq. (13) that reintroduces an explicitly time-dependent
  effective potential (via derivatives of `g(R)` and the energy-derivative of
  the width/shift functions), converging to a renormalized time-independent
  potential only in the long-time limit. **Not implemented anywhere in
  qscat** — see "Not used here".

## Entry/exit amplitude and the `sqrt(Gamma/2pi)` structure — a genuine but partial match

`qscat.core.lcp.lcp_da_cross_section`'s doorway is
`d = sqrt(Gamma/(2*pi)) * chi_{v_init}`, i.e. an entry amplitude built from
`sqrt(Gamma(R)/2pi)` where `Gamma(R)` is the **R-dependent** electronic
autodetachment width from the fixed-`R` resonance pole
(`qscat.core.lcp.local_complex_potential`).

This paper's Eq. (15)/(27a) give the analogous relation for its own
formalism: the discrete-continuum coupling is separable,
`V_dE(R) = f(E) g(R)`, and `f(E)` is fixed by the *energy-dependent* width
via `Gamma(E) = 2*pi*|f(E)|^2`, i.e. `f(E) = sqrt(Gamma(E)/2pi)` up to a
phase. `f(E)` is exactly the "entrance amplitude" that appears in the
initial condition (Eq. 5, `Psi_d(R,t=0) = V_{d,k_i}(R) chi_vi(R) = f(E_i)
g(R) chi_vi(R)`) and the "exit amplitude" in the T-matrix (Eq. 6). This is
the same `sqrt(Gamma/2pi)` algebraic structure the repo's doorway uses.

**The two `Gamma`s are not the same object, and this is worth stating
plainly so it is never conflated:**

| | This paper | `qscat.core.lcp` / Houfek 2006 |
|---|---|---|
| `Gamma` is a function of | electron energy `E` (via the coupling strength `f(E)`, R-independent) | nuclear geometry `R` (the width of the fixed-`R` electronic resonance pole) |
| `R`-dependence enters via | the separate factor `g(R)` in `V_dE(R) = f(E)g(R)` | `Gamma(R)` itself |
| Role in the doorway | `f(E_i)` scales the whole `R`-shape `g(R)` at the incident energy | `sqrt(Gamma(R)/2pi)` scales the R-shape `chi_{v_init}(R)` pointwise in `R` |

Both are legitimate, textbook forms of the LCP entry amplitude (this is the
kind of alternative-LCP-definition ambiguity the paper itself flags, p. 2927,
citing Berman et al 1983a: "the differences between alternative definitions
of the LCP become significant for broad resonances"), but they parametrize
the coupling along different variables. This paper does not itself contain
the `R`-dependent-`Gamma(R)` form the repo implements — that form is Houfek
et al. 2006's (see `reference/literature/houfek-2006-pra73-032721.md`,
Eq. 41, `eta_v(R) = sqrt(Gamma(R)/2pi)`), not this paper's.

## Parameters and numeric values

None. This paper cites its own model parameters (the N₂ ²Π_g harmonic model,
p. 2933–2934; the H₂ ²Σ_u⁺ Morse model, p. 2936–2937) only by reference to
other works not held in this collection — Estrada and Domcke (1989) for the
N₂ model, Berman et al (1985) and Gertitschke and Domcke (1991, 1993) for the
H₂ model — and tabulates no numeric constants of its own. Nothing here was
checked against `qscat.model` because there is nothing to check.

## Findings and limits

- **Standard LCP is quantitatively good for a narrow, not-too-near-threshold
  resonance.** For the model ²Π_g shape resonance in e⁻+N₂, "the differences
  between the results obtained with the standard LCP approximation, the
  FOCLCP approximation and the exact theory are rather small" — p. 2935. The
  main deficiency is in the *elastic* channel shape (p. 2935–2936), traced to
  the entrance/exit-amplitude approximation, not the LCP's memoryless
  short-time assumption per se.
- **Standard LCP fails seriously for a broad, near-threshold resonance.**
  For the model ²Σ_u⁺ shape resonance in e⁻+H₂ — "the prototypical example
  of a low-lying broad shape resonance" (p. 2936) — the standard LCP: (a)
  gives a survival probability `S(t)` with unphysical *finite* initial decay
  slope at `t=0` (exact `S(t)` has zero slope there), the qualitative
  signature of the memoryless assumption breaking down, p. 2937–2938; (b)
  overestimates the asymptotic (DA) survival probability by a factor of ~21,
  p. 2938; (c) overestimates the peak DA cross section by a factor of **14**
  (reduced to a factor of 4 by the paper's own FOCLCP correction), p. 2941;
  (d) gets the elastic cross-section shape wrong, p. 2941.
- **A temporary wavepacket splitting in the exact dynamics cannot be
  reproduced by *any* local theory** (standard LCP or FOCLCP), because it
  signifies "the complete breakdown of the Breit-Wigner single-pole
  approximation" underlying every local reduction — p. 2938–2939.
  Stated as an explicit limit on the whole LCP *concept*, not just the
  standard variant.
- **The paper's own scope**: it is a numerical test of one specific proposed
  improvement (FOCLCP) against the standard LCP and an exact non-local
  reference, for two illustrative 1-D nuclear models (N₂-like narrow/
  N₂-like-far-from-threshold vs H₂-like broad/near-threshold) — not a general
  theory paper (that role is filled by the Phys. Rep. 208, 97 (1991) review
  it repeatedly cites — `domcke-1991-physrep208-97.md`; see the scope caveat
  above).

## Terminology map

| Paper symbol | qModeling name / note |
|---|---|
| `Psi_d(R,t)` | doorway wavepacket concept behind `qscat.core.lcp.lcp_da_cross_section`'s `psi_sc` (TI resolvent form here, not this paper's TD propagation) |
| `V_d(R)` (fixed-nuclei discrete-state / resonance potential) | `qscat.core.lcp`'s `Vd` / `V_d(R)` — same role, same name |
| `Delta_L(R)`, `Gamma_L(R)` (LCP level-shift, width) | `qscat.core.lcp`'s `Gamma(R)` is this paper's `Gamma_L(R)` in spirit; the repo has no separate real level-shift `Delta_L(R)` term — `qscat.ecs.find_resonance_pole`'s `Re(E_pole)` already plays that role inside `Vd` |
| `f(E)`, `g(R)` (separable discrete-continuum coupling) | no direct qscat analogue by that split; see "Entry/exit amplitude" section — the repo's `sqrt(Gamma(R)/2pi)` plays an analogous but differently-parametrized role |
| `Gamma(E) = 2*pi*|f(E)|^2` | structurally the same relation as Houfek 2006 Eq. (41) `eta_v(R) = sqrt(Gamma(R)/2pi)`, but along `E` here vs `R` there/here — see table above |
| FOCLCP (first-order corrected LCP) | not implemented in qscat; see "Not used here" |
| `T_N`, `mu` | `qscat.dvr.kinetic`/`kinetic_sparse`, `model.mu` |

## Not used here

- The full time-dependent PO derivation (§2, Eq. 1–32) beyond the specific
  facts pulled out above — the repo's TI resolvent form
  (`qscat.core.lcp.lcp_da_cross_section`) is not derived from this paper; it
  follows Houfek et al. 2006's `R`-dependent LCP directly (see that note).
- The **FOCLCP (first-order corrected LCP) approximation itself** (Eq. 29–32
  and its H₂/N₂ numerical demonstration, Fig. 1–7) — this is the paper's own
  novel contribution and is **not implemented anywhere in qscat**. It is a
  candidate follow-on if a documented LCP departure (e.g. F₂'s near-threshold
  DA under-prediction) is ever chased further, but nothing in the current
  codebase depends on it.
- The N₂ and H₂ model parameter sets (harmonic N₂ discrete-state model,
  Morse H₂ discrete-state model) — neither is tabulated in this paper itself
  (both deferred to other, uncollected references), and neither model is the
  Houfek-et-al. Morse+sigmoid+Gaussian form `qscat.model` implements.
  These are illustrative 1-D toy models for this paper's own numerical test,
  unrelated to `qscat.model.{N2,NO,F2}`'s parametrization.
- Fig. 1–7 (survival probability, mean position/momentum, wavepacket
  snapshots, cross-section comparison curves) — read for their qualitative
  conclusions only (summarized in Findings above); no pixel/curve data
  extracted.
- The reference list itself (citations to Domcke 1991, Estrada and Domcke
  1989, Berman et al 1983a/1985, Cederbaum and Domcke 1981, etc.) — useful
  as pointers to sources this repo does not hold, not as content in its own
  right.
