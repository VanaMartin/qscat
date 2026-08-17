# P. L. Gertitschke, W. Domcke, Phys. Rev. A 47, 1031 (1993) — Time-dependent wave-packet description of dissociative electron attachment

**Source:** `reference/literature/gertitschke-1993-pra47-1031.pdf` (gitignored) ·
[DOI 10.1103/PhysRevA.47.1031](https://doi.org/10.1103/PhysRevA.47.1031)
**Pagination:** per-journal-volume, printed pages 1031–1044. Extractor/rendered
page N = printed page `1030 + N`. Offset checked at both ends: rendered page 1
carries the printed footer `47 1031`, rendered page 14 the printed header
`1044`. No mid-document shift (spot-checked at rendered 9 → printed 1039 and
rendered 12 → printed 1042).
**OCR quality:** **this PDF is a scan with no text layer.** `pypdf` extracts 285
bytes across all 14 pages — effectively nothing. There is no usable `.txt` to
grep; do not go looking for one. This note was written from rendered page images
(`pdftoppm -r 150 -png`), and every locator below is the printed page number
visible in those images.

**Not to be confused with** `gertitschke-1993-jpb26-2927.md` — a *different*
Gertitschke & Domcke paper from the same year (J. Phys. B **26**, 2927, on the
systematically-improved / FOCLCP local approximation). This one is the
time-dependent wave-packet paper in Phys. Rev. A **47**, 1031.

## Why this repository cares

**Nothing in the repository cites this paper today.** It is ingested ahead of
its consumer, deliberately: it is the time-dependent formulation of exactly the
problem `qscat.core.nrm` solves time-independently — dissociative attachment
within the nonlocal projection-operator model — and the planned time-dependent
NRM sub-project will be built from it. Recording it now, while the PDF is in
hand and the pages have been read, is cheaper than reconstructing the locators
later.

What it will supply when that work starts: the time-dependent nonlocal equation
of motion with its memory kernel (Eq. 2.1-2.4), the Laplace-transform route from
a single energy-independent propagated wave packet to VE and DA amplitudes at
all energies (Eq. 2.17-2.19), and the LCP limit written in the time domain
(Eq. 2.11-2.15) — the time-dependent sibling of the frequency-domain LCP limit
in `domcke-1991-physrep208-97.md`.

Independent of that, the paper carries one result the repository's existing LCP
work can already lean on: a quantitative, mechanism-level account of *why* the
LCP fails, in a case where it fails by more than an order of magnitude.

## What this repository uses

| Fact | Locator | Used by |
|---|---|---|
| **Time-dependent nonlocal equation of motion** for the wave packet `Ψ_d(R,t)` on the discrete state, with a memory (time-nonlocal) kernel | p. 1032, Eq. (2.1) | forward-looking: the equation a TD-NRM sub-project would propagate |
| Memory kernel `F(R,R',t)` built from the Feynman propagator of the *target* nuclear motion | p. 1032, Eq. (2.3)-(2.4) | forward-looking; the time-domain counterpart of Domcke Eq. (4.11) |
| Initial condition `Ψ_d(R,0) = V_{dk_i}(R) χ_v(R)` | p. 1032, Eq. (2.5) | forward-looking: what a TD-NRM launch state is |
| VE and DA T matrices by Laplace transform of the propagated packet | p. 1032, Eq. (2.6), (2.8)-(2.9) | forward-looking |
| Cross-section prefactors `σ = ν(4π³/k_i²)|T|²` for both VE and DA | p. 1032, Eq. (2.7), (2.10) | prefactor parity check below |
| **The LCP limit in the time domain**: `F(R,R',t) = i[Δ_L(R) − (i/2)Γ_L(R)]δ(R−R')δ(t)` — i.e. LCP is the *Markovian* approximation to Eq. (2.1) | p. 1032, Eq. (2.11)-(2.15) | the sharpest one-line statement of what LCP throws away; complements Domcke's frequency-domain derivation |
| The `E_res(R)` defining equation `E_res(R) − V_d(R) + V_0(R) − Δ(E_res(R),R) = 0` | p. 1032, Eq. (2.14) | how `E_res(R)` is pinned down in this formalism |
| Wigner threshold parametrization `Γ(E) = A(E/B)^α exp(−E/B)`, `α = l + 1/2` | p. 1033, Eq. (2.21) | forward-looking: the model coupling form a TD-NRM toy model could adopt |
| One energy-independent wave packet suffices: rescaling by `[Γ(E_i)/2π]^{-1/2}` removes the incident-energy dependence from the initial condition | p. 1033, Eq. (2.17a)-(2.19) | forward-looking: propagate once, extract all energies — the same economy `qscat.core.td_extractors` already exploits |
| **Finding: LCP is quantitatively reliable for the d-wave model** (cross sections "differ only little from the exact ones") | p. 1036, Fig. 3; p. 1042 | the control case — LCP failure is not universal |
| **Finding: LCP fails badly for the ²Σ_u⁺ resonance in e+H₂** — DA peak too large by `≃14×`, energy-integrated DA too large by `≃23×`, and a broader profile | p. 1039 | the headline result; the quantitative anchor |
| **Finding: LCP also gets the H₂ VE cross sections wrong** in both absolute and relative magnitude (LCP curves are plotted rescaled by ×0.3 to fit), though the deeply-inelastic fine structure is qualitatively reproduced | p. 1039, Fig. 5 | independent confirmation that LCP's VE errors are not a fixed scale factor — the same conclusion qscat's F₂/NO runs reach from the other direction |
| **The mechanism**: the exact nonlocal packet undergoes a pronounced temporary *splitting* between `≃2` and `≃5` fs, delaying dissociation; the LCP, being a single complex potential curve, cannot reproduce it, so it dissociates too fast and over-predicts the asymptotic anion population | p. 1041; p. 1042 | the physical explanation behind the numbers |
| The failure is generic and cannot be repaired by redefining the LCP | p. 1042 | limits how far an "improved LCP" can be pushed |

## Equations

Transcribed from rendered page images. **The paper works in eV / Å / fs, not
atomic units** (see p. 1035 and the Fig. 2 caption on p. 1036), with `ħ = 1`
stated at Eq. (2.1). Anything ported must be converted to
`libs/qscat/qscat/units.py`'s atomic units.

Time-dependent nonlocal dynamics:

```
i d/dt Psi_d(R,t) = [T_N + V_d(R)] Psi_d(R,t)
    + (1/i) Int_0^t dt' Int dR' F(R,R',t-t') Psi_d(R',t')        p. 1032, Eq. (2.1)
T_N = -(2 mu)^-1 d^2/dR^2                                        p. 1032, Eq. (2.2)
F(R,R',t) = Int k dk dOmega_k V_dk(R) e^{-i(k^2/2)t} K(R,R',t) V_kd(R')
                                                                 p. 1032, Eq. (2.3)
K(R,R',t) = <R| e^{-i(T_N + V_0)t} |R'>                          p. 1032, Eq. (2.4)
Psi_d(R,0) = V_dk_i(R) chi_v(R)                                  p. 1032, Eq. (2.5)
```

Amplitudes and cross sections:

```
T_v'v(E)  = (1/i) Int_0^inf dt e^{iEt} Int dR chi_v'*(R) V_dk_f*(R) Psi_d(R,t)
                                                                 p. 1032, Eq. (2.6)
sigma_v'v(E) = nu (4 pi^3 / k_i^2) Int dOmega_k_f |T_v'v(E)|^2_bar
                                                                 p. 1032, Eq. (2.7)
T_DA = [K/(2 pi mu)]^{1/2} lim_{R->inf} e^{-iKR} Int_0^inf dt Psi_d(R,t) e^{iEt}
                                                                 p. 1032, Eq. (2.8)
T_DA = [mu/(2 pi K)]^{1/2} lim_{t->inf} e^{iEt} Int dR Psi_d(R,t) e^{-iKR}
                                                                 p. 1032, Eq. (2.9)
sigma_DA(E) = nu (4 pi^3 / k_i^2) |T_DA|^2_bar                   p. 1032, Eq. (2.10)
```

(the overbar denotes averaging `|T|²` over target orientation; `ν` is the
spatial degeneracy of the resonance state.)

**The LCP limit, in the time domain** — the memory kernel collapses to a double
delta in space *and time*, i.e. LCP is the Markovian approximation:

```
F(R,R',t) = i[Delta_L(R) - (i/2) Gamma_L(R)] delta(R-R') delta(t)
                                                                 p. 1032, Eq. (2.11)
Gamma_L(R) = Gamma(E_res(R), R)                                  p. 1032, Eq. (2.12a)
Delta_L(R) = Delta(E_res(R), R)                                  p. 1032, Eq. (2.12b)
Gamma(E,R) = 2 pi Int dOmega_k |V_dk(R)|^2                       p. 1032, Eq. (2.13a)
Delta(E,R) = P Int (dE'/2pi) Gamma(E',R)/(E - E')                p. 1032, Eq. (2.13b)
E_res(R) - V_d(R) + V_0(R) - Delta(E_res(R), R) = 0              p. 1032, Eq. (2.14)
i d/dt Psi_d(R,t) = [T_N + V_d(R) + Delta_L(R) - (i/2) Gamma_L(R)] Psi_d(R,t)
                                                                 p. 1032, Eq. (2.15)
```

Threshold parametrization and the propagate-once trick:

```
Gamma(E,R) = Gamma(E) g(R)^2       (separability assumed)        p. 1033, Eq. (2.16)
Gamma(E)   = A (E/B)^alpha exp(-E/B),  alpha = l + 1/2           p. 1033, Eq. (2.21)
Psi~_d(R,t) = [Gamma(E_i)/2pi]^{-1/2} Psi_d(R,t)                 p. 1033, Eq. (2.17a)
Psi~_d(R,0) = g(R) chi_v(R)          (energy-independent)        p. 1033, Eq. (2.17b)
gamma(t)   = Int (dE/2pi) Gamma(E) e^{-iEt}
           = (AB/2pi) Gamma(1+alpha) / (1+iBt)^{1+alpha}         p. 1033, Eq. (2.23)
```

Dynamical observables used to diagnose the LCP failure:

```
S(t)   = <Psi~_d(t)|Psi~_d(t)> = Int dR |Psi~_d(R,t)|^2          p. 1036, Eq. (4.4)
<R>_t  = Int R dR |Psi~_d(R,t)|^2 / S(t)                         p. 1036, Eq. (4.5)
<P>_t  = Int dR Psi~_d*(R,t) (-i d/dR) Psi~_d(R,t) / S(t)        p. 1036, Eq. (4.6)
```

## Parameters and numeric values

The paper's **first** application is a synthetic test model (DA via a `d`-wave
shape resonance with a repulsive discrete-state curve), whose constants are
published in full on p. 1035, in **eV and Å**:

| Parameter | Value | Locator |
|---|---|---|
| `D_0` | 5 eV | p. 1035, Eq. (4.1) |
| `D_d` | 3.2 eV | p. 1035, Eq. (4.2) |
| `Q_0` | 5 eV | p. 1035, Eq. (4.1) |
| `Q_d` | 1.2 eV | p. 1035, Eq. (4.2) |
| `t_d` | 0.2 | p. 1035, Eq. (4.2) |
| `α_0 = α_d` | 1.96 Å⁻¹ (printed as "1.96 Å") | p. 1035 |
| `μ` | 1 amu | p. 1035 |
| `α` (threshold exponent) | 2.5 | p. 1035 |
| `A` | 1.757 eV | p. 1035 |
| `B` | 1.667 eV | p. 1035 |
| `C` | 0.98 Å⁻² | p. 1035 |

with `V_0(R) = D_0(e^{−2α_0(R−R_0)} − 2e^{−α_0(R−R_0)}) + Q_0`,
`V_d(R) = D_d(e^{−2α_d(R−R_0)} − 2 t_d e^{−α_d(R−R_0)}) + Q_d` (Eq. 4.1-4.2) and
`g(R) = exp[−c(R + 1/(2√C))²]` (Eq. 4.3). Note: the paper writes both `c` and
`C` in Eq. (4.3) and gives a value only for `C`; whether the lowercase `c` is a
typo for `C` **could not be determined from the source** and is recorded here as
unresolved rather than guessed.

**Checked against the repo: none of these constants appear in qscat.**

```
grep -rnF -e "1.757" -e "1.667" -e "0.98" -e "1.96" libs/qscat/qscat/model/library.py
```

returns nothing (exit 1) — as expected. (Use `-F`: without it, `.` is a regex
wildcard and `0.98` spuriously matches `D0=0.05980`.) `qscat.model.library` carries the Houfek
N₂/NO/F₂/H₂⁺ models (`houfek-2006-pra73-032721.md`,
`houfek-2008-pra77-012710.md`), which are a different, unrelated
parametrization. This paper is not a parameter source for anything shipped.

The **second** application (the ²Σ_u⁺ resonance in e+H₂) is built on *ab initio*
scattering data taken from Berman, Mündel & Domcke, Phys. Rev. A **31**, 641
(1985) (the paper's Ref. [35], p. 1043); no parameter table for it is printed
here, so this note carries none.

The **cross-section prefactor** is a genuine parity item and was checked:
Eq. (2.7) and Eq. (2.10) both give `σ = ν(4π³/k_i²)|T|²`; with `ν = 1` and
`k_i² = 2E` that is `4π³|T|²/(2E)`.

```
grep -rn "np.pi\*\*3" libs/qscat/qscat/core/
```

`libs/qscat/qscat/core/dissociation.py:230` and `:352`,
`libs/qscat/qscat/core/lcp.py:491` and `libs/qscat/qscat/core/driven.py:105` all
compute `4.0 * np.pi**3 * abs(t)**2 / (2.0 * float(e))` — **matches Eq. (2.10)
exactly (verified 2026-08-17)**, and agrees with Domcke Eq. (4.26) and Houfek
et al. 2008 Eq. (13)-(14). Three independent published statements and the code
all use the same energy-normalized convention.

## Findings and limits

- **The LCP failure for the H₂⁻ ²Σ_u⁺ resonance is large and quantified.** The
  LCP DA cross section peaks `≃14×` above the exact nonlocal result and, once
  integrated over energy, exceeds it by `≃23×`; the LCP profile is also broader
  (p. 1039). The VE cross sections are wrong in both absolute and relative
  magnitude — Fig. 5 plots the LCP curves rescaled by ×0.3 simply to fit them in
  the frame (p. 1039).
- **The mechanism, which is the paper's actual contribution.** In the exact
  nonlocal (non-Markovian) treatment the wave packet undergoes a pronounced
  *temporary splitting* between `≃2` and `≃5` fs; this delays dissociation, so
  the fast autodetachment of the complex suppresses the asymptotic H₂⁻
  population and hence σ_DA. The LCP, having no splitting, dissociates too fast
  and over-predicts (p. 1041). The smaller asymptotic fragment momentum in the
  exact treatment also narrows the DA profile relative to the LCP (p. 1041).
  The paper interprets the splitting as motion driven by *multiple* complex
  effective potential curves — the analytically continued S matrix has several
  singularities associated with the ²Σ_u⁺ resonance, i.e. the Breit-Wigner
  single-pole picture breaks down; it notes that independent *ab initio*
  attempts to locate the pole have themselves returned multiple solutions
  (p. 1041-1042).
- **The paper's own limit on repairing the LCP:** "the complete failure of the
  LCP model for the H₂⁻ ²Σ_u⁺ resonance is a generic feature of the dynamics of
  this system and cannot simply be overcome by some suitable redefinition of the
  LCP" (p. 1042, quoted). That is a stronger statement than "this particular LCP
  variant is inaccurate", and it bounds what an improved local model can be
  expected to achieve.
- **The control case matters as much as the failure.** For the synthetic `d`-wave
  model the LCP is "quite satisfactory", giving cross sections that "differ only
  little from the exact ones" (p. 1036), and the paper's summary confirms LCP
  "to be quantitatively reliable for this particular example" (p. 1042). Its
  exact results also reproduce independent time-independent calculations of the
  same model (p. 1036), which is how the paper validates its own propagator.
  LCP failure is therefore system-dependent, not universal — the same posture
  this repository takes toward `qscat.core.lcp`.
- **A frictional effect, reproduced qualitatively even by LCP.** The packet
  centroid initially moves faster than in the decay-free case, then slows in the
  curve-crossing region, ending with a significantly lower asymptotic velocity
  (p. 1036, Fig. 2). Ehrenfest's theorem does **not** hold here — `⟨P⟩_t` is not
  proportional to `d⟨R⟩_t/dt` — because the dynamics is non-Markovian open-system
  dynamics (p. 1036). Worth knowing before anyone builds an Ehrenfest-style
  diagnostic on a nonlocal propagation.
- **The oscillations in H₂ VE are *not* boomerang oscillations.** They arise from
  the extremely rapid broadening of the packet (whose tail "leaks backwards"),
  not from quasiperiodic motion of a localized packet as in the ²Π_g resonance
  in e+N₂ (p. 1042). Directly relevant to how this repository interprets its own
  N₂ boomerang structure (`docs/physics/n2-2d-td-cross-section.md`): the two
  look alike and are not.
- **Scope the paper claims for itself:** a computational scheme, not a new
  theory; for the systems treated, the time-dependent approach brings no saving
  in computational effort over steady-state methods — its value is qualitative
  insight, with the expectation that it becomes computationally superior only
  for systems with more degrees of freedom (p. 1043).

## Terminology map

| Paper symbol | qModeling name | Note |
|---|---|---|
| `Ψ_d(R,t)` | no current equivalent | the nuclear packet on the discrete state; `qscat.core.time_dependent` propagates the *full* 2-D wavefunction instead, with no `Q`/`P` split |
| `Ψ~_d(R,t)` (Eq. 2.17a) | no current equivalent | the energy-independent rescaled packet |
| `F(R,R',t)` | no current equivalent | time-domain memory kernel; the Fourier partner of `domcke-1991-physrep208-97.md`'s `F(R,R';E)` |
| `Γ_L(R)`, `Δ_L(R)` (Eq. 2.12) | `Gamma(R)` in `qscat.core.lcp`; `Δ_L` never formed | qscat's `V_d(R) = Re E_pole(R)` corresponds to `V_d + Δ_L` combined — same mapping as in the Domcke note |
| `V_L(R) = V_d(R) + Δ_L(R)` (p. 1036) | `V_d(R)` in `qscat.core.lcp` | this paper names the combined real curve explicitly, which the Domcke review does not |
| `E_res(R)` (Eq. 2.14) | `E_res(R)` from `qscat.ecs.find_resonance_pole` | same physical quantity; qscat obtains it as a complex ECS pole rather than by bisecting Eq. (2.14) |
| `S(t)` (Eq. 4.4) | no current equivalent | survival probability of the resonance state |
| `ν` (spatial degeneracy) | not carried in qscat | `ν = 1` implied by the `4π³/2E` prefactor |
| eV, Å, fs | hartree, bohr, a.u. of time | **units differ** — see the Equations preamble |

## Not used here

- **The numerical method.** The paper propagates on an equidistant space-time
  grid with a fourth-order predictor-corrector finite-difference scheme
  (p. 1034, Eq. 3.7-3.8), started with fourth-order Runge-Kutta, using FFT for
  the action of `T_N` (Kosloff & Kosloff) and Simpson's rule for the memory
  integrals, at `Δt = 4.04 × 10⁻³` fs (p. 1034). Reflections are suppressed by
  multiplying the last 25 grid values by a Gaussian absorber,
  `e^{−c(M−25−m)²}` with `c = 1.6 × 10⁻⁵` (p. 1034, Eq. 3.9). None of this is
  qscat's approach: qscat propagates with the order-N diagonal Padé stepper
  (`qscat.evolution.make_pade_stepper`, `vandijk-2007-pre75-036707.md`) on a
  FEM-DVR grid and absorbs with exterior complex scaling, not a Gaussian mask.
  Recorded so a future TD-NRM port does not assume the paper's discretization
  travels with its physics.
- The Lanczos-basis construction of the target vibrational states `χ̄_v(R)`
  (p. 1034, Eq. 3.3, `v_max = 30` from a 50-function Lanczos basis) — a
  discretized representation of the target dissociation continuum. qscat's
  FEM-DVR grid already provides such a discretization; the Lanczos route is not
  ported.
- The analytic level-shift expression `Δ(E,R)` in terms of the confluent
  hypergeometric function `₁F₁` (p. 1035, Eq. 3.10), which is available only
  because the width was parametrized as Eq. (2.21). This is exactly the class of
  "assume an energy dependence so the Hilbert transform is analytic" move that
  Houfek et al. 2008 avoid — see the "Not used here" section of
  `domcke-1991-physrep208-97.md`.
- Figs. 1-9 are read for their stated conclusions and the two numeric ratios
  quoted above (`≃14`, `≃23`, and the ×0.3 / ×0.1 plot rescalings); no curve or
  pixel data was extracted, and no cross-section values were digitized.
- The comparison against the experimental data of Ehrhardt et al. and of Schulz
  & Asundi shown in Figs. 5-6 (p. 1039) — noted as present, not used. The paper
  reports its nonlocal DA cross section reproduces the shape of the experimental
  H₂ curve but is too large by a factor of 2 (p. 1039); this repository makes no
  claims against experiment.
- Section II's derivation of the projection-operator formalism is deliberately
  brief — the paper defers to Domcke, Phys. Rep. **208**, 97 (1991) (its Ref.
  [9], p. 1043) and to Estrada & Domcke, Phys. Rev. A **40**, 1262 (1989) (its
  Ref. [33]) for the details. For the formalism, cite
  `domcke-1991-physrep208-97.md`, not this note.
