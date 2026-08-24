# The time-dependent nonlocal resonance model

The time-domain sibling of `docs/physics/nonlocal-resonance-model.md`. Same
model, same ingredients, same grids — a second, independent route to the same
`Ψ_d`, and a set of observables the first route cannot produce.

Method: Gertitschke & Domcke, Phys. Rev. A **47**, 1031 (1993)
(`reference/literature/gertitschke-1993-pra47-1031.md`), applied to the model of
Houfek, Rescigno & McCurdy, Phys. Rev. A **77**, 012710 (2008)
(`reference/literature/houfek-2008-pra77-012710.md`). Propagator: the order-3
diagonal Padé of van Dijk & Toyama, Phys. Rev. E **75**, 036707 (2007), already
shipped as `qscat.evolution.make_pade_stepper`.

## 1. Why, given that it is slower

It is slower — measurably. One F₂ propagation costs ~0.3–0.8 s per step over
~6000 steps; the resolvent answers the same energies in a fraction of that. The
time-dependent route is not justified on cost and this note does not claim it is.

It is justified because a resolvent returns `Ψ_d(R;E)` and says nothing about how
the collision complex reached it. A propagation returns `Ψ_d(R,t)`, and with it
the survival, the centroid and the mean momentum — the quantities that turn a
cross-section feature into a mechanism. Váňa & Houfek, Phys. Rev. A **95**,
022714 (2017) is built on exactly that: the correspondence between individual
peaks and quasibound states, the several time-separated contributions behind an
asymmetric peak, and the formation-time argument for NO's long-lived structure
(pp. 022714-13, -15). Gertitschke & Domcke's own headline result is of the same
kind — their explanation of why the LCP fails for H₂⁻ is a *temporary packet
splitting between ≈2 and ≈5 fs*, invisible in any energy-domain calculation
(p. 1041).

PRA 95 proposes this work in its closing paragraph: similar time-dependent
calculations "could be performed within the LCP approximation or the nonlocal
resonance model and thus interpret the results in the same way", offered as a
future direction and not executed there (p. 022714-16).

## 2. The formulation

### 2.1 The memory integral, and why it is not solved as one

PRA 47 Eq. (2.1) is non-Markovian:

```
i ∂_t Ψ_d(R,t) = [T_N + V_d(R)] Ψ_d(R,t)
                 + (1/i) ∫_0^t dt' ∫ dR' F(R,R',t−t') Ψ_d(R',t')
```

Evaluated literally this is `O(N_t²)` and needs the kernel at every lag. It does
not have to be. The kernel this repository builds
(`nrm.nonlocal_potential.nonlocal_operator`, PRA 77 Eq. 60–61) is already a sum
of resolvents:

```
F(E) = Σ_n diag(V_dn) · (E − H_n)^{-1} · diag(V_dn)
H_n  = T_R + V_0(R) + E_n(R)
```

### 2.2 The resummation

Introduce one auxiliary nuclear packet `φ_n` per projected electronic state. The
pair `(Ψ_d, {φ_n})` obeys a **time-local** system:

```
i ∂_t Ψ_d = (T_N + V_d) Ψ_d + Σ_n V_dn φ_n
i ∂_t φ_n = H_n φ_n + V_dn Ψ_d
```

i.e. propagation under one block Hamiltonian on `(1 + n)·N_R` unknowns:

```
        ⎡ T_N + diag(V_d)   diag(V_d1)           diag(V_d2)          … ⎤
H_ext = ⎢ diag(V_d1)        T_N + diag(V_0+E_1)   0                  … ⎥
        ⎣ diag(V_d2)        0                     T_N + diag(V_0+E_2)  ⎦
```

sparse (every off-diagonal block is diagonal), complex **symmetric** (ECS, never
Hermitian), and **energy-independent**.

**This is Eq. (2.1) resummed, not approximated.** Eliminating the auxiliary
blocks from the time-independent version, `φ_n = (E−H_n)^{-1} V_dn Ψ_d`, returns
PRA 77 Eq. (52) with exactly the `F(E)` above. Verified numerically: the
`d`-block of `(E − H_ext)^{-1}` matches `(E − T_N − V_d − F(E))^{-1}` to
**4.4e-14** relative (F₂, `N_R = 584`, `cond ≈ 1.1e6`).

### 2.3 Back to the frequency domain

For `i ∂_tΨ = H_ext Ψ` with `Ψ(0)` given, integrating `∫_0^∞ dt e^{i(E+i0)t}` by
parts and killing the `t→∞` boundary term by ECS absorption,

```
(E − H_ext) Ψ̂(E) = i Ψ(0)      ⟹      Ψ̂(E) = i (E − H_ext)^{-1} Ψ(0)
```

With `Ψ(0) = (ξ, 0, 0, …)` and `ξ = V_dk_i(R) χ_vi(R)` — PRA 47 Eq. (2.5), and
byte-identical to the right-hand side the time-independent route already builds
for Eq. (52) — the `d`-block gives

```
Ψ_d^TI(R;E) = −i ∫_0^∞ dt e^{iEt} Ψ_d(R,t)
```

The `−i` is derived here and confirmed independently by PRA 47 Eq. (2.6), whose
`(1/i)` makes `T_{v'v}(E) = ⟨χ_v' V_dk_f | Ψ_d^TI⟩` — precisely the shipped
`vibrational_excitation.t_resonant`.

**The consequence is the design.** The time-dependent route produces *the same
nuclear vector* the time-independent route produces, so every downstream step is
existing, tested code: `dissociation.da_sigma_from_psi` for DA, `t_resonant` +
`t_background` for VE. No cross-section normalization is re-derived, and the
conjugation conventions settled for the TI route (PRA 77 Eq. 34/35 collapsing
`V^{−*}_{dk}` to a non-conjugated `V^+_{dk}`; Eq. 37's bra carrying `φ⁺` at the
final channel energy) are inherited rather than re-litigated.

### 2.4 One propagation for every energy

`ξ` depends on the incident energy through `V_dk_i`, so exact algebra says one
launch state per energy. PRA 47 Eq. (2.17) removes that by rescaling with
`[Γ(E_i)/2π]^{−1/2}`, under the separability `Γ(E,R) = Γ(E) g(R)²` of Eq. (2.16)
— which the Houfek models do not satisfy exactly.

They satisfy it *numerically*. Because `α_c` is a constant in
`V_int(r,R) = −λ(R)e^{−α_c r²}`, the whole `R`-dependence of the electronic
problem enters through the single scalar `λ(R)` — `H_el(R) = H_0 + λ(R)W` — so
`V_dk⁺(R;E) = g(λ(R), E)`, a smooth function of two scalars. The launch matrix
`M[R,E_j] = ξ(R;E_j)` is therefore numerically low rank, and the scheme is to
SVD it, propagate the left singular vectors, and reconstruct every energy by
linearity of the resolvent:

```
M = U Σ Vᴴ,        Ψ_d(:,j) = Σ_m (σ_m V*_{mj}) · Ψ̂_m(:,j)
```

exact given the truncation, because `H_ext` is energy-independent so the
superposition commutes with the propagator. **`r = 1` is Eq. (2.17) exactly**;
higher ranks are its controlled generalization, with a measured rather than
assumed error. See §4.3 for the measured ranks — and note the rank is a property
of the discrete-state choice, not of the model.

## 3. What is verified

| Claim | Measured | Where |
|---|---|---|
| Elimination identity ⇒ `nonlocal_operator`'s `F(E)` | **4.4e-14** relative | `test_nrm_extended.py` |
| Half-Fourier transform vs closed-form finite-`T` | **1.8e-13**, converging as `dt⁶` | `test_nrm_propagation.py` |
| `⟨P⟩_t` vs analytic Gaussian (`p = 1.7`) | **1.6999999999999842** | `test_nrm_propagation.py` |
| **`Ψ_d^TD(R;E)` vs `Ψ_d^TI(R;E)`, vector to vector** (N₂) | **1.73e-04** | `test_nrm_propagation.py` |
| **`σ_TD/σ_TI`, F₂ DA, E = 0.02/0.03/0.05 Ha** | **1.0097 / 1.0138 / 1.0102** | `test_nrm_td_cross_section.py` |

![TD vs TI cross section](figures/f2-da-nrm-td-vs-ti.png)

![Vector-to-vector agreement](figures/n2-nrm-td-vs-ti-vector.png)

The N₂ vector error decomposes exactly. Fitting the measured `rel(T)`:

```
truncation(T) = 0.40 · sqrt(S(T)/S(0))        propagation(dt=1) = 1.43e-4
rel = sqrt(truncation² + propagation²)         (in quadrature, not linearly)
```

reproduces the convergence table to within 3% **from `T = 2000` onward**
(measured/budget 1.03, 1.03, 1.02, 1.00), which is what makes "this is
convergence, not a bug" a demonstration rather than an argument.

It does **not** hold at short times — 2.2× at `T = 500`, 1.4× at `T = 1000`.
The budget is asymptotic: `truncation ∝ √(S(T)/S(0))` assumes the tail
`−i∫_T^∞ e^{iEt}Ψ(t)dt` is dominated by a decaying remainder, and while the
packet is still inside the box that premise does not apply. The shipped gate
sits at `T = 4000`, well inside the regime where the budget is valid.

![Convergence](figures/nrm-td-convergence.png)

## 4. What the time domain shows that the resolvent cannot

### 4.1 The packet leaves — and on a coarse grid it does not

![Packet dynamics](figures/f2-da-nrm-td-packet.png)

On F₂'s fine production nuclear deck (65 points/bohr over R = 3–10.7, largest
node spacing 0.023 bohr) the DA packet dissociates: `⟨R⟩_t` climbs monotonically
**2.66 → 9.33 bohr**, `⟨P⟩_t` rises toward `K_R ≈ 58`, and `S(t)` falls
**0.94 → 0.007** as the wave crosses the ECS absorber at R = 10.7.

`⟨R⟩_t` then **turns over** near `t ≈ 4000` and falls back to ~6.4 bohr while
`⟨P⟩_t` decays toward zero. That is not the packet returning: it is the fast
dissociating flux being absorbed at the boundary, leaving the bound remnant of
§4.2 to dominate what is left inside the box. The turnover and the `S(t)`
plateau are the same fact seen twice — which is also why `⟨P⟩_t` peaks near 48
rather than reaching `K_R = 58.6`, since it averages the outgoing wave together
with a slow bound population.

On a coarse deck it does none of this — the centroid creeps to 2.77 bohr and
oscillates, survival plateaus at 0.673. The outgoing wavelength is
**0.107 bohr** against ~15 points/bohr, so the flux cannot propagate to the
absorber and rattles in place. It *looks* bound and is not. Only the time-domain
diagnostics distinguish the two; a resolvent returns a plausible number in both
cases.

### 4.2 `S(t)` plateaus on F₂, and σ converges anyway

F₂'s `V_d(R)` has a well — minimum **−0.149264 Ha at R = 3.363** against the
F+F⁻ asymptote **−0.126931**, a depth of **0.0223 Ha** — supporting **≥24**
near-real modes with `|Im E| = 1.5e-7 … 7.7e-6`. The launch populates them with
**5.08e-3** of its real-region norm, so `S(t)` flattens at 0.006–0.009 from
`T ≈ 12000` and **no absolute survival floor is reachable on this molecule**;
removing their ≈2e-2 contribution to the transform would need `T ≳ 1e7`.

`σ_DA` converges regardless, because it reads the wavefunction *value* at the
outermost real node, where those well-localized modes have almost no amplitude.
**The observable converges while the norm does not**, so convergence here is
defined as `σ_DA` stationary in `T` — not as `S(T)/S(0)` below a threshold. A
convergence check written the obvious way would warn permanently on every F₂ run.

### 4.3 The launch basis, and what it says about the two discrete states

![Launch rank](figures/nrm-td-launch-rank.png)

| window | choice | σ₂/σ₁ | σ₃/σ₁ | σ₄/σ₁ | rank @1e-6 |
|---|---|---|---|---|---|
| F₂ DA, 0.010–0.050 Ha | B | 5.7e-3 | 2.4e-4 | 5.3e-7 | **3** |
| F₂ DA, 0.010–0.050 Ha | A | 3.3e-1 | 9.0e-2 | 1.9e-3 | **7** |
| N₂ VE, 0.060–0.160 Ha | B | 9.8e-4 | 1.2e-6 | 1.1e-8 | **3** |
| N₂ VE, 0.060–0.160 Ha | A | 1.5e-1 | 6.2e-3 | 3.3e-4 | **5** |

PRA 47's "single energy-independent wave packet" is essentially exact for choice
B and materially approximate for choice A. The reason is structural: an
R-independent `φ_d` lets the launch state's R-dependence enter only through
`λ(R)`, while choice A's `φ_d(·;R)` injects genuine two-variable structure. This
is an independent symptom of the same R-dependence that degrades choice A in the
TI campaign (F₂ DA ratio 0.29–0.90; VE 0.565–1.140) and that PRA 77's Sec. VI A
predicts.

## 5. The arm set cannot be truncated

Truncating the sum over projected electronic states makes `H_ext`
**non-dissipative**: modes acquire `Im(E) > 0` and the propagation diverges
exponentially in `t`. Dense `eig` gives `max Im(E) = +2.787e-3` at `n_states=3`
on the N₂ gate deck and `+5.84e-4` on the F₂ production fixture. That the
divergence is identical to five digits across `dt = 4/2/1/0.5/0.25` and Padé
order 3 and 4 proves it is the operator, not the propagator.

**The asymmetry with the time-independent route is the point.** Truncating `F(E)`
perturbs a **resolvent** by a bounded amount — which is why the TI campaign runs
happily at `_N_STATES = 100` — whereas a single eigenvalue crossing into the
upper half-plane grows without bound in a **propagator**. The same approximation
is benign in one route and fatal in the other, so the TI precedent does not
transfer and an `n_states` convergence study is unavailable here.

The effect is **non-monotonic** — on one small F₂ fixture `max Im(E)` runs
+1.60e-4 / +2.02e-4 / +1.96e-4 / +9.03e-6 at `n_states` = 1/2/3/8, then ≤0 from
12 onward, and +2.4e-12 for the complete set. So the rule is not "more arms is
safer": *any* truncation may leave growing modes, and establishing that a
particular one does not requires a dense spectral check. `n_states=None` is the
only value needing no check. `propagate_nrm` warns at runtime when a packet grows.

**Spectral claims here require a dense eigensolver.** These matrices are strongly
non-normal and ARPACK `eigs(which="LI")` under-reports the largest imaginary
part — `+4.639e-4` against dense `eig`'s `+2.787e-3` on the same
716-dimensional matrix, six times low. Sparse results are lower bounds.

**The complete set's own `max Im(E)` was an under-resolved ECS tail, and that is
now settled.** On three decks it sits at `+1.4e-12 … +5.0e-13` — marginal
stability, which is correct, since genuinely bound nuclear states do not decay
(that is what §4.2 is). A small 68-arm F₂ fixture returned `+6.98e-5` instead,
five orders larger, while the propagation launched from `ξ` on that same fixture
**decayed**. Neither of the two explanations first offered (a real but unpopulated
mode; dense-eigensolver backward error) was right. Measured 2026-08-24 by
refining that fixture's **nuclear ECS tail** and nothing else:

| nuclear tail elements (6 → 20 bohr) | complete set `max Im(E)` | `Re E` there | `n_states=3` `max Im(E)` | `Re E` there |
|---|---|---|---|---|
| 2 | +6.98e-5 | +958 Ha | +3.305e-4 | −0.0754 Ha |
| 4 | +4.15e-8 | +666 Ha | +3.305e-4 | −0.0754 Ha |
| 6 | +4.71e-8 | +4806 Ha | +3.305e-4 | −0.0754 Ha |

Two things separate at once. The complete set's large value **moves four orders
under refinement and converges**, and the eigenvalue carrying it sits at
`Re E` of hundreds to thousands of hartree — a high-energy mode of the rotated
continuum, at an energy the packet never occupies. It was a discretisation
artifact of a tail spanning 14 bohr in two elements. The truncated set's value
**does not move at all**, on any tail, and sits at `Re E = −0.0754 Ha`, inside
the packet's own energy range. Artifact versus real growing mode is exactly this
contrast, and it is why figure `nrm-td-truncation-diverges.png` classifies by a
threshold rather than by the sign of `Im E`: the line it draws is where a mode
would grow the packet by 1% over the propagation (`ln 1.01 / 2T`), which the
truncations clear by four orders and the complete set misses by four.

Even so, **a spectral argument is not what makes a propagation safe here — the
runtime survival guard is.** `propagate_nrm`'s check that the packet does not
grow tests the thing that actually matters (is a growing mode populated?) rather
than a sufficient condition on the operator, and it needs no dense `eig` at
production size. Truncation is ruled out on both counts: the growing modes are
real under refinement AND demonstrably populated — on this fixture `S(t)/S(0)`
reaches 10 by `T = 4000` at `n_states = 3` while the complete set falls to ~0.6.

![Truncation diverges](figures/nrm-td-truncation-diverges.png)


## 6. Limits — what this does not establish

- **The TD/TI agreement validates the propagation, not the model.** Both routes
  run on the same grids with the same ingredients, so shared discretisation error
  cancels. Absolute scale is anchored elsewhere (`validation/n2/exact2d.py`
  against Houfek at `GATED_RTOL=1e-3`).
- **DA is F₂ only.** NO's DA channel is a known open failure for the TI route
  (`docs/physics/nonlocal-resonance-model.md` §7.2), so a TD run there would be
  measuring an already-broken oracle.
- **The DA fixture uses a reduced 55-point electronic grid.** The TD/TI ratio is
  insensitive to it (3–4 significant figures), but the fixture's absolute `σ_TI`
  is **not** the converged F₂ cross section.
- **No production-electronic-deck data point exists** on any molecule yet.
- **Vibrational excitation is not yet implemented** in the time domain, nor the
  Markovian (LCP) limit — so the comparison Gertitschke & Domcke's paper is built
  around is not yet reproducible here.

## 7. Literature

- **PRA 47, 1031 (1993)** — `reference/literature/gertitschke-1993-pra47-1031.md`.
  The time-dependent nonlocal equation of motion Eq. (2.1)–(2.5), the amplitude
  transforms Eq. (2.6)/(2.8)/(2.9), the LCP limit in the time domain
  Eq. (2.11)–(2.15), and the packet diagnostics Eq. (4.4)–(4.6).
- **PRA 77, 012710 (2008)** — `reference/literature/houfek-2008-pra77-012710.md`.
  The model: Eq. (52), Eq. (60)–(61), Eq. (31)/(37)/(38), Eq. (54).
- **PRA 95, 022714 (2017)** — `reference/literature/vana-houfek-2017-pra95-022714.md`.
  What time-dependent formulations of these models are *for*, and the closing
  proposal this sub-project executes (p. 022714-16).
- **PRE 75, 036707 (2007)** — `reference/literature/vandijk-2007-pre75-036707.md`.
  The diagonal-Padé propagator.
