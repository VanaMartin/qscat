# Time-dependent nonlocal resonance model — design

**Date:** 2026-08-19
**Status:** design, approved for planning
**Predecessors:** `2026-08-17-nrm-ti-da-design.md` (NRM core + dissociative
attachment), `2026-08-18-nrm-ve-design.md` (vibrational excitation). Both are
implemented and merged; this spec extends the same model into the time domain.

## 1. Goal

Compute `σ_DA(E)` and `σ_{v→v'}(E)` for N₂, F₂ and NO in the nonlocal resonance
model by **propagating a nuclear wave packet**, and gate the result against the
time-independent implementation `qscat.core.nrm` already ships.

The two routes solve the *same* equations with the *same* ingredients on the
*same* grids. Their agreement is therefore an exact identity, not a physical
approximation — which makes it the sharpest available test of the propagation,
and cleanly separates propagation error from model error. The comparison
against the exact 2-D solver rides on top of that, unchanged from the TI work.

## 2. Reference basis

| Source | Note | Supplies |
|---|---|---|
| Gertitschke & Domcke, Phys. Rev. A **47**, 1031 (1993) | `reference/literature/gertitschke-1993-pra47-1031.md` | The time-dependent nonlocal equation of motion Eq. (2.1)–(2.5), the amplitude transforms Eq. (2.6)/(2.8)/(2.9), the cross-section prefactors Eq. (2.7)/(2.10), the LCP limit in the time domain Eq. (2.11)–(2.15), and the packet diagnostics Eq. (4.4)–(4.6) |
| Houfek, Rescigno & McCurdy, Phys. Rev. A **77**, 012710 (2008) | `reference/literature/houfek-2008-pra77-012710.md` | The model itself: Eq. (52) nuclear equation, Eq. (60)–(61) `F(E)`, Eq. (31)/(37)/(38) the VE T-matrix split, Eq. (54) σ_DA — all already implemented |
| van Dijk & Toyama, Phys. Rev. E **75**, 036707 (2007) | `reference/literature/vandijk-2007-pre75-036707.md` | The diagonal-Padé propagator, already `qscat.evolution.make_pade_stepper` |

The paper works in **eV / Å / fs**. Everything here is atomic units
(`qscat.units`); no constant is carried across from its parameter tables, which
belong to a different (synthetic and H₂) model than the ones qModeling
implements.

**`reference/eMoScat`'s `module_NRM.cpp` specifies nothing.** The NRM was never
delivered as a working capability there. Its `MultiStep` path does propagate a
multi-component nuclear vector under one Padé operator — structurally the same
idea as §3 — and that is corroboration of shape at most. No code in this
sub-project may be justified by "eMoScat does it this way", and no test may
assert agreement with it. Note in particular that the `TimeIndependentSolution`
in the same file has its `V_dn` multiplications commented out.

## 3. The memory integral, and why it is not implemented literally

Eq. (2.1) is a non-Markovian equation of motion:

```
i ∂_t Ψ_d(R,t) = [T_N + V_d(R)] Ψ_d(R,t)
                 + (1/i) ∫_0^t dt' ∫ dR' F(R,R',t−t') Ψ_d(R',t')
```

Evaluated verbatim this is `O(N_t²)` in time and needs the kernel at every lag.
It does not have to be. The kernel this repository actually builds
(`nonlocal_potential.nonlocal_operator`, PRA 77 Eq. 60–61) is a sum of
resolvents:

```
F(E) = Σ_n diag(V_dn) · (E − H_n)^{-1} · diag(V_dn),
H_n  = T_R + V_0(R) + E_n(R)
```

Introduce one auxiliary nuclear packet `φ_n` per projected electronic state.
The pair `(Ψ_d, {φ_n})` then obeys a **time-local** coupled system,

```
i ∂_t Ψ_d = (T_N + V_d) Ψ_d + Σ_n V_dn φ_n
i ∂_t φ_n = H_n φ_n + V_dn Ψ_d
```

i.e. propagation under a single block Hamiltonian on `(1 + n_states)·N_R`
unknowns:

```
        ⎡ T_N + diag(V_d)   diag(V_d1)          diag(V_d2)        … ⎤
H_ext = ⎢ diag(V_d1)        T_N + diag(V_0+E_1)  0                … ⎥
        ⎣ diag(V_d2)        0                    T_N + diag(V_0+E_2) ⎦
```

**This is not an approximation of Eq. (2.1); it is Eq. (2.1) resummed.**
Eliminating `φ_n` from the time-independent version, `φ_n = (E−H_n)^{-1} V_dn Ψ_d`,
returns PRA 77 Eq. (52) with exactly the `F(E)` above. That elimination is the
first test written (Task 1), against `nonlocal_operator` + `solve_nuclear`.

Structural properties that matter:

- **Sparse.** Every off-diagonal block is *diagonal*; the arms are mutually
  uncoupled. The matrix is arrow-shaped.
- **Complex symmetric**, as everything under ECS is — `SparseLU` detects and
  exploits it per Padé factor.
- **Energy-independent.** One factorization serves every incident energy. The
  TI route must rebuild and invert `F(E)` at each energy; this is where the
  time-dependent route can be cheaper for a dense sweep.

## 4. From the propagated packet back to the existing extraction

Transform `Ψ̂(E) = ∫_0^∞ dt e^{i(E+i0)t} Ψ(t)`. Integrating `i∂_tΨ = H_ext Ψ`
by parts, with the `t→∞` boundary term killed by ECS absorption,

```
(E − H_ext) Ψ̂(E) = i Ψ(0)        ⟹    Ψ̂(E) = i (E − H_ext)^{-1} Ψ(0)
```

With the launch state `Ψ(0) = (ξ, 0, 0, …)` and `ξ = V_dk_i(R) χ_vi(R)` — Eq.
(2.5), and **byte-identical to the right-hand side `_psi_d_for_energy` already
builds for Eq. (52)** — the `d`-block gives

```
Ψ_d^TI(R;E) = −i ∫_0^∞ dt e^{iEt} Ψ_d(R,t)
```

The `−i` is derived here, and independently confirmed by Eq. (2.6), whose
`(1/i)` prefactor makes `T_{v'v}(E) = ⟨χ_v' V_dk_f | Ψ_d^TI⟩` — precisely
`vibrational_excitation.t_resonant`.

The consequence for the implementation is the point of the whole design: the
time-dependent route produces **the same nuclear vector the time-independent
route produces**, and every downstream step is then existing, tested code.

- **DA:** `dissociation.da_sigma_from_psi(nuclear_grid, mu, Ψ_d, e_total, eps_e, e_kin)`,
  unchanged.
- **VE:** `vibrational_excitation.t_resonant(chi_f, v_dk_f, Ψ_d)` plus
  `t_background(...)`, unchanged. `T^bg` is energy-domain and static — it
  contains no `Ψ_d` — so the time-dependent route replaces `T^res` only.

No cross-section normalization is re-derived, and the conjugation conventions
settled in the TI work (PRA 77 Eq. 34/35 collapsing `V^{−*}_{dk}` to a
non-conjugated `V^+_{dk}`; Eq. 37's bra carrying `φ⁺` at the final channel
energy) are inherited rather than re-litigated.

## 5. Incident energies — one propagation, not one per energy

`ξ = V_dk_i(R) χ_vi(R)` depends on the incident energy through `V_dk⁺`, so
exact algebra says one launch state per energy. PRA 47 Eq. (2.17) removes that
by rescaling with `[Γ(E_i)/2π]^{−1/2}`, leaving `Ψ̃_d(R,0) = g(R)χ_vi(R)` and one
propagation for all energies — resting on the separability
`Γ(E,R) = Γ(E) g(R)²` of Eq. (2.16), which the Houfek models do not satisfy
exactly.

They satisfy it *numerically*, which is what actually matters. Because
`α_c` is a constant in `V_int(r,R) = −λ(R)e^{−α_c r²}`, the whole `R`-dependence
of the electronic problem enters through the single scalar `λ(R)` —
`H_el(R) = H_0 + λ(R)W` — so `V_dk⁺(R;E) = g(λ(R), E)`, a smooth function of two
scalars. Measured on the production windows (2026-08-19, choice B, singular
values of `M[R,E_j] = ξ(R;E_j)` normalized to the first):

| window | σ₂/σ₁ | σ₃/σ₁ | σ₄/σ₁ | shape overlap across the window |
|---|---|---|---|---|
| F₂ DA, 0.010–0.050 Ha, 9 energies | 5.7e-3 | 2.4e-4 | 5.3e-7 | 0.999605 |
| N₂ VE, 0.060–0.160 Ha, 9 energies | 9.8e-4 | 1.2e-6 | 1.1e-8 | 0.999994 |

**The scheme.** SVD the launch matrix, propagate the left singular vectors as
`r` columns, and reconstruct by linearity of the resolvent:

```
M = U Σ Vᴴ,        ξ(:,j) = Σ_m (σ_m V*_{mj}) u_m
Ψ_d(:,j) = Σ_m (σ_m V*_{mj}) · Ψ̂_m(:,j),   Ψ̂_m(:,E) = −i ∫dt e^{iEt} U_m(R,t)
```

This is exact given the truncation: `H_ext` is energy-independent, so one
propagated `u_m` is transformed at **every** `E_j` and the superposition
commutes with the resolvent. Three columns replace forty-one at a residual of
5e-7 (F₂) / 1e-8 (N₂), and **`r = 1` reproduces Eq. (2.17) exactly** — the
paper's scheme is the rank-1 special case, not something discarded.

The difference from Eq. (2.16) is worth keeping straight: that is a claim of
*physical* separability, exact at all energies; this is *numerical* low rank
over a finite window, with a truncation error that is measured and reported
rather than assumed. It also degrades gracefully — a wider window or a harder
molecule raises `r` instead of invalidating the method.

Two limits on the numbers above. The rank is **window-dependent**, so a sweep
wider than the production window re-measures it; and both rows are **choice B**
(`AsymptoticDiscreteState`, `R`-independent `φ_d`). `PhysicalDiscreteState`
carries an `R`-dependent `φ_d` and is not assumed to be as low-rank — the
implementation measures it and picks `r` from a tolerance, never from a
hardcoded constant.

The half-Fourier transform is accumulated **on the fly**: each propagated
column carries an `(N_R, n_E)` accumulator, so no snapshot history is stored.
Quadrature reuses the weighting `time_dependent._quadrature_weights` already
applies to the Tannor-Weeks correlation transform.

## 6. The Markovian (LCP) limit

Eq. (2.11) states that LCP *is* the Markovian limit of Eq. (2.1): the kernel
collapses to `i[Δ_L(R) − (i/2)Γ_L(R)] δ(R−R') δ(t)`. Operationally that is the
same propagation with the arms removed — the `d`-block alone, with `F` replaced
by a local complex term. It costs one variant constructor.

**One convention must be checked, not assumed.** `qscat.core.nrm`'s
`v_d_discrete` is PRA 77 Eq. (20); `qscat.core.lcp`'s `Vd` is `E_res(R)`; and
Gertitschke Eq. (2.14) relates them through `Δ`. §4 of
`docs/physics/nonlocal-resonance-model.md` already records that these
"almost coincide" and differ measurably (0.0053 Ha on F₂, 0.0229 Ha on NO).
The spec does **not** decide which enters Eq. (2.15). The implementation
determines it by requiring the Markovian propagation to reproduce
`qscat.core.lcp.lcp_da_cross_section`, and the plan carries that as a gated
test with the discrepancy reported either way.

## 7. Diagnostics

Eq. (4.4)–(4.6), on the `d`-block, recorded per step:

```
S(t)   = ∫ dR |Ψ_d(R,t)|²                      survival
⟨R⟩_t  = ∫ R dR |Ψ_d|² / S(t)                   centroid
⟨P⟩_t  = ∫ dR Ψ_d* (−i ∂_R) Ψ_d / S(t)          mean momentum
```

These are what the time-independent route structurally cannot produce, and they
are the instrument Gertitschke uses to explain *why* LCP fails: the exact
nonlocal packet splits temporarily, delaying dissociation, while the local one
cannot. Two cautions are recorded in the physics note rather than discovered
later — Ehrenfest's theorem does **not** hold for this open-system dynamics
(`⟨P⟩_t` is not `μ d⟨R⟩_t/dt`; the paper states this at p. 1036), and oscillations from rapid packet broadening
look like N₂-style boomerang oscillations without being them.

`S(t)` does double duty as the convergence criterion of §9: the transform is
truncated, not converged, until the packet has left the real region.

## 8. API surface

Three new modules inside the existing `libs/qscat/qscat/core/nrm/` package:

**`extended.py`**

```python
extended_hamiltonian(ing, nuclear_grid, model, *, n_states=None) -> sp.csr_matrix
lcp_limit_hamiltonian(nuclear_grid, model, v_res) -> sp.csr_matrix

@dataclass(frozen=True)
class LaunchBasis:
    vectors: NDArray     # (N_ext, r) left singular vectors; arm blocks zero
    coeffs: NDArray      # (r, n_E) -- sigma_m * conj(V[m, j])
    energies: NDArray    # (n_E,) incident kinetic energies
    e_total: NDArray     # (n_E,) transform frequencies
    residual: float      # sigma_{r+1}/sigma_1 -- the truncation error, reported

initial_packet(nuclear_grid, elec_grid, model, phi_d, ing, eps, chi, v_init,
               energies, *, n_states=None, rank_tol=1e-6) -> LaunchBasis
```

`r` comes from `rank_tol`, never from a constant. `r == 1` is PRA 47
Eq. (2.17); `residual` is what the propagation reports so a reader can see
the economy's cost.

**`propagation.py`**

```python
@dataclass(frozen=True)
class TdNrmResult:
    psi_d: NDArray          # (N_R, n_E) -- the reconstructed Psi_d(R;E)
    time: NDArray           # (n_steps+1,)
    survival: NDArray       # (n_steps+1, n_E)
    centroid: NDArray
    momentum: NDArray
    residual: NDArray       # (n_E,) survival at t_max
    rank: int               # columns actually propagated

propagate_nrm(h_ext, launch, nuclear_grid, *, dt, n_steps, order=3) -> TdNrmResult
```

`r` columns are stepped; the per-energy packet used for the diagnostics and
the transform is `Sum_m coeffs[m,j] U_m(t)`, reconstructed per step.

**`td_cross_section.py`** (named to avoid colliding with `qscat.core.time_dependent`)
— signature parity with the TI entry points, so the
differential comparison is a one-identifier substitution:

```python
td_nrm_da_cross_section(nuclear_grid, elec_grid, model, phi_d, eps, chi,
                        v_init, E, *, ingredients=None, n_states=None,
                        dt, n_steps, order=3, rank_tol=1e-6,
                        markovian=False) -> NDArray
td_nrm_ve_cross_section(nuclear_grid, elec_grid, model, phi_d, eps, chi,
                        v_init, vprimes, E, *, ingredients=None, n_states=None,
                        dt, n_steps, order=3, rank_tol=1e-6,
                        include_background=True, markovian=False) -> NDArray
```

`markovian=True` selects the Eq. (2.15) propagation of §6.

## 9. Validation

Four gates, in strength order.

1. **The elimination identity (Task 1).** `(E − H_ext)^{-1}` restricted to the
   `d`-block equals `(E − T_N − V_d − F(E))^{-1}`, at machine precision, on a
   small deck. Tests the block construction with no propagation involved.
2. **The transform identity (Task 3) — the load-bearing gate.**
   `Ψ_d^TD(R;E)` against `Ψ_d^TI(R;E)` **vector to vector** on a small F₂ deck.
   Same ingredients, same grids, exact identity; a discrepancy is propagation
   or transform error and nothing else.
3. **Cross sections against the TI route**, F₂ DA and N₂/F₂ VE, at the anchors
   the TI work already gates. Tolerance is set by the convergence study of §10,
   not chosen in advance.
4. **Cross sections against the exact 2-D solver**, inherited: the TI route
   reproduces it to 0.06–0.33% (F₂ DA, choice B) and sub-0.7% (N₂/F₂ VE), so
   the time-dependent route must land in the same band.

The Markovian variant is gated separately against `lcp_da_cross_section` (§6).

**NO is a probe, not a promised fix.** The open collapse
(`docs/physics/nonlocal-resonance-model.md` §7.2) lives in the NRM-vs-exact
comparison. Since TD-NRM must reproduce TI-NRM by gate 2, it will show the
*shape* of the failure in time — packet splitting, survival, where the
dissociating flux goes — without by itself explaining it. That is the best
instrument currently available on the question and it is worth running; it is
not advertised as a resolution, and the spec commits to reporting whatever it
shows including "no new information".

## 10. Convergence knobs and defaults

Swept and documented per molecule, not guessed: `dt`, `n_steps` (through
`t_max`), Padé `order`, and `n_states`. `n_states` has a TI precedent
(`_N_STATES = 100`) and should transfer, which is itself a check. `t_max` is
determined by the residual criterion of §7 — the packet must be absorbed before
the transform is truncated — and gets a runtime warning, not a silent default.

Cost risk: the arrow LU may fill in through the head at `n_states ≈ 100`. If
measured to be the bottleneck, the fix is a Schur complement onto the `d`-block
— each arm is banded and solves independently — but it is implemented **only if
the measurement demands it**, not preemptively.

## 11. qscat-run

A new method, `nrm-td`, beside `ti` / `td` / `lcp` / `nrm`. This keeps artifact
keys flat (`nrm-td:da:v0`) and, decisively, lets one config request `[ti, nrm,
nrm-td]` so the gate of §9.3 is expressible as a single run on a single
discretisation. It reuses `presets.resolve_nrm_grids` unchanged and needs the
same `td:` evolution block the `td` method already parses.

## 12. Out of scope

- Dissociative recombination (H₂⁺) through the nonlocal model.
- The literal `O(N_t²)` memory convolution. §3's elimination is exact and gate 2
  tests the propagation end to end; a verbatim convolution would test a weaker
  claim at the cost of a task.
- Gertitschke's synthetic `d`-wave and H₂ models. Their parameters are recorded
  in the reference note and deliberately not implemented.
- Rust kernels. The propagation is a candidate hot path, but the lifecycle
  requires validation first.
