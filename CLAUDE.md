# CLAUDE.md — qModeling operating manual

## Purpose

qModeling is a Python-first quantum-mechanics (QM) research monorepo, home to
**QSCAT**. It resuscitates and extends prior C++/CUDA work (`eMoScat`,
`libXcuda`) as a maintainable, CPU-first Python codebase with Rust for the hot
paths that actually need it. GPU/CUDA and AWS deployment are recovered as
reference material but deliberately deferred — not part of this phase.

## The lifecycle

Every method or numerical capability moves through five stages:

1. **Design** — work out the math/algorithm (see `docs/physics/`), decide the
   API shape.
2. **Toy model** — a small, readable Python implementation under `projects/`.
3. **Validate** — check it against analytic benchmarks, conservation laws,
   convergence studies, or a reference implementation (`validation/`).
4. **Optimize in Rust** — once validated and proven a hot path, move it to a
   PyO3/maturin kernel under `native/`, keeping the Python version as the
   differential oracle.
5. **Promote to `qscat`** — validated, reusable code graduates into
   `libs/qscat`, the standard library other projects depend on.

**Invariants that hold at every stage:** the code is CPU-runnable on a laptop
(no GPU required), and it is containerizable (builds and runs under
`docker/`). See the `qm-method-lifecycle` skill for the enforced workflow.

## Repo map

```
libs/       qscat — the standard library: validated, reusable QM code
            (units, linalg, dvr, ecs, special, evolution, core, model
            submodules)
            - qscat.linalg: dimension-general sparse linear algebra --
              `kron_sum` (Kronecker sum over arbitrary D), `SparseLU` (cached
              factorization with fill-in/memory diagnostics), and `c_product`
              (the bilinear, non-conjugated ECS inner product) -- see
              docs/physics/nd-tensor-hamiltonian.md. `SparseLU` dispatches
              between SuperLU (`backend="scipy"`, numpy/scipy-only, the fallback
              AND the differential oracle) and a complex-symmetric MUMPS `SYM=2`
              backend (`backend="mumps"`, the optional `qscat[mumps]` extra):
              `backend="auto"` (default) picks MUMPS if present else SuperLU,
              and `set_default_backend`/`default_backend()` override what
              `"auto"` resolves to process-wide. On the ECS-complex-symmetric
              N₂ matrices MUMPS beats SuperLU by 72× in factor time / 9× in peak
              RSS at the 143k production deck (3.6 s / 0.8 GB vs 260 s / 7.4 GB)
              -- see docs/physics/mumps-sparse-backend.md. `SparseLU.refactor(
              A_new)` reuses the symbolic analysis (same sparsity pattern, e.g. a
              diagonal shift `E_tot·I − H` across an energy sweep): MUMPS
              `factor(reuse_analysis=True)` skips the SCOTCH ordering; scipy
              re-runs `splu` (correct, no reuse); a pattern guard raises on a
              structure mismatch. On the N₂ working grid a reuse sweep is ~80% /
              ~5× faster than fresh-per-energy (the analysis dominates the cheap
              numeric factor there; the fraction shrinks for larger decks) — see
              docs/physics/ti-energy-sweep-reuse.md. `ShiftInvertEigs` is the
              eigenvalue side of that same trick: the `k` eigenpairs nearest a
              complex shift, via shift-invert Arnoldi with `SparseLU` as the inner
              solve. Resonances are INTERIOR eigenvalues, which a plain Krylov
              iteration cannot reach, and `A − σ·I` keeps one sparsity pattern for
              every σ, so a sweep of shifts reuses the analysis exactly as an
              energy sweep does. Validated in 1-D against the dense
              `qscat.dvr.eigen`; the conventions that bite (the `A − σ·I` sign,
              nearest-shift ordering, Euclidean vs c-product normalization) and
              the measured working range are in
              docs/physics/shift-invert-eigensolver.md. Its first consumer is
              `qscat.core.exact_resonance_states` (below).
            - qscat.dvr: FEM-DVR-ECS radial grid (`FemDvrEcsGrid`), kinetic-
              energy assembly (`kinetic`), and diagonal-potential Hamiltonian
              + eigensolver helpers (`hamiltonian`, `eigen`) — see
              docs/physics/femdvr-ecs.md.
              Also `kinetic_sparse` and the N-dimensional tensor layer
              (`TensorGrid`, `kinetic_nd`, `potential_nd`, `hamiltonian_nd`):
              H = sum_d I x .. T_d .. x I + diag(V) for any D, sparse (CSR),
              validated on analytic box/oscillator benchmarks at D = 1, 2, 3.
            - qscat.ecs: exterior-complex-scaling coordinate map (`ecs_map`),
              the single source of the `z(x) = x` / `R0 + (x-R0)e^{i theta}`
              transform used by qscat.dvr's complex tail; also
              `find_resonance_pole(eigs_a, eigs_b, window)`, the general
              two-spectrum resonance-pole matcher (promoted from the N2
              resonance project) — see docs/physics/n2-resonance.md.
            - qscat.evolution: `make_cn_stepper(H, dt)`, a general Crank-
              Nicolson time propagator for `d/dt psi = -i H psi` (any complex,
              possibly non-Hermitian `H`) — promoted from the N2
              time-dependent cross-section project — see
              docs/physics/n2-td-cross-section.md. Also
              `make_sparse_cn_stepper(H, dt)`, the sparse sibling (factors
              once with `SparseLU`, matches the dense stepper to round-off) —
              promoted from the N2 2-D time-dependent cross-section project —
              see docs/physics/n2-2d-td-cross-section.md. Also
              `make_pade_stepper(H, dt, order)` (+ `pade_roots`), the order-N
              diagonal-Padé generalization of the sparse CN stepper (order 1 ==
              CN; `O(dt^(2N+1))` per step). Order-1 CN under-converges badly
              over a long propagation (~100% accumulated error at dt=0.5-1.0 vs
              `expm`); order-3 Padé is what makes the TD cross section converge
              to the TI oracle to ~1-2% — see docs/physics/n2-2d-td-cross-section.md.
            - qscat.core: the model-INDEPENDENT electron–diatomic VE-scattering
              engine, promoted from the N2 sub-projects. `driven.ve_cross_section`
              (exact TI driven Lippmann-Schwinger, `SparseLU.refactor` sweep,
              σ=π|S−δ|²/2E) and `time_dependent.td_ve_cross_section` (order-3
              Padé propagation + a `method`-selectable energy-extraction
              transform + elastic free-reference subtraction) — both take a
              `model`. `method="tw"` (default, Tannor-Weeks — a propagated
              Gaussian test packet) is joined by two siblings sharing the SAME
              propagate-once `Extractor` protocol (`td_extractors.py`):
              `method="delta"` (`Dirac`, eMoScat `DiracTestFunction2d` — a
              fixed-point line projection, needs `position`) and
              `method="flow"` (`Flux`, eMoScat `FluxTestFunction2d` — a
              fixed-surface Wronskian flux, needs `surface`; built on the new
              `qscat.dvr.dvr_first_derivative_at_node` DVR-derivative-at-a-node
              primitive). `td_ve_cross_sections_all` runs ONE shared propagation
              driving all three and returns `{"tw":, "delta":, "flow":}` — the
              honest, identical-dynamics three-way comparison. On N2, all three
              converge to the same `driven.ve_cross_section` TI oracle to ~3%
              at a converged grid (delta 0.971, flow 0.970 at E=0.10 Ha);
              cross-method spread at an under-converged grid (~20-25%) is a
              convergence diagnostic, not a disagreement — see
              docs/physics/td-extractors.md. All three extractors also
              implement `axis="nuclear"` (`Flux`/`Dirac`/`TannorWeeks`,
              `td_extractors.py`) — the DISSOCIATIVE ATTACHMENT (DA)
              generalization: the outgoing side moves from the electronic
              to the nuclear coordinate, projecting onto `n_channels` anion
              electronic bound states (`anion_electronic_states`) instead of
              neutral vibrational levels (no elastic free-reference
              subtraction — DA has no `v'==v_init` diagonal). Wired as
              `time_dependent.td_da_cross_section(method="flow"|"delta"|
              "tw")` (default `"flow"`, the natural DA extractor) and
              `td_da_cross_sections_all` (ONE shared propagation, all
              three), the TD sibling of `dissociation.da_cross_section`
              (below) — σ_DA uses `C_DA=π` (not the TI oracle's literal
              `4π³`; the two reconcile via `S=1−2πiT`). **Key finding: TD-DA
              needs a LARGE electronic launch-box grid (incident well
              inside `r_max`), NOT the small TI `da_grid` — an off-box
              incident diverges ~1e6×, a coarse nuclear grid reads σ≈0** (the
              fine per-molecule nuclear deck is unchanged/reused). F2/NO
              three-way validation converges to `da_cross_section` to a
              flow/delta plateau ~0.86-0.97, tw converges to order ~1 but
              oscillates ~0.55-1.42 (the noisiest, most test-packet-sensitive) — see
              docs/physics/td-da.md. The Coulomb-generalized
              `riccati_hankel_en_mass` (`qscat.special.radial`, the
              mass-generalized outgoing Hankel half already used by the
              nuclear extractors' `eta_outgoing`) is validated (μ=1
              byte-identical reduction + Wronskian) and is what SP3 (TD-DR,
              H₂⁺) will drive at nonzero charge. Also
              `dissociation.da_cross_section` (exact TI **dissociative
              attachment**: the same driven Ψ₊ solve projected onto the nuclear
              dissociation channel with the rearrangement interaction
              `V_DR = V_int + v0 − V_int(r,R→∞)`, σ_DA=4π³|T|²/2E; plus
              `anion_electronic_states`, `v_dr_diag`) — the DA magnitude needs a
              per-molecule NUCLEAR grid (fast K_R~58 exit wave), built by
              `grids.segmented_grid` from eMoScat's per-molecule decks — see
              docs/physics/diatomic-ve-cross-sections.md. Also
              `dissociation.dr_cross_section` (**dissociative recombination**
              for the IONIC H₂⁺: `da_cross_section` generalized to a Coulomb
              incident (`channel_vector(charge=−1)`) + a LOOP over the Rydberg
              exit series; the first non-laptop model, ~1.15M unknowns at full
              size — Docker/MUMPS) — see docs/physics/h2plus-dr.md. Also `lcp` — the
              **LCP (local-complex-potential) approximation** OF the DA (the
              approximation under test vs the exact `da_cross_section` oracle):
              `local_complex_potential` (model-independent R-dependent
              `V_d(R)/Γ(R)` via `qscat.ecs.find_resonance_pole`, seeded from the
              anion state) + `lcp_da_cross_section` (1-D TI resolvent, boundary-
              wavefunction-VALUE flux `ψ(X)=ψ_coeff[b]/√w_b`, on the fine
              per-molecule nuclear grid). Also `lcp_ve_cross_section` — the
              VIBRATIONAL-EXCITATION sibling on the same doorway/driven-equation
              formula (sparse `A(E) = E_tot·I − H_res` with `SparseLU.refactor`
              sweep reuse), graduated from the `n2_ti_cross_section` toy model
              — see docs/physics/n2-cross-section.md. The LCP's error is SYSTEMATIC and
              ENERGY-DEPENDENT, not a fixed percentage: on F₂ the ratio
              LCP/exact sweeps 0.263 → 1.736 across 0.010–0.050 Ha, crossing
              unity near E≈0.032 (dense 41-energy sweep, 2026-08-17), so it
              under-predicts below ~0.03 and over-predicts above. On NO it
              fails outright away from threshold: the exact σ_DA decays 13
              orders of magnitude while the LCP stays flat, reaching a ratio of
              1.8e5-2.2e7 (corrected 2026-08-24). Documented departures also in VE-elastic (LCP misses the
              non-resonant background) — see
              docs/physics/diatomic-ve-cross-sections.md. Also
              `resonance_levels`/`lcp_resonance_levels` (+ `ResonanceLevels`): the
              BORN-OPPENHEIMER approximation to the resonance energies -- the nuclear
              eigenvalue problem IN the complex curve, `H_N = T(mu) + V_d -
              i*Gamma/2` on the nuclear FEM-DVR-ECS grid, giving complex quasi-bound
              levels `E_v - i*Gamma_v/2` (the thesis's `omega_j`, promoted from
              eMoScat's real-part-only levels). Physical levels are picked by two
              NUCLEAR ECS angles (`qscat.ecs.match_angle_stable`, the multi-state
              sibling of `find_resonance_pole`); the electronic pole walk runs ONCE
              since `E_res(R)` at real `R` is angle-independent. A golden-rule
              comparator (`Gamma=0` levels + `<chi|Gamma|chi>`) rides along and
              reproduces what eMoScat/the thesis computed -- its divergence from the
              complex result is the non-perturbative signal. NOT Siegert
              pseudostates (see docs/physics/lcp-resonance-levels.md). Plus
              `exact_resonance_states` (+ `ExactResonanceStates`) is the
              approximation-free counterpart of those BO levels: poles of the FULL
              2-D S-matrix, found by generalizing ECS angle stability to TWO angles
              (three spectra — base, electronic-angle moved, nuclear-angle moved —
              with a state accepted only if it survives both, and both residuals
              reported). Seeds are passed in, so the exact solver never calls the
              approximation it measures. On N2 the exact pole lies BELOW the BO/LCP
              level in both position and width at every level; converged only for
              v=0 (0.22 meV in position, 0.30 meV in width, over an electronic box
              grown 24→72 bohr — order converges at 8, the real-region extent is the
              limiting knob, and widths converge slower than positions) — see
              docs/physics/exact-2d-resonances.md -- but read that claim with the
              overlap check below: on the finer 46k deck the POSITION difference
              turns positive at v>=3, while the width difference stays negative.
              Plus `bo` and `assignment`, the **verification layer** that says
              whether a pole is a resonance at all. Angle stability is NECESSARY
              AND NOT SUFFICIENT: on H2+ four of 57 angle-stable poles scored
              overlaps of 6e-4..7e-3 against a Born-Oppenheimer basis where genuine
              states score 0.87-0.99. `bo` builds the reference states
              `phi_j(r;R) chi_v(R)` -- `electronic_curves` for the BOUND (ion /
              Rydberg) case, `resonance_curve` for the RESONANT (neutral / anion)
              one, `bo_basis`/`bo_basis_from_levels` putting a vibrational ladder
              in either, all phase-aligned across R (without which the product
              flips sign at random R and every overlap integrates to ~0).
              `assignment` pairs a pole to a level BY OVERLAP (c-product, which is
              BILINEAR -- values above 1 are legitimate, inflated by
              `1/sqrt(rho_a rho_b)` with `rho = |c(psi,psi)|/||psi||^2`: N2's broad
              resonances score 1.02-1.19. Do NOT "fix" this by dividing by the
              Euclidean norms -- that denominator reweights the rotated tail the
              c-product exists to cancel, and it re-ranks H2+'s diffuse states onto
              the wrong partner; measured, not argued) and returns one of SEVEN
              verdicts. `real_weight` is the check the overlap CANNOT make: the
              c-product cancels the rotated tail by construction, so a state 97%
              outside the box still pairs at 0.99 with the BO product it genuinely
              is -- on H2+ that blindness hid 18 of 57 poles whose Rydberg orbitals
              are larger than the 300-bohr box (`box-limited`; it moved the measured
              BO regime split from 0.264/3.375 meV over 40 rows to 0.457/3.702 over
              24). N2's poles sit at `real_weight` 0.96 and are untouched.
              Separately, `spurious` vs `basis-limited` is separated by the
              CLOSED-CHANNEL energy constraint (`admissible_levels`/`basis_covers`,
              `n_eff = 1/sqrt(2*binding)` to the nearest threshold ABOVE, so a
              higher vibrational level admits only a LOWER Rydberg index and the
              admissible set is finite and computable). Also `pair_one_to_one`
              (Hungarian bijection, a cross-check) and `peak_positions`/
              `peak_alignment` (distance to observed cross-section peaks IN UNITS
              OF A RESONANCE WIDTH -- the only scale on which "lands on the peak"
              means anything). N2's own poles are now overlap-verified: 6/6 clean,
              0/6 pairing disagreements -- see
              docs/physics/h2plus-resonance-states.md.
              Also `nrm` — the **NONLOCAL
              RESONANCE MODEL** (Houfek/Rescigno/McCurdy, PRA 77, 012710
              (2008)): the rung between `lcp` and the exact solver, keeping the
              energy dependence and the NONLOCALITY the LCP discards. Solves the
              same 1-D nuclear equation but with a complex, energy-dependent
              kernel `F(E,R,R')` (Eq. 55-61, a P-projected resolvent expanded in
              the fixed-R electronic eigenbasis) in place of `−iΓ(R)/2`:
              `scattering` (`φ_k⁺` at real energy), `discrete_state` (the two
              `φ_d` choices), `coupling` (`V_dk⁺`, `Γ=2π|V_dk⁺|²`), `ingredients`
              (`E_n(R)`, `V_dn(R)`, `V_d(R)`), `nonlocal_potential`, and
              `nrm_da_cross_section` — plus `vibrational_excitation`
              (`j_dk`, `t_resonant`, `t_background`, `nrm_ve_cross_section`), the
              VIBRATIONAL-EXCITATION route on the same `F(E)`/ingredients/`Ψ_d⁺`:
              `σ = 4π³|T^res + T^bg|²/k_i²` from PRA 77's two-potential
              decomposition (Eq. 28/31/37/38), with `include_background=False`
              giving the paper's bare "nonlocal" curve against its "nonlocal +
              bg" one. Eq. (34)/(35) (radial case, real `φ_d`, `H_el` Hermitian)
              is what lets both T-matrix terms use the NON-conjugated `V_dk⁺`;
              that is a separate argument from p. 012710-6's complex-symmetric
              ECS c-product, and one citation must not do both jobs.
              `scattering_state_minus` (`φ⁻`, gated by a Hankel decomposition
              against analytic asymptotics) is the identity behind that choice,
              not a consumer of it — Eq. (37)'s bra carries `(φ⁻)* = φ⁺` at the
              FINAL channel energy. The method is completely determined by the
              discrete state, and BOTH of PRA 77's implementable choices are
              provided: `PhysicalDiscreteState` (A, the R-dependent scattering
              function at `Re E_res(R)`) and `AsymptoticDiscreteState` (B, the
              R-independent bound state); choice C ("compact") is not
              implemented. NOTE the naming collision: this package's
              `v_d_discrete` is PRA 77's Eq. (20) `V_d = V_0 + <φ_d|H_el|φ_d>`,
              NOT `qscat.core.lcp`'s `Vd` (which is `E_res(R)`) — they only
              "almost coincide", and differ by 0.0053 Ha (F₂) / 0.0229 Ha (NO) at
              the doorway. MEASURED: on F₂ choice B reproduces the exact
              `da_cross_section` oracle to 0.06-0.33% at four of five anchors and
              1.9% at the lowest (E=0.010, nearest threshold), against an oracle
              floor of 5.4e-5 — beating the LCP by 39-758× at every anchor —
              while choice A is degraded (ratio 0.29-0.90), the Born-Oppenheimer
              breakdown PRA 77 predicts. On NO, which is BEYOND the paper's
              tested range (PRA 77 publishes no DA cross section for NO at all —
              NO's DA channel opens at +0.1719 Ha against the paper's plotted
              0.01-0.08 Ha window, so it is energetically shut there), choice B
              reproduces the exact oracle to 1.5-1.9% — the SAME quality as F2.
              The "5-8 order collapse" this entry used to record was the exact
              2-D oracle's own error: `da_cross_section` used a post-form VOLUME
              T-matrix whose required cancellation is integrand/answer (x2.7 on
              F2, x1.5e6 on NO), so every box edge leaked. FIXED 2026-08-24 by
              reading the outgoing flux instead; gated against Vana 2017
              Fig. 3.14. `dr_cross_section` (H2+) still uses the volume form and
              is FLAGGED, not audited. MEASURED
              FOR VE, the channel PRA 77 plots for EVERY molecule in its study:
              choice B + background reproduces the exact `driven.ve_cross_section`
              oracle to better than 0.7% on BOTH N₂ (11 energies, 0.06-0.16 Ha)
              and F₂ (gated at 2 anchors, 0.02 and 0.04 Ha; recorded sweep:
              6 energies, 0.02-0.09 Ha), elastic and first-inelastic alike
              (0.99623-1.00692 worst over all four molecule/transition pairs),
              while choice A degrades to 0.565-1.140 — and the reason B is that
              good is PHYSICS, not luck: an R-INDEPENDENT `φ_d` carries no `∂_R`
              derivative couplings, so the model is FORMALLY EXACT and the
              residual is discretization error. The comparison is DIFFERENTIAL
              (both routes on the same grids), so it validates the model
              reduction, NOT the grid — absolute normalization is anchored by
              validation/n2/exact2d.py against Houfek. DA on F₂/NO, VE on
              N₂/F₂ (NO VE not run, though the paper publishes it — the natural
              follow-on); the NO DA question is CLOSED — it was the oracle, not the
              model — see
              docs/physics/nonlocal-resonance-model.md. The nrm package also has
              a TIME-DEPENDENT half (`extended`, `propagation`,
              `td_cross_section`): Gertitschke & Domcke PRA 47's memory integral
              is RESUMMED — one auxiliary nuclear packet per projected electronic
              state turns it into time-local propagation under a sparse arrow
              block Hamiltonian `H_ext`, whose auxiliary blocks eliminate back to
              PRA 77 Eq. (52) exactly (gated at 4.4e-14). The half-Fourier
              transform of the propagated packet IS the TI `Ψ_d(R;E)`
              (`Ψ_d^TI = −i∫₀^∞ e^{iEt}Ψ_d(t)dt`), so `da_sigma_from_psi` /
              `t_resonant` / `t_background` are reused unchanged: N₂ vector-to-
              vector 1.7e-4, F₂ σ_DA 0.986–1.014. **`n_states=None` is MANDATORY
              here** — truncating the arms makes `H_ext` non-dissipative and the
              propagation diverges exponentially, non-monotonically, which is
              benign for the TI RESOLVENT but fatal for a PROPAGATOR;
              `propagate_nrm` warns when a packet grows. Energies are propagated
              ONCE via an SVD of the launch matrix (PRA 47 Eq. 2.17's rank-1
              claim, generalized). The route is SLOWER than the resolvent (773 s
              vs 246 s on F₂'s deck) and is justified by `S(t)`/`⟨R⟩_t`/`⟨P⟩_t`,
              not cost. `td_nrm_ve_cross_section` is the VE sibling (N₂ 2.7e-4,
              F₂ 9.8e-5 against the TI route; `T^bg` is energy-domain so only
              `T^res` changes route), and `markovian=True` on either entry point
              propagates PRA 47 Eq. (2.15)'s LOCAL limit — the same solve with
              the arms removed, so `N_R` square instead of `(1+n_states)·N_R`,
              seconds instead of ~30 min. **Eq. (2.15) takes `qscat.core.lcp`'s
              `Vd` (= `E_res`), NOT `v_d_discrete`** — measured, 1.0002 against
              0.346/0.419/7.14, because Eq. (2.14) makes `V_d + Δ_L = E_res + V_0`
              and Eq. (20)'s `V_d` is short by the level shift; the wrong choice
              is locked out by a test. `markovian` substitutes the local doorway
              at BOTH ends (that is what reproduces the LCP rather than a hybrid)
              and REFUSES `include_background=True`, since Eq. (37)'s background
              needs a `φ_d` the local model does not have. Measured on F₂, the nonlocal
              and local packets are nearly identical (`⟨R⟩` within 0.01 bohr,
              both unimodal) — **no packet splitting**, unlike PRA 47's H₂⁻; the
              shipped LCP's packet differs because of its DOORWAY, not its
              kernel — see docs/physics/nrm-time-dependent.md. Plus
              `channels`, `grids`
              (parameterized FEM-DVR-ECS builders + `segmented_grid` for
              eMoScat's `(n_elem, endpoint)` deck format, plus `ecs_angle_family`
              -- the three-grid family `exact_resonance_states` needs, which now
              VALIDATES that each partner moves exactly one ECS angle and shares
              every real node; an identical partner grid used to be accepted and
              would pass the whole rotated continuum), `vibrational` (`v0`
              passed in), `wavepacket`, `correlation`, `plot`
              (`plot_resonance_levels` takes an explicit `pairing` -- its default
              sorted-index pairing is only correct when both level sets are
              complete and ordered alike). **`qscat.core` never
              imports `qscat.model`/`projects` at runtime** (depends only on the
              `ResonanceModel` protocol; enforced by
              `test_core_no_model_import.py`) — see
              docs/physics/qscat-core-scattering.md.
            - qscat.model: everything tied to a specific model — the
              `ResonanceModel` protocol (the contract `qscat.core` depends on;
              carries a `charge` attribute — 0 neutral, −1 for a cation),
              `DiatomicResonanceModel` (the shared Morse+sigmoid+Gaussian
              NEUTRAL form) + the `N2`/`NO`/`F2` registry, and
              `IonicResonanceModel` (the H₂⁺ Morse + σ-capture + `−1/r` Coulomb
              form) + the `H2P` registry entry (the first ION), and
              `FlexibleDiatomicModel` (`qscat.model.flexible`: EMO neutral +
              Gaussian well with `lam(R)`/`alpha(R)` as `SmoothR` sigmoids or
              long-range-correct `TailR` forms + optional shell; `from_diatomic`
              embeds N2/NO/F2 exactly) + the `O2` registry entry — the first
              FITTED model (the potential factory's fit to Alt & Houfek 2021's
              curves, its constants locked to the committed report by
              `validation/factory/test_o2_report.py`) and its spin–orbit
              components `O2_SO12`/`O2_SO32` (⅓ each). `qscat.model.N2`
              is the single source of truth for the N2 model (the N2 projects
              consume it via thin shims). Adding a molecule = a registry entry +
              validation, never solver code — see
              docs/physics/qscat-core-scattering.md. The Coulomb channel/special
              functions for ions live in `qscat.special.coulomb`
              (`coulomb_f_en`/`g`/`h1_en`, mpmath, the charge-z generalization of
              `riccati_bessel_en`).
            - qscat.tuning: the automatic FEM-DVR-ECS **discretisation tuner** —
              deterministic primitives that compute the minimal-DVR-point grid at
              a target precision from the potential + energy range, replacing the
              human "good eye" for element lengths. `analyze` (potential →
              local-wavenumber profile), `mesh` (adaptive equidistribution
              elements — `∫k dx` ~const per element — + the h/p quadrature sweep),
              `ecs` (the double-ECS-capped angle + exp-growth absorbing tail),
              `probes` (decoupled 1-D convergence: nuclear/electronic + the
              `channel_representation` probe that catches the K≈58-under-resolution
              failures), `metrics`+`propose` (`propose_grid` a-priori assembler +
              cost model), `incident` (`IncidentSpec`/`tw_analysis` — the TW
              wavepacket/test-function placement). Driven by the
              `discretisation-tuner` skill (the supervised loop). The mesh's
              de-Broglie phase constant `C` (`qscat.tuning.mesh._PHASE_COEFF`)
              is CALIBRATED (`validation/tuning/calibrate.py`) against F2's
              genuinely-open dissociative-attachment channel — the tuner
              reproduces-and-beats that eMoScat deck on the 1-D probes (37%
              fewer points, clean rtol=1e-3 convergence on the K~78 DA wave)
              and its cheapest probe correctly flags the coarse shared
              N2-style grid as under-resolved for that same wave
              (`validation/tuning/test_emoscat_decks.py`); H2+'s proxy nuclear
              deck is likewise a clean reproduce-and-beat; N2/NO's proposed
              nuclear grids cost more points than their decks (traced to a
              fixed real-region extent default, not to `C`) — a documented,
              reported limitation. The gate's `@slow` 2-D spot-check on F2 is
              a genuine, load-bearing finding, NOT a rubber stamp: the 1-D
              probes pass on the reproduce-and-beat grid, but the actual
              sigma_DA is NOT 2-D-converged there (one nuclear h-refinement
              changes it ~5x, toward the eMoScat deck's own value) — the
              a-priori mesh, built only from `v0`'s classical k(x) profile,
              cannot see the narrow R~2.5-2.7 bohr interaction feature
              eMoScat's deck hand-resolves; the 1-D probes are necessary but
              NOT sufficient for this observable. That gap is now CLOSED: a
              resonance-aware DA nuclear path (`propose_grid(..., channel=
              "dissociation")` -- exit-wave DVR order sized off the adiabatic
              resonance curve + a local crossing super-refine, with
              `refine_to_2d_convergence` as the general model-agnostic
              fallback) converges F2's sigma_DA on the FIRST a-priori pass
              (1.6562, matching the eMoScat deck) at deck-parity size
              (1000 vs 974 pts, 1.027x) and gives H2+'s resonant grid ~4%
              under its proxy deck (489 vs 510 pts) -- see
              docs/physics/discretisation-tuning.md.
apps/       qscat-run — THE single execution surface: one YAML config runs
            every observable (VE/DA/DR cross sections, wavefunction snapshots,
            vibrational eigenstates, resonance states and BO/LCP resonance
            levels) across the TI/TD/LCP methods, and writes csv/npz/png
            artifacts. Per-molecule curve drivers were retired into it, so a
            new figure is a config, not a script. See apps/qscat-run/README.md
            for the observables matrix and the config schema.
native/     Rust kernels (qscat-kernels crate) built with PyO3/maturin,
            mirroring validated Python APIs for hot paths
projects/   per-problem research and toy models — lifecycle stages 1-2
            - `n2_ti_cross_section`: time-independent (resolvent/driven-
              equation) N₂ vibrational-excitation cross-section inputs
              (`nuclear_grid.py`/`vibrational.py`/`vres.py`), built on
              `qscat.dvr`/`qscat.ecs` and the N₂ resonance pole finder; the
              solver itself graduated to `qscat.core.lcp.lcp_ve_cross_section`
              — see docs/physics/n2-cross-section.md.
            - `n2_td_cross_section`: time-dependent (Crank-Nicolson
              propagation + energy transform) route to the same N₂
              vibrational-excitation cross section (`propagator.py` — thin
              re-export of `qscat.evolution.make_cn_stepper`;
              `td_cross_section.py` — doorway wavepacket propagation under
              `H_res`, correlation function, energy transform), validated
              against the TI solver as an exact differential oracle — see
              docs/physics/n2-td-cross-section.md.
            - `n2_2d_cross_section`: the exact 2-D (electronic r × nuclear R)
              driven Lippmann-Schwinger solver for the same N₂
              vibrational-excitation cross section — no local-complex-
              potential reduction (`electronic_grid.py`/`channels.py`/
              `hamiltonian2d.py`/`cross_section_2d.py`/`convergence.py`/
              `nuclear_density.py`), validated standalone (free-particle and
              first-Born limits, S-matrix reciprocity/unitarity) and then
              gated against Houfek's data as an independent implementation of
              the same model/method; once gated, it is the ORACLE the
              1-D LCP solver is compared against — see
              docs/physics/n2-2d-cross-section.md. `ve_cross_section_2d` sweeps
              energies analyze-once/refactor-per-energy (via
              `SparseLU.refactor`) — same σ, cheaper sweep. The generic,
              experiment-agnostic `cross_section_plot.plot_cross_sections`
              (no physics, reference passed as an argument — reusable for
              F₂/NO) renders the dense σ_{0→v'}(E) curves; the committed N₂
              curve vs Houfek is docs/physics/figures/n2-2d-ti-cross-section.png
              — see docs/physics/ti-energy-sweep-reuse.md.
            - `n2_2d_td_cross_section`: time-dependent (sparse Crank-Nicolson
              propagation + Tannor-Weeks energy transform) route to the SAME
              exact 2-D N₂ vibrational-excitation cross section as
              `n2_2d_cross_section` — an incident Gaussian wavepacket
              `g(r) chi_0(R)` propagated under `H_2D`
              (`wavepacket.py`/`td_propagation.py`/`correlation.py`/
              `td_cross_section.py`/`convergence.py`),
              validated against the exact 2-D solver as an exact differential
              oracle (σ_TD/σ_TI = 0.973 at E=0.10, 0.988 at E=0.15 — measured
              2026-08-17 by `validation.n2.experiment` group F1). The
              ELASTIC (v'=v_init) channel subtracts a free-particle (V_int=0)
              reference S_free(E) instead of a literal 1 — the transform's
              outgoing normalization makes S_free≈2π²≠1, so |S−1|² left a ~500×
              spurious elastic background; `td_ve_cross_section_2d(...,
              subtract_free_reference=True)` (default) runs the reference
              propagation. The propagator is the order-3 diagonal Padé
              (`make_pade_stepper`, dt=1.0, eMoScat's setting) — order-1 CN
              under-converged, capping TD-vs-TI at ~10-15%; with order-3 the
              TD cross section matches the exact TI (and Houfek) to ~1-2%
              median across 0.04-0.18 Ha for elastic + first excitations,
              boomerang oscillations resolved point-by-point — see
              docs/physics/n2-2d-td-cross-section.md and the
              td-elastic-wavepacket-normalization note.
            - `potential_factory`: the toy-stage POTENTIAL FACTORY — fits a
              `FlexibleDiatomicModel` (EMO `v0` + Gaussian well with `lam(R)`
              AND `alpha(R)` + optional shell; embeds N2/NO/F2 exactly) to a
              tiered `Target` (T0 neutral curve, T1 pole curves, `asymptote`,
              T3 the published energy-dependent width) in stages that
              stop-and-report. `lam(R)` is either Houfek's sigmoid `SmoothR`
              (the published models) or the long-range-correct `TailR`
              (`f_inf + (1 − y_q)·P(y_p)`, dies as `R^−q`; O₂ needs it — the
              sigmoid could hold the table OR the asymptote, not both). The
              ASYMPTOTE is theory, not the figure: `V_0(∞) = 0`, `V_ion(∞) =
              −EA` (atomic) through a declared tail (`polarisation_tail`,
              `−α_d/2R⁴`), valid from a PER-MOLECULE `ResonanceTarget.R_inf`
              the operator sets from the literature; pinned at `R_inf`,
              2.2 `R_inf`, 5 `R_inf` in the polish and gated as its own tier
              (one node at 10 bohr used to leave V_ion 0.2 eV off at 20).
              Proven by round-tripping the existing models' OWN calculated
              curves (`extract_target`) back to their constants/curves. Nothing
              is fitted to experiment — see docs/physics/potential-factory.md.
validation/ analytic benchmarks, golden datasets, convergence studies
            - `validation/n2/`: N₂ electron-scattering harness; its C5 group
              anchors this solver's σ_{0→v'}(E) against Karel Houfek's
              independent `CSVE.V00.J00` data (documented cross-model
              tolerance, not exact agreement); its D1 group cross-checks the
              TD solver against both the TI solver (rtol<=0.10) and the same
              Houfek data at the GATED anchors (`td_check.py`); its E1 group
              (`exact2d.py`) reruns the same 6 anchors through the exact 2-D
              solver, GATED at `GATED_RTOL=1e-3` against Houfek (a tight
              differential-oracle bound, not the LCP's cross-model factor-3
              band) with the two DOCUMENTED-LIMITED LCP anchors reported as
              NOTE rows showing how the exact model closes that gap; its F1
              group (`td_exact2d.py`) reports the time-dependent 2-D solver's
              σ_TD-vs-σ_TI agreement at the two validated anchors as
              recorded, cited NOTE rows rather than a live in-harness
              propagation — a full run costs ~210–250s (measured), far over
              the harness's per-group budget, so the genuine PASS/FAIL gate
              lives in `n2_2d_td_cross_section/test_td_cross_section.py`'s
              `@slow` tests instead. Run the harness with
              `uv run python -m validation.n2.experiment`. `ti_curve.py` is the
              N₂-specific driver behind the dense exact-2D σ(E) figure: it reads
              Houfek (`loader`) and calls the generic `plot_cross_sections`
              (validation may import projects; projects must not import
              validation), and `test_ti_curve.py` gates the dense curve against
              Houfek at the anchors — see docs/physics/ti-energy-sweep-reuse.md.
              `pole_verification.py` (`python -m validation.n2.pole_verification`)
              answers the question N₂'s exact 2-D poles were published without:
              are they resonances, and is the sorted-index exact/BO pairing right?
              Both yes — 6/6 clean overlaps, 0/6 pairing disagreements — using the
              NEUTRAL basis path (`resonance_curve` + `resonance_levels` +
              `bo_basis_from_levels`), the counterpart of H₂⁺'s Rydberg path.
            - `validation/h2plus/`: the H₂⁺ DR exact-resonance campaign —
              `rydberg_levels` (a thin shim over `qscat.core.bo` now), `exact_poles`
              (the seeded pole campaign + `grid_family`), `bo_overlap` (the verdict
              report), `dr_levels_figure` and `resonance_state_figures` (the five
              committed figures), `reference_levels`/`reference_sweep` (the
              published ω_i^j table as a gated oracle, and the σ_DR sweep). The
              solver machinery moved into `qscat.core.bo`/`assignment`; what stays
              is the campaign — which curves, which windows, which seeds — see
              docs/physics/h2plus-resonance-states.md.
            - `validation/factory/`: the potential factory's base experiments
              (`python -m validation.factory.base_experiments --molecule N2|NO|F2
              --stage curves|fit|xs|all`): the published models' resonant
              curves `E_res(R)`/`Γ(R)`/`V_ion(R)` on a four-grid electronic
              ladder (converged to 1e-9..1e-7 Ha; the crossing node is gated
              out by design), the factory's refit from those curves, and the
              exact 2-D VE (+ DA for NO/F₂) cross sections on the `emoscat`
              decks for the published vs the refitted model — agree to 1e-9
              (NO's 1e-19-bohr² DA tail to 4e-7). Results under `results/`,
              figures `docs/physics/figures/{n2,no,f2}-factory-*.png`. Run it
              in the MUMPS container with `OMP_NUM_THREADS=1` — a 32-thread
              OpenBLAS is ~400× SLOWER on the tiny electronic eigenproblems.
              Also the O₂ IMAGE MATCH (`fit_o2.py`, `targets/o2.py`,
              `extract_fig2.py`): Alt & Houfek 2021's Fig. 2 is vector-extracted
              from the PDF (precision ~0.02 eV, no digitising), and the factory
              fits it over the FULL 1.85–6 bohr range to T0 met / T1 MET
              (E_res rms 20 meV = the extraction floor, Γ 8 %/14 %) / asymptote
              met (R_inf = 14 bohr, α_d(O) = 5.3 a.u.) / crossing 2.289 — in
              87 s on the laptop — with T3 open (not discrete-state-consistent:
              choice B vs the paper's Breit–Wigner pole width). `o2_levels.py`
              is the SPECTRAL CHECK, the metric that predicts the VE figure:
              the anion's quasi-bound levels in the fitted vs the extracted
              curve — peak positions within ±7 meV over v = 0..29 (0–2.6 eV),
              widths within ~10 % (`results/o2-anion-levels.csv`). The fitted
              model is `qscat.model.O2` (locked to the report by
              `test_o2_report.py`, including the y_p FRAME radius the report
              does not record); `o2_grids.py` builds its decks with the
              discretisation tuner (electronic as proposed; nuclear cut at 8
              bohr — DA is closed, the tuner's 18-bohr default is empty
              space — then h-REFINED ONCE: the 2-D spot check moved σ(0→1)
              at 1.36 eV by 69 % on one refinement, converged < 2 % after;
              a comb of meV peaks needs levels far tighter than the probe's
              1e-3; 324 × 549 = 178k unknowns; `test_o2_grids.py` locks the
              `O2:tuner` preset to it) and `o2_ve_energies.py` writes the
              LEVEL-AWARE energy mesh (background grid + 15 points across ±5
              widths of each level) into `apps/qscat-run/examples/o2-ve.yaml`
              — a uniform sweep walks past 0.01–8 meV peaks. `extract_fig5.py`
              vector-extracts the paper's OWN NRM/LCP VE curves (Fig. 5, six
              panels; tick labels are glyph outlines so the ranges are fixed
              in `PANELS`; upper-envelope centreline keeps meV peaks; the
              legend key samples are stroked INSIDE the curve paths and are
              masked by position — unmasked, every panel's "maximum" was the
              key line) and `o2_ve_figure.py` overlays the exact 2-D result
              on them — theory against theory, and only spin–orbit resolved
              (Fig. 5's peaks are doublets of two ⅓-weight components; an
              unsplit-vs-Fig. 5 comparison means nothing). Sweeps cost 0.38
              s/energy with MUMPS on sadaharu (3343 energies × six channels
              in ~1290 s; 46 s/energy with SuperLU on the laptop), and the
              mesh must DRAW the peaks — at Γ/1.5 spacing every height read
              ×0.69, a Lorentzian missed by Γ/3; the committed mesh is Γ/10,
              121 points across ±6Γ.
              SPIN–ORBIT: `extract_fig1.py` pulls the paper's Δ_SO(R) (Fig.
              1), `o2_target(so=±1)` moves the anion curve ∓Δ_SO/2, and
              `fit_o2_so.py` makes `qscat.model.O2_SO12`/`O2_SO32` from `O2`
              by POLISH ONLY (`refine_resonance` — the full re-track/re-smooth
              pipeline against a 10 meV shift fell into a wrong basin, E_res
              rms 31 mHa; the polish keeps the parent's 0.74 mHa); each
              component gets its own mesh and `o2_ve_figure.py --so12 --so32`
              sums them at ⅓ each against Fig. 5's doublets (`results/
              o2-so{12,32}-ve/`, figure `o2-2d-ti-ve-spin-orbit-vs-alt-
              houfek.png`): both members within 1–8 meV, heights 0.9–1.1,
              doublet separation 19.0–19.3 meV (paper's model 17.8, Allan's
              measurement 19.6 ± 1.0 — noted, not claimed). Nothing is
              compared with experiment — see
              docs/physics/potential-factory.md.
            - `validation/diatomic/`: the NO and F₂ exact-2D VE/DA/LCP cross
              sections — the model port, the first consumers of sub-project A
              beyond N₂. The per-molecule *curve/figure drivers* were RETIRED in
              the qscat-run consolidation (docs/superpowers/plans/2026-08-15-
              unified-experiment-observables.md): the exact-2D TI VE σ(E), TI
              σ_DA(E), and the LCP-vs-exact σ_DA overlay are now produced from
              config through `apps/qscat-run` (e.g. `apps/qscat-run/examples/
              f2-da-lcp-vs-exact.yaml`, `methods: [ti, lcp]`) — no dedicated
              solver code. The committed figures (`docs/physics/figures/{no,f2}-
              2d-ti-cross-section.png`, `{f2,no}-2d-ti-da-cross-section.png`,
              `{f2,no}-2d-da-lcp-vs-exact.png`) remain the sub-project-A/B
              deliverables (the LCP's error is energy-dependent, not a fixed
              percentage — see the qscat.core.lcp entry above). `ve_nrm.py` +
              `test_ve_nrm.py` are the NRM-vs-LCP-vs-exact VIBRATIONAL-EXCITATION
              gate (N₂ and F₂, four routes on one deck), and `ve_nrm_figure.py`
              is one of THREE surviving figure drivers — deliberate exceptions
              to the retirement above (the others are `td_nrm_figures.py`, whose
              panels need packet diagnostics no config exposes, and
              `da_figure.py`, whose LCP/exact ratio panel and published-reference
              overlay on a DA observable are likewise unreachable from
              `qscat_run.artifacts`) — because its figure needs both
              `include_background` settings AND the LCP's VE route, and neither
              is reachable from a single qscat-run config (`apps/qscat-run/
              examples/n2-ve-nrm-vs-exact.yaml` is the config form of the rest of
              it). It writes `docs/physics/figures/n2-ve-nrm-vs-exact.png`, the
              2x2 panel laid out to be compared by eye against PRA 77 Fig. 4
              (choice A) and Fig. 8 (choice B). `config.py`
              remains, trimmed to its one
              surviving job: the eMoScat per-molecule NUCLEAR deck definitions
              (`MoleculeConfig.da_grid()`) the **discretisation tuner** reads as
              its reproduce-and-beat reference; `test_da_grid.py::
              test_diatomic_decks_match_presets` locks those decks byte-identical
              to `qscat_run.presets`' copy (layering keeps the two separate).
              No independent golden data exists for NO/F₂ (only N₂ has Houfek's),
              so the exact solver IS the oracle — see
              docs/physics/diatomic-ve-cross-sections.md.
            - `validation/tuning/`: the `qscat.tuning` discretisation tuner's
              own calibration + gate (sub-project #8/final). `calibrate.py`
              (`uv run python -m validation.tuning.calibrate`) sweeps the
              mesh's de-Broglie phase constant `C` and picks the smallest
              value making `propose_grid`'s F₂ nuclear grid
              reproduce-and-beat the eMoScat F₂ DA deck (the molecule with a
              genuinely open dissociative-attachment channel in its tested
              range) — calibrated to `C = 0.10`; H₂⁺'s proxy nuclear deck is
              swept alongside N₂/NO/F₂ for reporting (a clean reproduce-and-
              beat, unlike N₂/NO — see docs/physics/discretisation-tuning.md).
              `test_emoscat_decks.py` gates `propose_grid` for N₂/NO/F₂/H₂⁺
              against their committed/proxy decks (F₂ and H₂⁺ strictly;
              N₂/NO comparatively, since their channel-representation floor
              wavenumber is a conservative bound their OWN decks don't clear
              either — a reported finding, not a loosened gate) and confirms
              the cheap `probe_channel_representation` probe still flags the
              coarse shared N₂-style grid as under-resolved for F₂ DA's
              K≈58-78 wave. Its `@pytest.mark.slow` 2-D spot-check
              (`test_f2_2d_da_cross_section_spot_check`) is the design
              spec's "final 2-D solve confirming the tensor-product grid" —
              and a genuine finding: F₂'s reproduce-and-beat nuclear grid
              (609 pts) passes both 1-D probes but gives an UNCONVERGED
              sigma_DA (one nuclear h-refinement moves it ~5x, toward the
              eMoScat deck's own value); the refined-grid FAMILY converges
              (refine¹ vs refine² agree to 2%), gated as such rather than
              forcing a false base-vs-refined match — see
              docs/physics/discretisation-tuning.md.

`projects/` and `validation/` (and their sub-project directories) are real
Python packages (`__init__.py` present at every level); all intra-repo
imports are package-absolute (e.g. `from projects.n2_resonance.potential
import v0`, `from validation.n2 import loader`) — no bare intra-dir imports,
`sys.path.insert`, or `importlib.util.spec_from_file_location` workarounds.
Root `pyproject.toml` sets `pythonpath = ["."]` under
`[tool.pytest.ini_options]` so pytest resolves these packages from the repo
root without any path hacking. A module meant to be run directly (not just
imported), such as `validation/n2/experiment.py`, is invoked as
`python -m validation.n2.experiment`, not by file path.
benchmarks/ standalone measurement scripts (a real package, run via
            `python -m benchmarks.<name>`; imports `projects`, never imported by
            `projects`/`qscat`) -- `mumps_vs_superlu` measures the MUMPS vs
            SuperLU factor/solve/RSS/fill on the real N₂ 2-D matrices;
            `sweep_reuse` measures the energy-sweep symbolic-reuse speedup
            (`SparseLU.refactor` vs fresh-per-energy); both run in the Docker
            `test` image (need system MUMPS + `qscat[mumps]`).
reference/  read-only oracles: eMoScat (C++/CUDA snapshot), libXcuda
            (CUDA submodule) — for porting reference only, never imported.
            reference/literature/ holds the published sources: the PDFs are
            gitignored, but a TRACKED *.md reference note per source carries
            every published fact the repo relies on, each anchored to a page
            plus equation/table/figure. Written with the mastering-references
            skill; start at reference/literature/README.md.
docs/       specs/plans (docs/superpowers), physics notes (docs/physics),
            and ADRs (docs/adr)
docker/     layered CPU images: base (architecture/vendor) + app (build/
            test/runtime)
.claude/    skills and agents (below) plus repo permissions settings
```

## Tech decisions

- **Python >=3.12**, managed with `uv`; the exact interpreter is pinned in
  `.python-version` (`3.12`) for reproducibility and Docker parity.
- **uv workspace**: `pyproject.toml` declares `members = ["libs/*", "native/*"]`
  with `package = false` at the root. The canonical setup command is
  `uv sync --all-packages` — a plain `uv sync` prunes the `qscat` and
  `qscat_kernels` workspace members.
- **Rust + PyO3/maturin** for compiled kernels (`native/qscat-kernels`),
  built and installed into the active environment via
  `uv run maturin develop --manifest-path native/qscat-kernels/Cargo.toml`.
- **Testing**: pytest, hypothesis for property-based tests, pytest-benchmark
  for kernel benchmarks; every Rust kernel gets a Python differential test
  against its oracle.
- **Quality tooling**: ruff (lint), mypy --strict (types), cargo clippy (Rust
  lint). Pre-commit (`.pre-commit-config.yaml`) currently enforces only
  **ruff**, **ruff-format**, and **cargo fmt --check**; mypy and clippy are
  run manually / in CI, not (yet) pre-commit hooks.
- **Units**: atomic units throughout (`libs/qscat/qscat/units.py`); no ad hoc
  unit conversions scattered through method code.
- **Docker is layered, not a single Dockerfile**:
  - `docker/base.Dockerfile` builds `qmodeling-base` — the architecture /
    BLAS-FFT vendor layer: OpenBLAS + LAPACK(E) + FFTW3 dev libs, pkg-config,
    a Rust toolchain, and uv/Python 3.12. Code targets the standard CBLAS /
    LAPACKE / FFTW3 ABIs, so this layer is swappable. Default vendor choice
    (OpenBLAS/FFTW3) keeps ARM/Graviton open; an MKL x86-64 variant is a
    planned future alternative, not implemented yet. It also provisions system
    **MUMPS** (`libmumps-seq-dev` + `libscotch-dev`, the sequential build) for
    `qscat.linalg`'s optional MUMPS backend, and **synthesizes the pkg-config
    `.pc` files** Debian omits (the conda-forge `{d,z,c,s}mumps_seq` names
    python-mumps looks for) so the extra builds against the system library. It
    also installs **ffmpeg** — the matplotlib `FFMpegWriter` backend for
    `qscat.viz`'s `.mp4` animation output (the `.gif` path needs only pillow) —
    which flows to the `build`/`test-deps`/`test` stages so the ffmpeg-gated
    `.mp4` viz test renders rather than skips.
  - `docker/Dockerfile` is `FROM ${BASE_IMAGE}` and layers `build` →
    `test-deps` → `test` → `runtime` on top, using `uv sync --all-packages`
    for setup. The `build` stage adds `--extra plot` (matplotlib) so
    `qscat.viz` animation works in `runtime`, which copies `build`'s venv
    verbatim (no toolchain there to re-sync/rebuild the Rust kernel);
    `runtime` also installs **ffmpeg** for the `.mp4` writer. `test-deps`
    (`FROM build`) additionally adds `--extra mumps` so the MUMPS backend is
    available but runs no tests; `test` is `FROM test-deps` and adds BOTH
    tiers, run differently (docs/adr/0005): the fast tier in parallel
    (`-m "not slow" -n auto --dist loadfile`) and the slow tier SERIAL, since
    those decks are sized in gigabytes and concurrency would OOM a container
    before it saved wall-clock. Measured on a 32-core x86 host: the fast tier
    31 s in parallel, the slow tier 13m31s serial. `test-deps`/`test` are
    deliberately separate targets so that pulling in MUMPS/plot for a compute
    run (`docker/run.sh`, which builds `test-deps`) never also pays for, or is
    blocked by, a full suite run — only `docker/build.sh test` reaches `test`. The
    `mumps` extra stays test-only (`runtime` omits it, keeping python-mumps
    out of the production image). Both MUMPS and the `.mp4` viz test run in
    the container and `@skipif`/`@skip`-absent on a bare Mac (no system
    MUMPS, no ffmpeg), so the Mac suite stays green while the same tests run
    + pass in Docker — see docs/physics/mumps-sparse-backend.md.
  - **Where `qscat.viz` is actually tested.** GitHub CI runs
    `uv sync --all-packages`, which does NOT install the `plot` extra, so every
    test guarded by `pytest.importorskip("matplotlib")` — the whole `qscat.viz`
    suite — SKIPS there. Those tests are exercised only by a local run with the
    extra installed and by the Docker `test` image. A green CI run is therefore
    not evidence that a viz change works; run `libs/qscat/tests/test_viz.py`
    locally or in the container before believing one.
  - `docker/build.sh [test|test-deps|runtime]` builds the base image then the
    requested app target. A `test` build runs both tiers and fails the build if
    either does; `runtime`'s `CMD` prints `qscat <version> ready` (currently
    `0.1.0.dev0` — the version is read from the package, not pinned here).
- **GPU/CUDA and AWS deployment are deferred.** `reference/libXcuda` is kept
  as future-GPU-kernel reference only.

## Skills & agents

| Name | Kind | Use when |
|---|---|---|
| `qm-method-lifecycle` | skill | Adding or porting any QM method/capability — enforces design → toy → validate → optimize-in-Rust → promote-to-qscat. |
| `numerical-validation` | skill | Validating quantum/numerical code where exact equality doesn't apply: analytic benchmarks, convergence studies, conservation checks, differential testing. |
| `python-to-rust-kernel` | skill | A validated Python method has a proven hot path — scaffolding a PyO3/maturin crate under `native/`, mirroring the API, benchmarking, keeping Python as the oracle. |
| `containerize-and-run` | skill | Packaging a capability to run in Docker locally or prep it for AWS — reproducible multi-stage `uv` + `maturin` builds, CPU-only. |
| `qscat-conventions` | skill | Unsure how the project names or measures things — atomic units, FEM-DVR-ECS notation, tolerance defaults, standard-library layout. |
| `discretisation-tuner` | skill | Setting up (or distrusting) a FEM-DVR-ECS grid — supervises the `qscat.tuning` loop (analyze the potential → adaptive equidistribution mesh + h/p + double-ECS-safe tail → convergence probes at the energy extremes → 2-D spot-check → minimal-cost grid at a target precision), instead of hand-picking element lengths. |
| `mastering-github` | skill | Preparing a branch for review, or deciding whether a file may cite a spec/plan/issue/PR. Holds the rule that **main must stand alone** — a reader with only the clone must understand every shipped file — and the two procedures built on it: `/review-ready` (dissolve working-file content into permanent homes, prune references that don't travel with a clone, self-audit, tidy, flip draft → ready) and `/tidy-history` (rewrite a fix-on-fix branch into logical commits without changing the end state). |
| `mastering-references` | skill | Adding a paper to `reference/literature/`, writing or updating its reference note, or citing a published result from code or docs. The source PDFs are gitignored, so the TRACKED `*.md` notes are what the repository cites: every extracted fact carries a locator, and every locator carries a **page**. Enforces a page-offset check before any locator is written, and requires parity checks against the code to be run rather than asserted — that discipline has already caught three real discrepancies. |
| `release-review` | skill | A harsh whole-repo or release-readiness review, or a consistency audit after a multi-commit campaign — parallel read-only review tracks plus the drift checks (prose-vs-code, N-places, all-branches) that single-commit review misses. |
| `port-scout` | agent | Before porting anything from `reference/eMoScat` or `reference/libXcuda` — read-only archaeologist that extracts the math/algorithm, not the C++. |
| `physics-reviewer` | agent | Before promoting a method into `qscat` — reviews for physical/numerical correctness (units, conservation, boundary conditions, ECS handling, convergence), not style. |
| `rust-kernel-engineer` | agent | During the optimize-in-Rust stage — builds PyO3/Rust kernels in `native/` mirroring a validated Python API, with benchmarks and differential tests. |

For general engineering workflow — brainstorming, TDD, systematic debugging,
writing/executing plans, requesting code review — use the **Superpowers**
skills (`superpowers:brainstorming`, `superpowers:test-driven-development`,
`superpowers:systematic-debugging`, `superpowers:writing-plans`,
`superpowers:executing-plans`, `superpowers:requesting-code-review`, etc.).
These are general-purpose and complement the qModeling-specific skills above.

## Reference oracles

`reference/` holds read-only prior codebases: `reference/eMoScat` (C++/CUDA
FEM-DVR-ECS electron-molecule scattering code, snapshot copy) and
`reference/libXcuda` (recovered CUDA layer, git submodule). Never edit,
build, or import these as dependencies — they exist only to be read. Always
use the `port-scout` agent before porting a method out of them; it extracts
the underlying math/algorithm so it can be reimplemented cleanly in Python.

## Common commands

```bash
# Setup — installs qscat and builds the Rust qscat_kernels in one step.
# Plain `uv sync` is NOT sufficient: it prunes the qscat/qscat_kernels members.
uv sync --all-packages

# The fast tier — the CI gate, and the one to run while iterating.
# See docs/adr/0005-test-tiers-fast-and-slow.md for what the tiers mean.
uv run pytest -m "not slow" -n auto --dist loadfile

# The production-scale tier: real decks, converged grids, published-value
# comparisons. SERIAL, deliberately — these are measured in minutes and
# gigabytes (the H2+ DR example peaks at ~19 GB on its own), so running them
# concurrently OOMs a laptop long before it saves wall-clock. CI does not run
# this tier by default -- it runs on demand, when a reviewer puts a `validate:*`
# label on the PR or dispatches .github/workflows/validation.yml by hand. Run it
# locally before merging anything that touches the solvers rather than relying
# on that.
uv run pytest -m slow

# The whole suite. Note there is no `-n` here: a bare `pytest -n 8` (no
# `-m` filter) schedules the multi-GB slow decks against each other with
# nothing bounding the total, which is how a routine parallel run reaches
# tens of GB. Parallelise the fast tier or nothing.
uv run pytest

# On `--dist loadfile` rather than the default `--dist load`: memory is the
# binding constraint, not CPU. Several modules build their grid at module
# scope, and per-test distribution makes every worker touching the module
# build its own copy; keeping a file on one worker bounds that to one copy.
# BLAS threads: pinning (OMP_NUM_THREADS=1 etc.) buys nothing ON THIS MAC --
# measured 136s pinned vs 133s unpinned at `-n 8`, and 76.1s pinned vs 74.6s
# unpinned at a saturating `-n 12`, because scipy here links Accelerate. That
# result does NOT transfer to Linux: on the 4-vCPU GitHub runner, OpenBLAS starts one thread per core in EACH
# xdist worker, and the unpinned parallel run went SLOWER than the serial one
# it replaced (>22 min vs 24.7 min; 186s once pinned). CI pins all three
# variables — see docs/adr/0005 point 6.

# Right after a maturin build in the same step, skip the (already-satisfied)
# sync to avoid a redundant resolve:
uv run --no-sync pytest

# Rebuild the Rust kernel extension in place during kernel iteration
uv run maturin develop --manifest-path native/qscat-kernels/Cargo.toml

# Lint / type-check
uv run ruff check .
# mypy is type-clean over the qscat PACKAGE. Point it at `libs/qscat/qscat`,
# not `libs/qscat` -- the latter also collects `libs/qscat/tests`, where strict
# mode reports a few hundred pre-existing errors (untyped test helpers), so it
# looks alarmingly broken while the shipped code is fine. Type stubs for the
# Rust qscat_kernels extension are pending, so repo-wide strict mypy isn't
# claimed to pass yet either:
uv run mypy libs/qscat/qscat

# Build and test the CPU Docker images (base, then test target)
docker/build.sh test
```
