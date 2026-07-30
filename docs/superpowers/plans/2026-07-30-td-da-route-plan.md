# TD dissociative-attachment route (σ_DA) via nuclear-axis extractors — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generalize the three TD energy extractors (`Flux`/`Dirac`/`TannorWeeks`) to the NUCLEAR axis and build `td_da_cross_section` (+ the one-propagation all-three helper), extracting σ_DA from the outgoing nuclear dissociation flux — validated three-way + against the TI `da_cross_section` oracle for F₂/NO.

**Architecture:** Each extractor gains an `axis="electronic"|"nuclear"` selector. The nuclear axis swaps: the surface/point/test-packet coordinate (`r`→`R`, `grids[0]`→`grids[1]`), the projected bound state (vibrational χ_{v'}(R) → **electronic anion state φ(r)** via `anion_electronic_states`), the outgoing wave (μ_e=1,l=`model.ell` → **μ_R=`model.mu`, l=0**), the threshold (eps[v'] → eps_e), and the σ form (VE `π|S−δ|²/2E` → **DA T-matrix `4π³|T|²/2E`**, no elastic subtraction). The transform *math* per method is otherwise unchanged. The electronic (VE) path stays BYTE-IDENTICAL. The DA propagation REUSES the VE propagation.

**Tech Stack:** Python ≥3.12, NumPy/SciPy, `qscat.core` (`td_extractors`, `time_dependent`, `correlation`, `dissociation` — `anion_electronic_states`/`v_dr_diag`/`da_cross_section` the TI oracle, `driven`), `qscat.dvr` (`dvr_first_derivative_at_node`), `qscat.special` (`riccati_bessel_en_mass` + a new `riccati_hankel_en_mass`). Models: **F₂/NO** (open DA channel + TI oracle). pytest; mypy --strict; ruff.

## Global Constraints

- **Atomic units.**
- **`qscat.core` never imports `qscat.model`/`projects` at runtime** (`test_core_no_model_import.py`). Nuclear extractors use `anion_electronic_states` (in `qscat.core.dissociation`, same layer — fine) + the `ResonanceModel` protocol; no model import.
- **The electronic (VE) path is BYTE-IDENTICAL:** `axis="electronic"` (the default) reproduces today's `Flux`/`Dirac`/`TannorWeeks` + `td_ve_cross_section(method=...)` bit-for-bit (SP1's golden + differential tests + the byte-identical-TW golden must still pass).
- **The σ_DA normalization is pinned by the TI oracle:** the ground-truth check is that each nuclear extractor's σ_DA CONVERGES to `qscat.core.dissociation.da_cross_section` for F₂/NO. Port-scout eMoScat's DA driver for the exact flux-S→σ_DA prefactor; the TI-convergence test is the arbiter (as SP1's TI-VE convergence pinned delta/flux).
- **Ported from eMoScat** (`reference/eMoScat/source/Model2d/{Dirac,Flux,Test}Function2d.cpp`, `axis_=='y'` branches; the DA driver in `source/simpleCoupling.cpp` / `Model2d/TimeDependentModel2d`). Run `port-scout` to confirm the nuclear `<<`/`contribution` + the σ_DA assembly before implementing each. Never build/import eMoScat.
- mypy --strict over `libs/qscat/qscat` clean; ruff clean.

## Extracted references (the reconciliation target)

TI `da_cross_section` (`dissociation.py:165-185`), per open channel `n`: `eps_e, phi = anion_electronic_states(g_r, model, R_inf=g_R.R0, n_states)`; `e_dr = e_tot − eps_e[n]`; `k_r = √(2μ e_dr)`; channel function `Φ_n = φ[n] ⊗ (riccati_bessel_en_mass(R, k_r, 0, μ)·√w_R)` masked; `T_n = c_product(Φ_n, V_DR·Ψ₊)`; `σ_n = 4π³|T_n|²/(2E)`. eMoScat nuclear flux (`FluxTestFunction2d.cpp`, `axis=='y'`): records `line_projection(bound_state=φ, 'x', position)` (project onto the electronic anion φ at nuclear node `position`) + its derivative; outgoing `sphHankel1En(yz(position), K_R, μ_R, l=0)/2`. So the TD nuclear extractor projects onto the SAME anion φ, at a nuclear surface/point, with a mass-μ_R, l=0 outgoing wave.

## File Structure

- `libs/qscat/qscat/special/radial.py` (modify) — add `riccati_hankel_en_mass(r, k, l, mu)` (the outgoing/Hankel analog of the existing `riccati_bessel_en_mass`).
- `libs/qscat/qscat/core/correlation.py` (modify) — `hankel_point_value`/`outgoing_surface_wave` gain a `mass` param (default 1.0 = electronic; μ_R for nuclear); + nuclear `eta`/`outgoing_channel` analogs for TW (or an `axis` param on the existing ones).
- `libs/qscat/qscat/core/td_extractors.py` (modify) — `axis` param on `Flux`/`Dirac`/`TannorWeeks`; the nuclear record/transform; the DA σ form; `free=` lifted into the `Extractor` protocol.
- `libs/qscat/qscat/core/time_dependent.py` (modify) — `Extractor` protocol gains `free=`; `td_da_cross_section(..., method=...)` + `td_da_cross_sections_all`.
- Tests: `libs/qscat/tests/test_td_extractors.py`, `test_correlation.py`, `test_radial.py` (or the special test file); `validation/diatomic/test_td_da.py` + `validation/diatomic/td_da.py`.

---

### Task 1: Lift `free=` into the `Extractor` protocol + `axis` scaffolding (VE byte-identical)

**Files:** modify `time_dependent.py` (protocol), `td_extractors.py`; test `test_td_extractors.py`.

**Interfaces:**
- `Extractor` protocol `sigma(self, E, *, free: "Extractor | None" = None) -> NDArray` (lift the class-typed `free` into the protocol — SP1 whole-branch-review tech-debt).
- `Flux`/`Dirac`/`TannorWeeks` gain a keyword-only `axis: str = "electronic"`. For Task 1, `axis="electronic"` is the ONLY implemented path (unchanged behavior); `axis="nuclear"` raises `NotImplementedError` (filled in Tasks 2-4). Store `axis`; branch the coordinate selection (`grids[0]` vs `grids[1]`) behind a small helper so Tasks 2-4 slot in.

- [ ] **Step 1:** the SP1 byte-identical golden + differential tests still pass unchanged (they exercise `axis="electronic"` by default). Add a test that `axis="nuclear"` currently raises `NotImplementedError` on each extractor.
- [ ] **Step 2-5:** implement the protocol `free=` lift + the `axis` param/scaffolding; run SP1's `test_td_extractors.py` (all pass, byte-identical) + the new NotImplementedError test; mypy/ruff; commit `refactor(core): lift free= into Extractor protocol + axis scaffolding (VE byte-identical)`.

---

### Task 2: Nuclear-axis `Flux` + `riccati_hankel_en_mass` + the DA σ (the primary dissociation extractor)

**Files:** modify `special/radial.py` (`riccati_hankel_en_mass`), `correlation.py` (`outgoing_surface_wave` mass param), `td_extractors.py` (`Flux` nuclear path); tests `test_radial.py`, `test_correlation.py`, `test_td_extractors.py`.

**Interfaces:**
- `special.riccati_hankel_en_mass(r, k, l, mu)` — the mass-μ energy-normalized outgoing Riccati-Hankel (analog of `riccati_bessel_en_mass`); unit-tested: reduces to `riccati_hankel_en` at μ=1, and its real/imag parts are consistent with `riccati_bessel_en_mass` (regular) + the Wronskian.
- `outgoing_surface_wave(grid, z_surface, k, l, charge, *, mass=1.0)` — gains `mass`; nuclear uses `riccati_hankel_en_mass(..., mass)` + its derivative. (`hankel_point_value` similarly gains `mass` in Task 3.)
- `Flux(..., axis="nuclear")`: constructed with `(tgrid, model, eps, chi, v_init, channels, surface, *, wp_in, dt, axis="nuclear")` where `channels` is the number of anion dissociation channels (`n_channels`) and `surface` is a NUCLEAR node index. `__init__` computes `eps_e, phi = anion_electronic_states(tgrid.grids[0], model, R_inf=tgrid.grids[1].R0, n_states=channels)`. `record`: per channel `c`, `b_c = ⟨φ_c | Ψ(·, R=surface)⟩` (project onto the ELECTRONIC anion φ_c, at nuclear node `surface`, via a c_product on the electronic axis) and `d_c = ⟨φ_c | (dvr_first_derivative_at_node(g_R, surface) applied on the NUCLEAR axis) Ψ⟩`. `sigma(E)`: `K_R=√(2μ_R(E_tot−eps_e,c))`, `(φ_out,φ_out')=outgoing_surface_wave(g_R, R_surface, K_R, 0, model.charge, mass=μ_R)`, Wronskian `S_c = −i/(2μ_R·ifc)·Σ_j w_j[conj(φ_out)d_c − b_c conj(φ_out')]e^{iE_tot t}dt`, then the DA σ: `σ = C_DA·|Σ_c S_c|²/(2E)` with `C_DA` the prefactor pinned to the TI oracle (port-scout eMoScat's DA driver; expect the `4π³`-family constant reconciling the Hankel-flux S to the TI regular-Bessel T). `ifc = eta_incident` on the ELECTRONIC incident axis (the incident electron — unchanged).

**Design notes:** port-scout `FluxTestFunction2d.cpp` `axis=='y'` (record) + constructor (nuclear `phi_out`) + the DA driver (σ assembly). The mass-μ_R outgoing wave + l=0 is the nuclear dissociation wave. `dvr_first_derivative_at_node` already works on any real grid — apply it on `grids[1]` (nuclear) here (vs `grids[0]` for VE). The anion φ_c projection is a c_product on the electronic axis at fixed nuclear node (transpose of the VE case). The `C_DA` prefactor is the one genuinely uncertain constant — resolve it by the TI-convergence test.

- [ ] **Step 1a:** `test_radial.py` — `riccati_hankel_en_mass` reduces to `riccati_hankel_en` at μ=1 (rtol 1e-12) + Wronskian check vs `riccati_bessel_en_mass`. Run→fail→implement→pass.
- [ ] **Step 1b (the LOAD-BEARING gate, `@slow`):** nuclear-`Flux` σ_DA converges to the TI `da_cross_section` for F₂ at an anchor. Harness: F₂ `MoleculeConfig.da_grid()` (`validation/diatomic/config.py`), `vibrational_states` for `eps`/`chi`; build `Flux(axis="nuclear")`, propagate (REUSE the VE propagation — same `propagate(tgrid, psi0, [], ..., extractors=[flux])`), `sigma_flux_da = flux.sigma(E)`; assert within ~15% of `da_cross_section(tgrid, F2, eps, chi, v_init, E)` at that anchor. This pins `C_DA`. (~minutes; run foreground, patient, do NOT background/Monitor.)
- [ ] **Steps 2-5:** implement; run→pass (tune `C_DA` against the TI oracle if the port-scouted constant needs a documented reconciliation factor — REPORT it, don't fudge); confirm VE byte-identical unchanged; mypy/ruff; commit `feat(core): nuclear-axis Flux + riccati_hankel_en_mass + DA cross section`.

---

### Task 3: Nuclear-axis `Dirac` (delta)

**Files:** modify `correlation.py` (`hankel_point_value` mass param), `td_extractors.py` (`Dirac` nuclear path); test `test_td_extractors.py`.

**Interfaces:** `hankel_point_value(grid, z_position, k, l, charge, *, mass=1.0)`. `Dirac(..., axis="nuclear")`: records `⟨φ_c|Ψ(·, R=position)⟩` at a nuclear point; `sigma` = the point-Hankel DA transform with `hankel_point_value(g_R, R_position, K_R, 0, charge, mass=μ_R)` + the same `C_DA` DA σ as Task 2.

- [ ] **Steps:** port-scout `DiracTestFunction2d.cpp` `axis=='y'`; differential test — nuclear-`Dirac` σ_DA agrees with nuclear-`Flux` from ONE propagation (`extractors=[flux, dirac]`, documented band ~0.25) + converges to TI `da_cross_section` at an anchor (`@slow`); implement; VE byte-identical unchanged; mypy/ruff; commit `feat(core): nuclear-axis Dirac (delta) DA extractor`.

---

### Task 4: Nuclear-axis `TannorWeeks` + nuclear `eta`/`outgoing_channel`

**Files:** modify `correlation.py` (nuclear `eta_incident`/`eta_outgoing`/`outgoing_channel` analogs, or `axis` params), `td_extractors.py` (`TannorWeeks` nuclear path); tests `test_correlation.py`, `test_td_extractors.py`.

**Interfaces:** nuclear analogs projecting on `grids[1]`, mass μ_R, l=0, against the anion φ_c (the "nuclear channel function" `g_out(R)·φ_c(r)`): `outgoing_channel_nuclear`/`eta_incident`/`eta_outgoing` (or the existing ones gain `axis`/`mass`). `TannorWeeks(..., axis="nuclear")`: records `c_c(t)=c_product(g_out(R)φ_c(r), Ψ)`; `sigma` = the eta-deconvolution DA transform + `C_DA`.

**Design notes:** this is the most new correlation code (nuclear TW test packet + eta). The incident `eta_in` stays electronic; the outgoing `eta_out` is nuclear (mass μ_R, l=0, against φ_c). port-scout `TestFunction2d.cpp` `axis=='y'`.

- [ ] **Steps:** port-scout; differential test — nuclear-`TannorWeeks` σ_DA agrees with nuclear-`Flux`/`Dirac` from ONE propagation + converges to TI at an anchor (`@slow`); implement the nuclear correlation helpers + the nuclear TW; VE byte-identical unchanged; mypy/ruff; commit `feat(core): nuclear-axis TannorWeeks + nuclear eta/outgoing_channel`.

---

### Task 5: `td_da_cross_section(method=)` + all-three helper + F₂/NO validation + docs

**Files:** modify `time_dependent.py` (`td_da_cross_section`, `td_da_cross_sections_all`); create `validation/diatomic/td_da.py` + `validation/diatomic/test_td_da.py`; modify `docs/physics/td-extractors.md` (or a new `td-da.md`), `CLAUDE.md`; export from `qscat.core`.

**Interfaces:**
- `td_da_cross_section(tgrid, model, eps, chi, v_init, E, *, dt, n_steps, wp_in, method="tw"|"delta"|"flow", position=None, surface=None, n_channels=1, order=3) -> NDArray` — builds the nuclear extractor for `method`, propagates (reusing the VE propagation), returns its σ_DA. (No `wp_out`/free-reference — DA is pure rearrangement.)
- `td_da_cross_sections_all(...) -> dict[str, NDArray]` — ONE propagation, `{"tw":σ_DA, "delta":σ_DA, "flow":σ_DA}`.

**Design notes:** the incident `psi0 = initial_state(tgrid, chi[v_init], **wp_in)` + `propagate` are IDENTICAL to VE; only the extractors differ. Validation harness: `td_da_cross_sections_all` on F₂/NO (their `da_grid()`) at the anchors where `da_cross_section` is the oracle — three-way agreement (documented band) + each converges to the TI σ_DA. Fast gate: three-way + TI at one reduced anchor; `@slow`: the F₂/NO multi-anchor comparison + a committed figure (σ_DA-vs-TI per method). Honest: the cross-method spread at under-converged grids is a convergence diagnostic (as SP1).

- [ ] **Step 1:** fast gate — `td_da_cross_sections_all` at one reduced F₂ anchor: three methods mutually agree + within band of `da_cross_section`. Run→fail→wire→pass.
- [ ] **Step 2:** the `@slow` F₂/NO validation (`validation/diatomic/td_da.py` + `test_td_da.py`) + the committed figure `docs/physics/figures/{f2,no}-td-da-three-way.png` (or a documented partial if budget-limited, honestly). 
- [ ] **Step 3:** docs (`td-da.md` or a section: the nuclear-axis generalization, the three nuclear extractors, the anion-φ projection, the DA σ + `C_DA` reconciliation to TI, the three-way + ~3%-to-TI result, that SP3 = TD-DR reuses this) + `CLAUDE.md` (`qscat.core` entry: `td_da_cross_section`/`td_da_cross_sections_all`, `axis=` extractors, `riccati_hankel_en_mass`).
- [ ] **Step 4:** mypy/ruff; commit `feat(core): td_da_cross_section (method-selectable nuclear extractors) + F2/NO validation + docs`.

---

## Verification (whole sub-project)

- `uv run pytest -q -m "not slow"` passes; SP1's VE tests + the byte-identical-TW golden still pass (electronic path byte-identical); `test_core_no_model_import.py` passes.
- Each nuclear extractor's σ_DA converges to the TI `da_cross_section` for F₂/NO (the `C_DA` normalization pinned); three-way agreement documented; the comparison figure committed.
- `riccati_hankel_en_mass` validated (μ=1 reduction + Wronskian). `uv run mypy libs/qscat/qscat` 0; `uv run ruff check .` clean.

## Out of scope (this plan)

- **SP3 — TD-DR (H₂⁺ ion)**: Coulomb incident + Coulomb outgoing (`coulomb_h1_en` + its derivative — clears the SP1 charged-Coulomb tech-debt) + Rydberg exit loop; Docker/MUMPS; its own spec. (Keep the `charge` plumbing in the nuclear extractors so SP3 is additive.)
- **NO/F₂ DR**, other ions; **Rust optimization**.
