# Alternative TD extractors (delta + flow) for VE — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the **delta** (`DiracTestFunction2d`) and **flow** (`FluxTestFunction2d`) energy extractors as siblings of Tannor–Weeks (TW) for the TD VE cross section, via a recorder+transform `Extractor` protocol that lets one wavepacket propagation feed all three — then validate three-way agreement + convergence to the TI oracle + Houfek (N₂).

**Architecture:** Refactor `qscat.core.time_dependent` so the propagation engine and the extraction are separated by an `Extractor` protocol (`record(psi)` per step + `sigma(E)` transform), mirroring eMoScat's `operator<<`/`contribution`. `propagate(..., extractors=[tw, delta, flow])` runs ONE trajectory; each extractor accumulates its own per-`v'` series and transforms it. TW becomes one `Extractor` — behavior-preserving. Delta and flow share TW's half-Fourier machinery (same `E_tot=E+eps[v_init]` phase, Simpson quadrature `w_n`, incident factor `eta_in`, threshold gating `E−eps_shift>0`, `σ=π|S−δ|²/2E`); they differ only in what is recorded and the outgoing factor.

**Tech Stack:** Python ≥3.12, NumPy/SciPy, `qscat.core` (`time_dependent`, `correlation`, `channels`, `driven` the TI oracle), `qscat.dvr`, `qscat.special` (`riccati_bessel_en`, coulomb), `qscat.evolution` (Padé stepper). pytest; mypy --strict over `libs/qscat/qscat`; ruff. Model: **N₂** (has TI oracle + Houfek data).

## Global Constraints

- **Atomic units** throughout.
- **`qscat.core` never imports `qscat.model`/`projects` at runtime** (enforced by `test_core_no_model_import.py`) — depend only on the `ResonanceModel` protocol.
- **TW is preserved byte-identical:** `td_ve_cross_section(..., method="tw")` (the new default) returns σ bit-for-bit equal to today's `td_ve_cross_section`. The elastic free-reference subtraction (`subtract_free_reference`, the `td-elastic-wavepacket-normalization` fix) and the order-3 Padé propagation are unchanged. Gated by `projects/n2_2d_td_cross_section/test_td_cross_section.py` (the TD-vs-TI tests) + a golden-value regression.
- **The delta/flow math is ported from eMoScat** (`reference/eMoScat/source/Model2d/{DiracTestFunction2d,FluxTestFunction2d}.cpp`). Per `qm-method-lifecycle`, run the `port-scout` agent to confirm the line-level formulas before implementing each port (the formulas below are already extracted for the plan; port-scout verifies them). Never import/build eMoScat.
- mypy --strict over `libs/qscat/qscat` clean; ruff clean.

## Extracted eMoScat math (the three extractors)

All three (per open channel `i` = energy `E_i`, outgoing `v'`, `E_tot = E_i + eps[v_init]`, quadrature weights `w_j`, incident factor `ifc[i] = eta_incident(E_i)`):

- **TW** (`TestFunction2d`, current): record `c_{v'}(t_n) = c_product(outgoing_channel_{v'}, Ψ(t_n))`.
  `S_i = 1/(2π·conj(eta_out_i)·ifc_i) · Σ_j w_j·c_{v'}(t_j)·e^{iE_tot t_j}·dt`.
- **delta** (`DiracTestFunction2d`): record `b_{v'}(t_n) = ⟨χ_{v'} | Ψ(position, ·; t_n)⟩` (a LINE projection of Ψ onto the bound vibrational state `χ_{v'}` in the coordinate perpendicular to the channel axis, at the fixed grid index `position` — snapped to an element end). Outgoing factor is the outgoing-Hankel VALUE at the point: `f_i = conj( sphHankel1En(z(position), k_i, μ, l) / 2 )` (`coulomb::sH1_en` for a charged residual). `S_i = 1/(2π·conj(f_i)·ifc_i) · Σ_j w_j·b_{v'}(t_j)·e^{iE_tot t_j}·dt`. (Identical to TW with `eta_out → f`.) `k_i = √(2μ(E_i − eps_shift))`.
- **flow** (`FluxTestFunction2d`): record BOTH `b_{v'}(t_n) = ⟨χ_{v'}|Ψ(surface;t_n)⟩` AND `d_{v'}(t_n) = ⟨χ_{v'}|∂_axis Ψ(surface;t_n)⟩` (value and normal derivative at a dividing surface `position`, an element border). `S_i = −i/(2μ·ifc_i) · Σ_j w_j·[ conj(φ_out_i)·d_{v'}(t_j) − b_{v'}(t_j)·conj(φ'_out_i) ]·e^{iE_tot t_j}·dt`, where `φ_out_i = sphHankel1En(z(surface), k_i, μ, l)/2` and `φ'_out_i` its spatial derivative at the surface (the Wronskian flux).

`σ_{v→v'}(E_i) = π·|S_i − δ_{v,v'}·S_ref_i|² / (2 E_i)`, with `S_ref` the free-reference for the elastic channel (unchanged from TW).

## File Structure

- `libs/qscat/qscat/core/time_dependent.py` (modify) — propagation engine takes `extractors`; the `Extractor` protocol; `td_ve_cross_section(..., method="tw"|"delta"|"flow")`; `td_ve_cross_sections_all(...)` (one propagation → dict of three σ). `PropagationResult` keeps diagnostics (t, norm, snapshots).
- `libs/qscat/qscat/core/td_extractors.py` (create) — `TannorWeeks`, `Dirac`, `Flux` extractor classes.
- `libs/qscat/qscat/core/correlation.py` (modify) — add `hankel_point_value` (delta's `f_i`) and `outgoing_surface_wave` (flux's `φ_out`, `φ'_out`) helpers, reusing the same `riccati_bessel_en`/coulomb special functions as `eta_outgoing`.
- `libs/qscat/qscat/dvr/derivative.py` (create) — `dvr_first_derivative_at_node(grid, node_index)`: the FEM-DVR first-derivative row at a real node (the flux surface). Exported from `qscat.dvr`.
- Tests: `libs/qscat/tests/test_td_extractors.py`, `test_dvr_derivative.py`; extend `test_correlation.py`.
- `validation/n2/` — `td_extractors.py` (three-way comparison harness) + a committed accuracy/cost figure; a gated anchor check.

---

### Task 1: `Extractor` protocol + propagation refactor (TW preserved byte-identical)

**Files:**
- Modify: `libs/qscat/qscat/core/time_dependent.py`
- Create: `libs/qscat/qscat/core/td_extractors.py`
- Modify: `libs/qscat/qscat/core/__init__.py` (exports)
- Test: `libs/qscat/tests/test_td_extractors.py`; the gate `projects/n2_2d_td_cross_section/test_td_cross_section.py` must still pass.

**Interfaces:**
- `Extractor` (Protocol, in `time_dependent`): `record(self, psi: NDArray[complex128]) -> None` (accumulate this step's per-`v'` datum); `sigma(self, E: NDArray[float64]) -> NDArray[float64]` (transform accumulated series → `(len(E), len(vprimes))`).
- `propagate(psi0, H, *, dt, n_steps, extractors, order=3, sample_period=..., snapshot_times=...) -> PropagationResult` — runs the order-3 Padé trajectory, calls `ex.record(psi)` on every extractor each step; returns diagnostics (`t`, `norm`, `snapshots`). The extractors hold the accumulated series.
- `TannorWeeks(Extractor)` (in `td_extractors`): constructed with `(tgrid, model, eps, chi, v_init, vprimes, wp_out)`; `record` appends `c_product(outgoing_channel_{v'}, psi)` per `v'`; `sigma(E)` is the current `sigma_from_correlations`/`_sigma_one_energy` transform (deconvolve `eta_out`/`eta_in`, elastic free-reference, `σ=π|S−δ|²/2E`) verbatim.
- `td_ve_cross_section(..., method: str = "tw")` unchanged signature + a `method` kwarg; `method="tw"` builds a `TannorWeeks` extractor, propagates, returns its `sigma(E)` — byte-identical to today.

**Design notes:** move the per-step `c[v']` recording out of `_propagate` into `TannorWeeks.record`; the free-reference propagation (`free=True`) becomes a second `TannorWeeks` (or a flag) whose `S_ref` the elastic channel subtracts — preserve the exact current logic. Keep `PropagationResult` for `t`/`norm`/`snapshots` (diagnostics), drop its `c` field (now inside the extractor) OR keep it populated by TW for back-compat — choose the smaller diff that keeps the project test green.

- [ ] **Step 1: Golden-value regression test (captures TODAY's σ before refactor).**

Assemble the small N₂ config by copying the harness in `projects/n2_2d_td_cross_section/test_td_cross_section.py` (the `propagation` fixture + how it builds `tgrid`, `eps`/`chi` via `vibrational_states`, `wp_in`/`wp_out`, `dt`/`n_steps`) — but use a REDUCED grid / `n_steps` so the test is fast. `qscat.core.driven.ve_cross_section` on the same grid is the TI oracle. Build via `qscat.model.N2` (validation/tests MAY import model; core may not).
```python
# libs/qscat/tests/test_td_extractors.py
# small, fast N2 2-D config (reduced grid + n_steps), 2 energies, 2 vprimes.
def test_tw_method_matches_prerefactor_golden():
    sigma = td_ve_cross_section(tgrid, N2, eps, chi, v_init, vprimes, E,
                                dt=dt, n_steps=n_steps, wp_in=wp_in, wp_out=wp_out,
                                method="tw")
    np.testing.assert_allclose(sigma, _GOLDEN_TW, rtol=0, atol=1e-12)
```
Capture `_GOLDEN_TW` by running the CURRENT `td_ve_cross_section` (no `method` kwarg) on that exact config once and pasting the array into the test — it pins today's behavior across the refactor.
- [ ] **Step 2:** run → fails (`method` kwarg / `TannorWeeks` don't exist).
- [ ] **Step 3:** implement the `Extractor` protocol + `propagate(extractors=...)` + `TannorWeeks` (move the TW record/transform verbatim) + `td_ve_cross_section(method="tw")` delegating to it.
- [ ] **Step 4:** run the golden test + `uv run pytest projects/n2_2d_td_cross_section/test_td_cross_section.py -q -m "not slow"` (the TD-vs-TI contract) — BOTH pass. (The `@slow` TD-vs-TI tests are the deeper gate; note they're slow, run if feasible.)
- [ ] **Step 5:** mypy + ruff; commit `refactor(core): Extractor protocol + propagate-once engine; TW is one extractor (byte-identical)`.

---

### Task 2: `Dirac` (delta) extractor

**Files:**
- Modify: `libs/qscat/qscat/core/td_extractors.py` (`Dirac`), `libs/qscat/qscat/core/correlation.py` (`hankel_point_value`)
- Modify: `libs/qscat/qscat/core/__init__.py`
- Test: `libs/qscat/tests/test_td_extractors.py`, `test_correlation.py`

**Interfaces:**
- `correlation.hankel_point_value(grid, z_position, k, l, charge) -> complex` — `sphHankel1En(z_position, k, μ, l)/2` (neutral) or `coulomb.sH1_en(...)` (charged): the outgoing-Hankel VALUE at the analysis point. Reuses the SAME special functions as `_outgoing_coeffs`.
- `Dirac(Extractor)` — constructed with `(tgrid, model, eps, chi, v_init, vprimes, position)` where `position` is a real electronic-grid index (snapped to an element end). `record`: append `⟨χ_{v'} | psi(position, ·)⟩` per `v'` (line projection of psi onto `χ_{v'}` in R, at electronic index `position`). `sigma(E)`: the TW transform with `eta_out_i → hankel_point_value(...)` at each `E_i`; same `eta_in`, `E_tot` phase, Simpson `w_j`, threshold gate, `σ=π|S−δ|²/2E`, elastic free-reference.

**Design notes:** port-scout `DiracTestFunction2d.cpp` first (`operator<<` = `line_projection(bound_state, perp_axis, position)`; `contribution` = the loop with `Q=1/(2π·conj(f_i)·ifc_i)`). The elastic free-reference: a `Dirac` free (`V_int=0`) propagation gives `S_ref`, subtracted for `v'==v_init` — same pattern as TW. The analysis `position` default: an element end in the asymptotic electronic region (past the interaction, inside the real grid) — mirror TW's `wp_out` standoff; expose it as a parameter.

- [ ] **Step 1: Differential test — delta agrees with TW from ONE propagation.**
```python
def test_delta_agrees_with_tw_same_trajectory():
    # one propagate() with extractors=[TannorWeeks(...), Dirac(...)]
    # at a converged-enough small config, delta sigma ~ TW sigma per open channel.
    res = propagate(psi0, H, dt=..., n_steps=..., extractors=[tw, delta])
    s_tw, s_delta = tw.sigma(E), delta.sigma(E)
    np.testing.assert_allclose(s_delta, s_tw, rtol=_DELTA_TW_RTOL)   # documented cross-method band
```
- [ ] **Step 2:** run → fails. **Step 3:** implement `hankel_point_value` + `Dirac`. **Step 4:** run → passes; add an oracle check `delta.sigma` vs `qscat.core.driven.ve_cross_section` at one anchor (looser rtol, `@slow` if needed). **Step 5:** mypy+ruff; commit `feat(core): Dirac (delta) TD extractor + hankel_point_value`.

---

### Task 3: `Flux` (flow) extractor + DVR first-derivative-at-node

**Files:**
- Create: `libs/qscat/qscat/dvr/derivative.py` (+ export from `qscat/dvr/__init__.py`)
- Modify: `libs/qscat/qscat/core/td_extractors.py` (`Flux`), `libs/qscat/qscat/core/correlation.py` (`outgoing_surface_wave`)
- Test: `libs/qscat/tests/test_dvr_derivative.py`, `test_td_extractors.py`

**Interfaces:**
- `dvr.dvr_first_derivative_at_node(grid, node_index) -> NDArray[complex128]` — the row of the FEM-DVR first-derivative operator at real `node_index` (an element-interior/boundary node): a vector `d` such that `d · psi_coeffs ≈ ∂_x ψ(x_node)`. Built from the element's Gauss–Lobatto Lagrange-derivative matrix (the same basis `kinetic` uses), respecting the ECS Jacobian on the complex tail (real region here — surface is inside the real grid). Test against analytic derivatives of known functions (`sin(kx)`, a Gaussian) on a single-element and multi-element grid to `rtol=1e-8`.
- `correlation.outgoing_surface_wave(grid, z_surface, k, l, charge) -> tuple[complex, complex]` — `(φ_out, φ'_out)` = `sphHankel1En(z_surface,...)/2` and its spatial derivative at the surface (analytic derivative of the Riccati-Bessel/Hankel, or the recurrence).
- `Flux(Extractor)` — constructed with `(tgrid, model, eps, chi, v_init, vprimes, surface)`; `record`: append `b_{v'} = ⟨χ_{v'}|psi(surface)⟩` AND `d_{v'} = ⟨χ_{v'}| (dvr_first_derivative_at_node · psi)(surface)⟩` per `v'`. `sigma(E)`: `S_i = −i/(2μ·ifc_i)·Σ_j w_j·[conj(φ_out_i)·d_{v'}(t_j) − b_{v'}(t_j)·conj(φ'_out_i)]·e^{iE_tot t_j}·dt`; then `σ=π|S−δ|²/2E` + elastic free-reference.

**Design notes:** port-scout `FluxTestFunction2d.cpp` first (`operator<<` fills `projection_` over the element `[index_start_, index_end_]` then `.derivative(position_)`; `contribution` = the Wronskian `conj(φ_out)·d − b·conj(φ'_out)` scaled by `−i/(2μ·ifc)`). The `surface` is an element border in the asymptotic electronic region (where the outgoing electron wave is clean). The DVR derivative is the genuinely new numerical piece — build it against `qscat.dvr`'s element structure (Gauss–Lobatto nodes/weights); `μ` here is the ELECTRONIC reduced mass (=1 in a.u. for the electron), matching eMoScat's `reduced_mass()` on the electronic axis.

- [ ] **Step 1a:** `test_dvr_derivative.py` — `dvr_first_derivative_at_node` matches analytic `∂_x sin(kx)`/Gaussian derivatives to `rtol=1e-8` (single + multi element). Run→fail→implement→pass.
- [ ] **Step 1b:** differential test — `Flux` agrees with TW (and delta) from ONE propagation within the documented band; oracle check vs TI at one anchor.
- [ ] **Steps 2–5:** implement `outgoing_surface_wave` + `Flux`; run→pass; mypy+ruff; commit `feat(core): Flux (flow) TD extractor + DVR first-derivative-at-node`.

---

### Task 4: `method=` selection, all-three helper, N₂ validation + accuracy/cost figure + docs

**Files:**
- Modify: `libs/qscat/qscat/core/time_dependent.py` (`td_ve_cross_section(method=...)` wires delta/flow; `td_ve_cross_sections_all`)
- Create: `validation/n2/td_extractors.py` (three-way comparison + figure), `validation/n2/test_td_extractors.py`
- Modify: `docs/physics/n2-2d-td-cross-section.md` (or a new note), `CLAUDE.md`
- Test: the above.

**Interfaces:**
- `td_ve_cross_section(..., method="tw"|"delta"|"flow")` — full wiring (delta/flow need `position`/`surface` params, default to the asymptotic standoff; document them).
- `td_ve_cross_sections_all(...) -> dict[str, NDArray]` — ONE propagation, returns `{"tw":σ, "delta":σ, "flow":σ}` (the honest three-way comparison — identical dynamics).

**Design notes:** the validation harness runs `td_ve_cross_sections_all` on N₂ across the C5/D1 anchor energies, asserts the three agree within the documented cross-method band AND each converges to `qscat.core.driven.ve_cross_section` (the TI oracle) + matches Houfek at the gated anchors. The accuracy/cost figure: per method, σ-vs-TI error and transform wall-cost (delta expected cheapest; flux robust near threshold). Follow `validation/n2` F1's pattern — a full multi-anchor run is `@slow`/recorded-note; the fast gate is the three-way agreement + oracle at ONE-TWO anchors on a reduced grid.

- [ ] **Step 1:** fast test — `td_ve_cross_sections_all` at one reduced-grid anchor: `tw`/`delta`/`flow` mutually agree within band AND within band of `ve_cross_section` (TI). Run→fail→wire→pass.
- [ ] **Step 2:** the `@slow` multi-anchor N₂ harness + the committed accuracy/cost figure (`docs/physics/figures/n2-td-extractors-comparison.png`) + Houfek anchoring. Run foreground/patient (minutes); if a method departs in a documented regime, RECORD it (a real finding, not a failure).
- [ ] **Step 3:** docs — `docs/physics/n2-2d-td-cross-section.md` (or new `td-extractors.md`): the three extractors, the recorder+transform architecture, the extracted formulas, the three-way agreement + cost/accuracy result, and that flow is the SP2 dissociation workhorse. `CLAUDE.md`: `qscat.core.time_dependent` entry — note the `method=` extractors + `td_extractors`.
- [ ] **Step 4:** mypy+ruff; commit `feat(core): method-selectable TD extractors + N2 three-way validation + docs`.

---

## Verification (whole sub-project)

- `uv run pytest -q -m "not slow"` passes; the TD-vs-TI contract (`projects/n2_2d_td_cross_section`) still passes (TW byte-identical); the golden-TW regression passes.
- `uv run mypy libs/qscat/qscat` 0; `uv run ruff check .` clean; `test_core_no_model_import.py` passes (no model/projects import).
- Three-way agreement (delta/flow vs TW from ONE propagation) + each converges to the TI `ve_cross_section` + matches Houfek at the anchors; the accuracy/cost figure committed.
- `dvr_first_derivative_at_node` validated against analytic derivatives.

## Out of scope (this plan)

- **SP2 — the TD-dissociation route** (σ_DA/σ_DR via outgoing NUCLEAR flux) — separate spec, reuses this `Extractor` + `Flux` infrastructure.
- **NO/F₂ VE** with the new extractors — trivial follow-on once proven on N₂.
- **Rust optimization** of the transforms.
- **Binary-storage/streaming** of buffers (eMoScat's `save_bin`/`opened_` machinery) — not needed in-memory.
