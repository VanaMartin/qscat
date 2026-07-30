# TD dissociative-attachment (σ_DA) route via nuclear-axis extractors — Design Spec

**Date:** 2026-07-30
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — spec for review
**Lifecycle:** `qm-method-lifecycle` — extends the time-dependent extractor family
(`qscat.core.time_dependent` + `qscat.core.td_extractors`, merged 2026-07-30) from the electronic
(VE) axis to the nuclear (dissociation) axis, for **dissociative attachment (DA)**. Sub-project
**2 of 3** (SP1 = the VE extractors; **SP3 = TD-DR**, the H₂⁺ ion route, its own later spec —
split out because DR is non-laptop (Docker/MUMPS) and cannot be laptop-validated, whereas DA has a
clean laptop TI oracle). Proving the nuclear extractors on DA against ground truth de-risks SP3.

## Context

SP1 added three interchangeable TD energy extractors — Tannor–Weeks (`TannorWeeks`), delta (`Dirac`),
flow (`Flux`) — on a recorder+transform `Extractor` protocol, driven propagate-once/record-all, for
the vibrational-excitation (VE) cross section. All three are ELECTRONIC-axis: the analysis
surface/point/test-packet lives in the electronic coordinate `r` (`grids[0]`), projected onto the
vibrational state χ_{v'}(R), with an electronic outgoing wave (μ_e=1, l=`model.ell`). They converge
to the exact TI VE oracle to ~3%.

Dissociation (dissociative attachment DA, dissociative recombination DR) is the natural home of the
**flux** method: the observable is the outgoing NUCLEAR current — the nuclei separating on the anion
curve while the electron stays bound to a fragment. eMoScat's `Flux`/`Dirac`/`TestFunction2d` each
carry an axis flag (`'x'` electronic / `'y'` nuclear); qModeling has only the electronic axis so far.
The exact TI oracles already exist: `qscat.core.dissociation.da_cross_section` (neutral F₂/NO) and
`dr_cross_section` (the H₂⁺ ion — Coulomb incident + a Rydberg exit series).

## Goal

Extend the extractor family to the **nuclear axis** and build the TD dissociative-attachment route:
propagate the incident wavepacket under `H_2D` and extract σ_DA via the outgoing nuclear flux (and
delta, and TW, for three-way cross-validation), validated against the TI `da_cross_section` oracle
(F₂/NO). Establishes the TD dissociation capability + the nuclear-axis extractor infrastructure that
SP3 (TD-DR) reuses. (The H₂⁺ ion — Coulomb outgoing + Rydberg series + Docker/MUMPS — is SP3.)

## Architecture — nuclear-axis generalization of the Extractor family

Each extractor gains an `axis`/`coordinate` selector ("electronic" | "nuclear"). The transform MATH
is unchanged per method (Wronskian flux; point-Hankel delta; eta-deconvolution TW); the axis swaps
four quantities:

| quantity | electronic (VE, existing) | **nuclear (dissociation, new)** |
|---|---|---|
| surface/point/test-packet coordinate | `r` = `tgrid.grids[0]` | `R` = `tgrid.grids[1]` |
| projected bound state | vibrational χ_{v'}(R) | **electronic anion state φ(r)** (`anion_electronic_states`) |
| outgoing wave (mass, l) | μ_e = 1, l = `model.ell` | **μ_R = `model.mu`, l = 0** (s-wave dissociation) |
| threshold / channels | eps[v'] (vibrational) | **eps_e** / the anion dissociation channels |
| flux prefactor μ | 1 (electronic) | **μ_R = `model.mu`** |
| DVR derivative grid (flux) | `grids[0]` | `grids[1]` |

Implementation is via an `axis` parameter on the existing `Flux`/`Dirac`/`TannorWeeks` (mirrors
eMoScat's single-class-with-axis-flag; DRY), NOT parallel new classes. The electronic path stays
BYTE-IDENTICAL (guarded by SP1's VE tests + the byte-identical-TW golden test).

**No elastic free-reference.** DA/DR are pure rearrangement channels (σ = 4π³|T|²/2E, no `−δ_{v,v'}`
elastic subtraction), so the `free=`/`subtract_free_reference` machinery is not used on the nuclear
axis. (SP1's tech-debt item — lifting `free=` into the `Extractor` protocol — is done as part of
generalizing the family cleanly.)

**The propagation is REUSED.** For DA the incident state + propagation under `H_2D` is the SAME as
TD-VE (electron in → resonance → nuclei dissociate); only the extraction differs. So `propagate` is
unchanged: build the three extractors on the nuclear axis and read σ_DA off the same trajectory. For
DR the incident is a Coulomb electron wavepacket (charge=−1); propagate; extract nuclear flux with a
Coulomb outgoing wave.

## The nuclear extractors (what each records/transforms)

Per open dissociation channel `c` (anion state φ_c, threshold eps_e,c), `E_DR = E_tot − eps_e,c`,
`K_R = √(2 μ_R E_DR)`, `ifc = eta_incident` (the incident deconvolution, computed on the ELECTRONIC
axis for the incident electron — unchanged):

- **flow** (`Flux`, axis="nuclear"): records `b_c(t) = ⟨φ_c | Ψ(·, R=surface)⟩` (value) and
  `d_c(t) = ⟨φ_c | ∂_R Ψ(·, R=surface)⟩` (nuclear-normal derivative, via
  `dvr_first_derivative_at_node` on `grids[1]`), both projecting onto the electronic anion state φ_c.
  Transform: `S_c = −i/(2 μ_R ifc) Σ_j w_j [conj(φ_out) d_c − b_c conj(φ_out')] e^{iE_tot t} dt` with
  `(φ_out, φ_out') = outgoing_surface_wave(g_R, R_surface, K_R, l=0, charge)`.
- **delta** (`Dirac`, axis="nuclear"): records `⟨φ_c | Ψ(·, R=position)⟩` at a nuclear point;
  transform = the point-Hankel form with `hankel_point_value(g_R, R_position, K_R, l=0, charge)`.
- **TW** (`TannorWeeks`, axis="nuclear"): a nuclear outgoing test packet `g_out(R) φ_c(r)`;
  transform = eta-deconvolution with NUCLEAR `eta_incident`/`eta_outgoing`/`outgoing_channel`
  analogs (projecting on `grids[1]`, mass μ_R, l=0) — the one genuinely new correlation-helper set.

σ = 4π³ |Σ_c T_c|² / (2E) style dissociation cross section, matching the TI `da_cross_section`
normalization (the plan pins the exact prefactor from `dissociation.py` + eMoScat).

## TD-DA (F₂/NO, laptop)

`td_da_cross_section(tgrid, model, eps, chi, v_init, E, *, dt, n_steps, wp_in, method=..., ...)` +
`td_da_cross_sections_all(...)` (one propagation → {tw, delta, flow} σ_DA). Reuses the VE
propagation. Validate: three-way agreement + convergence to the TI `da_cross_section` oracle for
F₂/NO (the molecules whose DA channel is genuinely open, on their per-molecule nuclear grids —
`validation/diatomic` / `MoleculeConfig.da_grid()`).

## Deliverables

- **D1** the nuclear-axis generalization of `Flux` (axis param; project onto anion φ; surface in R;
  outgoing μ_R/l=0; DVR derivative on `grids[1]`) + `free=` lifted into the `Extractor` protocol
  (needed to extend the family cleanly).
- **D2** the nuclear-axis `Dirac` + `TannorWeeks` (the latter incl. the nuclear
  `eta`/`outgoing_channel` helpers).
- **D3** `td_da_cross_section(method="tw"|"delta"|"flow")` + `td_da_cross_sections_all` (one
  propagation → three σ_DA); the F₂/NO three-way + TI-oracle validation + docs.

## Validation

- Nuclear extractors: three-way agreement (flux/delta/TW-nuclear from ONE propagation) + each
  converges to the TI `da_cross_section` for F₂/NO; the ~15-25% cross-method spread at
  under-converged grids is the documented convergence diagnostic (as in SP1). F₂/NO use their
  per-molecule nuclear grids (`validation/diatomic` / `MoleculeConfig.da_grid()`, the resolved DA
  grids), since the exact-2D DA needs them (per the resonance-aware / per-molecule-discretisation
  findings).
- The electronic (VE) path stays byte-identical (SP1 tests + golden test still pass).
- `uv run pytest -q -m "not slow"` passes; `uv run mypy libs/qscat/qscat` 0; `uv run ruff check .`
  clean; `qscat.core` still imports no model/projects at runtime.

## Out of scope

- **SP3 — the TD-DR route** (H₂⁺ ion): a Coulomb incident wavepacket (charge=−1) + nuclear flux with
  the **Coulomb outgoing wave** (`coulomb_h1_en` + its derivative — clears the SP1 charged-Coulomb
  `dphi_out` tech-debt, tested there) + the Rydberg exit-series loop (mirroring TI `dr_cross_section`).
  ~1.15M unknowns → Docker/MUMPS, 2-D convergence Docker-deferred. Its own spec; reuses SP2's
  nuclear extractors. The SP2 nuclear `Flux`/`Dirac`/`TW` should keep the `charge` plumbing so SP3 is
  mostly the Coulomb outgoing wave + the Rydberg loop, not a re-architecture.
- **NO/F₂ DR** and **other ions**.
- **Rust optimization** of the transforms.

## Build order (sequenced)

nuclear Flux + `free=` protocol lift → nuclear Dirac + nuclear TW → `td_da_cross_section` +
`td_da_cross_sections_all` + the F₂/NO three-way + TI-oracle validation + docs. Each step ends with a
laptop-testable deliverable.
