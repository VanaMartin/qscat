# Alternative TD extractors (delta + flow) for the VE cross section — Design Spec

**Date:** 2026-07-29
**Author:** Martin (martin@qscat.com) with Claude
**Status:** Approved design — spec for review
**Lifecycle:** `qm-method-lifecycle` — adds two new energy-extraction methods to the validated
time-dependent VE route (`qscat.core.time_dependent`), as alternatives to Tannor–Weeks. Sub-project
**1 of 2** (SP2 = the TD-dissociation route, a separate spec that reuses this infrastructure).

## Context

The time-dependent VE cross section (`qscat.core.time_dependent.td_ve_cross_section`) propagates an
incident Gaussian wavepacket under `H` and extracts `σ_{v→v'}(E)` by the **Tannor–Weeks** (TW)
energy transform: record the c-product correlation `c_{v'}(t_n)=⟨Φ_{v'}|Ψ(t_n)⟩` every step, then
`S_{v→v'}(E) = [2π·conj(eta_out)·eta_in]^{-1} Σ_n w_n e^{iE_tot t_n} c_{v'}(t_n) dt`,
`σ = π|S − δ|²/2E`. This is the port of eMoScat's `TestFunction2d`.

eMoScat implements **three** interchangeable TD extractors on one shared interface
(`TestFunctionInterface2d`): `TestFunction2d` (TW, ported), `DiracTestFunction2d` (**"delta"**), and
`FluxTestFunction2d` (**"flow"/flux**). The interface is already a **recorder + transform** pair —
`operator<<(psi)` accumulates a per-step series, `contribution(S, …)` flushes it to the S-matrix —
which is exactly the abstraction this sub-project builds in Python.

## Goal

Add the **delta** and **flow** extractors as siblings of TW, sharing the SAME wavepacket
propagation, for three reasons (all selected):
1. **Cross-validate the TD route** — three independent energy transforms of ONE trajectory should
   agree and converge to the TI oracle; strong evidence the TD extraction hides no artifact.
2. **Regime-specific accuracy/cost** — delta is cheaper (no Gaussian deconvolution); flux is the
   natural extractor for outgoing current (and the SP2 dissociation workhorse); TW needs
   test-function placement. Quantify accuracy vs cost per method.
3. **Port fidelity** — reproduce any eMoScat TD result, not just the TW slice.

## Architecture — recorder + transform, propagate once

Refactor `qscat.core.time_dependent` so the **propagation engine** and the **extraction** are
separated by an `Extractor` protocol (mirroring eMoScat's `<<`/`contribution`):

```python
class Extractor(Protocol):
    def record(self, psi: NDArray[complex128], t: float) -> None: ...   # per-step accumulate
    def sigma(self, E: NDArray[float64]) -> NDArray[float64]: ...        # energy-transform -> σ(E)
```

`propagate(psi0, H, dt, n_steps, extractors=[...], order=3)` runs ONE trajectory under the order-3
diagonal-Padé stepper; each step calls `ex.record(psi, t)` on every extractor. Because
cross-validation demands identical dynamics, all three extractors see the SAME `Ψ(t_n)` — no
re-propagation, no 3× cost, no separate-trajectory ambiguity. Each extractor then transforms its own
accumulated series to `σ(E)`.

The existing TW code becomes one `Extractor` implementation — a **behavior-preserving refactor**,
gated by the current TD-vs-TI differential test (`projects/n2_2d_td_cross_section` /
`validation/n2` D1/F1) and the `td-elastic-wavepacket-normalization` fix (the elastic free-reference
subtraction) must be preserved exactly.

## The three extractors

All share the propagation; they differ ONLY in what `record` accumulates and how `sigma` transforms.

- **Tannor–Weeks** (`TannorWeeks`, refactor of the current code): records the c-product
  `c_{v'}(t)=⟨Φ_{v'}|Ψ⟩` (against `correlation.outgoing_channel`); `sigma` deconvolves with
  `eta_incident`/`eta_outgoing` and does the elastic free-reference subtraction. Byte-identical σ to
  today.
- **delta** (`Dirac`, ports `DiracTestFunction2d` — "modified Tannor & Weeks with a
  delta-distribution instead of the test-function wave packet"): records `Ψ` at a fixed analysis
  **position** (a δ in the electronic coordinate) projected onto the outgoing vibrational channel
  `χ_{v'}(R)`; `sigma` is the same half-Fourier sum but with the delta's normalization — the
  free/Riccati-Bessel value at the point replaces the Gaussian `eta` deconvolution. Cheaper (no
  test-packet overlap), one analysis point.
- **flow** (`Flux`, ports `FluxTestFunction2d` — "time-energy Fourier transform of the probability
  flux projected to the outgoing state"): records the outgoing probability **flux** through a
  dividing **surface** (`j = (1/2μ_e or 1/2)·Im[Ψ* ∂_x Ψ]`-type bilinear, evaluated at the surface
  from the stored `phi_out`/`dphi_out` outgoing waves and their derivatives, projected onto
  `χ_{v'}`); `sigma` is the time→energy FT of that flux. Natural for measuring outgoing current;
  the SP2 dissociation extractor.

The EXACT per-step quantities and normalizations are extracted from
`reference/eMoScat/source/Model2d/{DiracTestFunction2d,FluxTestFunction2d}.cpp` during the plan
(via the `port-scout` agent, per `qm-method-lifecycle`) — this spec fixes the architecture and the
validation, not the line-level formulas.

## Layout

- `qscat.core.time_dependent` — keeps the propagation engine + `PropagationResult` + the `Extractor`
  protocol; `td_ve_cross_section(..., method="tw"|"delta"|"flow")` selects one, and a helper (e.g.
  `td_ve_cross_sections_all`) returns all three from a single propagation for the comparison.
- `qscat.core.td_extractors` (new sibling) — the three `Extractor` implementations (`TannorWeeks`,
  `Dirac`, `Flux`), plus any shared analysis-point / surface helpers.
- `qscat.core.correlation` — unchanged (TW's `eta`/`outgoing_channel`); delta/flow add their own
  small helpers (analysis-point value, surface flux) in `td_extractors` or `correlation`.
- `qscat.core` must still not import `qscat.model`/`projects` at runtime (enforced by
  `test_core_no_model_import.py`).

## Validation — N₂ primary (it has every oracle)

Per the `numerical-validation` skill:
1. **Differential (three-way agreement):** delta and flow agree with TW — all three from ONE
   propagation — within a documented cross-method tolerance, per open channel, across the anchor
   energies.
2. **Oracle convergence:** each of the three converges to the TI `qscat.core.driven.ve_cross_section`
   as `dt→0` / `n_steps→∞` (a small convergence study, as the current TW route already has).
3. **Experimental anchor:** all three match **Houfek's N₂ `CSVE.V00.J00`** data at the gated anchors
   (the same anchors `validation/n2` C5/D1 use), within the documented cross-model band.
4. **Accuracy/cost comparison (deliverable):** a committed table/figure — per method: accuracy vs the
   TI oracle, wall-cost of the transform (delta expected cheapest; flux expected robust near
   threshold), and any regime where one method is preferable. This is the concrete output of the
   "regime-specific" goal.

A full TD propagation is minutes-scale, so the live gate follows the existing pattern: the fast
suite asserts the three-way agreement + oracle convergence at ONE or TWO anchors on a reduced grid;
the full multi-anchor comparison is a `@slow`/documented run (like `validation/n2` F1).

## Deliverables

- **D1** the `Extractor` protocol + propagation-engine refactor (TW preserved byte-identical).
- **D2** the `Dirac` (delta) extractor + its transform, ported and validated.
- **D3** the `Flux` (flow) extractor + its transform, ported and validated.
- **D4** `td_ve_cross_section(method=...)` + the all-three helper; the N₂ three-way validation
  (differential + oracle + Houfek) and the accuracy/cost comparison figure/note.

## Out of scope (this sub-project)

- **SP2 — the TD-dissociation route** (σ_DA/σ_DR by propagating under `H` and extracting the
  outgoing NUCLEAR flux): its own spec, reuses this `Extractor` infrastructure; flux is its natural
  extractor. Starts only once delta/flow are proven here.
- **NO/F₂ VE with the new extractors** — a trivial follow-on once proven on N₂ (they lack Houfek
  data, so the exact TI solver is their only oracle); not needed to establish the methods.
- **Rust optimization** of the transforms — deferred until a hot path is proven.

## Validation summary

- `uv run pytest -q -m "not slow"` passes; the existing TD-vs-TI test still passes (TW refactor
  behavior-preserving). `uv run mypy libs/qscat/qscat` 0; `uv run ruff check .` clean.
- Three-way agreement (delta/flow vs TW) + oracle convergence gated at the anchor(s); the full
  comparison + Houfek anchoring committed (figure/note), the multi-anchor run `@slow`.
- `qscat.core` still imports no model/projects at runtime.
