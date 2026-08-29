# TD-NRM Memory Observables — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measure what the nonlocal kernel does during a propagation — where the
amplitude goes, when it comes back, and how badly a single decay rate describes
it — on N₂, F₂ and NO.

**Architecture:** Opt-in per-step recording inside `propagate_nrm`, whose
extended state already carries the auxiliary packets. The Markovian reference is
the local limit of `F(E)` itself, so both sides of every comparison come from one
object and no pole walk is involved.

**Tech Stack:** Python 3.12, numpy/scipy sparse, `qscat.core.nrm`, pytest.

**Spec:** `docs/superpowers/specs/2026-08-26-nrm-memory-observables-design.md`

## Global Constraints

- **Atomic units** (`qscat.units`).
- PRA 47, 1031 (1993) and PRA 77, 012710 (2008) are the sole specification;
  equation references in code name paper and equation. `reference/eMoScat`
  specifies nothing.
- **`n_states=None` for every nonlocal propagation** — truncation makes `H_ext`
  non-dissipative.
- ECS ⇒ the bilinear c-product for anything the model's algebra uses.
  **The observables here deliberately use the CONJUGATING product restricted to
  the real region**, because a population is not a c-product — every such use
  carries a comment saying so, and Task 1 establishes what it is worth.
- DVR coefficient space; Eq. (60) weights already absorbed. `√w` appears here
  only in the coefficient→value conversions of §3 and §4, each named.
- `uv run --no-sync ...`; `ruff check`, **`ruff format --check`** on touched
  files, and `mypy libs/qscat/qscat` clean.
- **The lint rule set has tightened since this plan was drafted.**
  `pyproject.toml` now selects `["E","F","I","UP","B","NPY","RUF","D","ANN"]` —
  `D` (pydocstyle) and `ANN` (annotations) are new. Every function this plan adds
  needs a docstring and complete type annotations; the code sketches below are
  illustrative and do NOT carry them. Check `[tool.ruff.lint] ignore` for the
  handful that are disabled rather than assuming a rule applies.
- Tests over ~10 s get `@pytest.mark.slow`.
- Commit trailers:
  `Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>`
  `Claude-Session: https://claude.ai/code/session_01NxwNdBLUXBampDLrusGrAS`

## File Structure

| File | Responsibility |
|---|---|
| `libs/qscat/qscat/core/nrm/memory.py` | `MemorySpec`, `local_width` (the `√w` row sum), and the per-step recorder. No propagation. |
| `libs/qscat/qscat/core/nrm/propagation.py` | `memory=` argument; the recorder call in the existing loop; the new `TdNrmResult` fields. |
| `libs/qscat/tests/test_nrm_memory.py` | The §4 identity, `local_width` against `Γ`, and the recorder's arithmetic. |
| `validation/diatomic/memory_observables.py` | The three-molecule campaign + figure. |
| `docs/physics/nrm-time-dependent.md` | A new section: what the memory says. |

---

### Task 1: `local_width`, and the identity that decides what may be claimed

**Files:**
- Create: `libs/qscat/qscat/core/nrm/memory.py`
- Test: `libs/qscat/tests/test_nrm_memory.py`

**Interfaces:**
- Produces `local_width(f_matrix, nuclear_grid) -> NDArray[float64]`, the local
  limit `Γ_loc(R)` of `F(E)`.

**This task decides what §2.1 of the spec is allowed to say.** Do the
measurement before the recorder exists, because a null result here changes the
next task's output contract.

- [ ] **Step 1: Write the failing test for `local_width`**

```python
def test_local_width_reproduces_the_eq_68_width_at_the_local_energy(deck):
    """F's local limit IS a width, and it matches the independently computed one.

    The oracle is **Eq. (68)**, `2*pi*|V_dk+(R; eps_loc)|^2` via
    `coupling.v_dk_plus` / `gamma_from_coupling`, at the same LOCAL electron
    energy -- which is what `nonlocal-resonance-model.md` §9's 0.977 / 1.011
    ratios are measured against.

    NOT `local_complex_potential`'s `Gamma`: that is the width at the RESONANCE
    energy, a different quantity, and on the N2 gate deck the ratio between them
    sweeps 0.12 -> 8.9 as `eps_loc` crosses zero at R = 2.43 while `eps_res`
    stays flat. No choice of E makes them comparable.

    Needs a bigger electronic box than the gate deck's: measured medians are
    1.218 at r_max=11 (the gate deck, sized for the algebraic transform identity
    which does not care), 0.996 at r_max=16, 0.997 at 24 and 30. Use r_max=16.
    """
```

Read `nonlocal-resonance-model.md` §9 for the exact definition it records, and
`resonance_pole_walk`'s warning for how to find the frozen radius. Implement:

```python
def local_width(f_matrix, nuclear_grid):
    """Local limit of `F(E)`: `Gamma_loc(R_i) = -2 Im[(F sqrt(w))_i / sqrt(w)_i]`.

    The `sqrt(w)` factors are the coefficient->value conversion. `diag F` alone
    is NOT the local limit -- the kernel spans ~10 nodes and `diag F` is 0.14x
    the local potential (§9).
    """
```

- [ ] **Step 2: Run, confirm it fails, implement, confirm it passes.**

- [ ] **Step 3: Measure the §4 imbalance — the load-bearing measurement**

```python
def test_the_coupling_exchange_balances_to_the_extent_V_dn_is_real(deck):
    """What the coupling removes from Psi_d it adds to the arms -- exactly only
    where `V_dn` is real, which under ECS it is not.

    Sum of the two rates is `4 sum_n Re[conj(Psi_d) phi_n] Im[V_dn]`. This test
    RECORDS that residual against the larger of the two rates; it does not
    assert a physics tolerance, because nobody has measured it before.
    """
```

Propagate a short run on the N₂ gate deck, and report, in the task report:
the imbalance as a fraction of `|exchange|`, its maximum over the run, and
whether it is concentrated in the interaction region or in the ECS tail.

- [ ] **Step 4: Turn the measurement into the gate.** Once measured, assert a
bound with headroom over what you measured — and **write the number into
`memory.py`'s module docstring**, because it is what licenses (or forbids)
calling `‖φ_n‖²` a population.

- [ ] **Step 5: Commit.**

---

### Task 2: the recorder, and `propagate_nrm`'s opt-in

**Files:**
- Modify: `libs/qscat/qscat/core/nrm/memory.py`, `propagation.py`
- Test: `libs/qscat/tests/test_nrm_memory.py`

**Interfaces:**
- `MemorySpec(gamma_local, n_channels=4)`; `propagate_nrm(..., memory=None)`;
  `TdNrmResult` gains `arm_norm`, `arm_norm_by_channel`, `exchange`,
  `exchange_local`, `imbalance`, each `None` when `memory` is None.

- [ ] **Step 1: Write the failing tests.** Three, and none of them is a
smoke test:

```python
def test_memory_off_changes_nothing():
    """The default path is byte-identical to before -- `psi_d` and every
    existing diagnostic, not just 'it still runs'. Every gate in the suite
    propagates, so a regression here is a regression everywhere."""

def test_exchange_matches_a_direct_two_state_calculation():
    """On a 2x2 `H_ext` with one arm, `2 Im<Psi_d|V phi>` is computable in
    closed form. The recorder must reproduce it to round-off -- an arithmetic
    check that does not depend on any physics being right."""

def test_the_markovian_exchange_is_never_positive():
    """`-<Psi_d|Gamma_loc|Psi_d>` with `Gamma_loc >= 0` cannot be a gain. This
    is what makes a POSITIVE nonlocal exchange meaningful rather than a sign
    convention, so it is asserted rather than assumed."""
```

- [ ] **Step 2: Run to confirm failure. Step 3: Implement. Step 4: Confirm.**

The recording block, for reference — the prototype's, which reproduced on the
N₂ deck (`arms peak 0.19 of S_d(0)` at t=25; exchange positive at 85 of 4001
steps from t=132):

```python
d = psi[:n_r]; arms = psi[n_r:].reshape(n_arm, n_r)
dr = d[real]
# CONJUGATING product, restricted to the real region: this is a probability,
# not a c-product. See memory.py's docstring for what Task 1 measured about
# how far that reading can be trusted.
arm_norm[m] = np.einsum("nr,nr->", arms[:, real].conj(), arms[:, real]).real
coup = np.einsum("nr,nr->r", V_dn[:, real], arms[:, real])
exchange[m] = 2.0 * np.vdot(dr, coup).imag
exchange_local[m] = -float(np.vdot(dr, gamma_local[real] * dr).real)
```

**Two properties of `Γ_loc` measured in Task 1 that bind this task:**

- **It is not zero where the channel is closed, and it is not safe to divide by.**
  `ε_loc > 0` only between R = 1.7426 and R = 2.4284 (31 of 153 real nodes,
  because `V_0` is a well). Outside that window `Γ_loc` decays smoothly rather
  than vanishing — 1.41e-4 one node out, ~1e-10 at R = 2.81, ~1e-12 by R = 4,
  and exactly zero at none of the 122 closed nodes. So `exchange_local =
  −⟨Ψ_d|Γ_loc|Ψ_d⟩` needs **no masking** and is well defined everywhere. But any
  RATIO against `Γ_loc` divides by ~1e-10 over most of the grid, which for a
  dissociating packet is where the packet spends its late time. **Compare the two
  exchange curves as a DIFFERENCE, or normalize both by `S_d` — never by each
  other.**
- **Its round-off negatives are deliberately not clamped.** 32 of 179 nodes carry
  a negative `Γ_loc`, worst −3.81e-9 (1.5e-7 of peak, median magnitude 7e-29),
  every one where the true width is zero. Clamping would hide the one case that
  matters — a negative entry meaning the ingredients are wrong. Their would-be
  gain contribution is 5e-44 of `|exchange_local|`, which stays strictly negative
  (max −3.47e-6), so `test_the_markovian_exchange_is_never_positive` is safe as
  written. **Do not add a clamp to make it safer.**

- [ ] **Step 5: Cost.** Measure the per-step overhead with `memory` on and off on
the N₂ gate deck and record both in the report. If it exceeds ~10 %, say so —
it is opt-in, so a real cost is acceptable, but it must be known.

- [ ] **Step 6: Re-point Task 1's imbalance test at the recorder.**
`test_nrm_memory.py` currently keeps its own per-step copy of the
`exch_d`/`exch_arm`/`residual` arithmetic. Once the recorder exists, that test
must call it — otherwise the gate measures the test's private implementation and
the shipped recorder is ungated, which is the same class of gap the Task 1 review
found in the algebraic assertion. Confirm the gate's numbers are unchanged by the
switch; if they move, the recorder and the test disagree and that is the finding.

- [ ] **Step 7: Commit.**

---

### Task 3: the three-molecule campaign

**Files:**
- Create: `validation/diatomic/memory_observables.py` (+ its test)

**Every quoted figure and every axis label states its normalization.** The
`exchange` field ships UNNORMALIZED: its raw positive maximum on the N₂ gate deck
is `+8.776e-7`, and the `+2.420e-4` quoted throughout this sub-project's history
is that number divided by `S_d(0)`. Two independent agents reproduced `2.420e-4`
because both normalized the same way to check against a figure that did not say
it was normalized — which is how a convention becomes invisible. Compare the two
exchange curves as a **difference**, or normalize **both** by `S_d`; never by
each other (see the `Γ_loc` note in Task 2).

**Check `arm_peak` per molecule rather than carrying N₂'s answer across.** On the
N₂ gate deck the first four blocks *are* the four largest (93.2 % of the peak arm
norm). That is a measurement of that deck. If it does not hold on F₂ or NO, that
is a result about where the flux goes, not a nuisance.

**The campaign decks need a converged electronic box.** `local_width` is 22 %
high at `r_max = 11` (the N₂ gate deck, which was sized for the algebraic
transform identity, which does not care) and converges by 14–16: measured medians
1.218 / 0.998 / 0.996 / 0.997 at r_max 11 / 14 / 16 / 24. `exchange_local` would
be 22 % wrong on an r_max=11 deck. Check each molecule's deck before running, and
report the box used.

**O₂ is deliberately excluded** even though `qscat.model.library` now registers
`O2`, `O2_SO12` and `O2_SO32`. Those are the potential factory's FIT, not a
published parameter set, so there is no external LCP comparison for them and the
comparative question below is not defined. N₂/F₂/NO only.

- [ ] **Step 1** — run all three molecules at their existing gate decks and
energies, recording for each: arm-norm peak and its time; exchange min/max, the
number and timing of net-return steps; the decay-law comparison at
`S/S₀ = 0.5 / 0.1 / 0.01`; and the Task 1 imbalance.
- [ ] **Step 2** — the figure: three panels as in the prototype, one column per
molecule. Follow `validation/diatomic/da_figure.py`'s structure and its
`--outdir` convention. **F₂ and NO are production decks — run them on sadaharu,
not locally.**
- [ ] **Step 3** — assert only what is cheap and stable in the test file. The
sign structure of §2.2 (Markovian never positive; nonlocal positive somewhere,
if it is) is the durable claim; exact peak times are not.
- [ ] **Step 4: Commit.**

**The comparative question this exists to answer, stated so the result is not
just three sets of curves:** in the energy domain the three molecules' LCP
failures are ordered N₂ (mild) < F₂ (sweeps through unity) < NO (undetermined).
**Does the return flux reproduce that ordering?** If it does, the exchange rate
is a predictor of LCP failure and that is the result. If it does not, the
ordering has a different cause and that is a more interesting one. Report which,
and do not round either into the other.

---

### Task 4: the physics note

**Files:** `docs/physics/nrm-time-dependent.md`, `CLAUDE.md`

- [ ] **Step 1** — a new section: the three observables, what each measures, the
Task 1 imbalance and what it licenses, and the three-molecule comparison with the
Task 3 verdict.
- [ ] **Step 2** — state the limits from spec §7 explicitly: these are
diagnostics of a model already gated, not evidence for it; NO/A stays unquotable;
the frozen pole walk is routed around, not fixed.
- [ ] **Step 3** — `CLAUDE.md`'s `qscat.core.nrm` entry gains the memory
observables in one clause.
- [ ] **Step 4: Commit.**

## Self-Review

**Spec coverage.** §2.1/§2.2/§2.3 → Tasks 2–3. §3 (`Γ_loc` from `F`) → Task 1.
§4 (the imbalance) → Task 1 Steps 3–4. §5 (API) → Task 2. §6 (three molecules)
→ Task 3. §7 (limits) → Task 4 Step 2.

**Placeholders.** One deliberate: Task 1 Step 4's bound is written only after
Step 3 measures it. Setting it first would be the defect this branch has already
been caught making three times.

**Type consistency.** `gamma_local` is `NDArray[float64]` on the full nuclear
grid throughout; the recorded series are `(n_steps+1, n_E)` except
`arm_norm_by_channel`, which is `(n_steps+1, k, n_E)`; every new
`TdNrmResult` field is `| None`.
